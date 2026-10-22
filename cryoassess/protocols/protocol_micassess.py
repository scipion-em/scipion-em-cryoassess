# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import os
import math
import time
import traceback
from enum import Enum
from emtable import Table

from pyworkflow.constants import PROD
import pyworkflow.protocol.constants as pwcts
from pyworkflow.object import Set, RELATION_SOURCE
from pyworkflow.protocol import params, ProtStreamingBase
from pyworkflow.utils.path import copyTree, cleanPath, makePath, createLink
from pwem.protocols import ProtPreprocessMicrographs
from pwem.objects import SetOfMicrographs

from .. import Plugin
from ..constants import CRYOASSESS_MODELS

OUTPUT = 'outputMicrographs'
OUTPUT_DISCARDED = 'discardedMicrographs'
# Folders micassess fills with the images of each label once its prediction is over
MICASSESS_LABELS = ['0Great', '1Decent', '2Contamination_Aggregate_Crack_Breaking_Drifting', '3Empty_no_ice',
                    '4Crystalline_ice', '5Empty_ice_no_particles_but_vitreous_ice']
BATCH_INPUT_STAR = 'input_micrographs.star'
BATCH_OUTPUT_DIR = 'output'
# Ids per filtered read
ID_BATCH_SIZE = 500
# Polls a closed input may show the same incomplete view, with nothing in flight, before giving up on it
TERMINAL_STALL_POLLS = 10


class outputs(Enum):
    outputMicrographs = SetOfMicrographs
    discardedMicrographs = SetOfMicrographs


def batchMicName(micId, fileName):
    """ Name of a micrograph inside a batch: micassess matches its results by basename, so two micrographs
    with the same basename would be mixed up without the id. """
    return '%06d__%s' % (micId, os.path.basename(fileName))


class CryoassessProtMics(ProtPreprocessMicrographs, ProtStreamingBase):
    """
    Protocol to assess micrographs from K2 or K3 cameras.
    """
    _label = 'assess micrographs'
    _devStatus = PROD
    _possibleOutputs = outputs
    stepsExecutionMode = pwcts.STEPS_PARALLEL

    def __init__(self, **kwargs):
        ProtPreprocessMicrographs.__init__(self, **kwargs)

    @classmethod
    def worksInStreaming(cls):
        return True

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label="Input micrographs", important=True)
        form.addParam('threshold', params.FloatParam, default=0.1,
                      label='Threshold 1',
                      help='Threshold for good/bad classification (step 1). Default is 0.1. '
                           'Higher number will cause more good micrographs '
                           'being classified as bad.')
        form.addParam('threshold2', params.FloatParam, default=0.1,
                      label='Threshold 2',
                      help='Threshold for great/decent classification (step 2). Default is 0.1. '
                           'Higher number will cause more good micrographs '
                           'being classified as bad.')
        form.addParam('batchSize', params.IntParam, default=32,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Batch size',
                      help='Batch size used in prediction. Default is 32. '
                           'Increasing this number will result in faster '
                           'prediction, if your GPU memory allows. '
                           'If memory error/warning appears, you should '
                           'lower this number.')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " Micassess can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")
        form.addParallelSection(threads=3, mpi=0)

        self._defineStreamingParams(form)
        form.getParam('streamingBatchSize').setDefault(5)
        form.getParam('streamingSleepOnWait').setDefault(5)

    # --------------------------- STEPS functions -----------------------------
    def stepsGeneratorStep(self):
        """ Discover the input micrographs as they arrive and assess them in batches. """
        self._prepareStreaming()
        batchSteps = []
        while True:
            if self._streamingMustStop():
                self.info('The run is over, no more micrographs will be assessed.')
                return
            inMics = self._loadInputSet()
            # Read before discovering: once the producer is closed, this poll sees everything it holds
            producerClosed = not inMics.isStreamOpen()
            with self._lock:
                complete = self._discoverNewMics(inMics, producerClosed)
                batchSteps.extend(self._insertBatchSteps(complete))
            inMics.close()
            if complete:
                self._insertFunctionStep(self.closeOutputStep, prerequisites=batchSteps, needsGPU=False)
                return
            self._waitForInput()

    def _waitForInput(self):
        if self._getStreamingSleepOnWait() > 0:
            self._streamingSleepOnWait()
        else:
            time.sleep(1)  # A generator must yield the CPU while it waits

    def assessBatchStep(self, micIds):
        """ Assess the micrographs of micIds and publish them as good or discarded.

        A batch that fails is not raised: the run would fail with it. Its micrographs stay out of the outputs,
        so Continue assesses them again, and the run fails only when every micrograph failed.
        """
        try:
            mics = self._loadMicsByIds(micIds)
            goodIds = self._runMicAssess(mics)
            with self._lock:
                self._publish([m for m in mics if m.getObjId() in goodIds], [m for m in mics if m.getObjId() not in goodIds])
        except Exception as e:
            self.error('Micrographs %s could not be assessed: %s' % (micIds, e))
            self.error(traceback.format_exc())
            with self._lock:
                self._failedIds.update(micIds)
        finally:
            with self._lock:
                self._activeBatches -= 1

    def closeOutputStep(self):
        self._closeOutputSet()
        self._ensureSourceRelations()
        if self._failedIds:
            published = sum(getattr(self, name).getSize() for name in (OUTPUT, OUTPUT_DISCARDED) if hasattr(self, name))
            if not published:
                raise RuntimeError('Micassess failed for every one of the %d micrographs; see the errors above. '
                                   'Continue assesses them again.' % len(self._failedIds))
            self.warning('%d micrographs could not be assessed and were left out; Continue assesses them again.'
                         % len(self._failedIds))

    # --------------------------- STREAMING functions -------------------------
    def _prepareStreaming(self):
        # Everything here is rebuilt on every execution: the durable outputs say what was already assessed
        self._watermark = 0
        self._seenIds = set()
        self._pendingIds = set()
        self._doneIds = set()
        for name in (OUTPUT, OUTPUT_DISCARDED):
            if hasattr(self, name):
                self._doneIds.update(getattr(self, name).getIdSet())
        self._failedIds = set()
        self._activeBatches = 0
        self._stallSignature, self._stallPolls = None, 0
        self._initTotalStars()

    def _streamingMustStop(self):
        """ Once a run is FAILED or ABORTED the executor joins its threads, the generator among them:
        a generator that keeps polling is never joined and the run hangs. """
        return self.getStatus() in (pwcts.STATUS_FAILED, pwcts.STATUS_ABORTED)

    def _loadInputSet(self):
        inMics = self._getInputMicrographs()
        inMics.close()
        inMics.load()
        inMics.loadAllProperties()
        return inMics

    def _discoverNewMics(self, inMics, producerClosed):
        """ Add to the pending ids the input ids not seen yet.

        A poll reads above the watermark only, so it costs what arrived and not the whole history. Once the
        producer is closed, ids that became visible below the watermark are listed and recovered.

        :return: whether the input is complete: its producer is closed and all it declares has been seen.
        """
        self._registerIds(self._idsAfter(inMics, self._watermark))
        if not producerClosed:
            return False
        declared = inMics.getSize()
        if len(self._seenIds) < declared:
            self._registerIds(set(inMics.getUniqueValues('id')) - self._seenIds)
        return self._isTerminalComplete(declared)

    @staticmethod
    def _idsAfter(inMics, watermark):
        try:
            return inMics.getUniqueValues('id', where='id > %d' % watermark)
        except NotImplementedError:
            return [micId for micId in inMics.getUniqueValues('id') if micId > watermark]

    def _registerIds(self, micIds):
        for micId in micIds:
            self._seenIds.add(micId)
            self._watermark = max(self._watermark, micId)
            if micId not in self._doneIds:
                self._pendingIds.add(micId)

    def _isTerminalComplete(self, declared):
        """ A producer can close while some of its rows are not visible yet: closing the output then would
        leave them out for good. When the view stays exactly the same poll after poll with nothing in flight,
        the rows are not coming, and waiting would leave the run RUNNING forever. """
        visible = len(self._seenIds)
        if visible >= declared or self._activeBatches:
            # Work in flight is progress, however long it takes
            self._stallSignature, self._stallPolls = None, 0
            return visible >= declared
        signature = (declared, visible)
        self._stallPolls = self._stallPolls + 1 if signature == self._stallSignature else 1
        self._stallSignature = signature
        if self._stallPolls >= TERMINAL_STALL_POLLS:
            raise RuntimeError('The input micrographs closed declaring %d items but only %d are visible, and that '
                               'has not changed in %d polls. Refusing to wait for items that are not coming.'
                               % (declared, visible, self._stallPolls))
        return False

    def _insertBatchSteps(self, complete):
        """ One step per full batch; once the input is complete, the remainder goes in a last smaller one.
        Batch size 0 puts every pending micrograph in the same step. """
        pendingIds = sorted(self._pendingIds)
        batchSize = self.streamingBatchSize.get() or len(pendingIds)
        steps = []
        while pendingIds and (len(pendingIds) >= batchSize or complete):
            batchIds, pendingIds = pendingIds[:batchSize], pendingIds[batchSize:]
            self._pendingIds.difference_update(batchIds)
            self._activeBatches += 1
            # Only persistable arguments: the step reads its micrographs again by id
            steps.append(self._insertFunctionStep(self.assessBatchStep, batchIds, prerequisites=[], needsGPU=True))
        return steps

    # --------------------------- UTILS functions -----------------------------
    def _loadMicsByIds(self, micIds):
        """ Clones of the input micrographs of micIds; the batch is assessed whole or not at all. """
        inMics = self._getInputMicrographs()
        mics = {}
        for start in range(0, len(micIds), ID_BATCH_SIZE):
            chunk = micIds[start:start + ID_BATCH_SIZE]
            try:
                for mic in inMics.iterItems(where='id IN (%s)' % ','.join(str(micId) for micId in chunk)):
                    mics[mic.getObjId()] = mic.clone()
            except NotImplementedError:
                for micId in chunk:
                    if micId in inMics:
                        mics[micId] = inMics[micId].clone()
        missing = sorted(set(micIds) - mics.keys())
        if missing:
            raise RuntimeError('Input micrographs %s are not visible' % missing)
        return [mics[micId] for micId in micIds]

    def _runMicAssess(self, mics):
        """ Run micassess on mics in a folder of their own and return the ids of the good ones. """
        workDir = self._getBatchDir(mics)
        cleanPath(workDir)
        makePath(workDir)
        byName = {}
        micsTable = Table(columns=['rlnMicrographName'])
        for mic in mics:
            name = batchMicName(mic.getObjId(), mic.getFileName())
            createLink(os.path.abspath(mic.getFileName()), os.path.join(workDir, name))
            byName[name] = mic.getObjId()
            micsTable.addRow(name)
        with open(os.path.join(workDir, BATCH_INPUT_STAR), 'w') as f:
            f.write("# Star file generated with Scipion\n")
            micsTable.writeStar(f, tableName='')

        self.runJob(Plugin.getProgram('micassess'), ' '.join(self._getArgs()), env=Plugin.getEnviron(),
                    cwd=workDir, numberOfThreads=1)

        outDir = os.path.join(workDir, BATCH_OUTPUT_DIR)
        # micassess creates the label folders once its prediction is over; no good star then just means
        # that no micrograph was good
        missingLabels = [label for label in MICASSESS_LABELS if not os.path.isdir(os.path.join(outDir, label))]
        if missingLabels:
            raise RuntimeError('micassess did not finish: %s has no %s' % (outDir, ', '.join(missingLabels)))
        goodNames = self._readGoodMicNames(os.path.join(workDir, self._goodStarName(BATCH_INPUT_STAR)))
        unknown = sorted(goodNames - byName.keys())
        if unknown:
            raise RuntimeError('micassess reported micrographs that are not in the batch: %s' % unknown)
        with self._lock:
            for label in MICASSESS_LABELS:
                copyTree(os.path.join(outDir, label), self._getExtraPath('MicAssess', label))
        return {byName[name] for name in goodNames}

    def _getBatchDir(self, mics):
        # The ids of a batch are not in any other batch of the same execution
        return self._getTmpPath('batch_%06d' % mics[0].getObjId())

    @staticmethod
    def _goodStarName(inputStar):
        return os.path.splitext(inputStar)[0] + '_good.star'

    @staticmethod
    def _readGoodMicNames(goodStar):
        if not os.path.exists(goodStar):
            return set()
        table = Table(fileName=goodStar, tableName='')
        return {os.path.basename(name) for name in table.getColumnValues('rlnMicrographName')}

    def _publish(self, goodMics, badMics):
        """ Append to each output the micrographs it does not hold yet: a batch published before a failure
        can be published again without duplicates. """
        for name, mics, suffix, label in ((OUTPUT, goodMics, '', 'good micrographs'),
                                          (OUTPUT_DISCARDED, badMics, '_discarded', 'discarded micrographs')):
            if not mics:
                continue
            existing = getattr(self, name, None)
            outMics = self._loadOutputSet(existing, suffix, label)
            newMics = [mic for mic in mics if existing is None or mic.getObjId() not in outMics]
            for mic in newMics:
                outMics.append(mic)
            if newMics:
                self._updateOutputSet(name, outMics, Set.STREAM_OPEN)
                self._appendTotalStars(newMics, name == OUTPUT)
            if existing is None:
                self._ensureSourceRelations()
            self._doneIds.update(mic.getObjId() for mic in mics)

    def _loadOutputSet(self, outMics, suffix, label):
        if outMics is not None:
            outMics.enableAppend()
            return outMics
        outMics = self._createSetOfMicrographs(suffix=suffix)
        outMics.setStreamState(Set.STREAM_OPEN)
        outMics.copyInfo(self._getInputMicrographs())
        # After copyInfo, which copies the label of the input
        outMics.setObjLabel(label)
        return outMics

    def _ensureSourceRelations(self):
        """ Define the source relation of every output that lacks it: it may be missing if the run stopped
        between publishing an output and relating it. This run relates its outputs to its only input, and
        stores the parent as the producer of the input, so the child is enough to know. """
        related = set()
        for relation in self.getRelations():
            try:
                if relation['name'] == RELATION_SOURCE:
                    related.add(relation['object_child_id'])
            except (KeyError, TypeError, IndexError):
                continue
        for name in (OUTPUT, OUTPUT_DISCARDED):
            output = getattr(self, name, None)
            if output is not None and output.getObjId() not in related:
                self._defineSourceRelation(self.inputMicrographs, output)

    def _initTotalStars(self):
        """ Star files in extra with every micrograph assessed, and the good ones. """
        header = 'data_\n\nloop_\n_rlnMicrographName \n'
        for fileName, comment in ((self.getInputFilename(), "# Star file generated with Scipion\n\n"),
                                  (self.getOutputFilename(), '')):
            if not os.path.exists(fileName):
                with open(fileName, 'w') as f:
                    f.write(comment + header)

    def _appendTotalStars(self, mics, good):
        fileNames = [self.getInputFilename()] + ([self.getOutputFilename()] if good else [])
        for fileName in fileNames:
            with open(fileName, 'a') as f:
                for mic in mics:
                    f.write(' ' + os.path.abspath(mic.getFileName()) + '\n')

    def getInputFilename(self):
        return self._getExtraPath('input_micrographs.star')

    def getOutputFilename(self):
        return self._getExtraPath('input_micrographs_good.star')

    def _getArgs(self):
        """ Return the list of args for the command. """
        args = ['-i %s' % BATCH_INPUT_STAR,
                '-o %s' % BATCH_OUTPUT_DIR,
                '-m %s' % Plugin.getVar(CRYOASSESS_MODELS),
                '-b %d' % self.batchSize.get(),
                '--t1 %0.2f' % self.threshold.get(),
                '--t2 %0.2f' % self.threshold2.get(),
                # One thread is the executor's and another the generator's
                '--threads %d' % max(1, self.numberOfThreads.get() - 2),
                '--gpus %(GPU)s']

        if self._getCameraType() is not None:
            args.append('-d %s' % self._getCameraType())

        return args

    def _getInputMicrographs(self):
        return self.inputMicrographs.get()

    def _getCameraType(self):
        """ Get camera type based on input mic size.
        :return string or None """
        micsizeX, micsizeY, _ = self._getInputMicrographs().getDim()
        x = max(micsizeX, micsizeY)
        y = min(micsizeX, micsizeY)
        if math.isclose(x / y, 1.0345, abs_tol=0.001):
            return 'K2'
        elif math.isclose(x / y, 1.4076, abs_tol=0.001):
            return 'K3'
        else:
            return None

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if not os.path.isdir(Plugin.getVar(CRYOASSESS_MODELS)):
            errors.append("Directory with models not found: %s" % CRYOASSESS_MODELS)
        if self.numberOfThreads.get() < 3:
            errors.append('At least 3 threads are needed: one is reserved by the executor and another one by '
                          'the step that watches the input micrographs.')

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, OUTPUT):
            summary.append("Output not ready or no good micrographs found")
        else:
            summary.append("Sorted micrographs into good and bad ones.")

        return summary

    def _warnings(self):
        warnings = []

        if self._getCameraType() is None:
            warnings.append("Micassess model was trained only on data from "
                            "Gatan K2 and K3 cameras.")

        return warnings

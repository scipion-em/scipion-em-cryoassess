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
import math
from emtable import Table

from pyworkflow.constants import BETA
import pyworkflow.protocol.constants as pwcts
from pyworkflow.protocol import params, STATUS_NEW
from pyworkflow.utils.path import copyTree
from pwem.protocols import ProtPreprocessMicrographs
from pwem.objects import SetOfMicrographs, Set

from .. import Plugin
from ..constants import CRYOASSESS_MODEL_MIC


class CryoassessProtMics(ProtPreprocessMicrographs):
    """
    Protocol to assess micrographs from K2 or K3 cameras.

    Find more information at https://github.com/cianfrocco-lab/Automatic-cryoEM-preprocessing
    """
    _label = 'assess micrographs'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtPreprocessMicrographs.__init__(self, **kwargs)
        self.stepsExecutionMode = pwcts.STEPS_PARALLEL

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label="Input micrographs", important=True)
        form.addParam('threshold', params.FloatParam, default=0.1,
                      label='Threshold',
                      help='Threshold for classification. Default is 0.1. '
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
        form.addParallelSection(threads=1, mpi=1)

        self._defineStreamingParams(form)
        form.getParam('streamingBatchSize').setDefault(5)
        form.getParam('streamingSleepOnWait').setDefault(5)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
      self._insertFunctionStep("initializeStep")
      self.closeSet = self._insertFunctionStep('closeSetStep', wait=True)

    def _stepsCheck(self):
      if not self.ended:
        closeStep = self._getFirstJoinStep()
        newMics = self._getNewInput()
        if len(newMics) >= self._getStreamingBatchSize():
          self.addDoneMicFns(newMics)
          numPass = self.asPass
          self.asPass += 1
          newDeps = self._insertNewMicsSteps(newMics, numPass)
          closeStep.addPrerequisites(*newDeps)

        if self.checkIfParentFinished():
          if len(newMics)==0:
            closeStep.setStatus(STATUS_NEW)
          else:
            self.lastRound = True
        self.updateSteps()

    def _insertNewMicsSteps(self, newMics, numPass):
      newSteps = []
      newSteps.append(self._insertFunctionStep('convertInputStep', newMics, numPass, prerequisites=[]))
      newSteps.append(self._insertFunctionStep('runMicAssessStep', numPass, prerequisites=newSteps[-1:]))
      newSteps.append(self._insertFunctionStep('createOutputStep', newMics, numPass, prerequisites=newSteps[-1:]))
      return newSteps


    # --------------------------- STEPS functions -----------------------------
    def initializeStep(self):
      '''Creates all the final output directories where each batch will be appended'''
      self.doneMicFns = set([])
      self.lastRound = False
      self.ended = False
      self.asPass = 1

      self.initTotalStars()

    def convertInputStep(self, newMics, numPass):
        """ Create a star file as expected by cryoassess."""
        micsTable = Table(columns=['rlnMicrographName'])
        for mic in newMics:
            micsTable.addRow(os.path.abspath(mic.getFileName()))
        with open(self.getInputFilename(numPass), 'w') as f:
            f.write("# Star file generated with Scipion\n")
            micsTable.writeStar(f, tableName='')
        self.appendTotalInputStar(numPass)

    def runMicAssessStep(self, numPass):
        """ Call cryoassess with the appropriate parameters. """
        params = ' '.join(self._getArgs(numPass))
        program = Plugin.getProgram('micassess')
        self.runJob(program, params, env=Plugin.getEnviron(),
                    cwd=self._getTmpPath(), numberOfThreads=1)
        self.appendTotalOutputStar(numPass)
        self.copyMicAssessOutput()

    def createOutputStep(self, newMics, numPass):
        outputName = "outputMicrographs"
        outMics = self._loadOutputSet(SetOfMicrographs, outputName+'.sqlite')

        # Parse output file and find good mics
        goodMicNames = self._getGoodMicFns(numPass)
        if len(goodMicNames):
            self.curGoodList = goodMicNames
            outMics.copyItems(newMics, updateItemCallback=self._addGoodMic)
            self._updateOutputSet(outputName, outMics)

    def closeSetStep(self):
      outputName = "outputMicrographs"
      outMics = self._loadOutputSet(SetOfMicrographs, outputName+'.sqlite')
      self._updateOutputSet(outputName, outMics, state=Set.STREAM_CLOSED)

      self._defineSourceRelation(self._getInputMicrographs(), self.outputMicrographs)
      self.ended = True


    # --------------------------- UTILS functions -----------------------------
    def _getStreamingBatchSize(self):
      if self.lastRound:
        return 1
      else:
        return self.streamingBatchSize.get()

    def getGoodInputMics(self):
      inpMics = self._getInputMicrographs()
      goodMics = SetOfMicrographs()
      goodMicFns = self._getGoodMicFns('')
      for inpMic in inpMics:
        if inpMic.getFileName() in goodMicFns:
          goodMics.append(inpMic)
      return goodMics
      
    def copyMicAssessOutput(self):
      copyTree(self._getTmpPath('MicAssess'), self._getExtraPath('MicAssess'))

    def initTotalStars(self):
      totalInputStarFn, totalOutputStar = self.getInputFilename(''), self.getOutputFilename('')
      sameTxt = 'data_\n\nloop_\n_rlnMicrographName \n'
      f1, f2 = open(totalInputStarFn, 'w'), open(totalOutputStar, 'w')
      f1.write("# Star file generated with Scipion\n\n")
      f1.write(sameTxt), f2.write(sameTxt)
      f1.close(), f2.close()

    def appendTotalInputStar(self, numPass):
      totalStarFn = self.getInputFilename('')
      newMicNames = self._getInputMicFns(numPass)
      if os.path.exists(totalStarFn):
        with open(totalStarFn, 'a') as f:
          for micName in newMicNames:
            f.write(' '+micName+'\n')

    def appendTotalOutputStar(self, numPass):
      totalStarFn = self.getOutputFilename('')
      newMicNames = self._getGoodMicFns(numPass)
      if os.path.exists(totalStarFn):
        with open(totalStarFn, 'a') as f:
          for micName in newMicNames:
            f.write(' '+micName+'\n')

    def addDoneMicFns(self, newMics):
      for newMic in newMics:
        self.doneMicFns.add(newMic.getFileName())

    def checkIfParentFinished(self):
        inpMics = self._getInputMicrographs()
        inpMics.loadAllProperties()
        if not inpMics.isStreamOpen():
          return True
        return False

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed
        return 'closeSetStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _updateOutputSet(self, outputName, outputSet, state=Set.STREAM_OPEN):
      outputSet.setStreamState(state)
      if self.hasAttribute(outputName):
        outputSet.write()  # Write to commit changes
        outputAttr = getattr(self, outputName)
        # Copy the properties to the object contained in the protcol
        outputAttr.copy(outputSet, copyId=False)
        # Persist changes
        self._store(outputAttr)
      else:
        # Here the defineOutputs function will call the write() method
        self._defineOutputs(**{outputName: outputSet})
        self._store(outputSet)

      # Close set databaset to avoid locking it
      outputSet.close()

    def _loadOutputSet(self, SetClass, baseName):
      """
      Load the output set if it exists or create a new one.
      """
      setFile = self._getPath(baseName)
      if os.path.exists(setFile) and os.path.getsize(setFile) > 0:
        outputSet = SetClass(filename=setFile)
        outputSet.loadAllProperties()
        outputSet.enableAppend()
      else:
        outputSet = SetClass(filename=setFile)
        outputSet.setStreamState(outputSet.STREAM_OPEN)
        outputSet.setObjLabel('good micrographs')
        outputSet.copyInfo(self._getInputMicrographs())

      return outputSet

    def _getArgs(self, numPass):
        """ Return the list of args for the command. """
        args = ['-i %s ' % os.path.basename(self.getInputFilename(numPass)),
                '-o %s ' % os.path.basename(self.getOutputFilename(numPass)),
                '-m %s' % Plugin.getVar(CRYOASSESS_MODEL_MIC),
                '-b %d' % self.batchSize.get(),
                '-t %0.2f' % self.threshold.get(),
                '--threads %d' % self.numberOfThreads.get(),
                '--gpus %s' % self.gpuList.get().strip().replace(" ", ",")]

        if self._getCameraType() is not None:
            args.append('-d %s' % self._getCameraType())

        return args

    def _getInputMicrographs(self):
        return self.inputMicrographs.get()

    def _getNewInput(self):
        inputMics = self._getInputMicrographs()
        newMics = []
        for mic in inputMics:
          if mic.getFileName() not in self.doneMicFns:
            newMic = mic.clone()
            newMics.append(newMic)
        return newMics

    def getInputFilename(self, numPass):
        if numPass == '':
            return self._getExtraPath('input_micrographs{}.star'.format(numPass))
        else:
            return self._getTmpPath('input_micrographs{}.star'.format(numPass))

    def getOutputFilename(self, numPass):
        if numPass == '':
            return self._getExtraPath('good_micrographs{}.star'.format(numPass))
        else:
          return self._getTmpPath('good_micrographs{}.star'.format(numPass))


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

    def _getRelPath(self, fn):
        """ Return relative path from cwd=extra. """
        return os.path.relpath(fn, self._getExtraPath())

    def _getInputMicFns(self, numPass):
      """ Parse input star file and get a list of mics. """
      table = Table(fileName=self.getInputFilename(numPass), tableName='')
      micNames = table.getColumnValues('rlnMicrographName')
      return micNames

    def _getGoodMicFns(self, numPass):
        """ Parse output star file and get a list of good mics. """
        table = Table(fileName=self.getOutputFilename(numPass), tableName='')
        micNames = table.getColumnValues('rlnMicrographName')
        return micNames

    def _addGoodMic(self, item, row):
        """ Callback function to append only good items. """
        if self._getRelPath(item.getFileName()) not in self.curGoodList:
            setattr(item, "_appendItem", False)


    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputMicrographs'):
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

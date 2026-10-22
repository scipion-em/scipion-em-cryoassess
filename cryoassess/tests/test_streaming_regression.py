# **************************************************************************
# *
# * Regression tests for the streaming and Continue behaviour of cryoassess.
# *
# * They run the protocols on a real project with real Sets; only the
# * cryoassess programs, the step scheduling and the poll sleep are
# * replaced, so neither cryoassess nor any dataset is needed.
# *
# **************************************************************************

import json
import os
from os.path import join
from unittest import mock

import mrcfile
import numpy as np
from emtable import Table

from pyworkflow.object import Pointer, Set, RELATION_SOURCE
from pyworkflow.protocol import FunctionStep
import pyworkflow.protocol.constants as pwcts
from pyworkflow.tests import BaseTest, setupTestProject
from pwem.objects import SetOfMicrographs, Micrograph, SetOfAverages, Particle
from pwem.protocols import EMProtocol

from ..protocols import CryoassessProtMics, CryoassessProt2D
from ..protocols.protocol_micassess import (OUTPUT, OUTPUT_DISCARDED, MICASSESS_LABELS, BATCH_INPUT_STAR,
                                            BATCH_OUTPUT_DIR)

MAX_POLLS = 20
# K2 frame proportions, so the camera type is detected
MIC_SHAPE = (30, 31)


class PollCap:
    """Stands for the poll wait of the generator: counts polls, fails past MAX_POLLS and lets the test act
    (e.g. produce new input) after a given poll."""

    def __init__(self, onPoll=None):
        self.polls = 0
        self.onPoll = onPoll

    def __call__(self):
        self.polls += 1
        if self.polls > MAX_POLLS:
            raise AssertionError("The generator polled %d times without leaving its loop." % self.polls)
        if self.onPoll:
            self.onPoll(self.polls)


class IdQueries:
    """Record the id queries of inputSet; hidden ids are missing from every id query until show(), and
    hideOnce ids from the first one only, as rows that become visible late."""

    def __init__(self, inputSet, hidden=(), hideOnce=()):
        self.wheres = []
        self.hidden, self.hideOnce = set(hidden), set(hideOnce)
        realValues = inputSet.getUniqueValues

        def getUniqueValues(attributes, where=None):
            values = realValues(attributes, where=where)
            if attributes != 'id':
                return values
            hide = self.hidden | (self.hideOnce if not self.wheres else set())
            self.wheres.append(where)
            return [value for value in values if value not in hide]

        self._patch = mock.patch.object(inputSet, 'getUniqueValues', side_effect=getUniqueValues)

    def show(self):
        self.hidden = set()

    def __enter__(self):
        self._patch.__enter__()
        return self

    def __exit__(self, *exc):
        self._patch.__exit__(*exc)


class CryoassessStreamingTest(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def setUp(self):
        self.dataDir = self.getOutputPath(self._testMethodName)
        os.makedirs(self.dataDir, exist_ok=True)
        # Runs resolve their paths from the project folder, as when launched
        cwd = os.getcwd()
        os.chdir(self.proj.getPath())
        self.addCleanup(os.chdir, cwd)

    def _publish(self, outputSet):
        """Make outputSet the output of a stored run, so protocols can point to it."""
        producer = self.proj.newProtocol(EMProtocol)
        self.proj.saveProtocol(producer)
        producer._defineOutputs(outputSet=outputSet)
        producer._store()
        return Pointer(producer, extended='outputSet')

    def _newSetFile(self, name):
        setFn = join(self.dataDir, '%s.db' % name)
        if os.path.exists(setFn):
            os.remove(setFn)  # Left by a previous run of the test
        return setFn


class TestCryoassessMicsStreaming(CryoassessStreamingTest):
    """Streaming contracts of assess micrographs."""

    def _newProtocol(self, names, closed=True, batchSize=2):
        """A saved run over a published Set of micrographs; names are their file names, relative to the
        data folder (a subfolder gives two micrographs the same basename)."""
        self.mics = SetOfMicrographs(filename=self._newSetFile('mics'))
        self.mics.setSamplingRate(1.0)
        self.mics.write()  # A producer stores the Set properties even before its first item
        self._appendMics(names, close=closed)
        prot = self.proj.newProtocol(CryoassessProtMics, streamingBatchSize=batchSize, numberOfThreads=3)
        prot.inputMicrographs.copy(self._publish(self.mics))
        self.proj.saveProtocol(prot)
        prot.makePathsAndClean()
        return prot

    def _appendMics(self, names, close=True):
        """Add micrographs to the input Set, as its producer would, and set its stream."""
        self.mics.close()
        self.mics.enableAppend()
        for name in names:
            micFn = join(self.dataDir, name)
            os.makedirs(os.path.dirname(micFn), exist_ok=True)
            with mrcfile.new(micFn, overwrite=True) as mrc:
                mrc.set_data(np.zeros(MIC_SHAPE, dtype=np.float32))
            mic = Micrograph(location=micFn)
            mic.setSamplingRate(1.0)
            self.mics.append(mic)
        self.mics.setStreamState(Set.STREAM_CLOSED if close else Set.STREAM_OPEN)
        self.mics.write()

    def _runGenerator(self, prot, onPoll=None):
        """Run the generator with its scheduling and poll wait intercepted.

        :return: (the micIds of every batch scheduled, whether the close step was scheduled); the batches
            are also kept in self.batchSteps as (function, args), ready to be executed.
        """
        self.batchSteps, closed = [], []

        def insert(func, *args, **kwargs):
            if func.__name__ == 'assessBatchStep':
                self.batchSteps.append((func, args))
            elif func.__name__ == 'closeOutputStep':
                closed.append(kwargs.get('prerequisites'))
            return len(self.batchSteps) + len(closed)

        with mock.patch.object(prot, '_insertFunctionStep', side_effect=insert), \
                mock.patch.object(prot, '_waitForInput', PollCap(onPoll)):
            prot.stepsGeneratorStep()
        return [list(args[0]) for _, args in self.batchSteps], bool(closed)

    def _fakeMicAssess(self, goodIds=(), fail=False, finish=True):
        """runJob standing for micassess: the micrographs of goodIds are good, the rest bad."""
        def runJob(program, args, cwd=None, **kwargs):
            if fail:
                raise RuntimeError('micassess crashed')
            if not finish:
                return
            for label in MICASSESS_LABELS:
                os.makedirs(join(cwd, BATCH_OUTPUT_DIR, label), exist_ok=True)
            names = Table(fileName=join(cwd, BATCH_INPUT_STAR), tableName='').getColumnValues('rlnMicrographName')
            good = [name for name in names if int(name.split('__')[0]) in goodIds]
            if good:  # micassess writes no good star when nothing is good
                table = Table(columns=['rlnMicrographName'])
                for name in good:
                    table.addRow(name)
                table.write(join(cwd, BATCH_INPUT_STAR.replace('.star', '_good.star')), tableName='')
        return runJob

    def _runBatches(self, prot, goodIds=(), **kwargs):
        with mock.patch.object(prot, 'runJob', side_effect=self._fakeMicAssess(goodIds, **kwargs)):
            for func, args in self.batchSteps:
                func(*args)  # A failed batch is logged and left out, never raised

    @staticmethod
    def _ids(outputSet):
        return sorted(item.getObjId() for item in outputSet) if outputSet is not None else []

    def _outputs(self, prot):
        return self._ids(getattr(prot, OUTPUT, None)), self._ids(getattr(prot, OUTPUT_DISCARDED, None))

    def testEveryBatchIsKeptInTheOutputs(self):
        prot = self._newProtocol(['mic%d.mrc' % i for i in range(1, 6)])
        batches, closed = self._runGenerator(prot)
        self.assertEqual(batches, [[1, 2], [3, 4], [5]])
        self.assertTrue(closed)

        self._runBatches(prot, goodIds={1, 3, 5})
        prot.closeOutputStep()

        self.assertEqual(self._outputs(prot), ([1, 3, 5], [2, 4]),
                         "A batch replaced the micrographs the previous batches had published.")
        self.assertEqual(getattr(prot, OUTPUT).getObjLabel(), 'good micrographs')
        self.assertEqual(getattr(prot, OUTPUT_DISCARDED).getObjLabel(), 'discarded micrographs')
        for name in (OUTPUT, OUTPUT_DISCARDED):
            self.assertFalse(getattr(prot, name).isStreamOpen())

    def testMicrographsWithTheSameBasenameAreAssessedApart(self):
        prot = self._newProtocol(['a/mic.mrc', 'b/mic.mrc'])
        self._runGenerator(prot)

        self._runBatches(prot, goodIds={2})

        self.assertEqual(self._outputs(prot), ([2], [1]),
                         "Two micrographs with the same basename got the verdict of each other.")

    def testPublishingABatchAgainAddsNoDuplicates(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'])
        self._runGenerator(prot)

        # Published, then the run stopped before the step was recorded as done: it runs again
        self._runBatches(prot, goodIds={1})
        self._runBatches(prot, goodIds={1})

        self.assertEqual(self._outputs(prot), ([1], [2]))

    def testContinueDoesNotAssessPublishedMicrographsAgain(self):
        prot = self._newProtocol(['mic%d.mrc' % i for i in range(1, 5)])
        self._runGenerator(prot)
        func, args = self.batchSteps[0]
        with mock.patch.object(prot, 'runJob', side_effect=self._fakeMicAssess({1})):
            func(*args)

        # What Continue rebuilds: the same run with every in-memory state lost
        batches, closed = self._runGenerator(prot)

        self.assertEqual(batches, [[3, 4]])
        self.assertTrue(closed)

    def testABatchRunsFromItsPersistedArguments(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'])
        self._runGenerator(prot)
        func, args = self.batchSteps[0]

        # Continue only has what the step stored, not the objects it was created with
        persisted = json.loads(FunctionStep(func, func.__name__, *args).argsStr.get())
        with mock.patch.object(prot, 'runJob', side_effect=self._fakeMicAssess({1, 2})):
            func(*persisted)

        self.assertEqual(self._outputs(prot), ([1, 2], []))

    def testAFailedBatchIsLeftOutAndTheRestFlows(self):
        prot = self._newProtocol(['mic%d.mrc' % i for i in range(1, 5)])
        self._runGenerator(prot)
        (failFunc, failArgs), (okFunc, okArgs) = self.batchSteps

        with mock.patch.object(prot, 'runJob', side_effect=self._fakeMicAssess(fail=True)):
            failFunc(*failArgs)
        with mock.patch.object(prot, 'runJob', side_effect=self._fakeMicAssess({3})):
            okFunc(*okArgs)
        prot.closeOutputStep()  # Some micrographs were assessed: the run does not fail

        self.assertEqual(self._outputs(prot), ([3], [4]),
                         "The micrographs of a failed batch must stay out of the outputs, so Continue assesses them.")
        batches, _ = self._runGenerator(prot)
        self.assertEqual(batches, [[1, 2]], "Continue did not assess the failed micrographs again.")

    def testTheRunFailsWhenEveryMicrographFailed(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'])
        self._runGenerator(prot)

        self._runBatches(prot, fail=True)

        with self.assertRaisesRegex(RuntimeError, 'every one of the 2 micrographs'):
            prot.closeOutputStep()

    def testAnUnfinishedMicassessFailsTheBatch(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'])
        self._runGenerator(prot)

        # micassess exits without its label folders: it never got to assess anything
        self._runBatches(prot, finish=False)

        self.assertEqual(self._outputs(prot), ([], []),
                         "An unfinished micassess run discarded micrographs that were never assessed.")
        self.assertEqual(prot._failedIds, {1, 2})

    def testNoGoodStarMeansEveryMicrographIsBad(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'])
        self._runGenerator(prot)

        self._runBatches(prot, goodIds=())

        self.assertEqual(self._outputs(prot), ([], [1, 2]))

    def testEachOutputGetsItsSourceRelationOnce(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc', 'mic3.mrc', 'mic4.mrc'])
        self._runGenerator(prot)
        self._runBatches(prot, goodIds={1, 3})
        prot.closeOutputStep()

        sources = [relation['object_child_id'] for relation in prot.getRelations()
                   if relation['name'] == RELATION_SOURCE]
        self.assertEqual(sorted(sources), sorted([getattr(prot, OUTPUT).getObjId(),
                                                  getattr(prot, OUTPUT_DISCARDED).getObjId()]))

    def testPollingOnlyReadsIdsAboveTheWatermark(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'], closed=False)

        def produce(polls):
            if polls == 3:  # More micrographs arrive after a few idle polls and the producer closes
                self._appendMics(['mic3.mrc', 'mic4.mrc'])

        with IdQueries(prot.inputMicrographs.get()) as queries:
            batches, closed = self._runGenerator(prot, onPoll=produce)

        self.assertEqual(batches, [[1, 2], [3, 4]])
        self.assertTrue(closed)
        self.assertNotIn(None, queries.wheres, "A poll listed every input id instead of only the new ones.")
        self.assertEqual(sorted(set(queries.wheres)), ['id > 0', 'id > 2'])

    def testAMicrographVisibleAfterLaterOnesIsStillAssessed(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc', 'mic3.mrc'], closed=False, batchSize=1)

        def close(polls):
            if polls == 1:
                self._appendMics([])

        # Micrograph 2 becomes visible only after micrograph 3 was seen
        with IdQueries(prot.inputMicrographs.get(), hideOnce=[2]):
            batches, closed = self._runGenerator(prot, onPoll=close)

        self.assertEqual(sorted(batches), [[1], [2], [3]], "The late micrograph was never assessed.")
        self.assertTrue(closed)

    def testAClosedInputIsWaitedForUntilAllItsItemsAreVisible(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'], batchSize=1)

        # The producer closed declaring two micrographs, but the second becomes visible two polls later
        with IdQueries(prot.inputMicrographs.get(), hidden=[2]) as queries:
            batches, closed = self._runGenerator(prot, onPoll=lambda polls: polls == 2 and queries.show())

        self.assertEqual(batches, [[1], [2]], "The output was closed before micrograph 2 was visible.")
        self.assertTrue(closed)

    def testAClosedInputThatNeverShowsAllItsItemsFails(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'], batchSize=2)

        with IdQueries(prot.inputMicrographs.get(), hidden=[2]):
            with self.assertRaisesRegex(RuntimeError, 'not coming', msg="Micrograph 2 was silently left out."):
                self._runGenerator(prot)

    def testTheWaitIsNotCountedWhileABatchIsInFlight(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc'], batchSize=1)

        # Micrograph 1 is being assessed (its step never ends here) while micrograph 2 stays hidden
        with IdQueries(prot.inputMicrographs.get(), hidden=[2]):
            with self.assertRaisesRegex(AssertionError, 'without leaving its loop'):
                self._runGenerator(prot)

    def testTheRemainderWaitsForAFullBatchUntilTheInputIsComplete(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc', 'mic3.mrc'], closed=False, batchSize=2)

        def produce(polls):
            if polls == 2:
                self.assertEqual(len(self.batchSteps), 1, "The remainder did not wait for a full batch.")
                self._appendMics([])

        batches, closed = self._runGenerator(prot, onPoll=produce)

        self.assertEqual(batches, [[1, 2], [3]])
        self.assertTrue(closed)

    def testBatchSizeZeroAssessesEveryPendingMicrographTogether(self):
        prot = self._newProtocol(['mic1.mrc', 'mic2.mrc', 'mic3.mrc'], batchSize=0)

        batches, closed = self._runGenerator(prot)

        self.assertEqual(batches, [[1, 2, 3]])
        self.assertTrue(closed)

    def testAnEmptyClosedInputCloses(self):
        prot = self._newProtocol([])
        # Not every storage can read a Set without items: use it as its producer stored it, with no ids
        inMics = prot.inputMicrographs.get()

        with mock.patch.object(prot, '_loadInputSet', return_value=inMics), mock.patch.object(inMics, 'close'), \
                mock.patch.object(inMics, 'getUniqueValues', return_value=[]):
            batches, closed = self._runGenerator(prot)
        prot.closeOutputStep()

        self.assertEqual((batches, closed), ([], True))

    def testAFailedRunStopsTheGenerator(self):
        prot = self._newProtocol(['mic1.mrc'], closed=False)
        prot.setStatus(pwcts.STATUS_FAILED)

        batches, closed = self._runGenerator(prot)

        self.assertEqual((batches, closed), ([], False),
                         "The executor is already joining the generator: it must not schedule more work.")


class TestCryoassess2DOrder(CryoassessStreamingTest):
    """assess 2D classes gives each verdict of 2dassess to the reference it belongs to."""

    def _newProtocol(self, ids):
        avgs = SetOfAverages(filename=self._newSetFile('averages'))
        avgs.setSamplingRate(1.0)
        stackFn = join(self.dataDir, 'averages.mrcs')
        with mrcfile.new(stackFn, overwrite=True) as mrc:
            mrc.set_data(np.stack([np.full((8, 8), i + 1, dtype=np.float32) for i in range(len(ids))]))
        for index, objId in enumerate(ids, start=1):
            avg = Particle(location=(index, stackFn))
            avg.setSamplingRate(1.0)
            avg.setObjId(objId)
            avgs.append(avg)
        avgs.write()
        prot = self.proj.newProtocol(CryoassessProt2D)
        prot.inputRefs.copy(self._publish(avgs))
        self.proj.saveProtocol(prot)
        prot.makePathsAndClean()
        prot._createFilenameTemplates()
        return prot

    def _fake2DAssess(self, prot, verdicts):
        """2dassess: verdicts is the label of each image of the stack, in order; None for a blank image."""
        for index, label in enumerate(verdicts, start=1):
            if label is None:
                continue
            labelDir = prot._getExtraPath('2DAssess', label)
            os.makedirs(labelDir, exist_ok=True)
            open(join(labelDir, 'particle_%d.jpg' % index), 'w').close()

    def testEachVerdictGoesToItsReference(self):
        # Ids whose Python set order is not their order: {1000, 2, 37}
        prot = self._newProtocol([2, 37, 1000])
        prot.convertInputStep()
        self._fake2DAssess(prot, ['Good', 'Noise', 'Noise'])

        prot.createOutputStep()

        self.assertEqual([avg.getObjId() for avg in prot.outputAverages], [2],
                         "The verdict of the first image of the stack went to another reference.")
        self.assertEqual(sorted(avg.getObjId() for avg in prot.outputAverages_discarded), [37, 1000])

    def testContinueStartingAtTheOutputStepKeepsTheOrder(self):
        prot = self._newProtocol([2, 37, 1000])
        prot.convertInputStep()
        self._fake2DAssess(prot, ['Noise', 'Good', 'Noise'])

        # What Continue rebuilds: a new instance of the run, with nothing from convertInputStep in memory
        prot = self.proj.getProtocol(prot.getObjId())
        prot._createFilenameTemplates()
        prot.createOutputStep()

        self.assertEqual([avg.getObjId() for avg in prot.outputAverages], [37])

    def testABlankImageIsDiscarded(self):
        prot = self._newProtocol([2, 37, 1000])
        prot.convertInputStep()
        # 2dassess writes no image for a blank reference: the second one gets no verdict
        self._fake2DAssess(prot, ['Good', None, 'Noise'])

        prot.createOutputStep()

        self.assertEqual([avg.getObjId() for avg in prot.outputAverages], [2])
        self.assertEqual(sorted(avg.getObjId() for avg in prot.outputAverages_discarded), [37, 1000],
                         "A blank reference was left out of both outputs.")

    def testOutputsAreNotSharedBetweenRuns(self):
        self.assertIsNot(CryoassessProt2D().outputDict, CryoassessProt2D().outputDict)

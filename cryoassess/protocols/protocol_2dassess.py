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

from glob import glob
import re
from enum import Enum

from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
from pwem.protocols import ProtProcessParticles
from pwem.objects import SetOfAverages, SetOfClasses2D

from .. import Plugin
from ..constants import CRYOASSESS_MODELS

REF_CLASSES = 0
REF_AVERAGES = 1
DISCARDED = '_discarded'

class outputs(Enum):
    outputAverages = SetOfAverages
    outputAverages_discarded = SetOfAverages
    outputClasses = SetOfClasses2D
    outputClasses_discarded = SetOfClasses2D


class CryoassessProt2D(ProtProcessParticles):
    """
    Protocol to assess 2D classes and 2D averages

    Find more information at https://github.com/cianfrocco-lab/Automatic-cryoEM-preprocessing
    """
    _label = 'assess 2D classes'
    _devStatus = PROD
    _possibleOutputs = outputs
    outputDict = {}

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {'input_cls': self._getExtraPath('input_classes.mrcs')}
        self._pathTemplate = self._getExtraPath('2DAssess/*/particle_*.jpg')
        self._goodClassesRegex = re.compile('Good/particle_(\d*)\.jpg')
        self._regex = re.compile('particle_(\d*)\.jpg')
        self.goodList = []
        self.badList = []

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputRefs', params.PointerParam,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      label="Input references", important=True)
        form.addParam('batchSize', params.IntParam, default=32,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Batch size',
                      help='Batch size used in prediction. Default is 32. '
                           'Increasing this number will result in faster '
                           'prediction. If memory error/warning appears, '
                           'you should lower this number.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('run2DAssessStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Create a mrcs stack as expected by cryoassess and store a lookup table
         to correctly associate the IDs."""
        imgSet = self.inputRefs.get()
        if isinstance(imgSet, SetOfClasses2D):
            self.useAsRef = REF_CLASSES
        else:
            self.useAsRef = REF_AVERAGES

        imgSet.writeStack(self._getFileName('input_cls'))
        self.createLookUpTable()

    def createLookUpTable(self):
        """ Create a lookup table to correctly associate the IDs."""
        classesIDs = self.inputRefs.get().getIdSet()
        self.lookUpTable = {}
        n = 1
        for classID in classesIDs:
            self.lookUpTable[classID] = n
            n += 1

    def run2DAssessStep(self):
        """ Call cryoassess with the appropriate parameters. """
        params = ' '.join(self._getArgs())
        program = Plugin.getProgram('2dassess')
        self.runJob(program, params, env=Plugin.getEnviron())

    def createOutputStep(self):
        inputRefs = self.inputRefs.get()
        if self.useAsRef == REF_CLASSES:
            outRefs = SetOfClasses2D.create(self._getPath())
            outBadRefs = SetOfClasses2D.create(self._getPath(), suffix=DISCARDED)
        elif self.useAsRef == REF_AVERAGES:
            outRefs = self._createSetOfAverages()
            outBadRefs = self._createSetOfAverages(suffix=DISCARDED)

        outRefs.copyInfo(inputRefs)
        outBadRefs.copyInfo(inputRefs)

        # Search through output files and find good and bad classes
        self._getReferences()
        if len(self.goodList):
            if self.useAsRef == REF_CLASSES:
                outRefs.setObjLabel('good 2D classes')
                outRefs.appendFromClasses(inputRefs, filterClassFunc=self._addGoodClass)
                self.outputDict[outputs.outputClasses.name] = outRefs
            else:
                outRefs.setObjLabel('good class averages')
                outRefs.copyItems(inputRefs, updateItemCallback=self._addGoodAvg)
                self.outputDict[outputs.outputAverages.name] = outRefs

        if len(self.badList):
            if self.useAsRef == REF_CLASSES:
                outBadRefs.setObjLabel('bad 2D classes')
                outBadRefs.appendFromClasses(inputRefs, filterClassFunc=self._addBadClass)
                self.outputDict[outputs.outputClasses_discarded.name + DISCARDED] = outBadRefs
            else:
                outBadRefs.setObjLabel('bad class averages')
                outBadRefs.copyItems(inputRefs, updateItemCallback=self._addBadAvg)
                self.outputDict[outputs.outputAverages_discarded.name + DISCARDED] = outBadRefs

        self._defineOutputs(**self.outputDict)
        if self.useAsRef == SetOfClasses2D:
            self._defineSourceRelation(self.inputRefs, outRefs)
            self._defineSourceRelation(self.inputRefs, outBadRefs)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputAverages') or not hasattr(self, 'outputClasses'):
            summary.append("Output not ready or no good references found")
        else:
            summary.append("Sorted class references into good and bad ones.")

        return summary

    # --------------------------- UTILS functions -----------------------------
    def _getArgs(self):
        """ Return the list of args for the command. """
        args = ['-i %s ' % self._getFileName('input_cls'),
                '-m %s/2dassess_062119.h5' % Plugin.getVar(CRYOASSESS_MODELS),
                '-b %d' % self.batchSize.get()]

        return args

    def _getReferences(self):
        """ Return the list of good and bad classes files. """
        files = sorted(glob(self._pathTemplate))
        for i in files:
            s = self._goodClassesRegex.search(i)
            if s:
                self.goodList.append(int(s.group(1)))
            else:
                s = self._regex.search(i)
                self.badList.append(int(s.group(1)))

    def _addGoodAvg(self, item, row):
        """ Callback function to append only good items. """
        if self.lookUpTable[item.getObjId()] not in self.goodList:
            setattr(item, "_appendItem", False)

    def _addBadAvg(self, item, row):
        """ Callback function to append only bad items. """
        if self.lookUpTable[item.getObjId()] not in self.badList:
            setattr(item, "_appendItem", False)

    def _addGoodClass(self, item):
        """ Callback function to append only good classes. """
        return False if self.lookUpTable[item.getObjId()] not in self.goodList else True

    def _addBadClass(self, item):
        """ Callback function to append only bad classes. """
        return False if self.lookUpTable[item.getObjId()] not in self.badList else True

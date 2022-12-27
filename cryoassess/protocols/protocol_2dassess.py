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

from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
from pwem.protocols import ProtProcessParticles
from pwem.objects import SetOfAverages, SetOfClasses2D

from .. import Plugin
from ..constants import *


class CryoassessProt2D(ProtProcessParticles):
    """
    Protocol to assess 2D classes and 2D averages
    """
    _label = 'assess 2D classes'
    _devStatus = PROD
    _possibleOutputs = {
        'outputAverages': SetOfAverages,
        'outputAverages_discarded': SetOfAverages,
        'outputClasses': SetOfClasses2D,
        'outputClasses_discarded': SetOfClasses2D
    }
    outputDict = {}

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {'input_cls': self._getExtraPath('input_classes.mrcs')}
        self._pathTemplate = self._getExtraPath('2DAssess/*/particle_*.jpg')
        self._goodClassesRegex = re.compile(r'Good/particle_(\d*)\.jpg')
        self._regex = re.compile(r'particle_(\d*)\.jpg')
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
        """ Create a mrcs stack as expected by cryoassess and
        store a list of objIds."""
        imgSet = self.inputRefs.get()
        imgSet.writeStack(self._getFileName('input_cls'))

        self.clsIds = list(self.inputRefs.get().getIdSet())

    def run2DAssessStep(self):
        """ Call cryoassess with the appropriate parameters. """
        params = ' '.join(self._getArgs())
        program = Plugin.getProgram('2dassess')
        self.runJob(program, params, env=Plugin.getEnviron())

    def createOutputStep(self):
        inputRefs = self.inputRefs.get()
        if self._getRefsType() == REF_CLASSES:
            outRefs = SetOfClasses2D.create(self._getPath())
            outBadRefs = SetOfClasses2D.create(self._getPath(),
                                               suffix=DISCARDED)
        else:
            outRefs = self._createSetOfAverages()
            outBadRefs = self._createSetOfAverages(suffix=DISCARDED)

        outRefs.copyInfo(inputRefs)
        outBadRefs.copyInfo(inputRefs)

        # Search through output files and find good and bad classes
        self._getReferences()
        if len(self.goodList):
            if self._getRefsType() == REF_CLASSES:
                outRefs.setObjLabel('good 2D classes')
                outRefs.appendFromClasses(inputRefs,
                                          filterClassFunc=self._addGoodClass)
                self.outputDict['outputClasses'] = outRefs
            else:
                outRefs.setObjLabel('good class averages')
                outRefs.copyItems(inputRefs,
                                  updateItemCallback=self._addGoodAvg)
                self.outputDict['outputAverages'] = outRefs

        if len(self.badList):
            if self._getRefsType() == REF_CLASSES:
                outBadRefs.setObjLabel('bad 2D classes')
                outBadRefs.appendFromClasses(inputRefs,
                                             filterClassFunc=self._addBadClass)
                self.outputDict['outputClasses_discarded'] = outBadRefs
            else:
                outBadRefs.setObjLabel('bad class averages')
                outBadRefs.copyItems(inputRefs,
                                     updateItemCallback=self._addBadAvg)
                self.outputDict['outputAverages_discarded'] = outBadRefs

        self._defineOutputs(**self.outputDict)
        if self._getRefsType() == REF_CLASSES:
            if len(self.goodList):
                self._defineSourceRelation(self.inputRefs, outRefs)
            if len(self.badList):
                self._defineSourceRelation(self.inputRefs, outBadRefs)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output not ready.")
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
        if self.clsIds.index(item.getObjId()) + 1 not in self.goodList:
            setattr(item, "_appendItem", False)

    def _addBadAvg(self, item, row):
        """ Callback function to append only bad items. """
        if self.clsIds.index(item.getObjId()) + 1 not in self.badList:
            setattr(item, "_appendItem", False)

    def _addGoodClass(self, item):
        """ Callback function to append only good classes. """
        return False if self.clsIds.index(item.getObjId()) + 1 not in self.goodList else True

    def _addBadClass(self, item):
        """ Callback function to append only bad classes. """
        return False if self.clsIds.index(item.getObjId()) + 1 not in self.badList else True

    def _getRefsType(self):
        imgSet = self.inputRefs.get()
        return REF_CLASSES if isinstance(imgSet, SetOfClasses2D) else REF_AVERAGES

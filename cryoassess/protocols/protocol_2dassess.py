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

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pwem.protocols import ProtProcessParticles

from .. import Plugin
from ..constants import CRYOASSESS_MODEL_2D


class CryoassessProt2D(ProtProcessParticles):
    """
    Protocol to assess 2D class averages.

    Find more information at https://github.com/cianfrocco-lab/Automatic-cryoEM-preprocessing
    """
    _label = 'assess 2D classes'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {'input_cls': self._getExtraPath('input_classes.mrcs')}
        self._goodTemplate = self._getExtraPath('2DAssess/Good/particle_*.jpg')
        self._regex = re.compile('particle_(\d*)\.jpg')
        self.goodList = []

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputAverages', params.PointerParam,
                      pointerClass='SetOfAverages',
                      label="Input averages", important=True)
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
        """ Create a mrcs stack as expected by cryoassess."""
        imgSet = self.inputAverages.get()
        imgSet.writeStack(self._getFileName('input_cls'))

    def run2DAssessStep(self):
        """ Call cryoassess with the appropriate parameters. """
        params = ' '.join(self._getArgs())
        program = Plugin.getProgram('2dassess')
        self.runJob(program, params, env=Plugin.getEnviron())

    def createOutputStep(self):
        inputAvg = self.inputAverages.get()
        outAvgs = self._createSetOfAverages()
        outAvgs.copyInfo(inputAvg)
        outAvgs.setObjLabel('good class averages')

        # Search through output files and find good classes
        self._getGoodAvgs()
        if len(self.goodList):
            outAvgs.copyItems(inputAvg, updateItemCallback=self._addGoodAvg)
            self._defineOutputs(outputAverages=outAvgs)
            self._defineSourceRelation(self.inputAverages, outAvgs)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputAverages'):
            summary.append("Output not ready or no good averages found")
        else:
            summary.append("Sorted class averages into good and bad ones.")

        return summary

    # --------------------------- UTILS functions -----------------------------
    def _getArgs(self):
        """ Return the list of args for the command. """
        args = ['-i %s ' % self._getFileName('input_cls'),
                '-m %s' % Plugin.getVar(CRYOASSESS_MODEL_2D),
                '-b %d' % self.batchSize.get()]

        return args

    def _getGoodAvgs(self):
        """ Return the list of good class files. """
        files = sorted(glob(self._goodTemplate))
        if files:
            for i in files:
                s = self._regex.search(i)
                self.goodList.append(int(s.group(1)))

    def _addGoodAvg(self, item, row):
        """ Callback function to append only good items. """
        if item.getObjId() not in self.goodList:
            setattr(item, "_appendItem", False)

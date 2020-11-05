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

from pyworkflow.constants import VERSION_3_0
import pyworkflow.protocol.params as params
from pwem.protocols import ProtPreprocessMicrographs

from .. import Plugin
from ..constants import CRYOASSESS_MODEL_MIC


class CryoassessProtMics(ProtPreprocessMicrographs):
    """
    Protocol to assess micrographs from K2 or K3 cameras.

    Find more information at https://github.com/cianfrocco-lab/Automatic-cryoEM-preprocessing
    """
    _lastUpdateVersion = VERSION_3_0
    _label = 'assess micrographs'

    def __init__(self, **kwargs):
        ProtPreprocessMicrographs.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'input_mics': self._getExtraPath('input_micrographs.star'),
            'output_mics': self._getExtraPath('good_micrographs.star')
        }
        self._goodList = []

        self._updateFilenamesDict(myDict)

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
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runMicAssessStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Create a star file as expected by cryoassess."""
        imgSet = self._getInputMicrographs()
        micsTable = Table(columns=['rlnMicrographName'])
        for img in imgSet:
            micsTable.addRow(self._getRelPath(img.getFileName()))
        with open(self._getFileName('input_mics'), 'w') as f:
            f.write("# Star file generated with Scipion\n")
            micsTable.writeStar(f, tableName='')

    def runMicAssessStep(self):
        """ Call cryoassess with the appropriate parameters. """
        params = ' '.join(self._getArgs())
        program = Plugin.getProgram('micassess')
        self.runJob(program, params, env=Plugin.getEnviron(),
                    cwd=self._getExtraPath(), numberOfThreads=1)

    def createOutputStep(self):
        inputMics = self._getInputMicrographs()
        outMics = self._createSetOfMicrographs()
        outMics.copyInfo(inputMics)
        outMics.setObjLabel('good micrographs')

        self._getGoodMics()
        if len(self._goodList):
            outMics.copyItems(inputMics, updateItemCallback=self._addGoodMic)
            self._defineOutputs(outputMicrographs=outMics)
            self._defineSourceRelation(self.inputMicrographs, outMics)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputMicrographs'):
            summary.append("Output not ready or no good micrographs found")
        else:
            summary.append("Sorted micrographs into good and bad ones.")

        return summary

    def _validate(self):
        errors = []

        if self._getCameraType() is None:
            errors.append("Wrong input dimensions!\n"
                          "This programs only supports K2 or K3 images!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getArgs(self):
        args = ['-i %s ' % os.path.basename(self._getFileName('input_mics')),
                '-o %s ' % os.path.basename(self._getFileName('output_mics')),
                '-m %s' % Plugin.getVar(CRYOASSESS_MODEL_MIC),
                '-b %d' % self.batchSize.get(),
                '-t %0.2f' % self.threshold.get(),
                '-d %s' % self._getCameraType(),
                '--threads %d' % self.numberOfThreads.get(),
                '--gpus %s' % self.gpuList.get().strip().replace(" ", ",")]

        return args

    def _getInputMicrographs(self):
        return self.inputMicrographs.get()

    def _getCameraType(self):
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

    def _getGoodMics(self):
        table = Table(fileName=self._getFileName('output_mics'), tableName='')
        micNames = table.getColumnValues('rlnMicrographName')
        print("FOUND:", micNames)
        self._goodList.extend(micNames)

    def _addGoodMic(self, item, row):
        if self._getRelPath(item.getFileName()) not in self._goodList:
            setattr(item, "_appendItem", False)

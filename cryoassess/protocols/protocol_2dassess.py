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

import pyworkflow.utils as pwutils
from pyworkflow.constants import VERSION_3_0
from pyworkflow.protocol.constants import STEPS_PARALLEL
import pyworkflow.protocol.params as params
from pwem.protocols import ProtProcessParticles

from cryoassess import Plugin


class CryoassessProt2D(ProtProcessParticles):
    """
    Protocol to assess 2D class averages.

    Find more information at https://github.com/cianfrocco-lab/Automatic-cryoEM-preprocessing
    """
    _lastUpdateVersion = VERSION_3_0
    _label = 'assess 2D classes'

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'input_cls': self._getExtraPath('input_classes.mrcs'),
            'output_cls': self._getExtraPath('good_classes.mrcs')
        }

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputAverages', params.PointerParam,
                      pointerClass='SetOfAverages',
                      label="Input averages", important=True)
        form.addParam('modelFile', params.FileParam,
                      label='Pre-trained model file',
                      help='Provide 2dassess model file.')
        form.addParam('batchSize', params.IntParam, default=32,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Batch size',
                      help='Batch size used in prediction. Default is 32. '
                           'Increasing this number will result in faster '
                           'prediction. If memory error/warning appears, '
                           'you should lower this number.')
        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('run2DAssessStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Create a mrcs file as expected by cryoassess."""
        imgSet = self.inputAverages.get()
        imgSet.writeStack(self._getFileName('input_cls'))

    def run2DAssessStep(self):
        """ Call cryoassess with the appropriate parameters. """
        params = ' '.join(self._getArgs())
        program = Plugin.getProgram('2dassess')
        self.runJob(program, params, env=Plugin.getEnviron())

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        self._createFilenameTemplates()
        if not pwutils.exists(self._getFileName("output_cls")):
            summary.append("Output not ready")
        else:
            summary.append("Sorted class averages into good and bad ones.")

        return summary

    def _validate(self):
        errors = []

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getArgs(self):
        args = ['-i %s ' % self._getFileName('input_cls'),
                '-m %s' % self.modelFile.get(),
                '-b %d' % self.batchSize.get()
                ]

        return args

    def _getRelPath(self, fn):
        """ Return relative path from cwd=extra. """
        return os.path.relpath(fn, self._getExtraPath())

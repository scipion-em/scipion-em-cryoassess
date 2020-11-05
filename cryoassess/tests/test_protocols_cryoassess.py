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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwem.protocols import ProtImportMicrographs, ProtImportAverages
from pyworkflow.utils import magentaStr

from ..protocols import *


class TestCryoassess(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.avgsFn = cls.ds.getFile('import/averages.mrcs')
        cls.micsFn = cls.ds.getFile('import/mics/*mrc')

        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImportMics = cls.newProtocol(
            ProtImportMicrographs,
            samplingRateMode=0,
            filesPath=cls.micsFn,
            samplingRate=3,
            magnification=50000,
            voltage=300,
            sphericalAberration=0.1)
        # cls.launchProtocol(cls.protImportMics)

        print(magentaStr("\n==> Importing data - averages:"))
        cls.protImportAvgs = cls.newProtocol(
            ProtImportAverages,
            filesPath=cls.avgsFn,
            samplingRate=3,
            checkStack=True)
        cls.launchProtocol(cls.protImportAvgs)

    def testMicAssess(self):
        print(magentaStr("\n==> Testing cryoassess - micassess:"))
        protAssessMics = self.newProtocol(
            CryoassessProtMics,
            inputMicrographs=self.protImportMics.outputMicrographs,
            modelFile="/home/gsharov/soft/scipion3/software/em/cryoassess-models/micassess_051419.h5"
        )
        self.launchProtocol(protAssessMics)
        # self._checkOutput(protAssessMics, 300, 400)

    def test2DAssess(self):
        print(magentaStr("\n==> Testing cryoassess - 2d assess:"))
        protAssess2D = self.newProtocol(
            CryoassessProt2D,
            inputAverages=self.protImportAvgs.outputAverages,
            modelFile="/home/gsharov/soft/scipion3/software/em/cryoassess-models/2dassess_062119.h5"
        )
        self.launchProtocol(protAssess2D)
        # self._checkOutput(protAssess2D, 300, 400)

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
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportMicrographs, ProtImportAverages

from ..protocols import CryoassessProt2D, CryoassessProtMics


class TestCryoassess(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.ds2 = DataSet.getDataSet('mda')
        cls.micsFn = cls.ds.getFile('micrographs/*mrc')
        cls.avgsFn = cls.ds2.getFile('averages/averages.stk')

        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImportMics = cls.newProtocol(
            ProtImportMicrographs,
            filesPath=cls.micsFn,
            samplingRate=7.08)
        cls.launchProtocol(cls.protImportMics)

        print(magentaStr("\n==> Importing data - averages:"))
        cls.protImportAvgs = cls.newProtocol(
            ProtImportAverages,
            filesPath=cls.avgsFn,
            samplingRate=5.04,
            checkStack=True)
        cls.launchProtocol(cls.protImportAvgs)

    def testMicAssess(self):
        print(magentaStr("\n==> Testing cryoassess - micassess:"))
        protAssessMics = self.newProtocol(
            CryoassessProtMics,
            inputMicrographs=self.protImportMics.outputMicrographs)
        self.launchProtocol(protAssessMics)
        micSet = getattr(protAssessMics, 'outputMicrographs', None)
        self.assertIsNotNone(micSet)

    def test2DAssess(self):
        print(magentaStr("\n==> Testing cryoassess - 2d assess:"))
        protAssess2D = self.newProtocol(
            CryoassessProt2D,
            inputAverages=self.protImportAvgs.outputAverages)
        self.launchProtocol(protAssess2D)
        avgSet = getattr(protAssess2D, 'outputAverages', None)
        self.assertIsNotNone(avgSet)

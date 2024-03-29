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


def getCryoAssessEnvName(version):
    return "cryoassess-%s" % version


V1_0_0 = "1.0.0"
VERSIONS = [V1_0_0]
CRYOASSESS_DEFAULT_VER_NUM = V1_0_0

DEFAULT_ENV_NAME = getCryoAssessEnvName(CRYOASSESS_DEFAULT_VER_NUM)
DEFAULT_ACTIVATION_CMD = 'conda activate ' + DEFAULT_ENV_NAME
CRYOASSESS_ENV_ACTIVATION = 'CRYOASSESS_ENV_ACTIVATION'
CRYOASSESS_MODELS = 'CRYOASSESS_MODELS'

REF_CLASSES = 0
REF_AVERAGES = 1
DISCARDED = '_discarded'

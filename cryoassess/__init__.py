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
import pwem
import pyworkflow.utils as pwutils
from pyworkflow import Config

from .constants import *


__version__ = '3.0.4'
_references = ['Li2020']
_logo = "cryoassess_logo.png"


class Plugin(pwem.Plugin):
    _url = "https://github.com/scipion-em/scipion-em-cryoassess"
    _supportedVersions = VERSIONS

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(CRYOASSESS_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)
        cls._defineEmVar(CRYOASSESS_MODEL_MIC,
                         'cryoassess-models/micassess_051419.h5')
        cls._defineEmVar(CRYOASSESS_MODEL_2D,
                         'cryoassess-models/2dassess_062119.h5')

    @classmethod
    def getCryoAssessEnvActivation(cls):
        """ Remove the scipion home and activate the conda environment. """
        activation = cls.getVar(CRYOASSESS_ENV_ACTIVATION)
        scipionHome = Config.SCIPION_HOME + os.path.sep

        return activation.replace(scipionHome, "", 1)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch cryoassess. """
        environ = pwutils.Environ(os.environ)
        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            del environ['PYTHONPATH']
        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
        activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def defineBinaries(cls, env):
        for ver in VERSIONS:
            cls.addCryoAssessPackage(env, ver,
                                     default=ver == CRYOASSESS_DEFAULT_VER_NUM)

    @classmethod
    def addCryoAssessPackage(cls, env, version, default=False):
        CRYOASSESS_INSTALLED = 'cryoassess_%s_installed' % version
        ENV_NAME = getCryoAssessEnvName(version)

        # try to get CONDA activation command
        installCmd = [cls.getCondaActivationCmd()]

        # Create the environment
        installCmd.append('conda create -y -n %s -c anaconda python=3.6 '
                          'pyqt=5 cudnn=7.1.2 intel-openmp=2019.4;' % ENV_NAME)

        # Activate the new environment
        installCmd.append('conda activate %s;' % ENV_NAME)

        # Install downloaded code
        url = "https://github.com/cianfrocco-lab/Automatic-cryoEM-preprocessing.git"
        installCmd.extend(['git clone %s cryoassess-master &&'
                           'cd cryoassess-master && git checkout master &&'
                           'cd .. && pip install -e cryoassess-master[gpu] &&' % url])

        # Flag installation finished
        installCmd.append('touch %s' % CRYOASSESS_INSTALLED)

        cryoassess_commands = [(" ".join(installCmd), CRYOASSESS_INSTALLED)]

        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        env.addPackage('cryoassess', version=version,
                       tar='void.tgz',
                       commands=cryoassess_commands,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def getProgram(cls, program):
        """ Create cryoAssess command line. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(),
                                       cls.getCryoAssessEnvActivation(),
                                       program)
        return fullProgram

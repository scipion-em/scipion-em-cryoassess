=================
cryoassess plugin
=================

This plugin provide a wrapper around `cryoassess <https://github.com/cianfrocco-lab/Automatic-cryoEM-preprocessing>`_ software: tools for automatic micrograph and 2D classes assessment.

.. image:: https://img.shields.io/pypi/v/scipion-em-cryoassess.svg
        :target: https://pypi.python.org/pypi/scipion-em-cryoassess
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-cryoassess.svg
        :target: https://pypi.python.org/pypi/scipion-em-cryoassess
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-cryoassess.svg
        :target: https://pypi.python.org/pypi/scipion-em-cryoassess
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-cryoassess?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-cryoassess
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-cryoassess
        :target: https://pypi.python.org/pypi/scipion-em-cryoassess
        :alt: Downloads


+--------------+----------------+
| prod: |prod| | devel: |devel| |
+--------------+----------------+

.. |prod| image:: http://scipion-test.cnb.csic.es:9980/badges/cryoassess_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/cryoassess_devel.svg


Installation
-------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-cryoassess

b) Developer's version

   * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-cryoassess.git

   * install

    .. code-block::

       scipion installp -p path_to_scipion-em-cryoassess --devel

cryoassess binaries will be installed automatically with the plugin at **software/em/cryoassess-0.1.0**, but you can also link an existing installation.

**Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program.

Configuration variables
-----------------------

*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen bellow but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"

*CRYOASSESS_HOME* (default = software/em/cryoassess-0.1.0):
Path where the cryoassess is installed.

*CRYOASSESS_ACTIVATION_CMD* (default = conda activate cryoassess-0.1.0):
Command to activate the cryoassess environment.

The deep-learning models can be downloaded from
`authors' website <https://cosmic-cryoem.org/software/cryo-assess/>`_ and set with the following variables:

*CRYOASSESS_MODEL_MIC* (default = software/em/cryoassess-models/micassess_051419.h5)

*CRYOASSESS_MODEL_2D* (default = software/em/cryoassess-models/2dassess_062119.h5)

Verifying
---------

To check the installation, simply run the following Scipion test:

``scipion test cryoassess.tests.test_protocols_cryoassess.TestCryoassess``

Supported versions
------------------

0.1.0

Protocols
----------

* assess micrographs
* assess 2D classes

References
-----------

1. High-Throughput Cryo-EM Enabled by User-Free Preprocessing Routines. Yilai Li, Jennifer N.Cash, John J.G. Tesmer, Michael A.Cianfrocco. Structure 2020, Volume 28 (7), Pages 858-869.e3

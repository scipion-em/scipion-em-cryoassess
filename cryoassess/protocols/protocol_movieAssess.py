# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo ( daniel.delhoyo.gomez@alumnos.upm.es )
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from pyworkflow.constants import BETA
import pyworkflow.protocol.constants as pwcts
from pyworkflow.protocol import params, STATUS_NEW
from pyworkflow.utils.path import copyTree
import pyworkflow.utils as pwutils

from pwem.protocols import ProtPreprocessMicrographs, ProtAlignMovies
from pwem.objects import SetOfMicrographs, Set, Micrograph, SetOfMovies

from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md
from pwem import emlib

from os.path import basename
import numpy as np
from .. import Plugin
from ..constants import CRYOASSESS_MODEL_MIC


def writeImageFromArray(array, fn):
  img = emlib.Image()
  img.setData(array)
  img.write(fn)


def readImage(fn):
  img = emlib.Image()
  img.read(fn)
  return img


def applyTransform(imag_array, M, shape):
  ''' Apply a transformation(M) to a np array(imag) and return it in a given shape
  '''
  imag = emlib.Image()
  imag.setData(imag_array)
  imag = imag.applyWarpAffine(list(M.flatten()), shape, True)
  return imag.getData()


def rotation(imag, angle, shape, P):
  '''Rotate a np.array and return also the transformation matrix
  #imag: np.array
  #angle: angle in degrees
  #shape: output shape
  #P: transform matrix (further transformation in addition to the rotation)'''
  (hsrc, wsrc) = imag.shape
  angle *= math.pi / 180
  T = np.asarray([[1, 0, -wsrc / 2], [0, 1, -hsrc / 2], [0, 0, 1]])
  R = np.asarray([[math.cos(angle), math.sin(angle), 0], [-math.sin(angle), math.cos(angle), 0], [0, 0, 1]])
  M = np.matmul(np.matmul(np.linalg.inv(T), np.matmul(R, T)), P)

  transformed = applyTransform(imag, M, shape)
  return transformed, M


def flipYImage(inFn, outFn=None, outDir=None):
  '''Flips an image in the Y axis'''
  if outFn == None:
    if not '_flipped' in basename(inFn):
      ext = pwutils.getExt(inFn)
      outFn = inFn.replace(ext, '_flipped' + ext)
    else:
      outFn = inFn.replace('_flipped', '')
  if outDir != None:
    outFn = outDir + '/' + basename(outFn)
  gainImg = readImage(inFn)
  imag_array = np.asarray(gainImg.getData(), dtype=np.float64)

  # Flipped Y matrix
  M, angle = np.asarray([[1, 0, 0], [0, -1, imag_array.shape[0]], [0, 0, 1]]), 0
  flipped_array, M = rotation(imag_array, angle, imag_array.shape, M)
  writeImageFromArray(flipped_array, outFn)
  return outFn

def writeMovieMd(movie, outXmd, f1, fN, useAlignment=False):
  movieMd = md.MetaData()
  frame = movie.clone()
  # get some info about the movie
  # problem is, that it can come from a movie set, and some
  # values might refer to different movie, e.g. no of frames :(
  firstFrame, _, frameIndex = movie.getFramesRange()  # info from the movie set
  _, _, lastFrame = movie.getDim()  # actual no. of frame in the current movie
  lastFrame += 1  # (convert no. of frames to index, one-initiated)
  if lastFrame == 0:
    # this condition is for old SetOfMovies, that has lastFrame = 0.
    frames = movie.getNumberOfFrames()
    if frames is not None:
      lastFrame = frames

  if f1 < firstFrame or fN > lastFrame:
    raise Exception("Frame range could not be greater"
                    " than the movie one.")

  ih = ImageHandler()

  if useAlignment:
    alignment = movie.getAlignment()
    if alignment is None:
      raise Exception("Can not write alignment for movie. ")
    a0, aN = alignment.getRange()
    if a0 < firstFrame or aN > lastFrame:
      raise Exception("Trying to write frames which have not been aligned.")
    shiftListX, shiftListY = alignment.getShifts()

  row = md.Row()
  stackIndex = frameIndex + (f1 - firstFrame)

  for i in range(f1, fN + 1):
    frame.setIndex(stackIndex)
    row.setValue(md.MDL_IMAGE, ih.locationToXmipp(frame))

    if useAlignment:
      shiftIndex = i - firstFrame
      row.setValue(emlib.MDL_SHIFT_X, shiftListX[shiftIndex])
      row.setValue(emlib.MDL_SHIFT_Y, shiftListY[shiftIndex])

    row.addToMd(movieMd)
    stackIndex += 1

  movieMd.write(outXmd)



CROP_NONE = 0
CROP_ALIGNMENT = 1
CROP_NEW = 2

class CryoassessProtMovies(ProtAlignMovies):
  """
  Protocol to assess movie averages from K2 or K3 cameras.
  """
  _label = 'assess movie averages'
  _devStatus = BETA

  INTERP_LINEAR = 0
  INTERP_CUBIC = 1

  # Map to xmipp interpolation values in command line
  INTERP_MAP = {INTERP_LINEAR: 1, INTERP_CUBIC: 3}

  outputAveragesName = 'allAverages'
  outputMicsName, outputMoviesName = 'outputGoodAverages', 'outputGoodMovies'

  def __init__(self, **kwargs):
    ProtPreprocessMicrographs.__init__(self, **kwargs)
    self.stepsExecutionMode = pwcts.STEPS_PARALLEL

  # --------------------------- DEFINE param functions ----------------------
  def _defineParams(self, form):
    form.addSection(label='Input')
    form.addParam('inputMovies', params.PointerParam,
                  pointerClass='SetOfMovies',
                  label="Input movies", important=True)
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

    form.addSection('Averaging')
    group = form.addGroup('Average')
    line = group.addLine('Frames to SUM',
                         help='Frames range to SUM on each movie. The '
                              'first frame is 1. If you set 0 in the final '
                              'frame to sum, it means that you will sum '
                              'until the last frame of the movie.')
    line.addParam('sumFrame0', params.IntParam, label='from')
    line.addParam('sumFrameN', params.IntParam, label='to')

    group.addParam('binFactor', params.FloatParam, default=1,
                   label='Binning factor',
                   help='Binning factor, it may be any floating number '
                        'Binning in Fourier is the first operation, so '
                        'that crop parameters are referred to the binned '
                        'images. ')

    group.addParam('cropRegion', params.EnumParam,
                   choices=['None', 'From Alignment', 'New'],
                   label="Define crop region", default=CROP_NONE,
                   help="Select if you want to crop the final micrograph. "
                        "If you select: \n"
                        "*None*, the final micrographs will have the same "
                        "dimensions as the input movies. The region of "
                        "interest, if the input movies have alignment, is "
                        "ignored. \n"
                        "*from Alignment*, if the movies has alignment, "
                        "the region of interest is defined and will be "
                        "apply; else, micrographs will have the same "
                        "dimensions as the input movies. \n"
                        "*New*, All crop parameters should be defined below."
                        "The region of interest, if the input movies have "
                        "alignment, is ignored. ")

    line = group.addLine('Crop offsets (px)',
                         condition='cropRegion==%d' % CROP_NEW)
    line.addParam('cropOffsetX', params.IntParam, default=0, label='X')
    line.addParam('cropOffsetY', params.IntParam, default=0, label='Y')

    line = group.addLine('Crop dimensions (px)',
                         condition='cropRegion==%d' % CROP_NEW,
                         help='How many pixels to crop from offset\n'
                              'If equal to 0, use maximum size.')
    line.addParam('cropDimX', params.IntParam, default=0, label='X')
    line.addParam('cropDimY', params.IntParam, default=0, label='Y')

    group.addParam('useAlignment', params.BooleanParam, default=True,
                   label="Use previous movie alignment to SUM frames?",
                   help="Input movies could have alignment information from"
                        "a previous protocol. If you select *Yes*, the "
                        "previous alignment will be taken into account.")

    form.addParam('splineOrder', params.EnumParam,
                  default=self.INTERP_CUBIC, choices=['linear', 'cubic'],
                  expertLevel=params.LEVEL_ADVANCED,
                  label='Interpolation',
                  help="linear (faster but lower quality), "
                       "cubic (slower but more accurate).")

    form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                   label="Choose GPU IDs",
                   help="GPU may have several cores. Set it to zero"
                        " if you do not know what we are talking about."
                        " First core index is 0, second 1 and so on."
                        " Micassess can use multiple GPUs - in that case"
                        " set to i.e. *0 1 2*.")
    form.addParallelSection(threads=4, mpi=1)

    self._defineStreamingParams(form)
    form.getParam('streamingBatchSize').setDefault(5)
    form.getParam('streamingSleepOnWait').setDefault(5)

  # --------------------------- INSERT steps functions ----------------------
  def _insertAllSteps(self):
    self._insertFunctionStep("initializeStep")
    self.closeSet = self._insertFunctionStep('closeSetStep', wait=True)

  def _stepsCheck(self):
    if not self.ended:
      closeStep = self._getFirstJoinStep()
      
      newMovies = self._getNewMovies()
      print('New movies: ', len(newMovies))
      if len(newMovies) > 0:
        self.addDoneMovieFns(newMovies)
        newDeps = self._insertNewMoviesSteps(newMovies)
        closeStep.addPrerequisites(*newDeps)
      
      newMics = self._getNewMics()
      print('New mics: ', len(newMics))
      if len(newMics) >= self._getStreamingBatchSize():
        self.addDoneMicFns(newMics)
        numPass = self.asPass
        self.asPass += 1
        newDeps = self._insertNewMicsSteps(newMics, numPass)
        closeStep.addPrerequisites(*newDeps)

      if self.checkIfParentFinished() and not self.checkIfMoviesPending(newMovies):
        if len(newMics) == 0:
          closeStep.setStatus(STATUS_NEW)
        else:
          self.lastRound = True
      self.updateSteps()

  def _insertNewMoviesSteps(self, newMovies):
    newSteps = []
    for movie in newMovies:
      # Prerequisites to save the Averages: movie is preprocessed and previous average is already saved (avoid concurrency)
      writePrer = newSteps[-1:] if len(newSteps) > 0 else []
      newSteps.append(self._insertFunctionStep('preprocessMovieStep', movie.clone(), prerequisites=writePrer))
    return newSteps

  def _insertNewMicsSteps(self, newMics, numPass):
    newSteps = []
    newSteps.append(self._insertFunctionStep('convertInputStep', newMics, numPass, prerequisites=[]))
    newSteps.append(self._insertFunctionStep('runMicAssessStep', numPass, prerequisites=newSteps[-1:]))
    newSteps.append(self._insertFunctionStep('createOutputStep', newMics, numPass, prerequisites=newSteps[-1:]))
    return newSteps

  # --------------------------- STEPS functions -----------------------------
  def initializeStep(self):
    '''Creates all the final output directories where each batch will be appended'''
    self.doneMovieFns, self.doneMicFns = set([]), set([])
    self.lastRound = False
    self.ended = False
    self.asPass = 1

    self.initTotalStars()
    
  def preprocessMovieStep(self, movie):
      avFn = self._processMovie(movie)
      newAverage = self.createMicrograph(avFn, movie)

      outMics = self._getAveragedMicrographs()
      outMics.append(newAverage)
      self._updateOutputSet(self.outputAveragesName, outMics)

  def convertInputStep(self, newMics, numPass):
    """ Create a star file as expected by cryoassess."""
    micsTable = Table(columns=['rlnMicrographName'])
    for mic in newMics:
      print('mic relPath: ', self._getRelPath(mic.getFileName()))
      micsTable.addRow(self._getRelPath(mic.getFileName()))
    with open(self.getInputFilename(numPass), 'w') as f:
      f.write("# Star file generated with Scipion\n")
      micsTable.writeStar(f, tableName='')
    self.appendTotalInputStar(numPass)

  def runMicAssessStep(self, numPass):
    """ Call cryoassess with the appropriate parameters. """
    params = ' '.join(self._getArgs(numPass))
    program = Plugin.getProgram('micassess')
    self.runJob(program, params, env=Plugin.getEnviron(),
                cwd=self._getTmpPath(), numberOfThreads=1)
    self.appendTotalOutputStar(numPass)
    self.appendProbsFile()
    self.copyMicAssessOutput()

  def createOutputStep(self, newMics, numPass):
    outMics = self._loadOutputSet(SetOfMicrographs, self.outputMicsName + '.sqlite', mode='mics')
    outMovies = self._loadOutputSet(SetOfMovies, self.outputMoviesName + '.sqlite', mode='movies')

    # Parse output file and find good mics
    goodMicNames = self._getGoodMicFns(numPass)
    if len(goodMicNames):
      self.curGoodList = goodMicNames
      outMics.copyItems(newMics, updateItemCallback=self._addGoodMic)
      outMovies.copyItems(self._getInputMovies(), updateItemCallback=self._addGoodMovie)

      self._updateOutputSet(self.outputMicsName, outMics)
      self._updateOutputSet(self.outputMoviesName, outMovies)

  def closeSetStep(self):
    outAvs = self._loadOutputSet(SetOfMicrographs, self.outputAveragesName + '.sqlite', mode='averages')
    outMics = self._loadOutputSet(SetOfMicrographs, self.outputMicsName + '.sqlite', mode='mics')
    outMovies = self._loadOutputSet(SetOfMovies, self.outputMoviesName + '.sqlite', mode='movies')

    self._updateOutputSet(self.outputAveragesName, outAvs, state=Set.STREAM_CLOSED)
    self._updateOutputSet(self.outputMicsName, outMics, state=Set.STREAM_CLOSED)
    self._updateOutputSet(self.outputMoviesName, outMovies, state=Set.STREAM_CLOSED)

    self._defineSourceRelation(self._getInputMovies(), getattr(self, self.outputAveragesName))
    self._defineSourceRelation(self._getInputMovies(), getattr(self, self.outputMicsName))
    self._defineSourceRelation(self._getInputMovies(), getattr(self, self.outputMoviesName))
    self.ended = True

  # --------------------------- UTILS functions -----------------------------
  def createMicrograph(self, micPath, movie):
    #print('Flipping Y output micrograph')
    newMcFn = flipYImage(micPath, outFn=micPath)
    movie.setMicName(newMcFn)

    mic = Micrograph(micPath)
    mic.copyObjId(movie)
    mic.setMicName(movie.getMicName())
    return mic

  def _processMovie(self, movie):
    movieFolder = self._getOutputMovieFolder(movie)
    os.mkdir(movieFolder)

    x, y, n = movie.getDim()
    s0, sN = self._getFrameRange(n, 'sum')

    inputMd = os.path.join(movieFolder, 'input_movie.xmd')
    writeMovieMd(movie, inputMd, s0, sN,
                 useAlignment=(movie.hasAlignment() and self.useAlignment))

    outputMicFn = self._getExtraPath(self._getOutputMicName(movie))

    if self.cropRegion == CROP_ALIGNMENT and movie.hasAlignment():
      roi = movie.getAlignment().getRoi()
    elif self.cropRegion == CROP_NEW:
      roi = [self.cropOffsetX.get(), self.cropOffsetY.get(),
             self.cropDimX.get(), self.cropDimY.get()]
    else:
      roi = None

    gainFn = self.inputMovies.get().getGain()
    ext = pwutils.getExt(self.inputMovies.get().getFirstItem().getFileName()).lower()
    if self.inputMovies.get().getGain() and ext in ['.tif', '.tiff']:
      self.flipY = True
      inGainFn = self.inputMovies.get().getGain()
      gainFn = flipYImage(inGainFn, outDir=self._getExtraPath())

    self.averageMovie(movie, inputMd, outputMicFn, self.binFactor.get(),
                      roi, self.inputMovies.get().getDark(), gainFn,
                      splineOrder=self.INTERP_MAP[self.splineOrder.get()])
    return outputMicFn

  def _getStreamingBatchSize(self):
    if self.lastRound:
      return 1
    else:
      return self.streamingBatchSize.get()

  def getGoodInputMics(self):
    inpMics = self._getAveragedMicrographs()
    goodMics = SetOfMicrographs()
    goodMicFns = self._getGoodMicFns('')
    for inpMic in inpMics:
      if inpMic.getFileName() in goodMicFns:
        goodMics.append(inpMic)
    return goodMics

  def copyMicAssessOutput(self):
    copyTree(self._getTmpPath('MicAssess'), self._getExtraPath('MicAssess'))

  def initTotalStars(self):
    totalInputStarFn, totalOutputStar = self.getInputFilename(''), self.getOutputFilename('')
    sameTxt = 'data_\n\nloop_\n_rlnMicrographName \n'
    f1, f2 = open(totalInputStarFn, 'w'), open(totalOutputStar, 'w')
    f1.write("# Star file generated with Scipion\n\n")
    f1.write(sameTxt), f2.write(sameTxt)
    f1.close(), f2.close()

  def appendTotalInputStar(self, numPass):
    totalStarFn = self.getInputFilename('')
    newMicNames = self._getInputMicFns(numPass)
    if os.path.exists(totalStarFn):
      with open(totalStarFn, 'a') as f:
        for micName in newMicNames:
          f.write(' ' + micName + '\n')

  def appendTotalOutputStar(self, numPass):
    totalStarFn = self.getOutputFilename('')
    newMicNames = self._getGoodMicFns(numPass)
    if os.path.exists(totalStarFn):
      with open(totalStarFn, 'a') as f:
        for micName in newMicNames:
          f.write(' ' + micName + '\n')

  def appendProbsFile(self):
    extraFile, tmpFile = self._getExtraPath('probs.tsv'), self._getTmpPath('probs.tsv')
    optionW = 'a' if os.path.exists(extraFile) else 'w'
    f = open(extraFile, optionW)
    with open(tmpFile, 'r') as filex:
      for line in filex:
        f.write(line)
    f.close()

  def addDoneMovieFns(self, newMovies):
    for newMovie in newMovies:
      self.doneMovieFns.add(newMovie.getFileName())

  def addDoneMicFns(self, newMics):
    for newMic in newMics:
      self.doneMicFns.add(newMic.getFileName())

  def checkIfParentFinished(self):
    inpMovies = self._getInputMovies()
    inpMovies.loadAllProperties()
    if not inpMovies.isStreamOpen():
      return True
    return False

  def checkIfMoviesPending(self, newMovies):
    if len(newMovies)==0 and len(self.doneMovieFns) == len(self._getAveragedMicrographs()):
      return False
    return True

  def _getFirstJoinStepName(self):
    # This function will be used for streaming, to check which is
    # the first function that need to wait for all micrographs
    # to have completed
    return 'closeSetStep'

  def _getFirstJoinStep(self):
    for s in self._steps:
      if s.funcName == self._getFirstJoinStepName():
        return s
    return None

  def _updateOutputSet(self, outputName, outputSet, state=Set.STREAM_OPEN):
    outputSet.setStreamState(state)
    if self.hasAttribute(outputName):
      outputSet.write()  # Write to commit changes
      outputAttr = getattr(self, outputName)
      # Copy the properties to the object contained in the protcol
      outputAttr.copy(outputSet, copyId=False)
      # Persist changes
      self._store(outputAttr)
    else:
      # Here the defineOutputs function will call the write() method
      self._defineOutputs(**{outputName: outputSet})
      self._store(outputSet)

    # Close set databaset to avoid locking it
    outputSet.close()

  def _loadOutputSet(self, SetClass, baseName, mode='mics'):
    """
    Load the output set if it exists or create a new one.
    """
    setFile = self._getPath(baseName)
    if os.path.exists(setFile) and os.path.getsize(setFile) > 0:
      outputSet = SetClass(filename=setFile)
      outputSet.loadAllProperties()
      outputSet.enableAppend()
    else:
      outputSet = SetClass(filename=setFile)
      outputSet.setStreamState(outputSet.STREAM_OPEN)
      if mode=='movies':
        label='good movies'
        originSet = self._getInputMovies()
      elif mode=='averages':
        label='averaged movies'
        originSet = self._getInputMovies()
      else:
        label = 'good micrographs'
        originSet = self._getAveragedMicrographs()
      outputSet.setObjLabel(label)
      outputSet.copyInfo(originSet)

    return outputSet

  def _getArgs(self, numPass):
    """ Return the list of args for the command. """
    args = ['-i %s ' % os.path.basename(self.getInputFilename(numPass)),
            '-o %s ' % os.path.basename(self.getOutputFilename(numPass)),
            '-m %s' % Plugin.getVar(CRYOASSESS_MODEL_MIC),
            '-b %d' % self.batchSize.get(),
            '-t %0.2f' % self.threshold.get(),
            '--threads %d' % self.numberOfThreads.get(),
            '--gpus %s' % self.gpuList.get().strip().replace(" ", ",")]

    if self._getCameraType() is not None:
      args.append('-d %s' % self._getCameraType())

    return args

  def _getInputMovies(self):
    return self.inputMovies.get()
  
  def _getAveragedMicrographs(self):
    return self._loadOutputSet(SetOfMicrographs, self.outputAveragesName + '.sqlite', mode='averages')

  def _getNewMovies(self):
    inputMovies = self._getInputMovies()
    newMovies = []
    for movie in inputMovies:
      if movie.getFileName() not in self.doneMovieFns:
        newMovie = movie.clone()
        newMovies.append(newMovie)
    return newMovies

  def _getNewMics(self):
    inputMics = self._getAveragedMicrographs()
    newMics = []
    for mic in inputMics:
      if mic.getFileName() not in self.doneMicFns:
        newMic = mic.clone()
        newMics.append(newMic)
    return newMics

  def getInputFilename(self, numPass):
    if numPass == '':
      return self._getExtraPath('input_micrographs{}.star'.format(numPass))
    else:
      return self._getTmpPath('input_micrographs{}.star'.format(numPass))

  def getOutputFilename(self, numPass):
    if numPass == '':
      return self._getExtraPath('good_micrographs{}.star'.format(numPass))
    else:
      return self._getTmpPath('good_micrographs{}.star'.format(numPass))

  def _getCameraType(self):
    """ Get camera type based on input mic size.
    :return string or None """
    micsizeX, micsizeY, _ = self._getAveragedMicrographs().getDim()
    x = max(micsizeX, micsizeY)
    y = min(micsizeX, micsizeY)
    if math.isclose(x / y, 1.0345, abs_tol=0.001):
      return 'K2'
    elif math.isclose(x / y, 1.4076, abs_tol=0.001):
      return 'K3'
    else:
      return None

  def _getRelPath(self, fn):
    """ Return relative path from cwd=tmp. """
    return os.path.relpath(fn, self._getTmpPath())

  def _getInputMicFns(self, numPass):
    """ Parse input star file and get a list of mics. """
    table = Table(fileName=self.getInputFilename(numPass), tableName='')
    micNames = table.getColumnValues('rlnMicrographName')
    return micNames

  def _getGoodMicFns(self, numPass):
    """ Parse output star file and get a list of good mics. """
    table = Table(fileName=self.getOutputFilename(numPass), tableName='')
    micNames = table.getColumnValues('rlnMicrographName')
    return micNames

  def _addGoodMic(self, item, row):
    """ Callback function to append only good items. """
    if self._getRelPath(item.getFileName()) not in self.curGoodList:
      setattr(item, "_appendItem", False)
      
  def _addGoodMovie(self, item, row):
    """ Callback function to append only good items. """
    if self._getRelPath(self._getExtraPath(self._getOutputMicName(item))) not in self.curGoodList:
      setattr(item, "_appendItem", False)

  # --------------------------- INFO functions ------------------------------
  def _summary(self):
    summary = []
    if not hasattr(self, 'outputMicrographs'):
      summary.append("Output not ready or no good micrographs found")
    else:
      summary.append("Sorted micrographs into good and bad ones.")
    return summary

  def _warnings(self):
    warnings = []
    return warnings

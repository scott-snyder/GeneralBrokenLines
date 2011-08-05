'''
Track fit with general broken lines.

Created on Jul 27, 2011

@author: kleinwrt
'''

import numpy as np
import math
from gblnum import BorderedBandMatrix
from mille import MilleRecord

class GblPoint(object):
  '''
  User supplied point on (initial) trajectory.
  
  Must have jacobians for propagation to previous or next point with offsets (first, last
  point, points with scatterer). May have:
  
    1. Measurement (1D or 2D)
    2. Scatterer (thin, 2D kinks)
    3. Additional local parameters (with derivatives). Fitted together with track parameters.
    4. Additional global parameters (with labels and derivatives). Not fitted, only passed
       on to (binary) file for fitting with Millepede-II. 
  '''
  def __init__(self):
    '''
    Create new point.
    '''
    self.__label = 0
    '''@ivar: label for referencing point (0,1,..,number(points)-1)
       @type: int'''
    self.__prev = 0
    '''@ivar: label of previous point with offset
       @type: int'''
    self.__offset = 0
    '''@ivar: >=0: offset number at point, <0: offset number at next point with offset
       @type: int'''
    self.__next = 0
    '''@ivar: label of next point with offset
       @type: int'''  
    self.__jacobians = []
    '''@ivar: jacobians for propagation to previous or next point with offsets
       @type: list(matrix(float))'''
    self.__measurement = []
    '''@ivar: measurement at point: projection (dm/du), residuals (to initial trajectory), precision
       @type: list(matrix(float))'''
    self.__scatterer = []
    '''@ivar: scatterer at point: initial kinks, precision (inverse covariance matrix)
       @type: list(matrix(float))'''
    self.__localDerivatives = []
    '''@ivar: local derivatives
       @type: list(matrix(float))'''    
    self.__globalLabels = []
    '''@ivar: global labels
       @type: list(matrix(int))'''  
    self.__localDerivatives = []
    '''@ivar: global derivatives
       @type: list(matrix(float))'''      
      
# for extension to retrieval of residuals, pulls    
#    self.__dataMeas = [0, 0]
# for extension to retrieval of residuals, pulls    
#    self.__dataScat = [0, 0]

  def addMeasurement(self, aMeasurement):
    '''
    Add a mesurement to a point.
   
    @param aMeasurement: measurement (projection residuals, precision)
    @type aMeasurement: list(matrix(float))
    '''
    self.__measurement = aMeasurement

  def hasMeasurement(self):
    '''
    Check point for a measurement.
    
    @rtype: bool 
    '''
    return (self.__measurement != [])

  def getMeasurement(self):
    '''
    Retrieve measurement of a point.
    
    @return: measurement (projection residuals, precision)
    @rtype: list(matrix(float))
    '''
    return self.__measurement
    
  def addScatterer(self, aScatterer):
    '''Add a (thin) scatterer to a point.
    
    @param aScatterer: scatterer (kinks, precision)
    @type aScatterer: list(matrix(float))
    '''
    self.__scatterer = aScatterer
  
  def hasScatterer(self):
    '''
    Check point for a scatterer.
    
    @rtype: bool     
    '''
    return (self.__scatterer != [])

  def getScatterer(self):
    '''
    Retrieve scatterer of a point.
    
    @return: scatterer (kinks, precision)
    @rtype: list(matrix(float))    
    '''
    return self.__scatterer
  
  def addLocals(self, derivatives):
    '''
    Add local derivatives.
    
    @param derivatives: local derivatives
    @type derivatives: list(matrix(float))    
    '''
    self.__localDerivatives = derivatives

  def addGlobals(self, labels, derivatives):
    '''
    Add global derivatives.
    
    @param labels: global labels
    @type labels: list(matrix(int))
    @param derivatives: global derivatives
    @type derivatives: list(matrix(float))
    '''
    self.__globalLabels = labels
    self.__globalDerivatives = derivatives

  def getNumLocals(self):
    '''
    Get number of local derivatives.
    
    @return: Number of local derivatives at point
    @rtype: int
    '''
    if (self.__localDerivatives != []):
      return self.__localDerivatives[0].shape[1]
    else:
      return 0
  
  def getLocalDerivatives(self):
    '''
    Get local derivatives.
    
    @return: local derivatives
    @rtype: matrix(float)
    '''
    if (self.__localDerivatives != []):
      return self.__localDerivatives[0]
    else:
      return [ [], [] ]

  def getGlobalLabels(self):
    '''
    Get global labels.
    
    @return: lglobal labels
    @rtype: matrix(int)   
    '''
    if (self.__globalLabels != []):
      return self.__globalLabels[0]
    else:
      return [ [], [] ]
    
  def getGlobalDerivatives(self):
    '''
    Get global derivatives.
    
    @return: global derivatives
    @rtype: matrix(float)   
    '''
    if (self.__globalDerivatives != []):
      return self.__globalDerivatives[0]
    else:
      return [ [], [] ]
    
  def setLabel(self, aLabel):
    '''
    Define label of a point.
    
    @param aLabel: label
    @type aLabel: int
    '''
    self.__label = aLabel

  def getLabel(self):
    '''
    Retrieve label of a point.
    
    @return: label
    @rtype: int
    '''
    return self.__label
 
  def setOffset(self, anOffset, aPrev, aNext):
    '''
    Define offset of a point and references to previous and next point with offset.
    
    @param anOffset: offset number (at point (>=0) or at next point with offset (<0))
    @type anOffset: int
    @param aPrev: label of previous point with offset
    @type aPrev: int
    @param aNext: label of next point with offset
    @type aNext: int
    '''
    self.__offset = anOffset
    self.__prev = aPrev
    self.__next = aNext

  def getOffset(self):
    '''
    Get offset of a point.
    
    @return: offset number (at point (>=0) or at next point with offset (<0))
    @rtype: int
    '''
    return self.__offset
  
#  def setDataMeas(self, aIndex, aData):
# for extension to retrieval of residuals, pulls    
#    self.__dataMeas[aIndex] = aData

#  def setDataScat(self, aIndex, aData):
# for extension to retrieval of residuals, pulls    
#    self.__dataScat[aIndex] = aData
        
  def queryJacobians(self):
    '''
    Query point for labels of enclosing offsets.
    
    @return: labels of previous and next point with offsets
    @rtype: pair(int)
    '''
    return [ self.__prev, self.__next ] 
    
  def addJacobians(self, twoJacobians):
    '''
    Add jacobians to enclosing offsets.
    
    @param twoJacobians: jacobians for propagation to previous and next point with offsets
    @type twoJacobians: pair(matrix(float))
    '''
    self.__jacobians = twoJacobians 
    
  def getDerivatives(self, index):
    '''
    Get derivatives for locally linearized track model (backward or forward propagation).
    
    @param index: 0 (previous) or 1 (next point with offsets)
    @type index: int
    @return: derivatives
    @rtype: list(matrix(float))
    '''
    aJacobian = self.__jacobians[index]
    matJ = aJacobian[3:5, 3:5] # J
    matS = aJacobian[3:5, 1:3] # S 
    vecD = aJacobian[3:5, 0:1] # d
    if (index < 1):
      matS = -matS
    matW = np.linalg.inv(matS) # W = +/- S^-1
    return [ matW, np.dot(matW, matJ), np.dot(matW, vecD) ] # W, W*J, W*d
    
  def printPoint(self):
    '''
    Print point.
    '''
    print " point ", self.__label, self.__offset, self.__prev, self.__next 

#------------------------------------------------------------------------------ 

class GblData(object):
  '''
  Data (block) containing value, precision and derivatives for measurements and kinks.
  
  Created from attributes of GblPoints, used to construct linear equation system for track fit.
  '''
  def __init__(self):
    '''
    Create new data.
    '''
    self.__label = 0
    '''@ivar: label of corresponding point
       @type: int'''
    self.__value = 0.
    '''@ivar: value (residual or kink)
       @type: float'''
    self.__precision = 0.
    '''@ivar: precision (diagonal element of inverse covariance matrix)
       @type: float'''
    self.__downWeight = 1.
    '''@ivar: down weighting factor (M-estimators)
       @type: float'''
    self.__prediction = 0.
    '''@ivar: prediction (for value from fit)
       @type: float'''
    self.__parameters = []
    '''@ivar: labels of fit parameters (with non zero derivative)
       @type: list(int)'''
    self.__derivatives = []
    '''@ivar: derivatives (prediction vs fit parameters)
       @type: list(float)'''
    self.__globalLabels = []
    '''@ivar: labels of global parameters
       @type: list(int)'''
    self.__globalDerivatives = []
    '''@ivar: derivatives (prediction vs global parameters)
       @type: list(float)'''

#    self.__predictionVar = 0.
               
  def addDerivatives(self, aLabel, aMeas, aPrec, parCurv, derCurv, \
               aDim, firstParBand, derBand, firstParLocal=0, \
               derLocal=[], labGlobal=[], derGlobal=[]):
    '''
    Add derivatives to data (block). Generate lists of labels.
    
    @param aLabel: label of corresponding point
    @type aLabel: int
    @param aMeas: value
    @type aMeas: float
    @param aPrec: precision
    @type aPrec: float
    @param parCurv: label for 'curvature' parameter (with (1) or without (0) 'curvature')
    @type parCurv: int
    @param derCurv: derivative vs 'curvature'
    @type derCurv: list(float)
    @param firstParBand: label of first band parameter (offset)
    @type firstParBand: int
    @param derBand: derivatives vs band parameters (offsets)
    @type derBand: list(float)
    @param firstParLocal: label of first local parameter
    @type firstParLocal: int
    @param labGlobal: labels of global parameters
    @type labGlobal: list(int)
    @param derGlobal: derivatives vs global parameters
    @type derGlobal: list(float)
    '''   
    self.__label = aLabel
    self.__value = aMeas
    self.__precision = aPrec
    if (parCurv > 0):
      if (derCurv[0] != 0.):                  # curvature derivatives
        self.__derivatives.append(derCurv[0])
        self.__parameters.append(parCurv)
    for i in range(len(derLocal)):           # local derivatives
      if (derLocal[i] != 0.):
        self.__derivatives.append(derLocal[i])
        self.__parameters.append(firstParLocal + i)
    iPar = firstParBand
    for der in derBand:                      # offset derivatives
      for i in aDim:
        if (der[i] != 0.):
          self.__derivatives.append(der[i])
          self.__parameters.append(iPar)
        iPar += 1  
    for i in range(len(labGlobal)):          # global derivatives
      if (derGlobal[i] != 0.):
        self.__globalLabels.append(labGlobal[i])
        self.__globalDerivatives.append(derGlobal[i])

  def getMatrices(self):  
    '''
    Calculate compressed matrix and right hand side from data.
    
    @return: indices, compressed right hand side and matrix
    @rtype: list
    '''
    aVector = np.array([ self.__derivatives ])
    aMatrix = np.dot(aVector.T, aVector)
    aValue = self.__value
    aWeight = self.__precision * self.__downWeight
    return [self.__parameters, aValue * aVector * aWeight, aMatrix * aWeight]
    
  def setPrediction(self, aVector):
    '''
    Calculate prediction for data from fit.
   
    @param aVector: values of fit parameters
    @type aVector: vector(float)
    '''
    self.__prediction = 0.
    for i in range(len(self.__parameters)):
      self.__prediction += self.__derivatives[i] * aVector[ self.__parameters[i] - 1 ]

#  def setPredictionVariance(self, aMatrix):
#    '''Calculate variance of prediction for data from fit.'''
# var(residual) = 1./precision**2 - var(prediction)
#    aBlockMatrix = aMatrix.getBlockMatrix(self.__parameters)
#    self.__predictionVar = np.dot(self.__derivatives.T, \
#                             np.dot(aBlockMatrix, self.__derivatives))
           
  def setDownWeighting(self, aMethod): 
    '''
    Outlier down weighting with M-estimators.
    
    @param aMethod: method (1=Tukey, 2=Huber, 3=Cauchy) 
    @type aMethod: int
    @return: weight (0..1)
    @rtype: float 
    '''
    scaledResidual = abs(self.__value - self.__prediction) * math.sqrt(self.__precision)   
    if (aMethod == 1):  # Tukey
      if (scaledResidual < 4.6851):
        aWeight = (1.0 - 0.045558 * scaledResidual ** 2) ** 2
      else:
        aWeight = 0.
    elif (aMethod == 2): # Huber
      if (scaledResidual < 1.345):
        aWeight = 1.
      else:
        aWeight = 1.345 / scaledResidual
    elif (aMethod == 3):  # Cauchy
      aWeight = 1.0 / (1.0 + (scaledResidual / 2.3849) ** 2)      
    self.__downWeight = aWeight
    return aWeight
        
  def getChi2(self):
    '''
    Calculate Chi2 (contribution) from data.
    
    @return: Chi2
    @rtype: float
    '''
    Chi2 = (self.__value - self.__prediction) ** 2 * self.__precision * self.__downWeight
    return Chi2
 
  def toRecord(self):
    '''
    Get data components (for copying to MP binaty record)
    
    @return: data components
    @rtype: list
    '''
    return [self.__value, self.__precision, self.__parameters, self.__derivatives, \
            self.__globalLabels, self.__globalDerivatives]

  def fromRecord(self, dataList):
    '''
    Set data components (from MP binaty record)
   
    @param dataList: data components
    @type dataList: list
    '''
    [self.__value, self.__precision, self.__parameters, self.__derivatives, \
            self.__globalLabels, self.__globalDerivatives] = dataList
            
  def analyzeData(self, maxBand):
    '''
    Analyze labels of fit parameters to determine number of parameters and
    border size with given maximal band width.
    
    @param maxBand: maximal band width
    @type maxBand: int   
    @return: number of parameters and border size (from this data)
    @rtype: pair(int) 
    '''
    maxPar = self.__parameters[-1]
    maxBor = 0
    for i in self.__parameters:
      if (i < maxPar - maxBand):
        maxBor = i
    return [maxPar, maxBor]
         
  def printData(self):
    '''
    Print data.
    '''
    print " meas. ", self.__label, self.__value, self.__precision
    print " param ", self.__parameters
    print " deriv ", self.__derivatives
 
#------------------------------------------------------------------------------ 
    
class GblTrajectory(object):
  '''
  General Broken Lines Trajectory.
  ================================
  
  For a track with an initial trajectory from a prefit of the
  measurements (internal seed) or an external prediction
  (external seed) the description of multiple scattering is 
  added by offsets in a local system. Along the initial
  trajectory points are defined with can describe a measurement
  or a (thin) scatterer or both. The refit provides corrections
  to the local track parameters (in the local system) and the 
  corresponding covariance matrix at any of those points.

  The broken lines trajectory is defined by (2D) offsets at the 
  first and last point and all points with a scatterer. The
  prediction for a measurement is obtained by interpolation of
  the enclosing offsets and for triplets of adjacent offsets
  kink angles are determined. This requires for all points the
  jacobians for propagation to the previous and next offset.

  Additional local or global parameters can be added and the
  trajectories can be written to special binary files for
  calibration and alignment with Millepede-II.
  (V. Blobel, NIM A, 566 (2006), pp. 5-13).

  The conventions for the coordinate systems follow:
  Derivation of Jacobians for the propagation of covariance
  matrices of track parameters in homogeneous magnetic fields
  A. Strandlie, W. Wittek, NIM A, 566 (2006) 687-698.
  
  Calling sequence:
  =================
    1. Create trajectory::
            traj = GblTrajectory()
    2. For all points on initial trajectory 
        - Create points and add appropiate attributes::
            point = GblPoint()
            point.addMeasurement(..)
            point.addScatterer(..)
            point.addLocals(..)
            point.addGlobals(..)
        - Add point (ordered by arc length) to trajectory, get label of point::  
            label = traj.addPoint(point)
    3. Define (points with) offsets::
            traj.defineOffsets()
    4. Optionally add external seed::
            traj.addExternalSeed(..)
    5. For all points on initial trajectory
        - Query for required jacobians::
            traj.queryJacobians(label)         
        - Add requested jacobians::
            traj.addJacobians(label,..)  
    6. Optionally write trajectory to MP binary file::
            traj.milleOut(..)
    7. Fit trajectory, bet Chi2, Ndf (and weight lost by M-estimators)::
            [..] = traj.fit()
    8. For any point on inital trajectory
        - Get corrections and covariance matrix for track parameters::
            [..] = traj.getResults(label) 
            
  Alternatively trajectories can by read from MP binary files and fitted. 
  As the points on the initial trajectory are not stored in this files results at
  points (corrections, covariance matrix) are not available.
  ''' 
  def __init__(self, hasCurv=True, aDim=[0, 1]):
    '''
    Create new trajectory.
    '''  
    self.__numPoints = 0
    '''@ivar: number of points on trajectory
       @type: int'''
    self.__numOffsets = 0
    '''@ivar: number of (points with) offsets on trajectory
       @type: int'''
    self.__numCurvature = (1 if hasCurv else 0)
    '''@ivar: 'curvature' is fit parameter (=1)
       @type: int'''
    self.__numParameters = 0
    '''@ivar: number fit parameters
       @type: int'''
    self.__numLocals = 0
    '''@ivar: number of local parameters
       @type: int'''    
    self.__externalPoint = 0
    '''@ivar: label of point with external seed
       @type: int'''    
    self.__externalNdf = 0
    '''@ivar: Ndf from external seed
       @type: int'''    
    self.__dimensions = aDim
    '''@ivar: active compoents of offsets (both ([0,1]) or single ([0] or [1])
       @type: list(int)'''
    self.__points = [] 
    '''@ivar: points on trajectory
       @type: list(GblPoint)'''
    self.__data = []
    '''@ivar: data (blocks) of trajectory
       @type: list(GblData)'''
    self.__externalSeed = []
    '''@ivar: external seed (for local, fit parameters)
       @type: list(matrix(float))'''
    self.__externalIndex = []
    '''@ivar: labels of fit parameters from external seed
       @type: list(int)'''

  def addPoint(self, point):
    '''
    Add point to trajectory. Points have to be ordered in arc length.
    
    @param point: point to be added
    @type point: GblPoint
    @return: label of added point
    @rtype: int    
    '''
    self.__numPoints += 1
    label = self.__numPoints
    point.setLabel(label)
    self.__points.append(point)
    self.__numLocals = max(self.__numLocals, point.getNumLocals())
    return label
    
  def getNumPoints(self):
    '''
    Get number of points on trajectory.
    
    @return: number of points
    @rtype: int
    '''  
    return self.__numPoints
  
  def addExternalSeed(self, aLabel, aSeed):
    '''
    Add external seed to trajectory.
    
    @param aLabel: label of point with external seed
    @type aLabel: int 
    @param aSeed: seed (covariance matrix of track parameters at point)
    @type aSeed: list(matrix(float)) 
    '''
    self.__externalPoint = aLabel
    self.__externalSeed = aSeed


  def queryJacobians(self, aLabel):
    '''
    Query point for adjacent scatterers.
    
    @param aLabel: label of point
    @type aLabel: int
    @return: labels of previous and next point with offset
    @rtype: pair(int)
    '''
    return self.__points[aLabel - 1].queryJacobians() 
    
  def addJacobians(self, aLabel, twoJacobians):
    '''
    Provide point for jacobians to adjacent scatterers.
    
    @param aLabel: label of point
    @type aLabel: int
    @param twoJacobians: jacobians for propagation to previous and next point with offset
    @type twoJacobians: pair(matrix(float))
    '''
    self.__points[aLabel - 1].addJacobians(twoJacobians) 
    
  def defineOffsets(self):
    '''
    Define offsets from list of points.
    '''

# built list of labels defining offsets
#   first point is offset    
    labels = [ self.__points[0].getLabel() ]
#   intermediate scatterers are offsets    
    for aPoint in self.__points[1:-1]:
      if (aPoint.hasScatterer()):
        labels.append(aPoint.getLabel()) 
#   last point is offset    
    labels.append(self.__points[-1].getLabel())
# set labels for previous/next offsets
#   first point is offset    
    self.__points[0].setOffset(0, labels[0], labels[1])
    nOffsets = 1
#   intermediate scatterers are offsets    
    for aPoint in self.__points[1:-1]:
      if (aPoint.hasScatterer()):
        aPoint.setOffset(nOffsets, labels[nOffsets - 1], labels[nOffsets + 1])
        nOffsets += 1
      else:  
        aPoint.setOffset(-nOffsets, labels[nOffsets - 1], labels[nOffsets])
#   last point is offset    
    self.__points[-1].setOffset(nOffsets, labels[nOffsets - 1], labels[nOffsets])
    self.__numOffsets = nOffsets + 1
    self.__numParameters = self.__numOffsets * len(self.__dimensions) \
                           + self.__numCurvature + self.__numLocals

  def dump(self):
    '''
    Dump trajectory.
    '''
    for p in self.__points:
      p.printPoint() 
      
  def milleOut(self, aFile):
    '''
    Write (data blocks of) trajectory to MP (binary) file.
    
    @param aFile: MP file
    @type aFile: file
    '''
    rec = MilleRecord()
#   data: measurements and kinks        
    for aData in self.__data:       
      rec.addData(aData.toRecord())
#   (compressed) external seed
    compressedIndex = []
    compressedSeed = []
    if (self.__externalIndex != []):
      aMatrix = self.__externalSeed[-1]    
      for i in range(len(self.__externalIndex)):
        iIndex = self.__externalIndex[i]
        for j in range(i + 1):
          jIndex = self.__externalIndex[j] 
          if (aMatrix[i, j] != 0.):
            ijIndex = iIndex * (iIndex - 1) / 2 + jIndex
            compressedIndex.append(ijIndex)
            compressedSeed.append(aMatrix[i, j])  
      rec.addSpecial(compressedIndex, compressedSeed, 1) 
                    
    rec.writeRecord(aFile)

  def milleIn(self, aFile):
    '''
    Read (data blocks of) trajectory from MP (binary) file.
    
    @param aFile: MP file
    @type aFile: file
    '''
    rec = MilleRecord()
    rec.readRecord(aFile)
    mPar = 0
    mBor = 0
    mBand = 3 * len(self.__dimensions) - 1  # max band width
    while (rec.moreData()):
      aTag = rec.specialDataTag()
      if (aTag < 0):
# get data
        aData = GblData()
        aData.fromRecord(rec.getData())
        self.__data.append(aData)
        [nPar, nBor] = aData.analyzeData(mBand)
        mPar = max(mPar, nPar)
        mBor = max(mBor, nBor)
      elif (aTag == 1):
# get seed        
        [compressedIndex, compressedSeed ] = rec.getSpecial()  
# uncompress index
        iIndex = 0
        iiIndex = 0
        rowList = []
        for ind in compressedIndex:
          while (ind > iiIndex):
            iIndex += 1
            iiIndex += iIndex
          rowList.append(iIndex)
          if (self.__externalIndex == []):
            self.__externalIndex.append(iIndex)
          elif (self.__externalIndex[-1] != iIndex):
            self.__externalIndex.append(iIndex)
# uncompress Matrix              
        nExtPar = len(self.__externalIndex)
        aMatrix = np.zeros((nExtPar, nExtPar))
        for i in range(len(compressedIndex)):
          iRow = rowList[i] - 1
          iCol = compressedIndex[i] - (iRow + 1) * iRow / 2 - 1
          aMatrix[iRow, iCol] = compressedSeed[i]
          if (iRow != iCol):
            aMatrix[iCol, iRow] = compressedSeed[i]
          else:
            self.__externalNdf += 1
        self.__externalSeed.append(aMatrix)
        
    self.__numParameters = mPar
    self.__numLocals = mBor - self.__numCurvature
   
  def __getJacobian(self, aLabel):
    '''
    Get jacobian for transformation from fit to track parameters at point.
    
    @param aLabel: label of point
    @type aLabel: int
    @return: labels of required fit parameters and jacobian
    @rtype: list
    '''
    def addMatrices(indexList, matList):
      '''
      Add matrices to jacobian.
      
      @param indexList: indices
      @type indexList: list
      @param matList: matrices
      @type matList: list
      '''
      for i in range(nDim):
        vecIndex.append(iPar + iOff + i + 1)
        for j in range(2):
          for k in range(len(indexList)):
            aJacobian[indexList[k] + j, iPar + i] = matList[k][j, aDim[i]] 

    def addVectors(indexList, vecList):
      '''
      Add vectors to jacobian.
      
      @param indexList: indices
      @type indexList: list
      @param vecList: vectors
      @type vecList: list     
      '''
      for j in range(2):
        for k in range(len(indexList)):
          aJacobian[indexList[k] + j, iPar] = vecList[k][j]

    aDim = self.__dimensions
    nDim = len(aDim)
    anIndex = abs(aLabel) - 1
#   check consistency of (index, direction)    
    if (aLabel > 0):
      nJacobian = 1
      if (anIndex >= self.__numPoints - 1):
        anIndex = self.__numPoints - 1
        nJacobian = 0
    else:
      nJacobian = 0
      if (anIndex <= 0):
        anIndex = 0
        nJacobian = 1
# Jacobian broken lines (q/p,..,u_i,u_i+1..) to local (q/p,u',u) parameters   
    nCurv = self.__numCurvature
    nLocals = self.__numLocals     
    nBorder = nCurv + nLocals
    nParBRL = nBorder + 2 * nDim
    nParLoc = nLocals + 5
    aJacobian = np.zeros((nParLoc, nParBRL))
    aPoint = self.__points[anIndex]
    nOffset = aPoint.getOffset() 
    vecIndex = []
    for i in range(nCurv):
      aJacobian[i, i] = 1.0
      vecIndex.append(i + 1)
    for i in range(nLocals):
      aJacobian[i + 5, i + nCurv] = 1.0
      vecIndex.append(i + 1 + nCurv)           
    if (nOffset < 0): # need interpolation
      [ prevW, prevWJ, prevWd ] = aPoint.getDerivatives(0) # W-, W- * J-, W- * d-
      [ nextW, nextWJ, nextWd ] = aPoint.getDerivatives(1) # W+, W+ * J+, W+ * d+
      sumWJ = prevWJ + nextWJ
      matN = np.linalg.inv(sumWJ)      # N = (W- * J- + W+ * J+)^-1 
#     derivatives for u_int      
      prevNW = np.dot(matN, prevW)     # N * W-
      nextNW = np.dot(matN, nextW)     # N * W+
      prevNd = np.dot(matN, prevWd)    # N * W- * d-
      nextNd = np.dot(matN, nextWd)    # N * W+ * d+
#     derivatives for u'_int
      prevWPN = np.dot(nextWJ, prevNW) # W+ * J+ * N * W-
      nextWPN = np.dot(prevWJ, nextNW) # W- * J- * N * W+
      prevWNd = np.dot(nextWJ, prevNd) # W+ * J+ * N * W- * d-
      nextWNd = np.dot(prevWJ, nextNd) # W- * J- * N * W+ * d+
      iPar = self.__numCurvature__
      iOff = nDim * (-nOffset - 1) # first offset ('i' in u_i)     
      iPar = 0
      if (self.__numCurvature > 0):
        addVectors([1, 3], [-nextWNd + prevWNd, -nextNd - prevNd]) # from curvature
        iPar += 1
      iPar += nLocals  
      addMatrices([1, 3], [-prevWPN, prevNW]) # from 1st offset
      iPar += nDim  
      addMatrices([1, 3], [ nextWPN, nextNW]) # from 2nd offset
    else:
      [ matW, matWJ, vecWd ] = aPoint.getDerivatives(nJacobian) # W, W * J, W * d
      mat1 = np.eye(2)
      iOff = nDim * (nOffset + nJacobian - 1) # first offset ('i' in u_i)     
      iPar = 0
      if (self.__numCurvature > 0):
        addVectors([1], [-vecWd])         # from curvature
        iPar += 1
      iPar += nLocals  
      addMatrices([1, 3], [-matWJ, mat1]) # from 1st offset
      iPar += nDim
      addMatrices([1], [matW])            # from 2nd offset 
    return [vecIndex, aJacobian]
 
  def getResults(self, aLabel):
    '''
    Get results (corrections, covarinace matrix) at point in forward or backward direction.
    
    @param aLabel: signed label of point (<0 backward, >0 forward)
    @type aLabel: int
    @return: correction vector, covarinace matrix for track parameters
    @rtype: list
    '''
    [ anIndex, aJacobian ] = self.__getJacobian(aLabel)
    nParBRL = len(anIndex)
    aVec = np.empty(nParBRL)
    for i in range(nParBRL):
      aVec[i] = self.__vector[anIndex[i]]         # compressed vector
    aMat = self.__matrix.getBlockMatrix(anIndex)  # compressed matrix
    locPar = np.dot(aJacobian, aVec) 
    locCov = np.dot(aJacobian, np.dot(aMat, aJacobian.T))
    return [ locPar, locCov ]
  
  def fit(self, optionList=""):
    '''
    Perform fit of trajectory.
    
    @param optionList: M-estimators to be used (one iteration per character)
    @type optionList: string
    @return: Chi2, Ndf, loss of weight from fit ([0., -1, 0.] if fit failed)
    @rtype: list
    '''
    def prepare():
      '''
      Prepare fit; generate data from points.
      '''
      aDim = self.__dimensions
      nCurv = self.__numCurvature
      nBorder = nCurv + self.__numLocals
# measurements
      for aPoint in self.__points:
        if (aPoint.hasMeasurement()):
          nOffset = aPoint.getOffset()
          nLabel = aPoint.getLabel()
          localDer = aPoint.getLocalDerivatives()
          globalLab = aPoint.getGlobalLabels()
          globalDer = aPoint.getGlobalDerivatives()
          [ matP, aMeas, aPrec ] = aPoint.getMeasurement()
          if (nOffset < 0): # need interpolation
            [ prevW, prevWJ, prevWd ] = aPoint.getDerivatives(0) # W-, W- * J-, W- * d-
            [ nextW, nextWJ, nextWd ] = aPoint.getDerivatives(1) # W+, W+ * J+, W+ * d+
            sumWJ = prevWJ + nextWJ         # W- * J- + W+ * J+
            sumWd = prevWd + nextWd         # W+ * d+ + W- * d-
            matN = np.linalg.inv(sumWJ)     # N = (W- * J- + W+ * J+)^-1
            matPN = np.dot(matP, matN)      # P * N
            prevPN = np.dot(matPN, prevW)   # P * N * W-
            nextPN = np.dot(matPN, nextW)   # P * N * W+
            vecPNd = np.dot(matPN, sumWd)   # P * N * (W+ * d+ + W- * d-)
            for i in range(2):
              iPar = (-nOffset - 1) * len(aDim) + nBorder + 1
              if (aPrec[i] > 0.):
                aData = GblData()
                aData.addDerivatives(nLabel, aMeas[i], aPrec[i], nCurv, -vecPNd[i], \
                            aDim, iPar, [prevPN[i], nextPN[i]], nCurv + 1, localDer[i], \
                            globalLab[i], globalDer[i])
                self.__data.append(aData)
#                aPoint.setDataMeas(i, len(self.__data)) 
          else:
            for i in range(2):
              iPar = nOffset * len(aDim) + nBorder + 1
              if (aPrec[i] > 0.):
                aData = GblData()
                aData.addDerivatives(nLabel, aMeas[i], aPrec[i], 0, [0.], \
                                     aDim, iPar, [matP[i]], nCurv + 1, localDer[i], \
                                     globalLab[i], globalDer[i])
                self.__data.append(aData)
#                aPoint.setDataMeas(i, len(self.__data)) 
# pseudo measurements from kinks
      for aPoint in self.__points[1:-1]:
        if (aPoint.hasScatterer()):
          nOffset = aPoint.getOffset()
          nLabel = aPoint.getLabel()        
#          print " kink ", nLabel, nOffset
          [ prevW, prevWJ, prevWd ] = aPoint.getDerivatives(0) # W-, W- * J-, W- * d-
          [ nextW, nextWJ, nextWd ] = aPoint.getDerivatives(1) # W+, W+ * J+, W+ * d+
          sumWJ = prevWJ + nextWJ         # W- * J- + W+ * J+
          sumWd = prevWd + nextWd         # W+ * d+ + W- * d-
          [ aMeas, aPrec ] = aPoint.getScatterer()
          for i in aDim:
            iPar = (nOffset - 1) * len(aDim) + nBorder + 1
            if (aPrec[i] > 0.):
              aData = GblData()
              aData.addDerivatives(nLabel, aMeas[i], aPrec[i], nCurv, -sumWd[i], \
                                   aDim, iPar, [prevW[i], -sumWJ[i], nextW[i]])            
              self.__data.append(aData)
#              aPoint.setDataScat(i, len(self.__data)) 
#     external seed
      if (self.__externalPoint != 0):
        [ index, aJacobian ] = self.__getJacobian(self.__externalPoint)
        self.__externalIndex = index
        aMatrix = np.dot(aJacobian.T, np.dot(self.__externalSeed[0] , aJacobian))
        self.__externalSeed.append(aMatrix)   # add transformed matrix
        self.__externalNdf = 0   # assume Ndf = number of non zero diagonal elements
        for i in range(len(index)):
          if (aMatrix[i, i] != 0.):
            self.__externalNdf += 1 

    def buildLinearEquationSystem():
      '''
      Build linear equation system from data.
      '''
      nBorder = self.__numCurvature + self.__numLocals
      self.__matrix = BorderedBandMatrix(self.__numParameters, nBorder)
      self.__vector = np.zeros(self.__numParameters)
      for aData in self.__data: 
        [ index, aVector, aMatrix ] = aData.getMatrices()
        for i in range(len(index)):
          self.__vector[ index[i] - 1 ] += aVector[0, i]        # update vector
        self.__matrix.addBlockMatrix(index, aMatrix)            # update matrix
      if (self.__externalSeed != []):                           # external seed
        index = self.__externalIndex
        aMatrix = self.__externalSeed[-1]
        self.__matrix.addBlockMatrix(index, aMatrix)            # update matrix

    def downWeight(aMethod):
      '''
      Down weight (data) outliers.
      
      @param aMethod: M-estimator
      @type: int
      @return: loss of weight (sum(1-down_weighting))
      @rtype: float
      '''
      aLoss = 0.
      for aData in self.__data: 
        aLoss += (1. - aData.setDownWeighting(aMethod))
      return aLoss
 
    def predict():
      '''
      Calculate predictions.
      '''
      for aData in self.__data: 
        aData.setPrediction(self.__vector)
        
    def externalChi2():
      '''
      Chi2 from external seed.
      
      @return: Chi2
      @rtype: float
      '''
      extChi2 = 0.
      nParSeed = len(self.__externalIndex)
      if (nParSeed > 0):
        aVector = np.empty(nParSeed)
        aMatrix = self.__externalSeed[-1]
        for i in range(nParSeed):
          aVector[i] = self.__vector[self.__externalIndex[i] - 1]
          extChi2 = np.dot(aVector.T, np.dot(aMatrix, aVector))
      return extChi2
                          
    prepare()
    buildLinearEquationSystem()
#
    try:
      aMethod = 0
      lostWeight = 0.
      self.__vector = self.__matrix.solveAndInvertBorderedBand(self.__vector)
      predict()
      
      for o in optionList:    # down weighting iterations    
        try:
          aMethod = "THC".index(o.upper()) + 1
          lostWeight = downWeight(aMethod)
          buildLinearEquationSystem()
          self.__vector = self.__matrix.solveAndInvertBorderedBand(self.__vector)
          predict()
        except ValueError:
          pass                  
             
      Ndf = len(self.__data) - self.__numParameters + self.__externalNdf
      Chi2 = externalChi2()
      for aData in self.__data: 
        Chi2 += aData.getChi2()
      Chi2 /= [1.0, 0.8737, 0.9326, 0.8228 ][aMethod]  
      return [Chi2, Ndf, lostWeight]
    
    except (ZeroDivisionError, np.linalg.linalg.LinAlgError):
      return [ 0., -1, 0.]
    

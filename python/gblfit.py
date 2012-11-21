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
  def __init__(self, aJacobian):
    '''
    Create new point.
    
    @param aJacobian: jacobian from previous point
    @type matrix(float)
    '''
    self.__label = 0
    '''@ivar: label for referencing point (0,1,..,number(points)-1)
       @type: int'''
    self.__offset = 0
    '''@ivar: >=0: offset number at point, <0: offset number at next point with offset
       @type: int'''
    self.__p2pJacobian = aJacobian
    '''@ivar: Point-to-point jacobian from previous point
       @type: matrix(float)'''
    self.__jacobians = [ [], [] ]
    '''@ivar: jacobians for propagation to previous or next point with offsets
       @type: pair(matrix(float))'''
    self.__measurement = None
    '''@ivar: measurement at point: projection (dm/du), residuals (to initial trajectory), precision
       @type: list(matrix(float))'''
    self.__measDim = 0
    '''@ivar: dimension of measurement (2D, 4D or 5D)
       @type: int''' 
    self.__measTransformation = None    
    '''@ivar: transformation (to eigen-vectors of precision matrix)
       @type: matrix(float))''' 
    self.__scatterer = None
    '''@ivar: scatterer at point: initial kinks, precision (inverse covariance matrix)
       @type: matrix(float)'''
    self.__localDerivatives = None
    '''@ivar: local derivatives
       @type: matrix(float))'''    
    self.__globalLabels = None
    '''@ivar: global labels
       @type: matrix(int)'''  
    self.__globalDerivatives = None
    '''@ivar: global derivatives
       @type: matrix(float)'''  
      
# for extension to retrieval of residuals, pulls    
#    self.__dataMeas = [0, 0]
# for extension to retrieval of residuals, pulls    
#    self.__dataScat = [0, 0]

  def addMeasurement(self, aMeasurement):
    '''
    Add a mesurement to a point.
   
    @param aMeasurement: measurement (projection (or None), residuals, precision
                         (diagonal of or full matrix))
    @type aMeasurement: list(matrix(float))
    '''
    self.__measurement = aMeasurement
    self.__measDim = aMeasurement[1].shape[0]
    if (aMeasurement[2].ndim == 2): # full precision matrix, need to diagonalize
      eigenVal, eigenVec = np.linalg.eigh(aMeasurement[2])
      self.__measTransformation = eigenVec.T
#     transform measurement
      if (aMeasurement[0] is None):
        self.__measurement[0] = self.__measTransformation
      else:
        self.__measurement[0] = np.dot(self.__measTransformation, aMeasurement[0])
      self.__measurement[1] = np.dot(self.__measTransformation, aMeasurement[1])
      self.__measurement[2] = eigenVal

  def hasMeasurement(self):
    '''
    Check point for a measurement.
    
    @rtype: bool 
    '''
    return (self.__measurement is not None)

  def getMeasurement(self):
    '''
    Retrieve measurement of a point.
    
    @return: measurement (projection residuals, precision)
    @rtype: list(matrix(float))
    '''
    return self.__measurement
 
  def getMeasDim(self):
    '''
    Retrieve dimension of measurement of a point.
    
    @return: measurement dimension (2, 4 or 5)
    @rtype: int
    '''
    return self.__measDim
     
  def addScatterer(self, aScatterer):
    '''
    Add a (thin) scatterer to a point.
    
    @param aScatterer: scatterer (kinks, precision)
    @type aScatterer: list(matrix(float))
    '''
    self.__scatterer = aScatterer
  
  def hasScatterer(self):
    '''
    Check point for a scatterer.
    
    @rtype: bool     
    '''
    return (self.__scatterer is not None)

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
    @type derivatives: matrix(float)    
    '''
    if (self.__measDim > 0):
      if (self.__measTransformation is None):
        self.__localDerivatives = derivatives
      else:
        self.__localDerivatives = np.dot(self.__measTransformation, derivatives)        

  def addGlobals(self, labels, derivatives):
    '''
    Add global derivatives.
    
    @param labels: global labels
    @type labels: matrix(int)
    @param derivatives: global derivatives
    @type derivatives: matrix(float)
    '''
    if (self.__measDim > 0):
      self.__globalLabels = labels
      if (self.__measTransformation is None):
        self.__globalDerivatives = derivatives
      else:
        self.__globalDerivatives = np.dot(self.__measTransformation, derivatives) 

  def getNumLocals(self):
    '''
    Get number of local derivatives.
    
    @return: Number of local derivatives at point
    @rtype: int
    '''
    if (self.__localDerivatives is not None):
      return self.__localDerivatives.shape[1]
    else:
      return 0
  
  def getLocalDerivatives(self):
    '''
    Get local derivatives.
    
    @return: local derivatives
    @rtype: matrix(float)
    '''
    return self.__localDerivatives

  def getGlobalLabels(self):
    '''
    Get global labels.
    
    @return: lglobal labels
    @rtype: matrix(int)   
    '''
    return self.__globalLabels
    
  def getGlobalDerivatives(self):
    '''
    Get global derivatives.
    
    @return: global derivatives
    @rtype: matrix(float)   
    '''
    return self.__globalDerivatives
    
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
 
  def setOffset(self, anOffset):
    '''
    Define offset of a point and references to previous and next point with offset.
    
    @param anOffset: offset number (at point (>=0) or at next point with offset (<0))
    @type anOffset: int
    '''
    self.__offset = anOffset

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
 
  def getP2pJacobian(self):
    '''
    Retrieve jacobian of a point.
    
    @return: point-to-point jacobian (from previous point)
    @rtype: matrix(float)
    '''
    return self.__p2pJacobian         
 
  def addPrevJacobian(self, aJacobian):
    '''
    Add jacobian to previous offset.
    
    @param aJacobian: jacobian for propagation to previous point with offsets
    @type aJacobian: matrix(float)
    '''
    self.__jacobians[0] = np.linalg.inv(aJacobian)

  def addNextJacobian(self, aJacobian):
    '''
    Add jacobian to next offset.
    
    @param aJacobian: jacobian for propagation to next point with offsets
    @type aJacobian: matrix(float)
    '''
    self.__jacobians[1] = aJacobian
    
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
    return matW, np.dot(matW, matJ), np.dot(matW, vecD)  # W, W*J, W*d
    
  def printPoint(self):
    '''
    Print point.
    '''
    print " point ", self.__label, self.__offset 

#------------------------------------------------------------------------------ 

class GblData(object):
  '''
  Data (block) containing value, precision and derivatives for measurements, kinks
  and external seed.
  
  Created from attributes of GblPoints, used to construct linear equation system for track fit.
  '''
  def __init__(self, aLabel=0, aValue=0., aPrec=0.):
    '''
    Create new data.
    
    @param aLabel: label of corresponding point
    @type aLabel: int
    @param aValue: value
    @type aValue: float
    @param aPrec: precision
    @type aPrec: float
    '''
    self.__label = aLabel
    '''@ivar: label of corresponding point
       @type: int'''
    self.__value = aValue
    '''@ivar: value (residual or kink)
       @type: float'''
    self.__precision = aPrec
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
               
  def addDerivatives(self, iRow, labDer, matDer, \
                     derLocal=None, labGlobal=None, derGlobal=None):
    '''
    Add derivatives to data (block) from measurements and kinks. Generate lists of labels.
    
    @param iRow: row of measurement (vector)
    @type iRow: int
    @param labDer: labels of derivatives vs curvature and band parameters (offsets)
    @type labDer: vector(int)
    @param matDer: derivatives vs curvature and band parameters (offsets)
    @type matDer: matrix(float)
    @param derLocal: derivatives vs local parameters
    @type derLocal: list(float)
    @param labGlobal: labels of global parameters
    @type labGlobal: list(int)
    @param derGlobal: derivatives vs global parameters
    @type derGlobal: list(float)
    '''   
    if (derLocal is not None):
      for i in range(derLocal.shape[1]):           # local derivatives
        if (derLocal[iRow, i] != 0.):
          self.__derivatives.append(derLocal[iRow, i])
          self.__parameters.append(i + 1)
        
    for i in range(len(labDer)):           # curvature, offset derivatives
      if (labDer[i] != 0 and matDer[iRow, i] != 0.):
        self.__derivatives.append(matDer[iRow , i])
        self.__parameters.append(labDer[i])

    if (derGlobal is not None):  
      for i in range(derGlobal.shape[1]):          # global derivatives
        if (derGlobal[iRow, i] != 0.):
          self.__globalLabels.append(labGlobal[iRow, i])
          self.__globalDerivatives.append(derGlobal[iRow, i])

  def addExtDerivatives(self, indexExt, derExt):
    '''
    Add derivatives to data (block) from external seed. Generate lists of labels.
    
    @param indexExt: labels from exteranl seed
    @type indexExt: list(int)
    @param derExt: derivatives from exteranl seed
    @type derExt: list(float)
    '''       
    for i in range(len(derExt)):           # external derivatives
      if (derExt[i] != 0.):
        self.__derivatives.append(derExt[i])
        self.__parameters.append(indexExt[i])  

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
    return self.__parameters, aValue * aVector * aWeight, aMatrix * aWeight
    
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
    return self.__value, self.__precision, self.__parameters, self.__derivatives, \
            self.__globalLabels, self.__globalDerivatives

  def fromRecord(self, dataList):
    '''
    Set data components (from MP binaty record)
   
    @param dataList: data components
    @type dataList: list
    '''
    self.__value, self.__precision, self.__parameters, self.__derivatives, \
            self.__globalLabels, self.__globalDerivatives = dataList
            
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
    return maxPar, maxBor
         
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
  (2D, 4D or 5D) measurements (internal seed) or an external 
  prediction(external seed) the description of multiple scattering
  is added by offsets in a local system. Along the initial
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
        - Create point (supply jacobian from previous point)::
            point = GblPoint(jacobian)
        - Optionally add measurement to point::    
            point.addMeasurement(..)
        - Optionally additional local or global parameters for measurement:: 
            point.addLocals(..)
            point.addGlobals(..)
        - Optionally add scatterer to point::    
            point.addScatterer(..)
        - Add point (ordered by arc length) to trajectory, get label of point::  
            label = traj.addPoint(point)
    3. Optionally add external seed::
            traj.addExternalSeed(..)
    4. Fit trajectory, bet Chi2, Ndf (and weight lost by M-estimators)::
            [..] = traj.fit()
    5. For any point on inital trajectory
        - Get corrections and covariance matrix for track parameters::
            [..] = traj.getResults(label) 
    6. Optionally write trajectory to MP binary file::
            traj.milleOut(..)
            
  Alternatively trajectories can by read from MP binary files and fitted. 
  As the points on the initial trajectory are not stored in this files results at
  points (corrections, covariance matrix) are not available.
  
  References:
  ===========  
    - V. Blobel, C. Kleinwort, F. Meier,
      Fast alignment of a complex tracking detector using advanced track models,
      Computer Phys. Communications (2011), doi:10.1016/j.cpc.2011.03.017
    - C. Kleinwort, General Broken Lines as advanced track fitting method,
      NIM A, 673 (2012), 107-110, doi:10.1016/j.nima.2012.01.024  
      
  ''' 
  def __init__(self, hasCurv=True, aDim=[0, 1]):
    '''
    Create new trajectory.
    
    @param hasCurv: flag for curvature
    @type hasCurv: bool
    @param aDim: active offset components (0:u_1, 1:u_2)
    @type aDim: list(int)
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
    self.__dimensions = aDim
    '''@ivar: active components of offsets (both ([0,1]) or single ([0] or [1])
       @type: list(int)'''
    self.__points = [] 
    '''@ivar: points on trajectory
       @type: list(GblPoint)'''
    self.__data = []
    '''@ivar: data (blocks) of trajectory
       @type: list(GblData)'''
    self.__externalSeed = None
    '''@ivar: external seed (for local, fit parameters)
       @type: matrix(float)'''

  def addPoint(self, point):
    '''
    Add point to trajectory. Points have to be ordered in arc length.
    
    @param point: point to be added
    @type point: GblPoint
    @return: label of added point (1..number(points))
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
    @type aSeed: matrix(float) 
    '''
    self.__externalPoint = aLabel
    self.__externalSeed = aSeed
    
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
        nPar, nBor = aData.analyzeData(mBand)
        mPar = max(mPar, nPar)
        mBor = max(mBor, nBor)
        
    self.__numParameters = mPar
    self.__numLocals = mBor - self.__numCurvature
   
  def __getJacobian(self, aLabel):
    '''
    Get jacobian for transformation from fit to track parameters at point.
    
    @param aLabel: (signed) label of point
    @type aLabel: int
    @return: labels of required fit parameters and jacobian
    @rtype: list
    '''
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
    anIndex = []
    labDer, matDer = self.__getFitToLocalJacobian(aPoint, 5, nJacobian)
#   from local parameters
    for i in range(nLocals):
      aJacobian[i + 5, i] = 1.0;
      anIndex.append(i + 1);
  
#   from trajectory parameters
    iCol = nLocals;
    for i in range(5):
      if (labDer[i] > 0):
        anIndex.append(labDer[i]);
        for j in range(5):
          aJacobian[j, iCol] = matDer[j, i];
        iCol += 1

    return anIndex, aJacobian
   
  def __getFitToLocalJacobian(self, aPoint, measDim, nJacobian=1):
    '''
    Get (part of) jacobian for transformation from (trajectory) fit to track parameters at point.

    Jacobian broken lines (q/p,..,u_i,u_i+1..) to local (q/p,u',u) parameters.
    @param aPoint: point to use
    @type aPoint: GblPoint
    @param measDim: dimension of 'measurement' (2, 4 or 5)
    @type measDim: int
    @param nJacobian: direction (0: to previous offset, 1: to next offset)
    @type nJacobian: int
    @return: labels for fit parameters with non zero derivatives, corresponding transformation matrix
    @rtype: list(vector(int), matrix(float)) 
    '''
    aDim = self.__dimensions
    nDim = len(aDim)
    nCurv = self.__numCurvature
    nLocals = self.__numLocals      
    nOffset = aPoint.getOffset() 
    anIndex = [0, 0, 0, 0, 0]
    aJacobian = np.zeros((measDim, 5))
    labOffset = measDim - 2
    labSlope = measDim - 4
    labCurv = measDim - 5
    
    if (nOffset < 0): # need interpolation
      prevW, prevWJ, prevWd = aPoint.getDerivatives(0) # W-, W- * J-, W- * d-
      nextW, nextWJ, nextWd = aPoint.getDerivatives(1) # W+, W+ * J+, W+ * d+
      sumWJ = prevWJ + nextWJ
      matN = np.linalg.inv(sumWJ)      # N = (W- * J- + W+ * J+)^-1 
#     local offset
      if (labOffset >= 0):
#       derivatives for u_int      
        prevNW = np.dot(matN, prevW)     # N * W-
        nextNW = np.dot(matN, nextW)     # N * W+
        prevNd = np.dot(matN, prevWd)    # N * W- * d-
        nextNd = np.dot(matN, nextWd)    # N * W+ * d+
        iOff = nDim * (-nOffset - 1) + nLocals + nCurv + 1 # first offset ('i' in u_i)
        if (nCurv > 0):
          aJacobian[labOffset:measDim, 0:1] = -prevNd - nextNd # from curvature
          anIndex[0] = nLocals + 1
        aJacobian[labOffset:measDim, 1:3] = prevNW
        aJacobian[labOffset:measDim, 3:5] = nextNW
        for i in range(nDim):
          anIndex[1 + aDim[i]] = iOff + i
          anIndex[3 + aDim[i]] = iOff + nDim + i
#     local slope 
      if (labSlope >= 0):
#       derivatives for u'_int
        prevWPN = np.dot(nextWJ, prevNW) # W+ * J+ * N * W-
        nextWPN = np.dot(prevWJ, nextNW) # W- * J- * N * W+
        prevWNd = np.dot(nextWJ, prevNd) # W+ * J+ * N * W- * d-
        nextWNd = np.dot(prevWJ, nextNd) # W- * J- * N * W+ * d+
        if (nCurv > 0):
          aJacobian[labSlope:labOffset, 0:1] = prevWNd - nextWNd # from curvature
        aJacobian[labSlope:labOffset, 1:3] = -prevWPN
        aJacobian[labSlope:labOffset, 3:5] = nextWPN

    else: # at point               
#     anIndex must be sorted
#     forward : iOff2 = iOff1 + nDim, index1 = 1, index2 = 3
#     backward: iOff2 = iOff1 - nDim, index1 = 3, index2 = 1
      iOff1 = nDim * nOffset + nCurv + nLocals + 1 # first offset ('i' in u_i)
      index1 = 3 - 2 * nJacobian # index of first offset
      iOff2 = iOff1 + nDim * (nJacobian * 2 - 1) # second offset ('i' in u_i)
      index2 = 1 + 2 * nJacobian # index of second offset
#     local offset
      if (labOffset >= 0):
        aJacobian[labOffset, index1] = 1.0 # from 1st Offset
        aJacobian[labOffset + 1, index1 + 1] = 1.0
        for i in range(nDim):
          anIndex[index1 + aDim[i]] = iOff1 + i
#     local slope and curvature
      if (labSlope >= 0):
        matW, matWJ, vecWd = aPoint.getDerivatives(nJacobian) # W, W * J, W * d
        if (nCurv > 0):
          aJacobian[labSlope:labOffset, 0:1] = -vecWd # from curvature
          anIndex[0] = nLocals + 1
        aJacobian[labSlope:labOffset, index1:index1 + 2] = -matWJ
        aJacobian[labSlope:labOffset, index2:index2 + 2] = matW
        for i in range(nDim):
          anIndex[index2 + aDim[i]] = iOff2 + i  
          
#   local curvature
    if (labCurv >= 0):
      if (nCurv > 0):
        aJacobian[labCurv, labCurv] = 1.0  
                    
    return anIndex, aJacobian    
  
  def __getFitToKinkJacobian(self, aPoint):              
    ''' 
    Get jacobian for transformation from (trajectory) fit to kink parameters at point.

    Jacobian broken lines (q/p,..,u_i-1,u_i,u_i+1..) to kink (du') parameters.
    @param aPoint: point to use
    @type aPoint: GblPoint
    @return: labels for fit parameters with non zero derivatives, corresponding transformation matrix
    @rtype: list(vector(int), matrix(float)) 
    '''
    aDim = self.__dimensions
    nDim = len(aDim)
    nCurv = self.__numCurvature
    nLocals = self.__numLocals        
    nOffset = aPoint.getOffset() 
    anIndex = [0, 0, 0, 0, 0, 0, 0]
    aJacobian = np.zeros((2, 7))    

    prevW, prevWJ, prevWd = aPoint.getDerivatives(0) # W-, W- * J-, W- * d-
    nextW, nextWJ, nextWd = aPoint.getDerivatives(1) # W+, W+ * J+, W+ * d+
    sumWJ = prevWJ + nextWJ         # W- * J- + W+ * J+
    sumWd = prevWd + nextWd         # W+ * d+ + W- * d-
    iOff = (nOffset - 1) * nDim + nCurv + nLocals + 1 # first offset ('i' in u_i)

#   local offset
    if (nCurv > 0):
      aJacobian[:, 0:1] = -sumWd # from curvature
      anIndex[0] = nLocals + 1      
    aJacobian[:, 1:3] = prevW # from 1st Offset
    aJacobian[:, 3:5] = -sumWJ # from 2nd Offset
    aJacobian[:, 5:7] = nextW # from 3rd Offset
    for i in range(nDim):
      anIndex[1 + aDim[i]] = iOff + i
      anIndex[3 + aDim[i]] = iOff + nDim + i
      anIndex[5 + aDim[i]] = iOff + nDim * 2 + i
        
    return anIndex, aJacobian    

  def getResults(self, aLabel):
    '''
    Get results (corrections, covarinace matrix) at point in forward or backward direction.
    
    @param aLabel: signed label of point (<0 backward, >0 forward)
    @type aLabel: int
    @return: correction vector, covarinace matrix for track parameters
    @rtype: list
    '''
    anIndex, aJacobian = self.__getJacobian(aLabel)
    nParBRL = len(anIndex)
    aVec = np.empty(nParBRL)
    for i in range(nParBRL):
      aVec[i] = self.__vector[anIndex[i] - 1]         # compressed vector
    aMat = self.__matrix.getBlockMatrix(anIndex)  # compressed matrix
    locPar = np.dot(aJacobian, aVec) 
    locCov = np.dot(aJacobian, np.dot(aMat, aJacobian.T))
    return  locPar, locCov 
  
  def fit(self, optionList=""):
    '''
    Perform fit of trajectory.
    
    @param optionList: M-estimators to be used (one iteration per character)
    @type optionList: string
    @return: Chi2, Ndf, loss of weight from fit ([0., -1, 0.] if fit failed)
    @rtype: list
    '''
 
    def defineOffsets():
      '''
      Define offsets from list of points.
      '''
# set labels for previous/next offsets
#     first point is offset    
      self.__points[0].setOffset(0)
      nOffsets = 1
#     intermediate scatterers are offsets    
      for aPoint in self.__points[1:-1]:
        if (aPoint.hasScatterer()):
          aPoint.setOffset(nOffsets)
          nOffsets += 1
        else:  
          aPoint.setOffset(-nOffsets)
#     last point is offset    
      self.__points[-1].setOffset(nOffsets)
      self.__numOffsets = nOffsets + 1
      self.__numParameters = self.__numOffsets * len(self.__dimensions) \
                           + self.__numCurvature + self.__numLocals
 
    def calcJacobians():
      '''
      Calculate Jacobians to previous/next scatterer from point to point ones.
      '''
      scatJacobian = np.empty((5, 5)) 
#     forward propagation (all)
      lastPoint = 0;
      numStep = 0;
      for iPoint in range(1, self.__numPoints):
        if (numStep == 0):
          scatJacobian = self.__points[iPoint].getP2pJacobian()
        else: 
          scatJacobian = np.dot(self.__points[iPoint].getP2pJacobian(), scatJacobian)
        numStep += 1
        self.__points[iPoint].addPrevJacobian(scatJacobian) # aPoint -> previous scatterer
        if (self.__points[iPoint].getOffset() >= 0):
          self.__points[lastPoint].addNextJacobian(scatJacobian) # lastPoint -> next scatterer
          numStep = 0;
          lastPoint = iPoint
#     backward propagation (without scatterers)
      for iPoint in range(self.__numPoints - 1, 0, -1):
        if (self.__points[iPoint].getOffset() >= 0):
          scatJacobian = self.__points[iPoint].getP2pJacobian()
          continue # skip offsets
        self.__points[iPoint].addNextJacobian(scatJacobian) # iPoint -> next scatterer
        scatJacobian = np.dot(scatJacobian, self.__points[iPoint].getP2pJacobian())
      
    def prepare():
      '''
      Prepare fit; generate data from points.
      '''
      aDim = self.__dimensions
# measurements
      for aPoint in self.__points:
        if (aPoint.hasMeasurement()):
          nLabel = aPoint.getLabel()
          measDim = aPoint.getMeasDim()
          localDer = aPoint.getLocalDerivatives()
          globalLab = aPoint.getGlobalLabels()
          globalDer = aPoint.getGlobalDerivatives()
          matP, aMeas, aPrec = aPoint.getMeasurement()
          labDer, matDer = self.__getFitToLocalJacobian(aPoint, measDim)
          matPDer = np.dot(matP, matDer)
          for i in range(measDim):
            if (aPrec[i] > 0.):
              aData = GblData(nLabel, aMeas[i], aPrec[i])
              aData.addDerivatives(i, labDer, matPDer, localDer, \
                                   globalLab, globalDer)
              self.__data.append(aData)
#                aPoint.setDataMeas(i, len(self.__data)) 
# pseudo measurements from kinks
      for aPoint in self.__points[1:-1]:
        if (aPoint.hasScatterer()):
          nLabel = aPoint.getLabel()        
          aMeas, aPrec = aPoint.getScatterer()
          labDer, matDer = self.__getFitToKinkJacobian(aPoint)
          for i in aDim:
            if (aPrec[i] > 0.):
              aData = GblData(nLabel, aMeas[i], aPrec[i])
              aData.addDerivatives(i, labDer, matDer)
              self.__data.append(aData)
#              aPoint.setDataScat(i, len(self.__data)) 
#     external seed
      if (self.__externalPoint != 0):
        externalIndex, aJacobian = self.__getJacobian(self.__externalPoint)
        eigenVal, eigenVec = np.linalg.eigh(self.__externalSeed)
        aMatrix = np.dot(eigenVec.T, aJacobian)
        for i in range(len(eigenVec)):
          if (eigenVal[i] > 0.):
            externalDerivatives = []
            for j in range(len(externalIndex)):
              externalDerivatives.append(aMatrix[i, j])
            aData = GblData(self.__externalPoint, 0., eigenVal[i])
            aData.addExtDerivatives(externalIndex, externalDerivatives)
            self.__data.append(aData)

    def buildLinearEquationSystem():
      '''
      Build linear equation system from data.
      '''
      nBorder = self.__numCurvature + self.__numLocals
      self.__matrix = BorderedBandMatrix(self.__numParameters, nBorder)
      self.__vector = np.zeros(self.__numParameters)
      for aData in self.__data: 
        index, aVector, aMatrix = aData.getMatrices()
        for i in range(len(index)):
          self.__vector[ index[i] - 1 ] += aVector[0, i]        # update vector
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
    
    if (self.__data == []): # generate data from points   
      defineOffsets()    
      calcJacobians()                      
      prepare()
    buildLinearEquationSystem() # create linear equations system from data
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
             
      Ndf = len(self.__data) - self.__numParameters 
      Chi2 = 0.
      for aData in self.__data: 
        Chi2 += aData.getChi2()
      Chi2 /= [1.0, 0.8737, 0.9326, 0.8228 ][aMethod]  
      return Chi2, Ndf, lostWeight
    
    except (ZeroDivisionError, np.linalg.linalg.LinAlgError):
      return  0., -1, 0.
    

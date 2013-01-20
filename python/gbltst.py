'''
Simple Test Program for General Broken Lines.

Created on Jul 27, 2011

@author: kleinwrt
'''

import numpy as np
import math
import time
from gblfit import GblPoint, GblTrajectory
#
def example1():
  '''
  Create points on initial trajectory, create trajectory from points,
  fit and write trajectory to MP-II binary file,
  get track parameter corrections and covariance matrix at points.
  
  Equidistant measurement layers and thin scatterers, propagation 
  with simple jacobian (quadratic in arc length differences).
  Curvilinear system (U,V,T) as local coordinate system.
  '''  
  def gblSimpleJacobian(ds, cosl, bfac):
    '''
    Simple jacobian: quadratic in arc length difference.
    
    @param ds: arc length difference
    @type ds: float
    @param cosl: cos(lambda)
    @type cosl: float
    @param bfac: Bz*c
    @type bfac: float
    @return: jacobian to move by 'ds' on trajectory
    @rtype: matrix(float)
    '''
    jac = np.eye(5)
    jac[1, 0] = -bfac * ds * cosl
    jac[3, 0] = -0.5 * bfac * ds * ds * cosl
    jac[3, 1] = ds
    jac[4, 2] = ds  
    return jac
#
  nTry = 10#00 #: number of tries
  nLayer = 5   #: number of detector layers
  print " Gbltst $Rev$ ", nTry, nLayer
  start = time.clock()
# track direction
  sinLambda = 0.3
  cosLambda = math.sqrt(1.0 - sinLambda ** 2)
  sinPhi = 0.
  cosPhi = math.sqrt(1.0 - sinPhi ** 2)
#  tDir = np.array([cosLambda * cosPhi, cosLambda * sinPhi, sinLambda])
# U = Z x T / |Z x T|, V = T x U
  uvDir = np.array([[-sinPhi, cosPhi, 0.], \
                  [-sinLambda * cosPhi, -sinLambda * sinPhi, cosLambda]])
# measurement resolution
  measErr = np.array([ 0.001, 0.001]) # 10 mu
  measPrec = 1.0 / measErr ** 2
# scattering error
  scatErr = np.array([ 0.001, 0.001]) # 1 mread
  scatPrec = 1.0 / scatErr ** 2
# RMS of CurviLinear track parameters (Q/P, slopes, offsets)
  clErr = np.array([0.001, -0.1, 0.2, -0.15, 0.25])
  clSeed = None
  seedLabel = 0 # label of point with seed
  if seedLabel != 0:
    print " external seed at label ", seedLabel
#
  bfac = 0.2998 # Bz*c for Bz=1
  step = 1.5 / cosLambda # constant steps in RPhi
#
  Chi2Sum = 0.
  NdfSum = 0
  LostSum = 0.
#
  binaryFile = open("milleBinaryISN.dat", "wb")
#
  for iTry in range(nTry):
# generate (CurviLinear) track parameters 
    clNorm = np.random.normal(0., 1., 5)  
    clPar = clErr * clNorm
    # covariance matrix
    clCov = np.eye(5)
    for i in range(5):
      clCov[i, i] = clErr[i] ** 2
# arclength
    s = 0.
    sPoint = []
# point-to-point jacobian (from previous point)    
    jacPointToPoint = np.eye(5)
# additional (local or global) derivatives    
    addDer = np.array([[1.0], [0.0]])
    labGlobal = np.array([[4711], [4711]])
# create trajectory
    traj = GblTrajectory(bfac != 0.)
    
    for iLayer in range(nLayer):
#     measurement directions   
      sinStereo = (0. if iLayer % 2 == 0 else 0.5) 
      cosStereo = math.sqrt(1.0 - sinStereo ** 2)    
      mDir = np.array([[sinStereo, cosStereo, 0.0], [0., 0, 1.]])
# projection measurement to local (curvilinear uv) directions (duv/dm)
      proM2l = np.dot(uvDir, mDir.T)
# projection local (uv) to measurement directions (dm/duv)
      proL2m = np.linalg.inv(proM2l)
# measurement - prediction in measurement system with error
      measNorm = np.random.normal(0., 1., 2)  
      meas = np.dot(proL2m, clPar[3:5]) + measErr * measNorm
# point with measurement
      point = GblPoint(jacPointToPoint)
      point.addMeasurement([proL2m, meas, measPrec])
# additional local parameters?
#      point.addLocals(addDer)
# additional global parameters?
      point.addGlobals(labGlobal, addDer)
      addDer = -addDer # locDer flips sign every measurement      
# add point to trajectory      
      iLabel = traj.addPoint(point)
      sPoint.append(s)
      if iLabel == abs(seedLabel):
        clSeed = np.linalg.inv(clCov)
# propagate to scatterer
      jacPointToPoint = gblSimpleJacobian(step, cosLambda, bfac)
      clPar = np.dot(jacPointToPoint, clPar)
      clCov = np.dot(jacPointToPoint, np.dot(clCov, jacPointToPoint.T))
      s += step
      if (iLayer < nLayer - 1):
        scat = np.array([0., 0.])
# point with scatterer
        point = GblPoint(jacPointToPoint)
        point.addScatterer([scat, scatPrec])
        iLabel = traj.addPoint(point)
        sPoint.append(s)
        if iLabel == abs(seedLabel):
          clSeed = np.linalg.inv(clCov)

# scatter a little    
        scatNorm = np.random.normal(0., 1., 2)  
        clPar[1:3] = clPar[1:3] + scatErr * scatNorm
# propagate to next measurement layer    
        clPar = np.dot(jacPointToPoint, clPar)
        clCov = np.dot(jacPointToPoint, np.dot(clCov, jacPointToPoint.T))
        s += step
 
# add external seed
    if clSeed is not None:    
      traj.addExternalSeed(seedLabel, clSeed)
# dump trajectory
#    traj.dump()
  
# fit trajectory
    Chi2, Ndf, Lost = traj.fit()
    print " Record, Chi2, Ndf, Lost", iTry, Chi2, Ndf, Lost
# write to MP binary file    
    traj.milleOut(binaryFile)
# sum up    
    Chi2Sum += Chi2
    NdfSum += Ndf
    LostSum += Lost
# get corrections and covariance matrix at points 
    if (iTry == 0):
      for i in range(1, 1):      
        locPar, locCov = traj.getResults(i)
        print " Point< ", i
        print " locPar ", locPar
        print " locCov ", locCov      
        locPar, locCov = traj.getResults(-i)
        print " Point> ", i
        print " locPar ", locPar
        print " locCov ", locCov  
#
  end = time.clock()
  print " Time [s] ", end - start
  print " Chi2Sum/NdfSum ", Chi2Sum / NdfSum
  print " LostSum/nTry ", LostSum / nTry
#
def example2():
  '''
  Read trajectory from MP-II binary file and refit.
  '''
#  
  binaryFile = open("milleBinaryISN.dat", "rb")
  nRec = 0
  maxRec = 10 #: maximum number of records to read
  Chi2Sum = 0.
  NdfSum = 0
  LostSum = 0.
  start = time.clock()
  
  try:
    while(nRec < maxRec):
# create trajectory
      traj = GblTrajectory(0)
# read from file      
      traj.milleIn(binaryFile) # get data blocks from file
      nRec += 1
# fit trajectory      
      Chi2, Ndf, Lost = traj.fit()
      print " Record, Chi2, Ndf, Lost", nRec, Chi2, Ndf, Lost
# sum up      
      Chi2Sum += Chi2
      NdfSum += Ndf
      LostSum += Lost 
         
  except EOFError:
    pass    
  
  print " records read ", nRec
  end = time.clock()
  print " Time [s] ", end - start
  print " Chi2Sum/NdfSum ", Chi2Sum / NdfSum
  print " LostSum/nTry ", LostSum / nRec
      
# create points on initial trajectory, create trajectory from points,
# fit and write trajectory to MP-II binary file
# get track parameter corrections and covariance matrix at points
example1()
# read trajectory from MP-II binary file and refit
example2()

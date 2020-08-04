'''
Created on 28 Sep 2018

@author: kleinwrt
'''

## \file
# silicon tracker example
#
# \author Claus Kleinwort, DESY, 2018 (Claus.Kleinwort@desy.de)
#
#  \copyright
#  Copyright (c) 2018-2020 Deutsches Elektronen-Synchroton,
#  Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
#  This library is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Library General Public License as
#  published by the Free Software Foundation; either version 2 of the
#  License, or (at your option) any later version. \n\n
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details. \n\n
#  You should have received a copy of the GNU Library General Public
#  License along with this program (see the file COPYING.LIB for more
#  details); if not, write to the Free Software Foundation, Inc.,
#  675 Mass Ave, Cambridge, MA 02139, USA.

import numpy as np
import math
import random
import time
from gblfit import GblPoint, GblTrajectory


## Simulate and reconstruct helical tracks in silicon pixel and (1D or 2D) strip detectors. 
#
#  Create points on initial trajectory, create trajectory from points,
#  fit and write trajectory to MP-II binary file (for rigid body alignment).
#
#  Setup:
#   - Beam (mainly) in X direction
#   - Constant magnetic field in Z direction
#   - Silicon sensors measuring in YZ plane, orthogonal (pixel) or non-orthogonal (stereo strips, double sided or composite) measurement systems
#   - Multiple scattering in sensors (air inbetween ignored)
#   - Curvilinear system (T,U,V) as local coordinate system and (q/p, slopes, offsets) as local track parameters
#
#  Local systems.
#  Up to three (different) local coordinate systems can be defined at each point:
#    - Track model linearization (propagation, fitting), e.g. curvilinear system
#    - Measurement, defined by two (optionally non-orthogonal) measurement directions,
#      normal to detector plane and detector position (offset)
#    - Alignment, defined by two orthogonal directions in detector plane, normal to that
#      and detector position (offset)
# 
# \remark To exercise (mis)alignment different sets of layers (with different geometry) for simulation and reconstruction can be used.
#
# Example steering file for Millepede-II (B=0):
# \code{.unparsed}
# Cfiles
# milleBinary.dat
# 
# method inversion 3 0.1
# chiscut 30. 6.
# printcounts
# ! fix first pixel and last stereo layer as reference
# parameter
#   1  0.  -1.
#   2  0.  -1.
#   3  0.  -1.
#   4  0.  -1.
#   5  0.  -1.
#   6  0.  -1.
#  61  0.  -1.
#  62  0.  -1.
#  63  0.  -1.
#  64  0.  -1.
#  65  0.  -1.
#  66  0.  -1.
# \endcode
#
def exampleSit():

  #
  random.seed(47117)
  
  # magnetic field
  bfac = 0.003  # B*c for 1 T
  #bfac = 0.  # B*c for 0 T
  # detector layers: name, position (x,y,z), thickness (X/X_0), (1 or 2) measurements (direction in YZ, resolution), spacing for composite layers  
  det = gblSiliconDet([ ['PIX1', (2.0, 0., 0.), 0.0033, [(0., 0.0010), (90., 0.0020)] ],  # pixel
                        ['PIX2', (3.0, 0., 0.), 0.0033, [(0., 0.0010), (90., 0.0020)] ],  # pixel
                        ['PIX3', (4.0, 0., 0.), 0.0033, [(0., 0.0010), (90., 0.0020)] ],  # pixel
                        ['S2D4', (6.0, 0., 0.), 0.0033, [(0., 0.0025), (5., 0.0025)], 0.1 ],  # strip 2D (composite), 5 deg stereo angle
                        ['S2D5', (8.0, 0., 0.), 0.0033, [(0., 0.0025), (-5., 0.0025)], 0.1 ],  # strip 2D (composite), -5 deg stereo angle
                        ['S2D6', (10.0, 0., 0.), 0.0033, [(0., 0.0025), (5., 0.0025)] ],  # strip 2D (double sided), 5 deg stereo angle
                        ['S2D7', (12.0, 0., 0.), 0.0033, [(0., 0.0025), (-5., 0.0025)] ],  # strip 2D (double sided), -5 deg stereo angle
                        ['S1D8', (15.0, 0., 0.), 0.0033, [(0., 0.0040)] ]  # strip 1D
                      ], bfac)
               
  nTry = 1000  #: number of tries
  qbyp = 0.2  # 5 GeV
  binaryFile = open("milleBinary.dat", "wb")
  #binaryFile = None
  #
  print " Gblsit $Rev: 151 $ ", nTry
  #
  start = time.clock()
  Chi2Sum = 0.
  NdfSum = 0
  LostSum = 0.
  #
  for iTry in range(nTry):
    # helix parameter for track generation
    dca = random.gauss(0., 0.1)
    z0 = random.gauss(0., 0.1)
    phi0 = 0.5 * random.uniform(-1., 1.)
    dzds = 0.3 * random.uniform(-1., 1.)
    curv = bfac * qbyp * math.sqrt(1. + dzds * dzds)
    genPar = [curv, phi0, dca, dzds, z0]
    #print " genPar ", iTry, genPar
    # generate hits
    genHits = det.generateHits(qbyp, genPar)
    # seed (with true parameters)
    seedPar = genPar
    seed = gblSimpleHelix(seedPar)
    sOld = 0.
    cosLambda = 1. / math.sqrt(1. + seedPar[3] ** 2)
    # construct GBL trajectory
    traj = GblTrajectory(bfac != 0.)
    # add GBL points
    for l, layer in enumerate(det.getLayers()):
      # prediction from seeding helix
      pred = layer.intersectWithHelix(seed)
      measPred = pred.getMeasPred() 
      sArc = pred.getArcLength()      
      # residuals
      res = np.array([genHits[l][0] - measPred[0], genHits[l][1] - measPred[1]])
      # measurement precision
      measPrec = np.array(layer.getPrecision())
      # Curvilinear system: track direction T, U = Z x T / |Z x T|, V = T x U
      # as local system
      curviDirs = pred.getCurvilinearDirs()
      # projection matrix (local to measurement)
      proL2m = np.linalg.inv(np.dot(curviDirs, np.linalg.inv(layer.getMeasSystemDirs())[:, :2]))
      # propagation
      jacPointToPoint = gblSimpleJacobian((sArc - sOld) / cosLambda, cosLambda, bfac)
      sOld = sArc
      # point with (independent) measurements (in measurement system)
      point = GblPoint(jacPointToPoint)
      # composite?
      if layer.isComposite():
        # 2nd prediction
        pred = layer.intersectWithHelix2(seed)
        measPred = pred.getMeasPred()
        # 4D measurement
        pro4D = np.zeros((4, 4)); pro4D[2:, 2:] = proL2m; pro4D[3, :2] = proL2m[1, :] * layer.getSpacing();
        res4D = np.array([0., 0., res[0], genHits[l][1] - measPred[1]])
        prec4D = np.array([0., 0., measPrec[0], measPrec[1]])
        point.addMeasurement([pro4D, res4D, prec4D]) 
      else:
        # 2D measurement  
        point.addMeasurement([proL2m, res, measPrec])
      # global parameters for rigid body alignment?
      if binaryFile is not None:
        # local (alignment) system, per layer
        labels = [l * 10 + 1, l * 10 + 2, l * 10 + 3, l * 10 + 4, l * 10 + 5, l * 10 + 6]
        labGlobal = np.array([labels, labels])
        derGlobal = layer.getRigidBodyDerLocal(pred.getPos(), pred.getDirection())
        # composite?
        if layer.isComposite():
          # 4D
          labG4D = np.array([labels, labels, labels, labels])
          derG4D = np.zeros((4, 6)); derG4D[2:] = derGlobal
          point.addGlobals(labG4D, derG4D)
        else:
          # 2D  
          point.addGlobals(labGlobal, derGlobal)
      # add scatterer to point
      radlen = layer.getRadiationLength() / abs(pred.getCosIncidence())
      scatErr = gblMultipleScatteringError(qbyp, radlen)  # simple model
      if scatErr > 0.:
        scat = np.array([0., 0.])
        scatP = np.array([1. / scatErr ** 2, 1. / scatErr ** 2])
        # composite?
        if layer.isComposite():
          # two similar sub layers
          scatP *= 0.5
        point.addScatterer([scat, scatP])
      # add point to trajectory      
      traj.addPoint(point)

    # fit trajectory
    Chi2, Ndf, Lost = traj.fit()
    print " Record, Chi2, Ndf, Lost", iTry, Chi2, Ndf, Lost
    # sum up    
    Chi2Sum += Chi2
    NdfSum += Ndf
    LostSum += Lost
    # write to binary file
    if binaryFile is not None:
      traj.milleOut(binaryFile)
    
  end = time.clock()
  print " Time [s] ", end - start
  print " Chi2Sum/NdfSum ", Chi2Sum / NdfSum
  print " LostSum/nTry ", LostSum / nTry


## Simple jacobian.
#  
#  Simple jacobian for (q/p, slopes, offsets) in curvilinear system,
#  constant magnetic field in Z direction, quadratic in arc length difference.
#
#  @param ds arc length difference; float
#  @param cosl cos(lambda); float
#  @param bfac Bz*c; float
#  @return jacobian to move by 'ds' on trajectory matrix(float)
#
def gblSimpleJacobian(ds, cosl, bfac):
  jac = np.eye(5)
  jac[1, 0] = -bfac * ds * cosl
  jac[3, 0] = -0.5 * bfac * ds * ds * cosl
  jac[3, 1] = ds
  jac[4, 2] = ds  
  return jac     


## Multiple scattering error
#
# Simple model (Rossi, Greisen)
#
# @param[in] qbyp    q/p [1/GeV]; float
# @param[in] xbyx0   thickness / radiation length; float
#
def gblMultipleScatteringError(qbyp, xbyx0):
  return 0.015 * abs(qbyp) * math.sqrt(xbyx0)
  

## Silicon layer
class gblSiliconLayer(object):

  ## Constructor
  #
  # @param[in]  layer   layer description; list
  #
  def __init__(self, layer):
    ## center
    self.__center = np.array(layer[1])
    ## radiation length
    self.__xbyx0 = layer[2] 
    # measurements (1D or 2D)
    meas = layer[3]
    ## resolution (for simulation)
    self.__resolution = (meas[0][1], meas[1][1] if len(meas) > 1 else 0.)
    ## precision (for reconstruction)
    self.__precision = (1. / meas[0][1] ** 2, 1. / meas[1][1] ** 2 if len(meas) > 1 else 0.)
    # measurement angles
    uPhi = meas[0][0] / 180. * math.pi
    vPhi = meas[1][0] / 180. * math.pi if len(meas) > 1 else uPhi + 0.5 * math.pi
    ## measurement direction u
    self.__uDir = np.array([0., math.cos(uPhi), math.sin(uPhi)])  
    ## measurement direction v
    self.__vDir = np.array([0., math.cos(vPhi), math.sin(vPhi)])
    ## normal to measurement plane
    self.__nDir = np.array([1., 0., 0.])
    ## measurement directions
    self.__measDirs = np.array([self.__uDir, self.__vDir, self.__nDir]) 
    ## local alignment system (IJK = YZX)
    self.__ijkDirs = np.array([[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]])
    ## spacing (for composite layers)
    self.__spacing = layer[4] if len(layer) > 4 else None
    
  ## Get radiation length
  def getRadiationLength(self):
    return self.__xbyx0
  
  ## Get resolution
  def getResolution(self):
    return self.__resolution

  ## Get precision
  def getPrecision(self):
    return self.__precision
  
  ## Get directions of measurement system
  def getMeasSystemDirs(self):
    return self.__measDirs
    
  ## Intersect with helix
  #
  # @param[in]  helix  helix
  # @return prediction
  #  
  def intersectWithHelix(self, helix):
    return helix.getPrediction(self.__center, self.__uDir, self.__vDir)  

  ## Intersect with helix (2nd sub layer)
  #
  # @param[in]  helix  helix
  # @return prediction
  #  
  def intersectWithHelix2(self, helix):
    return helix.getPrediction(self.__center + self.__spacing * self.__nDir, self.__uDir, self.__vDir)
  
  ## Is composite?
  def isComposite(self):
    return self.__spacing is not None  
  
  ## Get spacing
  def getSpacing(self):
    return self.__spacing
  
  ## Get rigid body derivatives in global frame
  #
  # @param[in] position   position (of prediction or measurement); vector
  # @param[in] trackDir   track direction; vector
  # @return global derivatives; matrix
  #
  def getRigidBodyDerGlobal(self, position, trackDir):
    # lever arms (for rotations)
    dist = position
    # dr/dm (residual vs measurement, 1-tdir*ndir^t/tdir*ndir)
    drdm = np.eye(3) - np.outer(trackDir, self.__nDir) / np.dot(trackDir, self.__nDir)
    # dm/dg (measurement vs 6 rigid body parameters)
    dmdg = np.zeros((3, 6))
    dmdg[0][0] = 1.; dmdg[0][4] = -dist[2]; dmdg[0][5] = dist[1]
    dmdg[1][1] = 1.; dmdg[1][3] = dist[2]; dmdg[1][5] = -dist[0]
    dmdg[2][2] = 1.; dmdg[2][3] = -dist[1]; dmdg[2][4] = dist[0]
    # drl/drg (local vs global residuals)
    drldrg = self.__measDirs    
    # drl/dg (local residuals vs rigid body parameters)
    drldg = np.dot(drldrg, np.dot(drdm, dmdg))
    return drldg

  ## Get rigid body derivatives in local (alignment) frame
  #
  # @param[in] position   position (of prediction or measurement); vector
  # @param[in] trackDir   track direction; vector
  # @return global derivatives
  #
  def getRigidBodyDerLocal(self, position, trackDir): 
    # track direction in local system 
    tLoc = np.dot(self.__ijkDirs, trackDir)
    # local slopes
    uSlope = tLoc[0] / tLoc[2]
    vSlope = tLoc[1] / tLoc[2]
    # (u,v) lever arms
    uPos, vPos = np.dot(self.__ijkDirs, position - self.__center)[:2]
    # wPos = 0 (in detector plane)
    # drl/dg (local residuals vs rigid body parameters)
    drldg = np.array([[1.0, 0.0, -uSlope, vPos * uSlope, -uPos * uSlope, vPos], \
                      [0.0, 1.0, -vSlope, vPos * vSlope, -uPos * vSlope, -uPos]])
    return drldg  

      
## Silicon detector
class gblSiliconDet(object):

  ## Constructor
  #
  # @param[in]  layers  layers; list
  # @param[in]  bfac    magnetic field (B*c); float
  #
  def __init__(self, layers, bfac):
    ## layers
    self.__layers = []  
    ## B*c
    self.__bfac = bfac
    
    for layer in layers:
      self.__layers.append(gblSiliconLayer(layer))
      
  ## Get layers
  def getLayers(self):
    return self.__layers    

  ## Generate hits on helix
  #
  # @param[in]  qbyp     q/p
  # @param[in]  genPar   helix parameters
  # @return list of hits
  #  
  def generateHits(self, qbyp, genPar): 

# list of hits
    hits = []
    localPar = genPar
    # print " track ", helix
    for layer in self.__layers:
      # local constant (Bfield) helix
      hlx = gblSimpleHelix(localPar)
      # prediction from local helix
      pred = layer.intersectWithHelix(hlx)
      meas = [pred.getMeasPred()] 
      # scatter at intersection point
      xPos, yPos = pred.getPos()[:2]
      radlen = layer.getRadiationLength() / abs(pred.getCosIncidence())
      errMs = gblMultipleScatteringError(qbyp, radlen)  # simple model
      cosLambda = 1. / math.sqrt(1. + localPar[3] ** 2)
      # move to intersection point
      newpar = hlx.moveTo((xPos, yPos))
      newpar[1] += random.gauss(0., errMs / cosLambda)  # phi
      newpar[3] += random.gauss(0., errMs / cosLambda ** 2)  # dzds
      newhlx = gblSimpleHelix(newpar)
      # move back
      localPar = newhlx.moveTo((-xPos, -yPos))
      # composite layer
      if layer.isComposite():
        # 2nd prediction from local helix
        pred = layer.intersectWithHelix2(hlx)
        meas.append(pred.getMeasPred())
        # scatter at intersection point
        xPos, yPos = pred.getPos()[:2] 
        cosLambda = 1. / math.sqrt(1. + localPar[3] ** 2)
        # move to intersection point
        newpar = hlx.moveTo((xPos, yPos))
        newpar[1] += random.gauss(0., errMs / cosLambda)  # phi
        newpar[3] += random.gauss(0., errMs / cosLambda ** 2)  # dzds
        newhlx = gblSimpleHelix(newpar)
        # move back
        localPar = newhlx.moveTo((-xPos, -yPos))
        
      # add (smeared) hit
      sigma = layer.getResolution()
      hits.append((random.gauss(meas[0][0], sigma[0]), random.gauss(meas[-1][1], sigma[1])))
    
    return hits


## Simple helix
#
# Assuming constant magnetic field in (positive) Z-direction
#
class gblSimpleHelix(object):  
  
  ## Constructor.
  #
  # @param[in] parameter helix parameter (curv, phi0, dca, dzds, z0); list
  #
  # For comparison: Generalized circle equation: 
  #  n_0 + x*n_1 + y*n_2 + (x*x+y*y)*n_3 = 0, 
  #  n_0 ~= -dca, (n_1, n_2) = -(cos(phi_0), sin(phi_0)), n_3 = 0.5*rinv      
  #
  def __init__(self, parameter): 
    ## curvature (in XY)
    self.__rinv = parameter[0]
    ## flight direction at point of closest approach (in XY)
    self.__phi0 = parameter[1]
    ## direction vector at point of closest approach (in XY)
    self.__dir0 = (math.cos(self.__phi0), math.sin(self.__phi0))
    ## distance of closest approach in (XY)
    self.__dca = parameter[2]
    ## dZ/ds
    self.__dzds = parameter[3]
    ## Z position at distance of closest approach
    self.__z0 = parameter[4]
    ## XY circle parameter: X position of center / R
    self.__xRelCenter = -(1. - self.__dca * self.__rinv) * self.__dir0[1]
    ## XY circle parameter: Y position of center / R
    self.__yRelCenter = (1. - self.__dca * self.__rinv) * self.__dir0[0]
    
  ## Get prediction
  #
  # Get prediction from intersection of helix with measurement plane.
  #
  # @param[in] refPos  reference position on detector plane; vector
  # @param[in] uDir    measurement direction 'u'; vector
  # @param[in] vDir    measurement direction 'v'; vector
  # @return prediction; class
  #
  def getPrediction(self, refPos, uDir, vDir):  
    # normal to (u,v) measurement plane
    nDir = np.cross(uDir, vDir); nDir /= np.linalg.norm(nDir)
    # ZS direction  
    cosLambda = 1. / math.sqrt(1. + self.__dzds * self.__dzds)
    sinLambda = self.__dzds * cosLambda
    # line (or helix)
    if self.__rinv == 0.:
      # track direction
      tDir = np.array([cosLambda * self.__dir0[0], cosLambda * self.__dir0[1], sinLambda])
      # distance (of point at dca to reference)
      pca = np.array([ self.__dca * self.__dir0[1] , -self.__dca * self.__dir0[0], self.__z0])
      dist = pca - refPos
      # arc-length
      sArc3D = -np.dot(dist, nDir) / np.dot(tDir, nDir); sArc2D = sArc3D * cosLambda
      # distance (of point at sArc to reference)
      pos = pca + sArc3D * tDir
      dist = pos - refPos
    else:
      # initial guess of 2D arc-length
      sArc2D = self.getArcLengthXY(refPos[0], refPos[1])
      nIter = 0
      while nIter < 10:
        nIter += 1
        # track direction
        dPhi = sArc2D * self.__rinv
        cosPhi = math.cos(self.__phi0 + dPhi); sinPhi = math.sin(self.__phi0 + dPhi)
        tDir = np.array([cosLambda * cosPhi, cosLambda * sinPhi, sinLambda])
        # distance (of point at sArc to reference)
        pos = np.array([(self.__xRelCenter + sinPhi) / self.__rinv, (self.__yRelCenter - cosPhi) / self.__rinv, self.__z0 + self.__dzds * sArc2D])
        dist = pos - refPos
        # arc-length correction (linearizing helix at sNew) 
        sCorr3D = -np.dot(dist, nDir) / np.dot(tDir, nDir)
        if abs(sCorr3D) > 0.0001:
          sArc2D += sCorr3D * cosLambda
        else:
          break
          
    # prediction in measurement directions
    pred = [np.dot(dist, uDir), np.dot(dist, vDir)]
    return gblHelixPrediction(sArc2D, pred, tDir, uDir, vDir, nDir, pos)  

  ## Get (2D) arc length for given point.
  #
  # Arc length from dca to point on circle on intersection with line
  # from circle center to given point
  #
  # @param[in] xPos   X Position; float
  # @param[in] yPos   Y Position; float
  # @return arc length from dca to point on circle; float
  #  
  def getArcLengthXY(self, xPos, yPos):
    # line
    if self.__rinv == 0:
      return self.__dir0[0] * xPos + self.__dir0[1] * yPos 
    # helix  
    dx = (xPos * self.__rinv - self.__xRelCenter)
    dy = (yPos * self.__rinv - self.__yRelCenter)
    dphi = math.atan2(dx, -dy) - self.__phi0
    if (abs(dphi) > math.pi):
      dphi -= cmp(dphi, 0.) * 2.0 * math.pi 
    return dphi / self.__rinv
  
  ## Change reference point
  #
  # @param[in] newRefPoint   new reference point (in XY); vector     
  # @return    new helix parameters; list      
  #
  # Based on V. Karimaki, NIM A305 (1991) 187-191, eqn (19)
  def moveTo(self, newRefPoint):
   
    rho = self.__rinv
    phi = self.__phi0
    dca = self.__dca
    dzds = self.__dzds
    z0 = self.__z0
    
    u = 1. - rho * dca
    dp = -newRefPoint[0] * self.__dir0[1] + newRefPoint[1] * self.__dir0[0] + dca
    dl = newRefPoint[0] * self.__dir0[0] + newRefPoint[1] * self.__dir0[1]
    sa = 2. * dp - rho * (dp * dp + dl * dl)
    sb = rho * newRefPoint[0] + u * self.__dir0[1]
    sc = -rho * newRefPoint[1] + u * self.__dir0[0]
    sd = math.sqrt(1. - rho * sa)
    # transformed parameters
    if rho == 0.:      
      dca = dp
      sArc = dl
      newPar = [rho, phi, dca]
    else:
      phi = math.atan2(sb, sc)
      dca = sa / (1. + sd)
      dphi = phi - self.__phi0 
      if abs(dphi) > math.pi: dphi -= cmp(dphi, 0.) * 2.0 * math.pi 
      sArc = dphi / rho
      newPar = [rho, phi, dca]  
    z0 += sArc * dzds 
    newPar += [dzds, z0]
    
    return newPar


## Prediction (from helix at measurement)
#
class gblHelixPrediction(object): 

  ## Constructor
  #
  # @param[in] sArc     arc length; float
  # @param[in] pred     prediction for measurement (u,v); list
  # @param[in] tDir     track direction at prediction; vector
  # @param[in] uDir     measurement direction for u; vector
  # @param[in] vDir     measurement direction for v; vector
  # @param[in] nDir     normal to measurement plane; vector
  # @param[in] pos      position at prediction; vector
  #
  def __init__(self, sArc, pred, tDir, uDir, vDir, nDir, pos):  
    ## arc-length
    self.__sarc = sArc
    ## prediction
    self.__pred = pred  
    ## track direction
    self.__tdir = tDir
    ## u direction
    self.__udir = uDir
    ## v direction
    self.__vdir = vDir
    ## normal to (u,v)
    self.__ndir = nDir
    ## position
    self.__pos = pos
    #
        
  ## Get arc-length   
  def getArcLength(self):
    return self.__sarc
  
  ## Get measurement prediction
  def getMeasPred(self):
    return self.__pred
  
  ## Get Position
  def getPos(self):
    return self.__pos
 
  ## Get (track) direction
  def getDirection(self):
    return self.__tdir
  
  ## Get cosine of incidence
  def getCosIncidence(self):
    return np.dot(self.__tdir, self.__ndir)
  
  ## Get curvilinear directions (U,V)
  #
  # Curvilinear system: track direction T, U = Z x T / |Z x T|, V = T x U
  #
  def getCurvilinearDirs(self):
    cosTheta = self.__tdir[2]; sinTheta = math.sqrt(self.__tdir[0] ** 2 + self.__tdir[1] ** 2)
    cosPhi = self.__tdir[0] / sinTheta; sinPhi = self.__tdir[1] / sinTheta
    return np.array([[-sinPhi, cosPhi, 0.], [-cosPhi * cosTheta, -sinPhi * cosTheta, sinTheta]])
  
      
if __name__ == '__main__':
  exampleSit()

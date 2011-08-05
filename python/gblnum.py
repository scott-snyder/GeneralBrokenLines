'''
Algebra for linear equation system with bordered band matrix.  

Created on Jul 27, 2011

@author: kleinwrt
'''

import numpy as np

class BorderedBandMatrix(object):
  '''
  (Symmetric) Bordered Band Matrix. 
    
  Separate storage of border, mixed and band parts. 
  Example for matrix size=8 with border size and band width of two::
  
     +-                                 -+
     |  B11 B12 M13 M14 M15 M16 M17 M18  |
     |  B12 B22 M23 M24 M25 M26 M27 M28  |
     |  M13 M23 C33 C34 C35  0.  0.  0.  |
     |  M14 M24 C34 C44 C45 C46  0.  0.  |
     |  M15 M25 C35 C45 C55 C56 C57  0.  |
     |  M16 M26  0. C46 C56 C66 C67 C68  |
     |  M17 M27  0.  0. C57 C67 C77 C78  |
     |  M18 M28  0.  0.  0. C68 C78 C88  |
     +-                                 -+
     
  Is stored as::
  
     +-         -+     +-                         -+
     |  B11 B12  |     |  M13 M14 M15 M16 M17 M18  | 
     |  B12 B22  |     |  M23 M24 M25 M26 M27 M28  |
     +-         -+     +-                         -+
     
                       +-                         -+
                       |  C33 C44 C55 C66 C77 C88  |
                       |  C34 C45 C56 C67 C78  0.  |
                       |  C35 C46 C57 C68  0.  0.  |
                       +-                         -+
  '''
  def __init__(self, nSize, nBorder=1, nBand=5):
    '''
    Create new BBmatrix.
    
    @param nSize: size of matrix
    @type nSize: int
    @param nBorder: size of border (default: 1, 'curvature')
    @type nBorder: int
    @param nBand: (maximal) band width (5) 
    @type nBand: int
    '''
    nSizeBand = nSize - nBorder
    self.__numSize = nSize
    '''@ivar: size of matrix
       @type: int'''
    self.__numBorder = nBorder
    '''@ivar: size of border
       @type: int'''
    self.__numBand = 0        # actual band width
    '''@ivar: (actual) band width
       @type: int'''
    self.__numCol = nSizeBand
    '''@ivar: size of band part of matrix
       @type: int'''
    self.__border = np.zeros((nBorder, nBorder))
    '''@ivar: border part B
       @type: matrix(float)'''
    self.__mixed = np.zeros((nBorder, nSizeBand))
    '''@ivar: mixed part M
       @type: matrix(float)'''
    self.__band = np.zeros((nBand + 1, nSizeBand))
    '''@ivar: band part C
       @type: matrix(float)'''
#    print " new BBM ", self.__border__.shape, self.__mixed__.shape, self.__band__.shape
    
# add symmetric block matrix (with 'expansion' according to index)  
  def addBlockMatrix(self, aIndex, aMatrix):
    '''
    Add (compressed) block to BBmatrix::
    
      BBmatrix(aIndex(i),aIndex(j)) += aMatrix(i,j)
    
    @param aIndex: list of indices
    @type aIndex: list(int)
    @param aMatrix: (compressed) matrix
    @type aMatrix: matrix(float)
    '''
    nBorder = self.__numBorder
    for i in range(len(aIndex)):
      iIndex = aIndex[i] - 1
      for j in range(i + 1):
        jIndex = aIndex[j] - 1
        if (iIndex < nBorder):
          self.__border[iIndex, jIndex] += aMatrix[i, j]
          if (iIndex != jIndex):
            self.__border[jIndex, iIndex] += aMatrix[i, j]            
        elif (jIndex < nBorder):
          self.__mixed[jIndex, iIndex - nBorder] += aMatrix[i, j]        
        else:
          nBand = iIndex - jIndex
          self.__band[nBand, jIndex - nBorder] += aMatrix[i, j] 

          self.__numBand = max(self.__numBand, nBand)
    return self.__numBand 
    
  def getBlockMatrix(self, aIndex):
    '''
    Retrieve (compressed) block from BBmatrix::
    
      aMatrix(i,j) = BBmatrix(aIndex(i),aIndex(j))
    
    @param aIndex: list of indices
    @type aIndex: list(int)
    @return: (compressed) matrix
    @rtype: matrix(float)
    '''
    nBorder = self.__numBorder
    nSize = len(aIndex)
    aMatrix = np.empty((nSize, nSize))  
    for i in range(nSize):
      iIndex = aIndex[i] - 1 
      for j in range(i + 1):
        jIndex = aIndex[j] - 1
        iMax = max(iIndex, jIndex)
        iMin = min(iIndex, jIndex)
        if (iMax < nBorder):
          aMatrix[i, j] = self.__border[iIndex, jIndex]
        elif (iMin < nBorder):
          aMatrix[i, j] = self.__mixed[iMin, iMax - nBorder]       
        else:
          nBand = iIndex - jIndex
          aMatrix[i, j] = self.__band[nBand, jIndex - nBorder]
        aMatrix[j, i] = aMatrix[i, j]
    return aMatrix    
          
  def printMatrix(self):
    '''
    Print BBmatrix.
    '''
    print " block part "
    nRow = self.__numBorder
    for i in range(nRow):  
      print " row ", i, self.__border[i]
    print " mixed part "
    for i in range(nRow):  
      print " row ", i, self.__mixed[i]
    nRow = self.__numBand + 1
    print " band part "
    for i in range(nRow):
      print " diagonal ", i, self.__band[i]
      
# solve BorderedBandMatrix * aSolution = aRightHandSide,
# calculate bordered band part of inverse of BorderedBandMatrix
  def solveAndInvertBorderedBand(self, aRightHandSide):
    '''
    Solve linear equation A*x=b system with BBmatrix A, calculate BB part of inverse of A.
    
    @param aRightHandSide: right hand side 'b' of linear equation system
    @type aRightHandSide: vector(float)
    @return: solution
    @rtype: vector(float)
    @raise ZeroDivisionError: Band matrix is not positive definite
    @note: BBmatrix is replaced by BB part of it's inverse
    '''
#============================================================================
## from Dbandmatrix.F (MillePede-II by V. Blobel, Univ. Hamburg) 
#============================================================================
    def decomposeBand():
      '''
      (root free) Cholesky decomposition of band part: C=LDL^T
      
      @note: band part (C) is replaced by its decomposition (D,L)
      '''
      nRow = self.__numBand + 1
      nCol = self.__numCol
      auxVec = np.copy(self.__band[0]) * 16.0 # save diagonal elements
      for i in range(nCol):
        if ((self.__band[0, i] + auxVec[i]) != self.__band[0, i]):
          self.__band[0, i] = 1.0 / self.__band[0, i]
        else:
          self.__band[0, i] = 0.0  # singular
        for j in range(min(nRow, nCol - i) - 1):
          rxw = self.__band[j + 1, i] * self.__band[0, i]
          for k in range(min(nRow, nCol - i) - j - 1):
            self.__band[k, i + j + 1] -= self.__band[k + j + 1, i] * rxw
          self.__band[j + 1, i] = rxw  

# solve band * aSolution = aRightHandSide
    def solveBand(aRightHandSide):
      '''
      Solve linear equation system for band part.
     
      @param aRightHandSide: right hand side
      @type aRightHandSide:vector(float)
      @return: solution
      @rtype: vector(float)
      '''
      nRow = self.__numBand + 1
      nCol = self.__numCol  
      aSolution = np.copy(aRightHandSide)
      for i in range(nCol):   # forward substitution
        for j in range(min(nRow, nCol - i) - 1):
          aSolution[j + i + 1] -= self.__band[j + 1, i] * aSolution[i]    
      for i in range(nCol - 1, -1, -1):   # backward substitution
        rxw = self.__band[0, i] * aSolution[i]    
        for j in range(min(nRow, nCol - i) - 1):
          rxw -= self.__band[j + 1, i] * aSolution[j + i + 1]    
        aSolution[i] = rxw       
      return aSolution

# invert band part
    def invertBand():
      '''
      Invert band part.
      
      @return: band part
      @rtype: matrix(float)
      '''
      nRow = self.__numBand + 1
      nCol = self.__numCol       
      inverseBand = np.zeros((nRow, nCol))
      for i in range(nCol - 1, -1, -1):   
        rxw = self.__band[0, i]
        for j in range(i, max(0, i - nRow + 1) - 1, -1):
          for k in range(j + 1, min(nCol, j + nRow)):
            rxw -= inverseBand[abs(i - k), min(i, k)] * self.__band[k - j, j]
          inverseBand[i - j, j] = rxw
          rxw = 0.0
      return inverseBand
    
# band part of: anArray * aSymArray * anArray.T
    def bandOfAVAT(anArray, aSymArray):
      '''
      Calculate band part of A*V*A^T.

      @param anArray: matrix A
      @type anArray: matrix(float)
      @param aSymArray: symmetric matrix V
      @type aSymArray: matrix(float)
      @return: band part
      @rtype: matrix(float)      
      '''
      nCol = self.__numCol
      nBand = self.__numBand
      aBand = np.empty((nBand + 1, nCol))
      for i in range(nCol):
        for j in range(max(0, i - nBand), i + 1):
          aBand[i - j, j] = np.dot(anArray[i], np.dot(aSymArray, anArray[j]))
      return aBand

# setup
    nBorder = self.__numBorder
#    nRow = self.__numBand +1
    nCol = self.__numCol
    aSolution = np.empty(nBorder + nCol)
# decompose band
    decomposeBand()
    if ((self.__band[0] <= 0.0).any()):
      raise ZeroDivisionError("Band matrix not positive definite")

    if (nBorder > 0):
      auxMat = np.empty((nBorder, nCol))
# solve for mixed part      
      for i in range(nBorder):    
        auxMat[i] = solveBand(self.__mixed[i])
      auxMatT = auxMat.T      	
# solve for border
      auxVec = aRightHandSide[:nBorder] - np.dot(auxMat, aRightHandSide[nBorder:])
      invBorder = np.linalg.inv(self.__border - np.dot(self.__mixed, auxMatT))
      aSolution[:nBorder] = np.dot(invBorder, auxVec)
# solve for band part 
      aSolution[nBorder:] = solveBand(aRightHandSide[nBorder:]) \
                          - np.dot(auxMatT, aSolution[:nBorder])
# parts of inverse
      self.__border = invBorder
      self.__mixed = np.dot(-invBorder, auxMat)
      self.__band = invertBand() + bandOfAVAT(auxMatT, invBorder)
    else:
# solve for band part (only)
      aSolution = solveBand(aRightHandSide)
      self.__band = invertBand()
    return aSolution

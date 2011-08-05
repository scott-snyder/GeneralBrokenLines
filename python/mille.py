'''
Input/output of MP-II binary records.

Created on Aug 1, 2011

@author: kleinwrt
'''

import array, math

class MilleRecord(object):
  '''
  Millepede-II (binary) record.
  
  Containing information for local (track) and global fit.
 
  The data blocks are collected in two arrays, a real array and
  an integer array, of same length.  The content of the arrays::
         real array              integer array    
     0   0.0                     error count (this record)  
     1   RMEAS, measured value   0                            __iMeas   -+
     2   local derivative        index of local derivative               |
     3   local derivative        index of local derivative               |
     4    ...                                                            | block
         SIGMA, error (>0)       0                            __ iErr    |
         global derivative       label of global derivative              |
         global derivative       label of global derivative             -+
         RMEAS, measured value   0                            __position
         local derivative        index of local derivative
         local derivative        index of local derivative
         ...
         SIGMA, error            0
         global derivative       label of global derivative
         global derivative       label of global derivative
         ...
         global derivative       label of global derivative   __recLen
  '''

  def __init__(self):
    '''
    Create MP-II binary record.
    '''
    self.__position = 1 
    '''@ivar: position in record, usually start of next data block (
       @type: int'''    
    self.__numData = 0  
    '''@ivar: number of data blocks in record 
       @type: int'''
    self.__recLen = 0 
    '''@ivar: record length 
       @type: int'''   
    self.__iMeas = 0    
    '''@ivar: position of value in current data block
       @type: int'''
    self.__iErr = 0     
    '''@ivar: position of error in current data block 
       @type: int'''
    self.__inder = array.array('i') 
    '''@ivar: array with markers (0) and labels 
       @type: array(int32)'''
    self.__glder = array.array('f') 
    '''@ivar: array with values, errors and derivatives 
       @type: (float32)'''
    
  def addData(self, dataList):
    '''
    Add data block to (end of) record.
    
    @param dataList: list with measurement, error, labels and derivatives
    @type dataList: list
    '''
    if (self.__numData == 0): # first word is error counter
      self.__inder.append(0)
      self.__glder.append(0.) 
    self.__numData += 1     
 
    [aMeas, aPrec, indLocal, derLocal, labGlobal, derGlobal] = dataList
    self.__inder.append(0)
    self.__glder.append(aMeas)
    self.__inder.fromlist(indLocal)
    self.__glder.fromlist(derLocal)    
    self.__inder.append(0)
    self.__glder.append(1.0 / math.sqrt(aPrec)) # convert to error
    self.__inder.fromlist(labGlobal)
    self.__glder.fromlist(derGlobal)
    
  def getData(self):
    '''
    Get data block from current position in record.
    
    @return: list with measurement, error, labels and derivatives
    @rtype: list
    '''
    aMeas = self.__glder[self.__iMeas]
    indLocal = []
    derLocal = []
    for i in range(self.__iMeas + 1, self.__iErr):
      indLocal.append(self.__inder[i])
      derLocal.append(self.__glder[i])
    aPrec = 1.0 / self.__glder[self.__iErr] ** 2 # convert to precision 
    indGlobal = []
    derGlobal = []
    for i in range(self.__iErr + 1, self.__position):
      indGlobal.append(self.__inder[i])
      derGlobal.append(self.__glder[i])    
    return [ aMeas, aPrec, indLocal, derLocal, indGlobal, derGlobal ]   
    
  def addSpecial(self, iSpecial, gSpecial, aTag=0):
    '''
    Add special data block to record::

         real array              integer array    
         0.0                     0    
         -float(NSP)-0.1*tag     0   ! indicates special data of length NSP
         following NSP floating and NSP integer data
    
    @param iSpecial: list of labels
    @type  iSpecial: list(int)
    @param gSpecial: list of values
    @type  gSpecial: list(float)
    @param aTag: tag (1: external seed (matrix))
    @type  aTag: int
    '''
    self.__inder.append(0)
    self.__glder.append(0.)
    self.__inder.append(0)
    self.__glder.append(-len(iSpecial) - 0.1 * aTag)
    self.__inder.fromlist(iSpecial)
    self.__glder.fromlist(gSpecial)
    
  def getSpecial(self):
    '''
    Get special data block from current position.
    
    @return: special data block (list of labels, list of values)
    @rtype: list
    '''
    iSpecial = []
    gSpecial = []
    for i in range(self.__iErr + 1, self.__position):
      iSpecial.append(self.__inder[i])
      gSpecial.append(self.__glder[i])    
    return [ iSpecial, gSpecial ]         
    
  def printRecord(self):
    '''
    Print record.
    '''
    print " MilleRecord, len: ", len(self.__inder)
    print self.__inder
    print self.__glder
    
  def writeRecord(self, aFile):
    '''
    Write record to file.
    
    @param aFile: (binary) file
    @type  aFile: file
    '''
    header = array.array('i') # header with number of words
    header.append(len(self.__inder) * 2)
    header.tofile(aFile)
    self.__glder.tofile(aFile)
    self.__inder.tofile(aFile)
    
  def readRecord(self, aFile):
    '''
    Read record from file.
    
    @param aFile: (binary) file
    @type  aFile: file
    '''
    header = array.array('i') # header with number of words
    header.fromfile(aFile, 1)
    self.__recLen = header[0] / 2
    self.__glder.fromfile(aFile, self.__recLen)
    self.__inder.fromfile(aFile, self.__recLen)
    
  def moreData(self):
    '''
    Locate next data block.
    
    @return: next block exists
    @rtype: bool
    '''
    if (self.__position < self.__recLen):
      while (self.__position < self.__recLen and self.__inder[self.__position] != 0):
        self.__position += 1
      self.__iMeas = self.__position 
      self.__position += 1
      while (self.__position < self.__recLen and self.__inder[self.__position] != 0):
        self.__position += 1
      self.__iErr = self.__position 
      self.__position += 1
      while (self.__position < self.__recLen and self.__inder[self.__position] != 0):
        self.__position += 1
      self.__numData += 1  
      return True
    else:
      return False
    
  def specialDataTag(self):
    '''
    Get special data tag from block.
    
    @return: tag or -1 for ordinary data block 
    @rtype: int
    '''
    aTag = -1
    if (self.__iMeas + 1 == self.__iErr and self.__glder[self.__iErr] < 0):
      aTag = int(-self.__glder[self.__iErr] * 10. + 0.5) % 10
    return aTag

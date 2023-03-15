#include "MilleBinary.h"
#include <iostream>

using namespace gbl;

extern "C" {

MilleBinary* MilleBinaryCtor(const char* fileName, int filenamesize, int doublePrecision, int keepZeros, int aSize) {
#ifdef JNA_DEBUG
  std::cout << "MilleBinaryCtor(" << fileName << ", " << filenamesize << ", " 
    << doublePrecision << ", " << keepZeros << ", " << aSize << ")" << std::endl;
#endif
	std::string binName(fileName,filenamesize);
	MilleBinary* mb = new MilleBinary(binName,doublePrecision!=0,keepZeros!=0,aSize);
#ifdef JNA_DEBUG
  std::cout << "MilleBinary created at " << mb << std::endl;
#endif
  return mb;
}

void MilleBinary_close(MilleBinary* self) {
#ifdef JNA_DEBUG
  std::cout << "MilleBinary_close(" << self << ")" << std::endl;
#endif
	self->Close();
}
}

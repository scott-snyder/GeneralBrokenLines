#include "MilleBinary.h"
#include <iostream>

using namespace gbl;

extern "C" {

/**
 * MilleBinaryCtor
 *
 * Non-mangled constructor name for creating a new MilleBinary file.
 *
 * We return the raw pointer to this object so that JNA can effectively access it.
 *
 * Unfortunately the translation of booleans to/from java is not very stable, so it
 * is safer to simply pass integers and check if they are 0 (for false) or nonzero (for true).
 */
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

/**
 * Closing a MilleBinary file is the same as destructing it
 *
 * Since the destructor of the MilleBinary class is what handles
 * performing the final write operations, we simply `delete` the 
 * object pointed to by the passed pointer.
 *
 * @note This means the object on the JNA side will be invalid
 * and will cause a program crash if it is accessed after using this
 * function on it!
 */
void MilleBinary_close(MilleBinary* self) {
#ifdef JNA_DEBUG
	std::cout << "MilleBinary_close(" << self << ")" << std::endl;
#endif
	if (self) delete self;
}
}

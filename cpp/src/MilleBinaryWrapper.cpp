#include "MilleBinary.h"
#include <iostream>

using namespace gbl;

extern "C" {
    
    MilleBinary* MilleBinaryCtor(const char* fileName, int filenamesize, int doublePrecision, int aSize) {
        
        std::string binName(fileName,filenamesize);
        
        return new MilleBinary(binName,doublePrecision,aSize);
    }
    
    void MilleBinary_close(MilleBinary* self) {
        
        self->Close();
    }
}

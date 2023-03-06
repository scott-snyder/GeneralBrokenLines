#include "GblPoint.h"
#include "Eigen/Core"


const int NROW = 5;
const int NCOL = 5;

using namespace gbl;
using namespace Eigen;

extern "C" { 
    GblPoint* GblPointCtor(double matrixArray[NROW*NCOL]) {
        
        Map<Matrix5d> jacobian(matrixArray,5,5);
        return new GblPoint(jacobian);
        
    }
    
    void GblPoint_printPoint(const GblPoint* self, unsigned int level) {
        self->printPoint(level);
    }
    
    unsigned int GblPoint_hasMeasurement(const GblPoint* self) {
        return self->hasMeasurement();
    }
    
    
    double GblPoint_getMeasPrecMin(const GblPoint* self) {
        return self->getMeasPrecMin();
    }
    
    //Only supporting:
    //2D position residual
    //2x2 projection matrix
    
    void GblPoint_addMeasurement2D(GblPoint* self, 
                                   double *projArray,
                                   double *resArray,
                                   double *precArray, 
                                   double minPrecision) { 
        
        Map<Matrix2d> aProjection(projArray,2,2);
        Map<Vector2d> aResiduals(resArray, 2);
        Map<Vector2d> aPrecision(precArray,2);
        
        self->addMeasurement(aProjection, aResiduals, aPrecision, minPrecision);
    }

    
    //Only support vector precision
    void GblPoint_addScatterer(GblPoint* self, double *resArray, double *precArray) {
        
        Map<Vector2d> aResiduals(resArray,2);
        Map<Vector2d> aPrecision(precArray,2);
        
        self->addScatterer(aResiduals,aPrecision);
    }
    
    
    void GblPoint_addGlobals(GblPoint* self, int *labels, int nlabels, double* derArray) {
        std::vector<int> aLabels;
        for (int i=0; i<nlabels; i++) {
            aLabels.push_back(labels[i]);
        }
        Map<Eigen::MatrixXd> derivatives(derArray,1,nlabels);
        self->addGlobals(aLabels, derivatives);
    }

    int GblPoint_getNumGlobals(GblPoint* self) {
        return self->getNumGlobals();
    }
    
    //TODO revisit these!

    //Should I add the number of nlabels?
    void GblPoint_getGlobalLabels(GblPoint* self, int* labels) {
        
        std::vector<int> glabels;
        self->getGlobalLabels(glabels);

        //std::cout<<"GblPointWrapper::glabels"<<std::endl;
        
        for (unsigned int il = 0 ; il < glabels.size(); il++) {
            labels[il] = glabels.at(il);
            //std::cout<<glabels.at(il)<<std::endl;
        }
    }  
    
    //Should I add the number of derivatives? -  Row/Col? CHECK CHECK CHECK
    void GblPoint_getGlobalDerivatives(GblPoint* self, double* gders) {
        
        int ngders = self->getNumGlobals();
        Eigen::MatrixXd e_gders(1,ngders);
        self->getGlobalDerivatives(e_gders);
        
        //std::cout<<"GblPointWrapper::e_ders"<<std::endl;
        //std::cout<<e_gders<<std::endl;

        Map<MatrixXd>(gders,1,ngders) = e_gders;
        
        
                
    }

}


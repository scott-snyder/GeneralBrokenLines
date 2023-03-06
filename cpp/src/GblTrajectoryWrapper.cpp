#include "GblTrajectory.h"

using namespace gbl;
using namespace Eigen;

extern "C" {
    
    GblTrajectory* GblTrajectoryCtor(int flagCurv, int flagU1dir, int flagU2dir) {
        
        return new GblTrajectory(flagCurv, flagU1dir, flagU2dir);
        
    }
        
    //Simple trajectory constructor wrapper
    GblTrajectory* GblTrajectoryCtorPtrArray(GblPoint* points[], int npoints, 
                                             int flagCurv, int flagU1dir, int flagU2dir) {
        
        
        std::vector<GblPoint> aPointList;
        
        for (int i=0; i<npoints; i++) {
            //get the point pointer
            GblPoint* gblpoint = points[i];
            
            //add it to the vector
            aPointList.push_back(*(gblpoint));
        }
        
        return new GblTrajectory(aPointList, flagCurv, flagU1dir, flagU2dir);
    }

    //Simple trajectory constructor with seed wrapper
    
    GblTrajectory* GblTrajectoryCtorPtrArraySeed(GblPoint* points[], int npoints,
                                                 int aLabel, double seedArray[],
                                                 int flagCurv, int flagU1dir, int flagU2dir) {
        
        std::vector<GblPoint> aPointList;
        
        for (int i=0; i<npoints; i++) {
            //get the point pointer
            GblPoint* gblpoint = points[i];
            
            //add it to the vector
            aPointList.push_back(*(gblpoint));
        }
        
        Map<Matrix5d> seed(seedArray,5,5);
        
        return new GblTrajectory(aPointList, aLabel, seed, flagCurv, flagU1dir, flagU2dir);
        
    }
    
    //Composed trajectory constructor for 2 body decay

    GblTrajectory* GblTrajectoryCtorPtrComposed(GblPoint* points_1[], int npoints_1, double trafo_1[],
                                                GblPoint* points_2[], int npoints_2, double trafo_2[]) {
        
        std::vector<std::pair<std::vector<GblPoint>, Eigen::MatrixXd> > pointsAndTransList;
        
        //first track
        std::vector<GblPoint> points_trk_1;
        
        for (int i=0; i<npoints_1; i++) {
            //get the point pointer
            GblPoint* gblpoint = points_1[i];
            //add it to the vector 
            points_trk_1.push_back(*(gblpoint));
        }
        
        MatrixXd inner_1(2,3);
        inner_1(0,0)=trafo_1[0];
        inner_1(0,1)=trafo_1[1];
        inner_1(0,2)=trafo_1[2];

        inner_1(1,0)=trafo_1[3];
        inner_1(1,1)=trafo_1[4];
        inner_1(1,2)=trafo_1[5];

        std::pair<std::vector<GblPoint>, MatrixXd> track_trafo_1 = std::make_pair(points_trk_1, inner_1);
        
        //second track
        std::vector<GblPoint> points_trk_2;
        
        for (int i=0; i<npoints_2; i++) {
            //get the point pointer
            GblPoint* gblpoint = points_2[i];
            //add it to the vector 
            points_trk_2.push_back(*(gblpoint));
        }

        MatrixXd inner_2(2,3);
        inner_2(0,0)=trafo_2[0];
        inner_2(0,1)=trafo_2[1];
        inner_2(0,2)=trafo_2[2];

        inner_2(1,0)=trafo_2[3];
        inner_2(1,1)=trafo_2[4];
        inner_2(1,2)=trafo_2[5];

        std::pair<std::vector<GblPoint>, MatrixXd> track_trafo_2 = std::make_pair(points_trk_2, inner_2);

        pointsAndTransList.push_back(track_trafo_1);
        pointsAndTransList.push_back(track_trafo_2);
        
        return new GblTrajectory(pointsAndTransList);

    }

    void GblTrajectory_fit(GblTrajectory* self, double* Chi2, int* Ndf, double* lostWeight, char* c_optionList, unsigned int aLabel) {
        
        std::string optionList(c_optionList);
        self->fit(*Chi2, *Ndf, *lostWeight, optionList,aLabel);
    }
  
    void GblTrajectory_delete(GblTrajectory* self) {
      if (self)
        delete self;
    }

    int GblTrajectory_isValid(GblTrajectory* self) {
        
        return (int) self->isValid();
    }

    int GblTrajectory_getNumPoints(GblTrajectory* self) {
      return (int) self->getNumPoints();
    }

    void GblTrajectory_printTrajectory(GblTrajectory* self, int level) {
        return self->printTrajectory();
    }

    void GblTrajectory_printData(GblTrajectory* self) {
        return self->printData();
    }
    
    void GblTrajectory_printPoints(GblTrajectory* self, int level) {
        return self->printPoints(level);
    }

    //Only 5-vector and 5x5 cov matrix.
    void GblTrajectory_getResults(GblTrajectory* self, int aSignedLabel, double* localPar, int* nLocalPar,
                                  double * localCov, int* sizeLocalCov) {

        Eigen::VectorXd e_localPar(5);
        Eigen::MatrixXd e_localCov(5,5);
        
        self->getResults(aSignedLabel, e_localPar, e_localCov);

        //std::cout<<"gblTrajectoryWrapper::getResults"<<std::endl;
        //std::cout<<e_localPar<<std::endl;
        //std::cout<<e_localCov<<std::endl;
        
        Map<Vector5d>(localPar, 5) = e_localPar;
        Map<Matrix5d>(localCov, 5, 5) = e_localCov; 
        *nLocalPar = 5;
        *sizeLocalCov = 5;
        
    }

    //Wrapper to get the residuals - Assume 2d residuals max

    void GblTrajectory_getMeasResults(GblTrajectory* self, int aLabel, int* numData, 
                                      double* aResiduals, double* aMeasErrors, double* aResErrors, 
                                      double* aDownWeights) {
        
        Eigen::VectorXd e_aResiduals(2);
        Eigen::VectorXd e_aMeasErrors(2);
        Eigen::VectorXd e_aResErrors(2);
        Eigen::VectorXd e_aDownWeights(2);
        unsigned int num_data = 0;
        
        unsigned int out = self->getMeasResults(aLabel, num_data, e_aResiduals, e_aMeasErrors,
                                                e_aResErrors, e_aDownWeights);
        
        *numData = num_data;
        
        for (unsigned int i = 0; i < num_data; i++) {
            aResiduals[i] = e_aResiduals(i);
            aMeasErrors[i] = e_aMeasErrors(i);
            aResErrors[i] = e_aResErrors(i);
            aDownWeights[i] = e_aDownWeights(i);
        }
        
    }
    
    
    void GblTrajectory_milleOut(GblTrajectory* self, MilleBinary* millebinary) {
        self->milleOut(*millebinary);
    }

    
    

}




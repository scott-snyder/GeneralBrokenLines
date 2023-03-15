#include "GblTrajectory.h"

using namespace gbl;
using namespace Eigen;

extern "C" {

/**
 * convert the pointer array of GblPoints into a vector holding the objects
 *
 * @note We *move* the data pointed to into the vector so the array *should
 * not* be accessed after this call is made.
 *
 * @param[in] points array of pointers to GblPoints to put into vector
 * @param[in] npoints number of points (size of array)
 * @return vector of GblPoints with same content as array
 */
std::vector<GblPoint> ptr_array_to_vector(GblPoint* points[], int npoints) {
	std::vector<GblPoint> points_vec;
	// since we already know the size of the vector, we 'reserve' the size
	// so that the vector doesn't need to waste time copying/moving the GblPoints
	// around as it grows in size
	// we do *not* use 'resize' since that would involve default-constructing
	// all the GblPoints
	points_vec.reserve(npoints);

	for (int i{0}; i < npoints; ++i) {
		// get the pointer
		GblPoint* gblpoint = points[i];
		// MOVE the data into the vector,
		//   this transfers ownership into the vector and invalidates the array
		points_vec.emplace_back(*(gblpoint));
	}

	return points_vec;
}

GblTrajectory* GblTrajectoryCtor(int flagCurv, int flagU1dir, int flagU2dir) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectoryCtor(" << flagCurv << ", " << flagU1dir << ", " << flagU2dir << ")" << std::endl;
#endif
	return new GblTrajectory(flagCurv, flagU1dir, flagU2dir);
	
}
	
//Simple trajectory constructor wrapper
GblTrajectory* GblTrajectoryCtorPtrArray(GblPoint* points[], int npoints, 
										 int flagCurv, int flagU1dir, int flagU2dir) {
#ifdef JNA_DEBUG
	std::cout << "GblTracjectoryCtorPtrArray("
		<< points << ", " << npoints << ", "
		<< flagCurv << ", " << flagU1dir << ", " << flagU2dir
		<< ")" << std::endl;
#endif
	
	return new GblTrajectory(ptr_array_to_vector(points, npoints), flagCurv, flagU1dir, flagU2dir);
}

//Simple trajectory constructor with seed wrapper

GblTrajectory* GblTrajectoryCtorPtrArraySeed(GblPoint* points[], int npoints,
											 int aLabel, double seedArray[],
											 int flagCurv, int flagU1dir, int flagU2dir) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectoryCtorPtrArraySeed("
		<< points << ", " << npoints << ", "
		<< aLabel << ", " << seedArray << ", "
		<< flagCurv << ", " << flagU1dir << ", " << flagU2dir
		<< ")" << std::endl;
#endif
	
	Map<Matrix5d> seed(seedArray,5,5);
	
	return new GblTrajectory(ptr_array_to_vector(points, npoints), aLabel, seed, flagCurv, flagU1dir, flagU2dir);
}

//Composed trajectory constructor for 2 body decay

GblTrajectory* GblTrajectoryCtorPtrComposed(GblPoint* points_1[], int npoints_1, double trafo_1[],
											GblPoint* points_2[], int npoints_2, double trafo_2[]) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectoryCtorPtrComposed("
		<< points_1 << ", " << npoints_1 << ", " << trafo_1 << ", "
		<< points_2 << ", " << npoints_2 << ", " << trafo_2 << ")"
		<< std::endl;
#endif
	
	
	MatrixXd inner_1(2,3);
	inner_1(0,0)=trafo_1[0];
	inner_1(0,1)=trafo_1[1];
	inner_1(0,2)=trafo_1[2];

	inner_1(1,0)=trafo_1[3];
	inner_1(1,1)=trafo_1[4];
	inner_1(1,2)=trafo_1[5];

	std::pair<std::vector<GblPoint>, MatrixXd> track_trafo_1 = std::make_pair(ptr_array_to_vector(points_1, npoints_1), inner_1);
	
	//second track

	MatrixXd inner_2(2,3);
	inner_2(0,0)=trafo_2[0];
	inner_2(0,1)=trafo_2[1];
	inner_2(0,2)=trafo_2[2];

	inner_2(1,0)=trafo_2[3];
	inner_2(1,1)=trafo_2[4];
	inner_2(1,2)=trafo_2[5];

	std::pair<std::vector<GblPoint>, MatrixXd> track_trafo_2 = std::make_pair(ptr_array_to_vector(points_1, npoints_1), inner_2);
	
	return new GblTrajectory({track_trafo_1, track_trafo_2});

}

void GblTrajectory_fit(GblTrajectory* self, double* Chi2, int* Ndf, double* lostWeight, char* c_optionList, unsigned int aLabel) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_fit("
		<< self << ", " << Chi2 << ", " << Ndf << ", " << lostWeight << ", " << c_optionList << ", " << aLabel << ")"
		<< std::endl;
#endif
	
	std::string optionList(c_optionList);
	self->fit(*Chi2, *Ndf, *lostWeight, optionList,aLabel);
}

void GblTrajectory_delete(GblTrajectory* self) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_delete(" << self << ")" << std::endl;
#endif
	if (self) delete self;
}

int GblTrajectory_isValid(GblTrajectory* self) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_isValid(" << self << ")" << std::endl;
#endif
	return (int) self->isValid();
}

int GblTrajectory_getNumPoints(GblTrajectory* self) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_getNumPoints(" << self << ")" << std::endl;
#endif
	return (int) self->getNumPoints();
}

void GblTrajectory_printTrajectory(GblTrajectory* self, int level) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_printTrajectory(" << self << ", " << level << ")" << std::endl;
#endif
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
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_getResults(" << self 
		<< ", " << aSignedLabel << ", " << localPar << ", " << nLocalPar
		<< ", " << localCov << ", " << sizeLocalCov << ")" << std::endl;
#endif

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
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_getMeasResults(" << self 
		<< ", " << aLabel << ", " << numData << ", " << aResiduals
		<< ", " << aMeasErrors << ", " << aResErrors << ", " << aDownWeights << ")" << std::endl;
#endif
	
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
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_milleOut(" << self 
		<< ", " << millebinary << ")" << std::endl;
#endif
	self->milleOut(*millebinary);
}

}




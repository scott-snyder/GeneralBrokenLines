#include "MilleBinary.h"
#include "GblPoint.h"
#include "GblTrajectory.h"
#include "GblUtilities.h" 

#include "Eigen/Core"

const int NROW = 5;
const int NCOL = 5;

using namespace gbl;
using namespace Eigen;

#ifdef JNA_DEBUG
#include <iostream>
int num_gbl_point = 0;
#endif

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
GblPoint* GblPointCtor(double matrixArray[NROW*NCOL]) {
	Map<Matrix5d> jacobian(matrixArray,5,5);
	GblPoint* self = new GblPoint(jacobian);
#ifdef JNA_DEBUG
	//std::cout << "GblPointCtor at " << self << " " << ++num_gbl_point << std::endl;
#endif
	return self;
}

void GblPoint_delete(GblPoint* self) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_delete(" << self << ") " << --num_gbl_point << std::endl;
#endif
	if (self) delete self;
}


void GblPoint_printPoint(const GblPoint* self, unsigned int level) {
#ifdef JNA_DEBUG
	//std::cout << "GblPoint_printPoint(" << self << ", " << level << ")" << std::endl;
#endif
	self->printPoint(level);
}

int GblPoint_getNumMeasurements(GblPoint* self) {
  return (self->getMeasEnd() - self->getMeasBegin());
}

//Only supporting:
//2D position residual
//2x2 projection matrix

void GblPoint_addMeasurement2D(GblPoint* self, 
								 double *projArray,
								 double *resArray,
								 double *precArray, 
								 double minPrecision) { 
#ifdef JNA_DEBUG
	std::cout << "GblPoint_addMeasurement2D("
		<< self << ", " << projArray << ", " << resArray << ", " << precArray << ", " << minPrecision
		<< ")" << std::endl;
#endif
	
	Map<Matrix2d> aProjection(projArray,2,2);
	Map<Vector2d> aResiduals(resArray, 2);
	Map<Vector2d> aPrecision(precArray,2);
	
	self->addMeasurement(aProjection, aResiduals, aPrecision, minPrecision);
}


//Only support vector precision
void GblPoint_addScatterer(GblPoint* self, double *resArray, double *precArray) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_addScatterer(" << self << ", " 
		<< resArray << ", " << precArray << ")" << std::endl;
#endif
	// chose to do the Vector2d addScatterer since
	// PF's original comment "only support vector precision"
	Eigen::Vector2d aResiduals(resArray);
	Eigen::Vector2d aPrecision(precArray);
	
	self->addScatterer(aResiduals,aPrecision);
}

void GblPoint_addGlobals(GblPoint* self, int *labels, int nlabels, double* derArray) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_addGlobals(" << self
		<< ", " << labels << ", " << nlabels << ", " << derArray << ")" << std::endl;
#endif
	std::vector<int> aLabels;
	for (int i=0; i<nlabels; i++) {
		aLabels.push_back(labels[i]);
	}
	Map<Eigen::MatrixXd> derivatives(derArray,1,nlabels);
	self->addGlobals(aLabels, derivatives);
}

void GblPoint_getGlobalLabelsAndDerivatives(GblPoint* self, int* nlabels, int** labels, double** ders) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_getGlobalLabelsAndDerivatives("
		<< self << ", " << labels << ", " << ders
		<< ")" << std::endl;
#endif
	std::vector<int> glabels;
	std::vector<double> gders;

#ifdef JNA_DEBUG
  std::cout << "  GblPoint has " << std::flush
    << self->getMeasEnd() - self->getMeasBegin()
    << " measurements." << std::endl;
#endif

	//Should I add the number of derivatives? -  Row/Col? CHECK CHECK CHECK
	self->getGlobalLabelsAndDerivatives(
			0 /* aMeas */, 0 /* aRow	*/,
			glabels, gders);

#ifdef JNA_DEBUG
	std::cout <<"  Aquired " << glabels.size() << " labels and " << gders.size() << " derivatives." <<std::endl;
#endif
	
	*nlabels = glabels.size();
	*labels = new int[*nlabels];
	*ders   = new double[*nlabels];

#ifdef JNA_DEBUG
	std::cout <<"  Allocated return arrays." << std::endl;
#endif
	
	
	for (std::size_t il{0}; il < *nlabels; ++il) {
		(*labels)[il] = glabels.at(il);
		//std::cout<<glabels.at(il)<<std::endl;
	}

	//std::cout<<"GblPointWrapper::gders"<<std::endl;
	
	for (std::size_t id{0}; id < *nlabels; ++id) {
		(*ders)[id] = gders.at(id);
		//std::cout<<gders.at(il)<<std::endl;
	}

	// set array using Eigen
	//Map<MatrixXd>(ders,1,gders.size()) = gders;
#ifdef JNA_DEBUG
  std::cout << "  Done with assignment and leaving." << std::endl;
#endif
}

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
		// COPY the data into the vector,
		points_vec.emplace_back(*(gblpoint));
#ifdef JNA_DEBUG
		std::cout << "COPY GblPoint " << gblpoint << " -> " << &(points_vec.back()) << std::endl;
#endif
	}

	return points_vec;
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

//Gbl Detector Layer representation
//aCenter 3-vector
//aResolution 2-vector
//aPrecision 2-vector
//measTrafo 3x3 matrix
//alignTrafo 3x3 matrix

GblDetectorLayer* GblDetectorLayerCtor(const char* aName, int aLayer, int aDim, double thickness,
																			 double aCenter[], double aResolution[], double aPrecision[],
																			 double measTrafo[], double alignTrafo[]) {
		
	// uses eigen's Map structure to decompose an array into our type of Vector/Matrix
	// need to then copy that object into a local variable since Map's only allow const
	// reference access
	Vector3d e_aCenter(Map<Vector3d>(aCenter,3));
	Vector2d e_aResolution(Map<Vector2d>(aResolution,2));
	Vector2d e_aPrecision(Map<Vector2d>(aPrecision,2));
	Matrix3d e_measTrafo(Map<Matrix3d>(measTrafo,3,3));
	Matrix3d e_alignTrafo(Map<Matrix3d>(alignTrafo,3,3));
	
	GblDetectorLayer* layer = new GblDetectorLayer(aName, aLayer, aDim, thickness, e_aCenter, e_aResolution, e_aPrecision, e_measTrafo, e_alignTrafo);
	return layer;
}

void GblDetectorLayer_delete(GblDetectorLayer* self) {
	if (self) delete self;
}

void GblDetectorLayer_print(GblDetectorLayer* self) {
	self->print();
}

double GblDetectorLayer_getRadiationLength(GblDetectorLayer* self) {
	return self->getRadiationLength();
}

//TODO Implement wrappers for
//getResolution
//getPrecision
//getCenter
//getMeasSystemDirs
//getAlignSystemDirs


//Helix prediction on layer

GblHelixPrediction* GblDetectorLayer_intersectWithHelix(GblDetectorLayer* self, GblSimpleHelix* hlx) {
	
	Vector3d center = self->getCenter();
	Vector3d udir		= (self->getMeasSystemDirs()).row(0);
	Vector3d vdir		= (self->getMeasSystemDirs()).row(1);
	
	return new GblHelixPrediction(hlx->getPrediction(center, udir, vdir));
}



//Simple Helix
GblSimpleHelix* GblSimpleHelixCtor(double aRinv, double aPhi0, double aDca, double aDzds, double aZ0) {
	return new GblSimpleHelix(aRinv, aPhi0, aDca, aDzds, aZ0);
}

void GblSimpleHelix_delete(GblSimpleHelix* self) {
	if (self) delete self;
}

double GblSimpleHelix_getPhi(GblSimpleHelix* self, double aRadius) {
	return self->getPhi(aRadius);
}

double GblSimpleHelix_getArcLengthR(GblSimpleHelix* self, double aRadius) {
	return self->getPhi(aRadius);
}

double GblSimpleHelix_getArcLengthXY(GblSimpleHelix* self, double xPos, double yPos) {
	return self->getArcLengthXY(xPos, yPos);
}


void GblSimpleHelix_moveToXY(GblSimpleHelix* self, double xPos, double yPos,
														 double* newPhi0, double* newDca, double* newZ0) {
	self->moveToXY(xPos, yPos,
								 *newPhi0, *newDca, *newZ0);
}

//refPos, uDir and vDir are 3-vectors
GblHelixPrediction* GblSimpleHelix_getPrediction(GblSimpleHelix* self, double refPos[], double uDir[], double vDir[]) {
	
	Map<Vector3d> e_refPos(refPos,3);
	Map<Vector3d> e_uDir(uDir,3);
	Map<Vector3d> e_vDir(vDir,3);
	GblHelixPrediction prediction = self->getPrediction(e_refPos,e_uDir,e_vDir);
	
	/*std::cout<<"Cross Check GBL predicted position!"<<std::endl;
	std::cout<<prediction->getPosition()<<std::endl;
	
	std::cout<<"Cross Check GBL meas predicted position!"<<std::endl;
	std::cout<<prediction->getMeasPred()<<std::endl;*/
	
	// JNA only deals with pointers so we need to dynamically create a new copy
	return new GblHelixPrediction(prediction);
}


//Helix Prediction
GblHelixPrediction* GblHelixPredictionCtor(double sArc, double aPred[], double tDir[], double uDir[], double vDir[],
																				 double nDir[], double aPos[]) {
	Map<Vector2d> e_aPred(aPred,2);
	Map<Vector3d> e_tDir(tDir,3);
	Map<Vector3d> e_uDir(uDir,3);
	Map<Vector3d> e_vDir(vDir,3);
	Map<Vector3d> e_nDir(nDir,3);
	Map<Vector3d> e_aPos(aPos,3);
	
	
	return new GblHelixPrediction(sArc, e_aPred, e_tDir, e_uDir, e_vDir, 
																e_nDir, e_aPos);
}

void GblHelixPrediction_delete(GblHelixPrediction* self) {
	if (self) delete self;
}

double GblHelixPrediction_getArcLength(GblHelixPrediction* self) {
	return self->getArcLength();
}

void GblHelixPrediction_getMeasPred(GblHelixPrediction* self, double* prediction) {
	Vector2d e_pred = self->getMeasPred();
	
	prediction[0] = e_pred(0);
	prediction[1] = e_pred(1);
}

void GblHelixPrediction_getPosition(GblHelixPrediction* self, double* position) {
	Vector3d e_pos = self->getPosition();
	
	position[0] = e_pos(0);
	position[1] = e_pos(1);
	position[2] = e_pos(2);
}

void GblHelixPrediction_getDirection(GblHelixPrediction* self, double direction[]) {
	Vector3d e_dir = self->getDirection();
	
	direction[0] = e_dir(0);
	direction[1] = e_dir(1);
	direction[2] = e_dir(2);
}

double GblHelixPrediction_getCosIncidence(GblHelixPrediction* self) {
	return self->getCosIncidence();
}

void GblHelixPrediction_getCurvilinearDirs(GblHelixPrediction* self, double curvilinear[]) {
	Matrix<double,2,3> curDirs = self->getCurvilinearDirs();
	
	//std::cout<<"Check curvilinear Directions" <<std::endl;
	//std::cout<<curDirs<<std::endl;
					
	curvilinear[0] = curDirs(0,0);
	curvilinear[1] = curDirs(0,1);
	curvilinear[2] = curDirs(0,2);
	curvilinear[3] = curDirs(1,0);
	curvilinear[4] = curDirs(1,1);
	curvilinear[5] = curDirs(1,2);
}
}


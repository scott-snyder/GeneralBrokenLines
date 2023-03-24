#include "GblUtilities.h" 

using namespace gbl;
using namespace Eigen;

extern "C" {
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

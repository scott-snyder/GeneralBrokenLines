#include "GblPoint.h"
#include "Eigen/Core"

const int NROW = 5;
const int NCOL = 5;

using namespace gbl;
using namespace Eigen;

#ifdef JNA_DEBUG
int num_gbl_point = 0;
#endif

extern "C" { 
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

void GblPoint_getGlobalLabelsAndDerivatives(GblPoint* self, int* nlabels, int labels[], double ders[]) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_getGlobalLabelsAndDerivatives("
		<< self << ", " << labels << ", " << ders
		<< ")" << std::endl;
#endif
	std::vector<int> glabels;
	std::vector<double> gders;

	//Should I add the number of derivatives? -  Row/Col? CHECK CHECK CHECK
	self->getGlobalLabelsAndDerivatives(
			0 /* aMeas */, 0 /* aRow	*/,
			glabels, gders);

	//std::cout<<"GblPointWrapper::glabels"<<std::endl;
	
	*nlabels = glabels.size();
	labels = new int[*nlabels];
	ders   = new double[*nlabels];
	
	for (std::size_t il{0}; il < *nlabels; ++il) {
		labels[il] = glabels.at(il);
		//std::cout<<glabels.at(il)<<std::endl;
	}

	//std::cout<<"GblPointWrapper::gders"<<std::endl;
	
	for (std::size_t id{0}; id < *nlabels; ++id) {
		ders[id] = gders.at(id);
		//std::cout<<gders.at(il)<<std::endl;
	}

	// set array using Eigen
	//Map<MatrixXd>(ders,1,gders.size()) = gders;
}

}


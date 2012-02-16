/*
 * GblPoint.cpp
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

#include "GblPoint.h"

/// Create a point.
/**
 * Create point on (initial) trajectory. Needs transformation jacobian from previous point.
 * \param [in] aJacobian Transformation jacobian from previous point
 */
GblPoint::GblPoint(const TMatrixD &aJacobian) :
		theLabel(0), theOffset(0), measDim(0), transFlag(false), measTransformation(), scatFlag(
				false), localDerivatives(), globalLabels(), globalDerivatives() {
	// TODO Auto-generated constructor stub
	for (unsigned int i = 0; i < 5; i++) {
		for (unsigned int j = 0; j < 5; j++) {
			p2pJacobian(i, j) = aJacobian(i, j);
		}
	}
}

GblPoint::~GblPoint() {
	// TODO Auto-generated destructor stub
}

/// Add a mesurement to a point.
/**
 * Add measurement with diagonal precision (inverse covariance) matrix.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aProjection Projection from measurement to local system
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (diagonal)
 */
void GblPoint::addMeasurement(const TMatrixD &aProjection,
		const TVectorD &aResiduals, const TVectorD &aPrecision) {
	measDim = aResiduals.GetNrows();
	unsigned int iOff = 5 - measDim;
	for (unsigned int i = 0; i < measDim; i++) {
		measResiduals(iOff + i) = aResiduals[i];
		measPrecision(iOff + i) = aPrecision[i];
		for (unsigned int j = 0; j < measDim; j++) {
			measProjection(iOff + i, iOff + j) = aProjection(i, j);
		}
	}
}

/// Add a mesurement to a point.
/**
 * Add measurement on local system with arbitrary precision (inverse covariance) matrix.
 * Will be diagonalized.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (matrix)
 */
void GblPoint::addMeasurement(const TVectorD &aResiduals,
		const TMatrixDSym &aPrecision) {
	measDim = aResiduals.GetNrows();
	TMatrixDSymEigen measEigen(aPrecision);
	measTransformation = measEigen.GetEigenVectors();
	measTransformation.T();
	transFlag = true;
	TVectorD transResiduals = measTransformation * aResiduals;
	TVectorD transPrecision = measEigen.GetEigenValues();
	unsigned int iOff = 5 - measDim;
	for (unsigned int i = 0; i < measDim; i++) {
		measResiduals(iOff + i) = transResiduals[i];
		measPrecision(iOff + i) = transPrecision[i];
		for (unsigned int j = 0; j < measDim; j++) {
			measProjection(iOff + i, iOff + j) = measTransformation(i, j);
		}
	}
}

/// Check for measurement at a point.
/**
 * Get dimension of measurement (0 = none).
 * \return measurement dimension
 */
unsigned int GblPoint::hasMeasurement() {
	return measDim;
}

/// Retrieve measurement of a point.
/**
 * \param [out] aProjection Projection from (diagonalized) measurement to local system
 * \param [out] aResiduals Measurement residuals
 * \param [out] aPrecision Measurement precision (diagonal)
 */
void GblPoint::getMeasurement(SMatrix55 &aProjection, SVector5 &aResiduals,
		SVector5 &aPrecision) {
	aProjection = measProjection;
	aResiduals = measResiduals;
	aPrecision = measPrecision;
}

/// Add a (thin) scatterer to a point.
/**
 * \param [in] aResiduals Scatterer residuals
 * \param [in] aPrecision Scatterer precision (diagonal of inverse covariance matrix)
 */
void GblPoint::addScatterer(const TVectorD &aResiduals,
		const TVectorD &aPrecision) {
	scatFlag = true;
	scatResiduals(0) = aResiduals[0];
	scatResiduals(1) = aResiduals[1];
	scatPrecision(0) = aPrecision[0];
	scatPrecision(1) = aPrecision[1];
}

/// Check for scatterer at a point.
bool GblPoint::hasScatterer() {
	return scatFlag;
}

/// Retrieve scatterer of a point.
/**
 * \param [out] aResiduals Scatterer residuals
 * \param [out] aPrecision Scatterer precision (diagonal)
 */
void GblPoint::getScatterer(SVector2 &aResiduals, SVector2 &aPrecision) {
	aResiduals = scatResiduals;
	aPrecision = scatPrecision;
}

/// Add local derivatives to a point.
/**
 * Point needs to have a measurement.
 * \param [in] aDerivatives Local derivatives (matrix)
 */
void GblPoint::addLocals(const TMatrixD &aDerivatives) {
	if (measDim) {
        	localDerivatives.ResizeTo(aDerivatives);
		if (transFlag) {
			localDerivatives = measTransformation * aDerivatives;
		} else {
			localDerivatives = aDerivatives;
		}
	}
}

/// Retrieve number of local derivatives from a point.
unsigned int GblPoint::getNumLocals() {
	return localDerivatives.GetNcols();
}

/// Retrieve local derivatives from a point.
TMatrixD GblPoint::getLocalDerivatives() {
	return localDerivatives;
}

/// Add global derivatives to a point.
/**
 * Point needs to have a measurement.
 * \param [in] aLabels Global derivatives labels
 * \param [in] aDerivatives Global derivatives (matrix)
 */
void GblPoint::addGlobals(const std::vector<int> &aLabels,
		const TMatrixD &aDerivatives) {
	if (measDim) {
		globalLabels = aLabels;
		globalDerivatives.ResizeTo(aDerivatives);
		if (transFlag) {
			globalDerivatives = measTransformation * aDerivatives;
		} else {
			globalDerivatives = aDerivatives;
		}
	}
}

/// Retrieve number of global derivatives from a point.
unsigned int GblPoint::getNumGlobals() {
	return globalDerivatives.GetNcols();
}

/// Retrieve global derivatives labels from a point.
std::vector<int> GblPoint::getGlobalLabels() {
	return globalLabels;
}

/// Retrieve global derivatives from a point.
TMatrixD GblPoint::getGlobalDerivatives() {
	return globalDerivatives;
}

/// Define label of point
/**
 * \param [in] aLabel Label identifying point
 */
void GblPoint::setLabel(unsigned int aLabel) {
	theLabel = aLabel;
}

/// Retrieve label of point
unsigned int GblPoint::getLabel() {
	return theLabel;
}

/// Define offset for point
/**
 * \param [in] anOffset Offset number
 */
void GblPoint::setOffset(int anOffset) {
	theOffset = anOffset;
}

/// Retrieve offset for point
int GblPoint::getOffset() {
	return theOffset;
}

/// Retrieve point-to-(previous)point jacobian
SMatrix55 GblPoint::getP2pJacobian() {
	return p2pJacobian;
}

/// Define jacobian to previous scatterer
/**
 * \param [in] aJac Jacobian
 */
void GblPoint::addPrevJacobian(const SMatrix55 aJac) {
	int ifail = 0;
// to optimize: need only two last rows of inverse
	prevJacobian = aJac.Inverse(ifail);
}

/// Define jacobian to next scatterer
/**
 * \param [in] aJac Jacobian
 */
void GblPoint::addNextJacobian(const SMatrix55 aJac) {
	nextJacobian = aJac;
}

/// Retrieve derivatives of local track model
/**
 * Linearized track model: F_u(q/p,u',u) = J*u + S*u' + d*q/p,
 * W is inverse of S, negated for backward propagation.
 * \param [in] aDirection Propagation direction (>0 forward, else backward)
 * \param [out] matW W
 * \param [out] matWJ W*J
 * \param [out] vecWd W*d
 */
void GblPoint::getDerivatives(int aDirection, SMatrix22 &matW, SMatrix22 &matWJ,
		SVector2 &vecWd) {

	if (aDirection < 1) {
		matWJ = prevJacobian.Sub<SMatrix22>(3, 3);
		matW = -prevJacobian.Sub<SMatrix22>(3, 1);
		vecWd = prevJacobian.SubCol<SVector2>(0, 3);
	} else {
		matWJ = nextJacobian.Sub<SMatrix22>(3, 3);
		matW = nextJacobian.Sub<SMatrix22>(3, 1);
		vecWd = nextJacobian.SubCol<SVector2>(0, 3);
	}
	matW.Invert();
	matWJ = matW * matWJ;
	vecWd = matW * vecWd;
}

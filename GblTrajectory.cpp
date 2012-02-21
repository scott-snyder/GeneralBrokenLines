/*
 * GblTrajectory.cpp
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

/** \mainpage General information
 *
 *  \section intro_sec Introduction
 *
 *  For a track with an initial trajectory from a prefit of the
 *  measurements (internal seed) or an external prediction
 *  (external seed) the description of multiple scattering is
 *  added by offsets in a local system. Along the initial
 *  trajectory points are defined with can describe a measurement
 *  or a (thin) scatterer or both. Measurements are arbitrary
 *  functions of the local track parameters at a point (e.g. 2D:
 *  position, 4D: slope+position). The refit provides corrections
 *  to the local track parameters (in the local system) and the
 *  corresponding covariance matrix at any of those points.
 *  Outliers can be down-weighted by use of M-estimators.
 *
 *  The broken lines trajectory is defined by (2D) offsets at the
 *  first and last point and all points with a scatterer. The
 *  prediction for a measurement is obtained by interpolation of
 *  the enclosing offsets and for triplets of adjacent offsets
 *  kink angles are determined. This requires for all points the
 *  jacobians for propagation to the previous and next offset.
 *  These are calculated from the point-to-point jacobians along
 *  the initial trajectory.
 *
 *  Additional local or global parameters can be added and the
 *  trajectories can be written to special binary files for
 *  calibration and alignment with Millepede-II.
 *  (V. Blobel, NIM A, 566 (2006), pp. 5-13).
 *
 *  The conventions for the coordinate systems follow:
 *  Derivation of Jacobians for the propagation of covariance
 *  matrices of track parameters in homogeneous magnetic fields
 *  A. Strandlie, W. Wittek, NIM A, 566 (2006) 687-698.
 *
 *  \section call_sec Calling sequence
 *
 *    -# Create trajectory:\n
 *            <tt>traj = GblTrajectory(..)</tt>
 *    -# For all points on initial trajectory:
 *        - Create points and add appropriate attributes:\n
 *           - <tt>point = GblPoint(..)</tt>
 *           - <tt>point.addMeasurement(..)</tt>
 *           - <tt>point.addScatterer(..)</tt>
 *           - <tt>point.addLocals(..)</tt>
 *           - <tt>point.addGlobals(..)</tt>
 *        - Add point (ordered by arc length) to trajectory, return label of point:\n
 *            <tt>label = traj.addPoint(point)</tt>
 *    -# Optionally add external seed:\n
 *            <tt>traj.addExternalSeed(..)</tt>
 *    -# Construct and fit trajectory, return error code,
 *       get Chi2, Ndf (and weight lost by M-estimators):\n
 *            <tt>ierr = traj.fit(..)</tt>
 *    -# For any point on initial trajectory:
 *        - Get corrections and covariance matrix for track parameters:\n
 *            <tt>[..] = traj.getResults(label)</tt>
 *    -# Optionally write trajectory to MP binary file:\n
 *            <tt>traj.milleOut(..)</tt>
 *
 *  \section impl_sec Implementation
 *
 *  Matrices are implemented with ROOT (root.cern.ch). User input or output is in the
 *  form of TMatrices. Internally SMatrices are used for fixes sized and simple matrices
 *  based on std::vector<> for variable sized matrices.
 *
 *  \section ref_sec References
 *    - V. Blobel, C. Kleinwort, F. Meier,
 *      Fast alignment of a complex tracking detector using advanced track models,
 *      Computer Phys. Communications (2011), doi:10.1016/j.cpc.2011.03.017
 *    - C. Kleinwort, General Broken Lines as advanced track fitting method,
 *      NIM A, 673 (2012), 107-110, doi:10.1016/j.nima.2012.01.024
 */

#include "GblTrajectory.h"

/// Create new trajectory.
/**
 * Curved trajectory in space (default) or without curvature (q/p) or in one
 * plane (u-direction) only.
 * \param [in] flagCurv Use q/p
 * \param [in] flagU1dir Use in u1 direction
 * \param [in] flagU2dir Use in u2 direction
 */
GblTrajectory::GblTrajectory(bool flagCurv, bool flagU1dir, bool flagU2dir) :
		numPoints(0), numOffsets(0), numCurvature(flagCurv ? 1 : 0), numParameters(
				0), numLocals(0), externalPoint(0), theDimension(0), thePoints(), theData(), externalIndex(), externalSeed() {

	if (flagU1dir)
		theDimension.push_back(0);
	if (flagU2dir)
		theDimension.push_back(1);
	thePoints.reserve(100); // reserve some space for points
}

GblTrajectory::~GblTrajectory() {
}

/// Add point to trajectory.
/**
 * Points have to be ordered in arc length.
 * \param [in] aPoint Point to add
 */
unsigned int GblTrajectory::addPoint(GblPoint aPoint) {
	numPoints++;
	aPoint.setLabel(numPoints);
	thePoints.push_back(aPoint);
	numLocals = std::max(numLocals, aPoint.getNumLocals());
	return numPoints;
}

/// Retrieve number of point from trajectory
unsigned int GblTrajectory::getNumPoints() const {
	return numPoints;
}

/// Add external seed to trajectory.
/**
 * \param [in] aLabel (Signed) label of point for external seed
 * (<0: in front, >0: after point, slope changes at scatterer!)
 * \param [in] aSeed Precision matrix of external seed
 */
void GblTrajectory::addExternalSeed(unsigned int aLabel,
		const TMatrixDSym &aSeed) {
	externalPoint = aLabel;
	externalSeed.ResizeTo(aSeed);
	externalSeed = aSeed;
}

/// Define offsets from list of points.
/**
 * Define offsets at points with scatterers and first and last point.
 * All other points need interpolation from adjacent points with offsets.
 */
void GblTrajectory::defineOffsets() {
//   first point is offset
	thePoints[0].setOffset(0);
	int nOffsets = 1;
//   intermediate scatterers are offsets
	for (unsigned int i = 1; i < numPoints - 1; i++) {
		if (thePoints[i].hasScatterer()) {
			thePoints[i].setOffset(nOffsets);
			nOffsets++;
		} else {
			thePoints[i].setOffset(-nOffsets);
		}
	}
//   last point is offset
	thePoints[numPoints - 1].setOffset(nOffsets);
	numOffsets = nOffsets + 1;
	numParameters = numOffsets * theDimension.size() + numCurvature + numLocals;
}

/// Calculate Jacobians to previous/next scatterer from point to point ones.
void GblTrajectory::calcJacobians() {

	SMatrix55 scatJacobian;
// forward propagation (all)
	unsigned int lastPoint = 0;
	unsigned int numStep = 0;
	for (unsigned int iPoint = 1; iPoint < numPoints; iPoint++) {
		if (numStep == 0) {
			scatJacobian = thePoints[iPoint].getP2pJacobian();
		} else {
			scatJacobian = thePoints[iPoint].getP2pJacobian() * scatJacobian;
		}
		numStep++;
		thePoints[iPoint].addPrevJacobian(scatJacobian); // iPoint -> previous scatterer
		if (thePoints[iPoint].getOffset() >= 0) {
			thePoints[lastPoint].addNextJacobian(scatJacobian); // lastPoint -> next scatterer
			numStep = 0;
			lastPoint = iPoint;
		}
	}
// backward propagation (without scatterers)
	numStep = 0;
	unsigned iPoint = numPoints - 1;
	for (unsigned int i = 1; i < numPoints - 1; i++) {
		iPoint--;
		if (thePoints[iPoint].getOffset() >= 0) {
			numStep = 0;
			continue; // skip offsets
		}
		if (numStep == 0) {
			scatJacobian = thePoints[iPoint].getP2pJacobian();
		} else {
			scatJacobian = scatJacobian * thePoints[iPoint].getP2pJacobian();
		}
		numStep++;
		thePoints[iPoint].addNextJacobian(scatJacobian); // iPoint -> next scatterer
	}
}

/// Get jacobian for transformation from fit to track parameters at point.
/**
 * Jacobian broken lines (q/p,..,u_i,u_i+1..) to track (q/p,u',u) parameters
 * including additional local parameters.
 * \param [in] aSignedLabel (Signed) label of point for external seed
 * (<0: in front, >0: after point, slope changes at scatterer!)
 * \return List of fit parameters with non zero derivatives and
 * corresponding transformation matrix
 */
std::pair<std::vector<unsigned int>, TMatrixD> GblTrajectory::getJacobian(
		int aSignedLabel) const {

	unsigned int nDim = theDimension.size();
	unsigned int nCurv = numCurvature;
	unsigned int nLocals = numLocals;
	unsigned int nBorder = nCurv + nLocals;
	unsigned int nParBRL = nBorder + 2 * nDim;
	unsigned int nParLoc = nLocals + 5;
	std::vector<unsigned int> anIndex;
	anIndex.reserve(nParBRL);
	TMatrixD aJacobian(nParLoc, nParBRL);
	aJacobian.Zero();

	unsigned int aLabel = abs(aSignedLabel);
	int nJacobian; // 0: prev, 1: next
	// check consistency of (index, direction)
	if (aSignedLabel > 0) {
		nJacobian = 1;
		if (aLabel >= numPoints) {
			aLabel = numPoints;
			nJacobian = 0;
		}
	} else {
		nJacobian = 0;
		if (aLabel <= 0) {
			aLabel = 1;
			nJacobian = 1;
		}
	}
	GblPoint aPoint = thePoints[aLabel - 1];
	std::vector<unsigned int> labDer(5);
	SMatrix55 matDer;
	getFitToLocalJacobian(labDer, matDer, aPoint, 5, nJacobian);

	// from local parameters
	for (unsigned int i = 0; i < nLocals; i++) {
		aJacobian(i + 5, i) = 1.0;
		anIndex.push_back(i + 1);
	}
	// from trajectory parameters
	unsigned int iCol = nLocals;
	for (unsigned int i = 0; i < 5; i++) {
		if (labDer[i] > 0) {
			anIndex.push_back(labDer[i]);
			for (unsigned int j = 0; j < 5; j++) {
				aJacobian(j, iCol) = matDer(j, i);
			}
			iCol++;
		}
	}
	return std::make_pair(anIndex, aJacobian);
}

/// Get (part of) jacobian for transformation from (trajectory) fit to track parameters at point.
/**
 * Jacobian broken lines (q/p,..,u_i,u_i+1..) to local (q/p,u',u) parameters.
 * \param [out] anIndex List of fit parameters with non zero derivatives
 * \param [out] aJacobian Corresponding transformation matrix
 * \param [in] aPoint Point to use
 * \param [in] measDim Dimension of 'measurement'
 * (<=2: calculate only offset part, >2: complete matrix)
 * \param [in] nJacobian Direction (0: to previous offset, 1: to next offset)
 */
void GblTrajectory::getFitToLocalJacobian(std::vector<unsigned int> &anIndex,
		SMatrix55 &aJacobian, GblPoint &aPoint, unsigned int measDim,
		unsigned int nJacobian) const {

	unsigned int nDim = theDimension.size();
	unsigned int nCurv = numCurvature;
	unsigned int nLocals = numLocals;

	int nOffset = aPoint.getOffset();

	if (nOffset < 0) // need interpolation
			{
		SMatrix22 prevW, prevWJ, nextW, nextWJ, sumWJ, matN, prevNW, nextNW;
		SVector2 prevWd, nextWd, prevNd, nextNd;
		int ierr;
		aPoint.getDerivatives(0, prevW, prevWJ, prevWd); // W-, W- * J-, W- * d-
		aPoint.getDerivatives(1, nextW, nextWJ, nextWd); // W-, W- * J-, W- * d-
		sumWJ = prevWJ + nextWJ;
		matN = sumWJ.Inverse(ierr); // N = (W- * J- + W+ * J+)^-1
		// derivatives for u_int
		prevNW = matN * prevW; // N * W-
		nextNW = matN * nextW; // N * W+
		prevNd = matN * prevWd; // N * W- * d-
		nextNd = matN * nextWd; // N * W+ * d+

		unsigned int iOff = nDim * (-nOffset - 1) + nLocals + nCurv + 1; // first offset ('i' in u_i)

		// local offset
		if (nCurv > 0) {
			aJacobian.Place_in_col(-prevNd - nextNd, 3, 0); // from curvature
			anIndex[0] = nLocals + 1;
		}
		aJacobian.Place_at(prevNW, 3, 1); // from 1st Offset
		aJacobian.Place_at(nextNW, 3, 3); // from 2nd Offset
		for (unsigned int i = 0; i < nDim; i++) {
			anIndex[1 + theDimension[i]] = iOff + i;
			anIndex[3 + theDimension[i]] = iOff + nDim + i;
		}

		// local slope and curvature
		if (measDim > 2) {
			SMatrix22 prevWPN, nextWPN;
			SVector2 prevWNd, nextWNd;
			// derivatives for u'_int
			prevWPN = nextWJ * prevNW; // W+ * J+ * N * W-
			nextWPN = prevWJ * nextNW; // W- * J- * N * W+
			prevWNd = nextWJ * prevNd; // W+ * J+ * N * W- * d-
			nextWNd = prevWJ * nextNd; // W- * J- * N * W+ * d+
			if (nCurv > 0) {
				aJacobian(0, 0) = 1.0;
				aJacobian.Place_in_col(prevWNd - nextWNd, 1, 0); // from curvature
			}
			aJacobian.Place_at(-prevWPN, 1, 1); // from 1st Offset
			aJacobian.Place_at(nextWPN, 1, 3); // from 2nd Offset
		}
	} else {
		unsigned int iOff = nDim * (nOffset + nJacobian - 1) + nCurv + nLocals
				+ 1; // first offset ('i' in u_i)

		// local offset
		aJacobian(3, 1) = 1.0; // from 1st Offset
		aJacobian(4, 2) = 1.0;
		for (unsigned int i = 0; i < nDim; i++) {
			anIndex[1 + theDimension[i]] = iOff + i;
		}

		// local slope and curvature
		if (measDim > 2) {
			SMatrix22 matW, matWJ;
			SVector2 vecWd;
			aPoint.getDerivatives(nJacobian, matW, matWJ, vecWd); // W, W * J, W * d
			if (nCurv > 0) {
				aJacobian(0, 0) = 1.0;
				aJacobian.Place_in_col(-vecWd, 1, 0); // from curvature
				anIndex[0] = nLocals + 1;
			}
			aJacobian.Place_at(-matWJ, 1, 1); // from 1st Offset
			aJacobian.Place_at(matW, 1, 3); // from 2nd Offset
			for (unsigned int i = 0; i < nDim; i++) {
				anIndex[3 + theDimension[i]] = iOff + nDim + i;
			}
		}
	}
}

/// Get jacobian for transformation from (trajectory) fit to kink parameters at point.
/**
 * Jacobian broken lines (q/p,..,u_i-1,u_i,u_i+1..) to kink (du') parameters.
 * \param [out] anIndex List of fit parameters with non zero derivatives
 * \param [out] aJacobian Corresponding transformation matrix
 * \param [in] aPoint Point to use
 */
void GblTrajectory::getFitToKinkJacobian(std::vector<unsigned int> &anIndex,
		SMatrix27 &aJacobian, GblPoint &aPoint) const {

	unsigned int nDim = theDimension.size();
	unsigned int nCurv = numCurvature;
	unsigned int nLocals = numLocals;

	int nOffset = aPoint.getOffset();

	SMatrix22 prevW, prevWJ, nextW, nextWJ, sumWJ;
	SVector2 prevWd, nextWd, sumWd;
	aPoint.getDerivatives(0, prevW, prevWJ, prevWd); // W-, W- * J-, W- * d-
	aPoint.getDerivatives(1, nextW, nextWJ, nextWd); // W-, W- * J-, W- * d-
	sumWJ = prevWJ + nextWJ; // W- * J- + W+ * J+
	sumWd = prevWd + nextWd; // W+ * d+ + W- * d-

	unsigned int iOff = (nOffset - 1) * nDim + nCurv + nLocals + 1; // first offset ('i' in u_i)

	// local offset
	if (nCurv > 0) {
		aJacobian.Place_in_col(-sumWd, 0, 0); // from curvature
		anIndex[0] = nLocals + 1;
	}
	aJacobian.Place_at(prevW, 0, 1); // from 1st Offset
	aJacobian.Place_at(-sumWJ, 0, 3); // from 2nd Offset
	aJacobian.Place_at(nextW, 0, 5); // from 1st Offset
	for (unsigned int i = 0; i < nDim; i++) {
		anIndex[1 + theDimension[i]] = iOff + i;
		anIndex[3 + theDimension[i]] = iOff + nDim + i;
		anIndex[5 + theDimension[i]] = iOff + nDim * 2 + i;
	}
}

/// Get fit results at point.
/**
 * Get corrections and covariance matrix for local track and additional parameters
 * in forward or backward direction.
 * \param [in] aSignedLabel (Signed) label of point for external seed
 * (<0: in front, >0: after point, slope changes at scatterer!)
 * \param [out] localPar Corrections for local parameters
 * \param [out] localCov Covariance for local parameters
 */
void GblTrajectory::getResults(int aSignedLabel, TVectorD &localPar,
		TMatrixDSym &localCov) const {
	std::pair<std::vector<unsigned int>, TMatrixD> indexAndJacobian =
			getJacobian(aSignedLabel);
	unsigned int nParBrl = indexAndJacobian.first.size();
	TVectorD aVec(nParBrl); // compressed vector
	for (unsigned int i = 0; i < nParBrl; i++) {
		aVec[i] = theVector(indexAndJacobian.first[i] - 1);
	}
	TMatrixDSym aMat = theMatrix.getBlockMatrix(indexAndJacobian.first); // compressed matrix
	localPar = indexAndJacobian.second * aVec;
	localCov = aMat.Similarity(indexAndJacobian.second);
}

/// Build linear equation system from data (blocks).
void GblTrajectory::buildLinearEquationSystem() {
	unsigned int nBorder = numCurvature + numLocals;
	theVector.resize(numParameters);
	theMatrix.resize(numParameters, nBorder);
	double aValue, aWeight;
	std::vector<unsigned int>* indLocal;
	std::vector<double>* derLocal;
	std::vector<GblData>::iterator itData;
	for (itData = theData.begin(); itData < theData.end(); itData++) {
		itData->getLocalData(aValue, aWeight, indLocal, derLocal);
		for (unsigned int j = 0; j < indLocal->size(); j++) {
			theVector((*indLocal)[j] - 1) += (*derLocal)[j] * aWeight * aValue;
		}
		theMatrix.addBlockMatrix(aWeight, indLocal, derLocal);
	}
}

/// Prepare fit
/**
 * Generate data (blocks) from measurements, kinks and external seed.
 */
void GblTrajectory::prepare() {
	unsigned int nDim = theDimension.size();
	unsigned int maxData = 2 * (numPoints + numOffsets - 2); // upper limit
	theData.reserve(maxData);
	// measurements
	SMatrix55 matP;
	std::vector<GblPoint>::iterator itPoint;
	for (itPoint = thePoints.begin(); itPoint < thePoints.end(); itPoint++) {
		SVector5 aMeas, aPrec;
		unsigned int measDim = itPoint->hasMeasurement();
		if (measDim) {
			unsigned int nLabel = itPoint->getLabel();
			TMatrixD localDer = itPoint->getLocalDerivatives();
			std::vector<int> globalLab = itPoint->getGlobalLabels();
			TMatrixD globalDer = itPoint->getGlobalDerivatives();
			itPoint->getMeasurement(matP, aMeas, aPrec);
			unsigned int iOff = 5 - measDim; // first active component
			std::vector<unsigned int> labDer(5);
			SMatrix55 matDer, matPDer;
			getFitToLocalJacobian(labDer, matDer, *itPoint, measDim);
			if (measDim > 2) {
				matPDer = matP * matDer;
			} else { // 'shortcut' for position measurements
				matPDer.Place_at(
						matP.Sub<SMatrix22>(3, 3) * matDer.Sub<SMatrix25>(3, 0),
						3, 0);
			}

			for (unsigned int i = iOff; i < 5; i++) {
				if (aPrec(i) > 0.) {
					GblData aData(nLabel, aMeas(i), aPrec(i));
					aData.addDerivatives(i, labDer, matPDer, iOff, localDer,
							globalLab, globalDer);
					theData.push_back(aData);
				}
			}

		}
	}
	// pseudo measurements from kinks
	for (itPoint = thePoints.begin() + 1; itPoint < thePoints.end() - 1;
			itPoint++) {
		SVector2 aMeas, aPrec;
		if (itPoint->hasScatterer()) {
			unsigned int nLabel = itPoint->getLabel();
			itPoint->getScatterer(aMeas, aPrec);
			std::vector<unsigned int> labDer(7);
			SMatrix27 matDer;
			getFitToKinkJacobian(labDer, matDer, *itPoint);
			for (unsigned int i = 0; i < nDim; i++) {
				unsigned int iDim = theDimension[i];
				if (aPrec(iDim) > 0.) {
					GblData aData(nLabel, aMeas(iDim), aPrec(iDim));
					aData.addDerivatives(iDim, labDer, matDer);
					theData.push_back(aData);
				}
			}
		}
	}
	// external seed
	if (externalPoint > 0) {
		std::pair<std::vector<unsigned int>, TMatrixD> indexAndJacobian =
				getJacobian(externalPoint);
		externalIndex = indexAndJacobian.first;
		std::vector<double> externalDerivatives(externalIndex.size());
		TMatrixDSymEigen externalEigen(externalSeed);
		TVectorD valEigen = externalEigen.GetEigenValues();
		TMatrixD vecEigen = externalEigen.GetEigenVectors();
		vecEigen = vecEigen.T() * indexAndJacobian.second;
		for (int i = 0; i < externalSeed.GetNrows(); i++) {
			if (valEigen(i) > 0.) {
				for (int j = 0; j < externalSeed.GetNcols(); j++) {
					externalDerivatives[j] = vecEigen(i, j);
				}
				GblData aData(externalPoint, 0., valEigen(i));
				aData.addDerivatives(externalIndex, externalDerivatives);
				theData.push_back(aData);
			}
		}
	}
}

/// Calculate predictions for all points.
void GblTrajectory::predict() {
	std::vector<GblData>::iterator itData;
	for (itData = theData.begin(); itData < theData.end(); itData++) {
		itData->setPrediction(theVector);
	}
}

/// Down-weight all points.
/**
 * \param [in] aMethod M-estimator (1: Tukey, 2:Huber, 3:Cauchy)
 */
double GblTrajectory::downWeight(unsigned int aMethod) {
	double aLoss = 0.;
	std::vector<GblData>::iterator itData;
	for (itData = theData.begin(); itData < theData.end(); itData++) {
		aLoss += (1. - itData->setDownWeighting(aMethod));
	}
	return aLoss;
}

/// Perform fit of trajectory.
/**
 * Optionally iterate for outlier down-weighting.
 * \param [out] Chi2 Chi2 sum (corrected for down-weighting)
 * \param [out] Ndf  Number of degrees of freedom
 * \param [out] lostWeight Sum of weights lost due to down-weighting
 * \param [in] optionList Iterations for down-weighting
 * (One character per iteration: t,h,c (or T,H,C) for Tukey, Huber or Cauchy function)
 * \return Error code (non zero value indicates failure of fit)
 */
unsigned int GblTrajectory::fit(double &Chi2, int &Ndf, double &lostWeight,
		std::string optionList) {
	const double normChi2[4] = { 1.0, 0.8737, 0.9326, 0.8228 };
	const std::string methodList = "TtHhCc";

	unsigned int aMethod = 0;

	defineOffsets();
	calcJacobians();
	prepare();
	buildLinearEquationSystem();
	lostWeight = 0.;
	unsigned int ierr = 0;
	try {

		theMatrix.solveAndInvertBorderedBand(theVector, theVector);
		predict();

		for (unsigned int i = 0; i < optionList.size(); i++) // down weighting iterations
				{
			size_t aPosition = methodList.find(optionList[i]);
			if (aPosition != std::string::npos) {
				aMethod = aPosition / 2 + 1;
				lostWeight = downWeight(aMethod);
				buildLinearEquationSystem();
				theMatrix.solveAndInvertBorderedBand(theVector, theVector);
				predict();
			}
		}
		Ndf = theData.size() - numParameters;
		Chi2 = 0.;
		for (unsigned int i = 0; i < theData.size(); i++) {
			Chi2 += theData[i].getChi2();
		}
		Chi2 /= normChi2[aMethod];

	} catch (int e) {
		std::cout << " GblTrajectory::fit exception " << e << std::endl;
		Chi2 = 0.;
		Ndf = -1;
		lostWeight = 0.;
		ierr = e;
	}
	return ierr;
}

/// Write trajectory to Millepede-II binary file.
void GblTrajectory::milleOut(MilleBinary &aMille) {
	float fValue;
	float fErr;
	std::vector<unsigned int>* indLocal;
	std::vector<double>* derLocal;
	std::vector<int>* labGlobal;
	std::vector<double>* derGlobal;

//   data: measurements, kinks and external seed
	std::vector<GblData>::iterator itData;
	for (itData = theData.begin(); itData < theData.end(); itData++) {
		itData->getAllData(fValue, fErr, indLocal, derLocal, labGlobal,
				derGlobal);
		aMille.addData(fValue, fErr, *indLocal, *derLocal, *labGlobal,
				*derGlobal);
	}
	aMille.writeRecord();
}

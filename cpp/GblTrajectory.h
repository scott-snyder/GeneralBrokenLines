/*
 * GblTrajectory.h
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

#ifndef GBLTRAJECTORY_H_
#define GBLTRAJECTORY_H_

#include "GblPoint.h"
#include "GblData.h"
#include "GblPoint.h"
#include "BorderedBandMatrix.h"
#include "MilleBinary.h"
#include "TMatrixDSymEigen.h"

/// GBL trajectory.
/**
 * List of GblPoints ordered by arc length.
 * Can be fitted and optionally written to MP-II binary file.
 */
class GblTrajectory {
public:
	GblTrajectory(unsigned int nReserve = 100, bool flagCurv = true,
			bool flagU1dir = true, bool flagU2dir = true);
	virtual ~GblTrajectory();
	unsigned int addPoint(GblPoint aPoint);
	unsigned int getNumPoints() const;
	GblPoint* getPoint(unsigned int aLabel);
	void addExternalSeed(unsigned int aLabel, const TMatrixDSym &aSeed);
	unsigned int getResults(int aSignedLabel, TVectorD &localPar,
			TMatrixDSym &localCov) const;
	unsigned int getMeasResults(unsigned int aLabel, unsigned int &numRes,
			TVectorD &aResiduals, TVectorD &aMeasErrors, TVectorD &aResErrors,
			TVectorD &aDownWeights);
	unsigned int getScatResults(unsigned int aLabel, unsigned int &numRes,
			TVectorD &aResiduals, TVectorD &aMeasErrors, TVectorD &aResErrors,
			TVectorD &aDownWeights);
	unsigned int fit(double &Chi2, int &Ndf, double &lostWeight,
			std::string optionList = "");
	void milleOut(MilleBinary &aMille);

private:
	unsigned int numPoints; ///< Number of point on trajectory
	unsigned int numOffsets; ///< Number of (points with) offsets on trajectory
	unsigned int numCurvature; ///< Number of curvature parameters (0 or 1)
	unsigned int numParameters; ///< Number of fit parameters
	unsigned int numLocals; ///< Total number of (additional) local parameters
	unsigned int externalPoint; ///< Label of external point (or 0)
	bool dataGenerated; ///< GblData prepared for fit or milleOut
	bool fitOK; ///< Trajectory has been successfully fitted (results are valid)
	std::vector<unsigned int> theDimension; ///< List of active dimensions (0=u1, 1=u2) in fit
	std::vector<GblPoint> thePoints; ///< List of points on trajectory
	std::vector<GblData> theData; ///< List of data blocks
	std::vector<unsigned int> measDataIndex; ///< mapping points to data blocks from measurements
	std::vector<unsigned int> scatDataIndex; ///< mapping points to data blocks from scatterers
	std::vector<unsigned int> externalIndex; ///< List of fit parameters used by external seed
	TMatrixDSym externalSeed; ///< Precision (inverse covariance matrix) of external seed
	VVector theVector; ///< Vector of linear equation system
	BorderedBandMatrix theMatrix; ///< (Bordered band) matrix of linear equation system

	std::pair<std::vector<unsigned int>, TMatrixD> getJacobian(
			int aSignedLabel) const;
	void getFitToLocalJacobian(std::vector<unsigned int> &anIndex,
			SMatrix55 &aJacobian, const GblPoint &aPoint, unsigned int measDim,
			unsigned int nJacobian = 1) const;
	void getFitToKinkJacobian(std::vector<unsigned int> &anIndex,
			SMatrix27 &aJacobian, const GblPoint &aPoint) const;
	void defineOffsets();
	void calcJacobians();
	void buildLinearEquationSystem();
	void prepare();
	void predict();
	double downWeight(unsigned int aMethod);
	void getResAndErr(unsigned int aData, double &aResidual,
			double &aMeadsError, double &aResError, double &aDownWeight);
};

#endif /* GBLTRAJECTORY_H_ */

/*
 * example4.cpp
 *
 *  Created on: Dec 16, 2022
 *      Author: kleinwrt
 */

/** \file
 *  Example application (with thick scatterers).
 *
 *  \author Claus Kleinwort, DESY, 2022 (Claus.Kleinwort@desy.de)
 *
 *  \copyright
 *  Copyright (c) 2022 - 2023 Deutsches Elektronen-Synchroton,
 *  Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
 *  This library is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version. \n\n
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details. \n\n
 *  You should have received a copy of the GNU Library General Public
 *  License along with this program (see the file COPYING.LIB for more
 *  details); if not, write to the Free Software Foundation, Inc.,
 *  675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <time.h>
#include "exampleUtil.h"
#include "GblTrajectory.h"

using namespace gbl;
using namespace Eigen;

void example4() {
	/// Simple technical example (curvilinear as local system, demonstrating **thick** GBL scatterers).
	/**
	 * Create points on initial trajectory, create trajectory from points,
	 * fit and write trajectory to MP-II binary file,
	 * get track parameter corrections and covariance matrix at points.
	 *
	 * Equidistant measurement layers and multiple thin scatterers, propagation
	 * with simple jacobian (quadratic in arc length differences).
	 * Curvilinear system (U,V,T) as local coordinate system.
	 *
	 * This example simulates and refits tracks in a system of planar detectors
	 * with 2D measurements in a constant magnet field in Z direction using
	 * the curvilinear system as local system and (Q/P, slopes, offsets) as
	 * local track parameters. The true track parameters are
	 * randomly smeared with respect to a (constant and straight) reference
	 * trajectory with direction (lambda, phi) and are used (only) for the
	 * on-the-fly simulation of the measurements and scatterers. The predictions
	 * from the reference trajectory are therefore always zero and the residuals
	 * needed (by addMeasurement) are equal to the measurements.
	 */

//MP	MilleBinary mille; // for producing MillePede-II binary file
	unsigned int nTry = 1000; //: number of tries
	unsigned int nLayer = 10; //: number of detector layers
	bool useThickScatterer = true; //: use thick (GBL) scatterers at measurements
	std::cout << " Gbltst-thickScat $Id$ " << nTry << ", " << nLayer
			<< ", useThickScat: " << useThickScatterer << std::endl;

	srand(4711);

	clock_t startTime = clock();
// track direction
	double sinLambda = 0.3;
	double cosLambda = sqrt(1.0 - sinLambda * sinLambda);
	double sinPhi = 0.;
	double cosPhi = sqrt(1.0 - sinPhi * sinPhi);
// tDir = (cosLambda * cosPhi, cosLambda * sinPhi, sinLambda)
// U = Z x T / |Z x T|, V = T x U
	Matrix<double, 2, 3> uvDir;
	uvDir(0, 0) = -sinPhi;
	uvDir(0, 1) = cosPhi;
	uvDir(0, 2) = 0.;
	uvDir(1, 0) = -sinLambda * cosPhi;
	uvDir(1, 1) = -sinLambda * sinPhi;
	uvDir(1, 2) = cosLambda;
// measurement resolution
	Vector2d measErr;
	measErr << 0.001, 0.001;
	Vector2d measPrec; // (independent) precisions
	measPrec << 1.0 / (measErr(0) * measErr(0)), 1.0 / (measErr(1) * measErr(1));
// scattering error
	Vector2d scatErr;
	scatErr << 0.001, 0.001;
	Vector2d scatPrec;
	Vector2d scat(0., 0.);
	scatPrec << 1.0 / (scatErr(0) * scatErr(0)), 1.0 / (scatErr(1) * scatErr(1));
	Matrix4d scatPrecThick;
	Vector4d scatThick(0., 0., 0., 0.);
// (RMS of) CurviLinear track parameters (Q/P, slopes, offsets)
	Vector5d clPar;
	Vector5d clErr;
	clErr << 0.001, -0.1, 0.2, -0.15, 0.25;
//  scattering cov matrix for thick scatterers
	Matrix4d scatCov;
//  number of thin scatterers in thick scatterer
	unsigned int numScat;
// additional parameters
	Vector2d addPar;
	addPar << 0.0025, -0.005;
	std::vector<int> globalLabels;
	globalLabels.push_back(4711);
	globalLabels.push_back(4712);
// global labels for MP
	/*MP	std::vector<int> globalLabels(2);
	 globalLabels[0] = 11;
	 globalLabels[1] = 12; */

	double bfac = 0.2998; // Bz*c for Bz=1
	double step = 1.5 / cosLambda; // constant steps in RPhi

	double Chi2Sum = 0.;
	int NdfSum = 0;
	double LostSum = 0.;
	int numFit = 0;

	for (unsigned int iTry = 1; iTry <= nTry; ++iTry) {
		// curvilinear track parameters
		for (unsigned int i = 0; i < 5; ++i) {
			clPar[i] = clErr[i] * unrm();
		}
//		std::cout << " Try " << iTry << ":" << clPar << std::endl;
		Matrix2d addDer;
		addDer.setZero();
		addDer(0, 0) = 1.;
		addDer(1, 1) = 1.;
// arclength
		//double s = 0.;
		Matrix5d jacPointToPoint, jacToNextPlane;
		jacPointToPoint.setIdentity();
		scatCov.setZero();
		numScat = 0;
// create list of points
		std::vector<GblPoint> listOfPoints;
		listOfPoints.reserve(2 * nLayer);

		for (unsigned int iLayer = 0; iLayer < nLayer; ++iLayer) {
//			std::cout << " Layer " << iLayer << ", " << s << std::endl;
//     measurement directions
			double sinStereo = (iLayer % 2 == 0) ? 0. : 0.1;
			double cosStereo = sqrt(1.0 - sinStereo * sinStereo);
			Matrix<double, 3, 2> mDirT;
			mDirT.setZero();
			mDirT(1, 0) = cosStereo;
			mDirT(2, 0) = sinStereo;
			mDirT(1, 1) = -sinStereo;
			mDirT(2, 1) = cosStereo;
// projection measurement to local (curvilinear uv) directions (duv/dm)
			Matrix2d proM2l = uvDir * mDirT;
// projection local (uv) to measurement directions (dm/duv)
			Matrix2d proL2m = proM2l.inverse();
			// point with (independent) measurements (in measurement system)
			GblPoint pointMeas(jacPointToPoint);
			// measurement - prediction in measurement system with error
			Vector2d meas = proL2m * clPar.tail(2);
			//MP			meas += addDer * addPar; // additional parameters
			for (unsigned int i = 0; i < 2; ++i) {
				meas[i] += measErr[i] * unrm();
			}
			pointMeas.addMeasurement(proL2m, meas, measPrec);

			// additional local parameters?
//			point.addLocals(addDer);
//MP			point.addGlobals(globalLabels, addDer);
			addDer *= -1.; // Der flips sign every measurement
// describe scattering with thick scatterer ?
			if (useThickScatterer) {
				if (numScat) {
					if (numScat == 1) {
						std::cout
								<< " Singular cov. matrix for thick scatterer at layer "
								<< iLayer << std::endl;
						std::cout << scatCov << std::endl;
						return;
					}
					scatPrecThick = scatCov.inverse();
					pointMeas.addThickScatterer(scatThick, scatPrecThick);
					scatCov.setZero();
					numScat = 0;
				}
			} else {
// describe scattering with two thin scatterers
				pointMeas.addScatterer(scat, scatPrec);
			}
// scatter at measurement too
			// scatter a little
			numScat++;
			for (unsigned int i = 0; i < 2; ++i) {
				clPar[i + 1] += scatErr[i] * unrm();
				scatCov(i, i) += scatErr[i] * scatErr[i];
			}
// add point to trajectory
			listOfPoints.push_back(pointMeas);
// no scattering after last measurement
			if (iLayer == nLayer - 1)
				break;
// propagate to scattering plane
			jacToNextPlane = gblSimpleJacobian(step, cosLambda, bfac);
			jacPointToPoint = jacToNextPlane;
			//jac2 = gblSimpleJacobian2(step, cosLambda, bfac);
			clPar = jacToNextPlane * clPar;
			scatCov = jacToNextPlane.bottomRightCorner<4, 4>() * scatCov
					* jacToNextPlane.bottomRightCorner<4, 4>().transpose();
			//s += step;

			if (useThickScatterer) {
// no intermediate GblPoint, update jacobian
				jacPointToPoint = jacToNextPlane * jacPointToPoint;
			} else {
// describe scattering with two thin scatterers, requires intermediate GblPoint
				GblPoint pointScat(jacToNextPlane);
				pointScat.addScatterer(scat, scatPrec);
				listOfPoints.push_back(pointScat);
			}
			// scatter a little
			numScat++;
			for (unsigned int i = 0; i < 2; ++i) {
				clPar[i + 1] += scatErr[i] * unrm();
				scatCov(i, i) += scatErr[i] * scatErr[i];
			}
// propagate to next measurement layer
			clPar = jacToNextPlane * clPar;
			scatCov = jacToNextPlane.bottomRightCorner<4, 4>() * scatCov
					* jacToNextPlane.bottomRightCorner<4, 4>().transpose();
			//s += step;
		}
//
		// create trajectory
		GblTrajectory traj(listOfPoints);
		//GblTrajectory traj(listOfPoints, seedLabel, clSeed); // with external seed
		//traj.printPoints();
		/*
		 if (not traj.isValid()) {
		 std::cout << " Invalid GblTrajectory -> skip" << std::endl;
		 continue;
		 }*/
// fit trajectory
		double Chi2;
		int Ndf;
		double lostWeight;
		traj.fit(Chi2, Ndf, lostWeight);
		//std::cout << " Fit: " << Chi2 << ", " << Ndf << ", " << lostWeight << std::endl;
		/* look at (track parameter) corrections
		 VectorXd aCorrection(5);
		 MatrixXd aCovariance(5, 5);
		 traj.getResults(1, aCorrection, aCovariance);
		 std::cout << " cor " << std::endl << aCorrection << std::endl;
		 std::cout << " cov " << std::endl << aCovariance << std::endl;
		 */
		/* look at residuals
		 for (unsigned int label = 1; label <= listOfPoints.size(); ++label) {
		 unsigned int numData = 0;
		 std::cout << " measResults, label " << label << std::endl;
		 VectorXd residuals(2), measErr(2), resErr(2), downWeights(2);
		 traj.getMeasResults(label, numData, residuals, measErr, resErr,
		 downWeights);
		 std::cout << " measResults, numData " << numData << std::endl;
		 for (unsigned int i = 0; i < numData; ++i) {
		 std::cout << " measResults " << label << " " << i << " "
		 << residuals[i] << " " << measErr[i] << " " << resErr[i]
		 << std::endl;
		 }
		 } */
// debug printout
		//traj.printTrajectory(1);
		//traj.printPoints();
		//traj.printData();
// write to MP binary file
//MP		traj.milleOut(mille);
		Chi2Sum += Chi2;
		NdfSum += Ndf;
		LostSum += lostWeight;
		numFit++;
	}

	clock_t endTime = clock();
	double diff = endTime - startTime;
	double cps = CLOCKS_PER_SEC;
	std::cout << " Time elapsed " << diff / cps << " s" << std::endl;
	std::cout << " Chi2/Ndf = " << Chi2Sum / NdfSum << std::endl;
	std::cout << " Tracks fitted " << numFit << std::endl;
	if (LostSum > 0.)
		std::cout << " Weight lost   " << LostSum << std::endl;
}


/*
 * exampleSit.cpp
 *
 *  Created on: 11 Oct 2018
 *      Author: kleinwrt
 */

/** \file
 *  Example silicon tracker application.
 *
 *  \author Claus Kleinwort, DESY, 2018 (Claus.Kleinwort@desy.de)
 *
 *  \copyright
 *  Copyright (c) 2018 Deutsches Elektronen-Synchroton,
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
#include "exampleSit.h"
#include "GblTrajectory.h"

using namespace gbl;
using namespace Eigen;

/// Silicon tracker example
/**
 * Simulate and reconstruct helical tracks in silicon pixel and (1D or 2D) strip detectors.
 *
 *  Create points on initial trajectory, create trajectory from points,
 *  fit and write trajectory to MP-II binary file (for rigid body alignment).
 *
 *  Setup:
 *   - Beam (mainly) in X direction
 *   - Constant magnetic field in Z direction
 *   - Silicon sensors measuring in YZ plane, orthogonal (pixel) or non-orthogonal (stereo strips) measurement systems
 *   - Multiple scattering in sensors (air in between ignored)
 *   - Curvilinear system (T,U,V) as local coordinate system and (q/p, slopes, offsets) as local track parameters
 *
 * \remark To exercise (mis)alignment different sets of layers (with different geometry)
 * for simulation and reconstruction can be used.
 *
 * Example steering file for Millepede-II (B=0):
 * \code{.unparsed}
 * Cfiles
 * milleBinaryISN.dat
 *
 * method inversion 3 0.1
 * chiscut 30. 6.
 * printcounts
 * ! fix first pixel and last stereo layer as reference
 * parameter
 *   1  0.  -1.
 *   2  0.  -1.
 *   3  0.  -1.
 *   4  0.  -1.
 *   5  0.  -1.
 *   6  0.  -1.
 *  61  0.  -1.
 *  62  0.  -1.
 *  63  0.  -1.
 *  64  0.  -1.
 *  65  0.  -1.
 *  66  0.  -1.
 * end
 * \endcode
 */
void exampleSit() {

	// detector layers (ordered in X):
	// name, position (x,y,z), thickness (X/X_0), (1 or 2) measurements (direction in YZ, resolution)
	std::vector<GblDetectorLayer> layers;
	layers.push_back(
			CreateLayerSit("PIX1", 0, 2.0, 0., 0., 0.0033, 0., 0.0010, 90.,
					0.0020)); // pixel
	layers.push_back(
			CreateLayerSit("PIX2", 1, 3.0, 0., 0., 0.0033, 0., 0.0010, 90.,
					0.0020)); // pixel
	layers.push_back(
			CreateLayerSit("PIX3", 2, 4.0, 0., 0., 0.0033, 0., 0.0010, 90.,
					0.0020)); // pixel
	layers.push_back(
			CreateLayerSit("S2D4", 3, 6.0, 0., 0., 0.0033, 0., 0.0025, 5.0,
					0.0025)); // strip 2D, +5 deg stereo angle
	layers.push_back(
			CreateLayerSit("S2D5", 4, 8.0, 0., 0., 0.0033, 0., 0.0025, -5.,
					0.0025)); // strip 2D, -5 deg stereo angle
	layers.push_back(
			CreateLayerSit("S2D6", 5, 10., 0., 0., 0.0033, 0., 0.0025, 5.0,
					0.0025)); // strip 2D, +5 deg stereo angle
	layers.push_back(
			CreateLayerSit("S2D7", 6, 12., 0., 0., 0.0033, 0., 0.0025, -5.,
					0.0025)); // strip 2D, -5 deg stereo angle
	layers.push_back(
			CreateLayerSit("S1D8", 7, 15., 0., 0., 0.0033, 0., 0.0040)); // strip 1D, no sensitivity to Z

	/* print layers
	 for (unsigned int iLayer = 0; iLayer < layers.size(); ++iLayer) {
	 layers[iLayer].print();
	 } */

	unsigned int nTry = 10000; //: number of tries
	std::cout << " GblSit $Rev$ " << nTry << ", " << layers.size()
			<< std::endl;
	srand(4711);
	clock_t startTime = clock();

	double qbyp = 0.2; // 5 GeV
	const double bfac = 0.003;  // B*c for 1 T
	// const double bfac = 0.;  // B*c for 0 T

	MilleBinary mille; // for producing MillePede-II binary file

	double Chi2Sum = 0.;
	int NdfSum = 0;
	double LostSum = 0.;
	int numFit = 0;

	for (unsigned int iTry = 0; iTry < nTry; ++iTry) {

		// helix parameter for track generation
		const double genDca = 0.1 * unrm(); // normal
		const double genZ0 = 0.1 * unrm(); // normal
		const double genPhi0 = 0.2 * (2. * unif() - 1.); // uniform
		const double genDzds = 0.3 * (2. * unif() - 1.); // uniform
		const double genCurv = bfac * qbyp * sqrt(1. + genDzds * genDzds);

		//
		// generate hits
		//
		std::vector<Vector2d> hits;
		double curv(genCurv), phi0(genPhi0), dca(genDca), dzds(genDzds), z0(
				genZ0);
		const double cosLambda = 1. / sqrt(1. + dzds * dzds);
		for (unsigned int iLayer = 0; iLayer < layers.size(); ++iLayer) {
			// local constant (Bfield) helix
			GblSimpleHelix hlx = GblSimpleHelix(curv, phi0, dca, dzds, z0);
			// prediction from local helix
			GblHelixPrediction pred = layers[iLayer].intersectWithHelix(hlx);
			// std::cout << " layer " << iLayer << " arc-length " << pred.getArcLength() << std::endl;
			Vector2d meas = pred.getMeasPred();
			// smear according to resolution
			Vector2d sigma = layers[iLayer].getResolution();
			meas[0] += sigma[0] * unrm();
			meas[1] += sigma[1] * unrm();
			// save hit
			hits.push_back(meas);
			// scatter at intersection point
			Vector3d measPos = pred.getPosition();
			double radlen = layers[iLayer].getRadiationLength()
					/ fabs(pred.getCosIncidence());
			double errMs = gblMultipleScatteringError(qbyp, radlen); // simple model
			// move to intersection point
			hlx.moveToXY(measPos[0], measPos[1], phi0, dca, z0); // update phi0, dca, z0
			phi0 += unrm() * errMs / cosLambda; // scattering for phi
			dzds += unrm() * errMs / (cosLambda * cosLambda); // scattering for dzds
			GblSimpleHelix newhlx = GblSimpleHelix(curv, phi0, dca, dzds, z0); // after scattering
			// move back
			newhlx.moveToXY(-measPos[0], -measPos[1], phi0, dca, z0); // update phi0, dca, z0
		}

		//
		// fit track with GBL
		//
		// seed (with true parameters)
		double seedCurv(genCurv), seedPhi0(genPhi0), seedDca(genDca), seedDzds(
				genDzds), seedZ0(genZ0);
		GblSimpleHelix seed = GblSimpleHelix(seedCurv, seedPhi0, seedDca,
				seedDzds, seedZ0);
		// (previous) arc-length
		double sOld = 0.;
		const double cosLambdaSeed = 1. / sqrt(1. + (seedDzds * seedDzds));
		// list of points on trajectory
		std::vector<GblPoint> listOfPoints;
		for (unsigned int iLayer = 0; iLayer < layers.size(); ++iLayer) {
			// std::cout << " hit " << iLayer << " " << hits[iLayer].transpose() << std::endl;
			GblDetectorLayer& layer = layers[iLayer];
			// prediction from seeding helix
			GblHelixPrediction pred = layer.intersectWithHelix(seed);
			double sArc = pred.getArcLength(); // arc-length
			Vector2d measPrediction = pred.getMeasPred(); // measurement prediction
			Vector2d measPrecision = layer.getPrecision(); // measurement precision
			// residuals
			Vector2d res(hits[iLayer][0] - measPrediction[0],
					hits[iLayer][1] - measPrediction[1]);
			// transformation global system to local (curvilinear) (u,v) (matrix from row vectors)
			Matrix<double, 2, 3> transG2l = pred.getCurvilinearDirs();
			// transformation measurement system to global system
			Matrix3d transM2g = layer.getMeasSystemDirs().inverse();
			// projection matrix (measurement plane to local (u,v))
			Matrix2d proM2l = transG2l * transM2g.block<3, 2>(0, 0); // skip measurement normal
			// projection matrix (local (u,v) to measurement plane)
			Matrix2d proL2m = proM2l.inverse();
			// propagation
			Matrix5d jacPointToPoint = gblSimpleJacobian(
					(sArc - sOld) / cosLambdaSeed, cosLambdaSeed, bfac);
			sOld = sArc;
			// point with (independent) measurements (in measurement system)
			GblPoint point(jacPointToPoint);
			point.addMeasurement(proL2m, res, measPrecision);
			// global labels and parameters for rigid body alignment
			std::vector<int> labGlobal(6);
			unsigned int layerID = layer.getLayerID();
			for (int p = 0; p < 6; p++)
				labGlobal[p] = layerID * 10 + p + 1;
			Vector3d pos = pred.getPosition();
			Vector3d dir = pred.getDirection();
			Matrix<double, 2, 6> derGlobal = layer.getRigidBodyDerLocal(pos,
					dir);
			point.addGlobals(labGlobal, derGlobal);
			// add scatterer to point
			double radlen = layer.getRadiationLength()
					/ fabs(pred.getCosIncidence());
			double errMs = gblMultipleScatteringError(qbyp, radlen); // simple model
			if (errMs > 0.) {
				Vector2d scat(0., 0.);
				Vector2d scatPrec(1. / (errMs * errMs), 1. / (errMs * errMs)); // scattering precision matrix is diagonal in curvilinear system
				point.addScatterer(scat, scatPrec);
			}
			// add point to trajectory
			listOfPoints.push_back(point);
		}
		// create trajectory
		GblTrajectory traj(listOfPoints, bfac != 0.);
		// fit trajectory
		double Chi2;
		int Ndf;
		double lostWeight;
		unsigned int ierr = traj.fit(Chi2, Ndf, lostWeight);
		// std::cout << " Fit " << iTry << ": "<< Chi2 << ", " << Ndf << ", " << lostWeight << std::endl;
		// successfully fitted?
		if (!ierr) {
			// write to MP binary file
			traj.milleOut(mille);
			// update statistics
			Chi2Sum += Chi2;
			NdfSum += Ndf;
			LostSum += lostWeight;
			numFit++;
		}
	}
	clock_t endTime = clock();
	double diff = endTime - startTime;
	double cps = CLOCKS_PER_SEC;
	std::cout << " Time elapsed " << diff / cps << " s" << std::endl;
	std::cout << " Chi2/Ndf = " << Chi2Sum / NdfSum << std::endl;
	std::cout << " Tracks fitted " << numFit << std::endl;
}

namespace gbl {

/// Create a silicon layer with 1D measurement.
/**
 * Create silicon layer with 1D measurement (u) at fixed X-position.
 *
 * \param [in] aName      name
 * \param [in] layer      layer ID
 * \param [in] xPos       X-position (of center)
 * \param [in] yPos       Y-position (of center)
 * \param [in] zPos       Z-position (of center)
 * \param [in] thickness  thickness / radiation_length
 * \param [in] uAngle     angle of u-direction in YZ plane
 * \param [in] uRes       resolution in u-direction
 */
GblDetectorLayer CreateLayerSit(const std::string aName, unsigned int layer,
		double xPos, double yPos, double zPos, double thickness, double uAngle,
		double uRes) {
	Vector3d aCenter(xPos, yPos, zPos);
	Vector2d aResolution(uRes, 0.);
	Vector2d aPrecision(1. / (uRes * uRes), 0.);
	Matrix3d measTrafo;
	const double cosU = cos(uAngle / 180. * M_PI);
	const double sinU = sin(uAngle / 180. * M_PI);
	measTrafo << 0., cosU, sinU, 0., -sinU, cosU, 1., 0., 0.; // U,V,N
	Matrix3d alignTrafo;
	alignTrafo << 0., 1., 0., 0., 0., 1., 1., 0., 0.; // Y,Z,X
	return GblDetectorLayer(aName, layer, 1, thickness, aCenter, aResolution,
			aPrecision, measTrafo, alignTrafo);
}

/// Create a silicon layer with 2D measurement.
/**
 * Create silicon layer with 2D measurement (u,v) at fixed X-position.
 * The measurement directions in the YZ plane can be orthogonal or non-orthogonal
 * (but must be different).
 *
 * \param [in] aName      name
 * \param [in] layer      layer ID
 * \param [in] xPos       X-position (of center)
 * \param [in] yPos       Y-position (of center)
 * \param [in] zPos       Z-position (of center)
 * \param [in] thickness  thickness / radiation_length
 * \param [in] uAngle     angle of u-direction in YZ plane
 * \param [in] uRes       resolution in u-direction
 * \param [in] vAngle     angle of v-direction in YZ plane
 * \param [in] vRes       resolution in v-direction
 */
GblDetectorLayer CreateLayerSit(const std::string aName, unsigned int layer,
		double xPos, double yPos, double zPos, double thickness, double uAngle,
		double uRes, double vAngle, double vRes) {
	Vector3d aCenter(xPos, yPos, zPos);
	Vector2d aResolution(uRes, vRes);
	Vector2d aPrecision(1. / (uRes * uRes), 1. / (vRes * vRes));
	Matrix3d measTrafo;
	const double cosU = cos(uAngle / 180. * M_PI);
	const double sinU = sin(uAngle / 180. * M_PI);
	const double cosV = cos(vAngle / 180. * M_PI);
	const double sinV = sin(vAngle / 180. * M_PI);
	measTrafo << 0., cosU, sinU, 0., cosV, sinV, 1., 0., 0.; // U,V,N
	Matrix3d alignTrafo;
	alignTrafo << 0., 1., 0., 0., 0., 1., 1., 0., 0.; // Y,Z,X
	return GblDetectorLayer(aName, layer, 2, thickness, aCenter, aResolution,
			aPrecision, measTrafo, alignTrafo);
}

}

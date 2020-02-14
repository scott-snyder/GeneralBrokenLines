/*
 * exampleDc.cpp
 *
 *  Created on: 6 Nov 2018
 *      Author: kleinwrt
 */

/** \file
 *  Example drift chamber application.
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
#include "exampleDc.h"
#include "GblTrajectory.h"

using namespace gbl;
using namespace Eigen;

/// Drift chamber example
/**
 * Simulate and reconstruct helical tracks in a sector of (forward) drift chambers.
 *
 *  Create points on initial trajectory, create trajectory from points,
 *  fit and write trajectory to MP-II binary file (for rigid body alignment).
 *
 *  Setup:
 *   - Beam forward (+Z) direction
 *   - Constant magnetic field in Z direction
 *   - Planar drift chambers with normal in XZ plane, center at Y=0, 1D measurement from (azimuthal) wires, cell size 2 cm.
 *   - Multiple scattering in detectors (gas, wires, walls) (air in between ignored)
 *   - Curvilinear system (T,U,V) as local coordinate system and (q/p, slopes, offsets) as local track parameters
 *
 * \remark To exercise (mis)alignment different sets of layers (with different geometry)
 * for simulation and reconstruction can be used.
 *
 * Example steering file for Millepede-II (B=0, chamber alignment):
 * \code{.unparsed}
 * Cfiles
 * milleBinaryISN.dat
 *
 * method inversion 3 0.1
 * chiscut 30. 6.
 * printcounts
 * ! fix first chamber as reference
 * parameter
 * 1001  0.  -1.
 * 1002  0.  -1.
 * 1003  0.  -1.
 * 1004  0.  -1.
 * 1005  0.  -1.
 * 1006  0.  -1.
 * end
 * \endcode
 */
void exampleDc() {

	// detector layers (ordered in Z):
	// name, position (x,y,z), thickness (X/X_0), xz-angle, stereo-angle, resolution
	double thickness[12];
	thickness[0] = 0.0025; // gas, wire, wall
	for (unsigned int iLayer = 1; iLayer < 11; ++iLayer)
		thickness[iLayer] = 0.0005; // gas, wire
	thickness[11] = 0.0025; // gas, wire, wall
	// list of layers
	std::vector<GblDetectorLayer> layers;
	const double theta(25.), cosTheta(cos(theta / 180. * M_PI)), sinTheta(
			sin(theta / 180. * M_PI));
	// 1. chamber at distance of 200 cm
	double dist = 200.;
	for (unsigned int iLayer = 0; iLayer < 6; ++iLayer) {
		layers.push_back(
				CreateLayerDc("CH1+", 1, dist * sinTheta, 0.0, dist * cosTheta,
						thickness[iLayer], theta, 6., 0.030)); // +6 deg stereo layers
		dist += 2.;
	}
	for (unsigned int iLayer = 6; iLayer < 12; ++iLayer) {
		layers.push_back(
				CreateLayerDc("CH1-", 1, dist * sinTheta, 0.0, dist * cosTheta,
						thickness[iLayer], theta, -6., 0.030)); // -6 deg stereo layers
		dist += 2.;
	}
	// 2. chamber at distance of 300 cm
	dist = 300.;
	for (unsigned int iLayer = 0; iLayer < 6; ++iLayer) {
		layers.push_back(
				CreateLayerDc("CH2+", 2, dist * sinTheta, 0.0, dist * cosTheta,
						thickness[iLayer], theta, 6., 0.030)); // +6 deg stereo layers
		dist += 2.;
	}
	for (unsigned int iLayer = 6; iLayer < 12; ++iLayer) {
		layers.push_back(
				CreateLayerDc("CH2-", 2, dist * sinTheta, 0.0, dist * cosTheta,
						thickness[iLayer], theta, -6., 0.030)); // -6 deg stereo layers
		dist += 2.;
	}
	// 3. chamber at distance of 400 cm
	dist = 400.;
	for (unsigned int iLayer = 0; iLayer < 6; ++iLayer) {
		layers.push_back(
				CreateLayerDc("CH3+", 3, dist * sinTheta, 0.0, dist * cosTheta,
						thickness[iLayer], theta, 6., 0.030)); // +6 deg stereo layers
		dist += 2.;
	}
	for (unsigned int iLayer = 6; iLayer < 12; ++iLayer) {
		layers.push_back(
				CreateLayerDc("CH3-", 3, dist * sinTheta, 0.0, dist * cosTheta,
						thickness[iLayer], theta, -6., 0.030)); // -6 deg stereo layers
		dist += 2.;
	}

	/* print layers
	 for (unsigned int iLayer = 0; iLayer < layers.size(); ++iLayer) {
	 layers[iLayer].print();
	 } */

	unsigned int nTry = 10000; //: number of tries
	std::cout << " GblDc $Rev$ " << nTry << ", " << layers.size()
			<< std::endl;
	srand(4711);
	clock_t startTime = clock();

	double qbyp = 0.2; // 5 GeV
	// const double bfac = 0.003;  // B*c for 1 T
	const double bfac = 0.;  // B*c for 0 T

	MilleBinary mille; // for producing MillePede-II binary file

	double Chi2Sum = 0.;
	int NdfSum = 0;
	double LostSum = 0.;
	int numFit = 0;

	for (unsigned int iTry = 0; iTry < nTry; ++iTry) {

		// helix parameter for track generation
		const double genDca = 0.1 * unrm(); // normal
		const double genZ0 = 0.1 * unrm(); // normal
		const double genPhi0 = 0.52 * (2. * unif() - 1.); // uniform, [-30..30] deg
		const double genDzds = 10. * unif() + 1.2; // uniform, lambda ~ [50..85] deg
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
			phi0 += unrm() * errMs / cosLambda;		       // scattering for phi
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
			double sArc = pred.getArcLength();	// arc-length
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
			Matrix2d proM2l = transG2l * transM2g.block<3, 2>(0, 0);// skip measurement normal
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
			Vector3d pos = pred.getPosition();
			Vector3d dir = pred.getDirection();
			/* Layer alignment in local (measurement) system
			 for (int p = 0; p < 6; p++)
			 labGlobal[p] = (iLayer + 1) * 10 + p + 1;
			 Matrix<double, 2, 6> derGlobal = layer.getRigidBodyDerLocal(pos,
			 dir); */
			// Chamber alignment in global system (as common system for both stereo orientations)
			unsigned int layerID = layer.getLayerID();
			for (int p = 0; p < 6; p++)
				labGlobal[p] = layerID * 1000 + p + 1;		// chamber alignment
			Matrix<double, 2, 6> derGlobal = layer.getRigidBodyDerGlobal(pos,
					dir).block<2, 6>(0, 0);
			point.addGlobals(labGlobal, derGlobal);
			// add scatterer to point
			double radlen = layer.getRadiationLength()
					/ fabs(pred.getCosIncidence());
			double errMs = gblMultipleScatteringError(qbyp, radlen);// simple model
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
		//std::cout << " Fit " << iTry << ": "<< Chi2 << ", " << Ndf << ", " << lostWeight << std::endl;
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

/// Create a drift chamber layer with 1D measurement.
/**
 * Create drift chamber layer with 1D measurement (u)
 * \param [in] aName       name
 * \param [in] layer       layer ID
 * \param [in] xPos        X-position (of center)
 * \param [in] yPos        Y-position (of center)
 * \param [in] zPos        Z-position (of center)
 * \param [in] thickness   thickness / radiation_length
 * \param [in] xzAngle     angle of normal in XZ plane
 * \param [in] stereoAngle stereo angle
 * \param [in] uRes        resolution in u-direction
 */
GblDetectorLayer CreateLayerDc(const std::string aName, unsigned int layer,
		double xPos, double yPos, double zPos, double thickness, double xzAngle,
		double stereoAngle, double uRes) {
	Vector3d aCenter(xPos, yPos, zPos);
	Vector2d aResolution(uRes, 0.);
	Vector2d aPrecision(1. / (uRes * uRes), 0.);
	Matrix3d measTrafo;
	const double cosXz = cos(xzAngle / 180. * M_PI);
	const double sinXz = sin(xzAngle / 180. * M_PI);
	const double cosSt = cos(stereoAngle / 180. * M_PI);
	const double sinSt = sin(stereoAngle / 180. * M_PI);
	measTrafo << cosSt * cosXz, sinSt, -cosSt * sinXz, -sinSt * cosXz, cosSt, sinSt
			* sinXz, sinXz, 0., cosXz; // U,V,N
	Matrix3d alignTrafo;
	alignTrafo << cosXz, 0., -sinXz, 0., 1., 0., sinXz, 0., cosXz; // I,J,K
	return GblDetectorLayer(aName, layer, 1, thickness, aCenter, aResolution,
			aPrecision, measTrafo, alignTrafo);
}

}

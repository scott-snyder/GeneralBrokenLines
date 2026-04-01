/*
 * exampleCdc.cpp
 *
 *  Created on: 23 Mar 2023
 *      Author: kleinwrt
 */

/** \file
 *  Example (cylindrical) drift chamber application (composed trajectory with geometric constraint).
 *
 *  \author Claus Kleinwort, DESY, 2023 (Claus.Kleinwort@desy.de)
 *
 *  \copyright
 *  Copyright (c) 2023 Deutsches Elektronen-Synchroton,
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

#include "exampleUtilCdc.h"
#include "GblTrajectory.h"
#include "Mille/MilleFactory.h"
//#include <fstream>

using namespace gbl;
using namespace Eigen;

/// Drift chamber example
/**
 * Simulate and reconstruct helical tracks in a (cylindrical) drift chamber.
 *
 *  Create points on initial trajectories, create \b COMPOSED trajectory (for two tracks from same vertex)
 *  from points, fit and write trajectory to MP-II binary file (for rigid body alignment).
 *
 *  Setup:
 *   - Collider, beams in Z direction
 *   - Constant magnetic field in Z direction
 *   - Cylindrical drift chamber
 *   - No multiple scattering in detectors (gas, wires, walls) (air in between ignored)
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
 * ! fix first layer as reference
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
int main() {

	// detector setup, ~ Belle-II CDC
	const unsigned int nSuper = 9; // number of super layers
	double nLayer[nSuper] = { 8, 6, 6, 6, 6, 6, 6, 6, 6 }; // number of layers per super layer
	double rInner[nSuper] = { 16.80, 25.70, 36.52, 47.69, 58.41, 69.53, 80.25,
			91.37, 102.00 }; // inner radius of super layer
	double rOuter[nSuper] = { 23.80, 34.80, 45.57, 56.69, 67.41, 78.53, 89.25,
			100.37, 111.14 }; // outer radius of super layer
	double stereo[nSuper] = { 0., 0.068, 0., -0.060, 0., 0.064, 0., -0.072, 0. }; // stereo angle of super layer
	double zStart[nSuper] = { -35.9, -51.4, -57.5, -59.6, -61.8, -63.9, -66.0,
			-68.2, -70.2 }; // -Z end of wires per super layer
	double zEnd[nSuper] = { 67.9, 98.6, 132.9, 144.7, 146.8, 148.9, 151.0,
			153.2, 155.3 }; // +Z end of wires per super layer

	unsigned int nTry = 1000; //: number of tries
	std::cout << " GblComposedGeo " << nTry << ", " << nSuper << std::endl;
	srand(4711);
	clock_t startTime = clock();

	const double bfac = 0.0045;  // B*c for 1.5 T

	double beamPos[] = { 0., 0., 0. };
	double beamSize[] = { 0.005, 0.005, 0.3 };
	const bool useBeamSpot = true;

	auto mille = Mille::spawnMilleRecord("milleBinary1.dat"); // for producing MillePede-II binary file

	double Chi2Sum = 0.;
	int NdfSum = 0;
	double LostSum = 0.;
	int numFit = 0;
	unsigned int iEvent = 0;

	// event loop
	for (unsigned int iTry = 0; iTry < nTry; ++iTry) {
		unsigned int nTrack;
		std::vector<std::array<double, 6>> allGenPar; // collect parameters for generation of tracks

		iEvent++;
		// vertex
		std::array<double, 3> vtx = { beamPos[0] + beamSize[0] * unrm(),
				beamPos[1] + beamSize[1] * unrm(), beamPos[2]
						+ beamSize[2] * unrm() };
		std::cout << " gen(vertex)  " << iEvent << ": " << vtx[0] << " "
				<< vtx[1] << " " << vtx[2] << " " << std::endl;
		// inline generation of track pair
		nTrack = 2;
		for (unsigned int iTrack = 0; iTrack < nTrack; ++iTrack) {

			// helix parameter for track generation
			const double qbyp = 0.2; // 5 GeV
			const double genPhi0 = M_PI * (2. * unif() - 1.); // uniform, [-30..30] deg
			const double genDzds = (2. * unif() - 1.); // uniform, lambda ~ [-45..45] deg
			const double genCurv = bfac * qbyp * sqrt(1. + genDzds * genDzds);
			// from common vertex
			const double genDca = sin(genPhi0) * vtx[0] - cos(genPhi0) * vtx[1];
			const double genZ0 = vtx[2];

			std::cout << " gen(inline)  " << iEvent << ": " << genCurv << " "
					<< genPhi0 << " " << genDca << " " << genDzds << " "
					<< genZ0 << std::endl;
			std::array<double, 6> par = { genCurv, genPhi0, genDca, genDzds,
					genZ0, 1. };
			allGenPar.push_back(par);
		}

		std::vector<std::vector<GblPoint>> listsOfPoints;
		std::vector<Matrix<double, 2, 3>> listOfTrafos;
		for (unsigned int iTrack = 0; iTrack < nTrack; ++iTrack) {
			//
			// generate hits: (virtual) wires as detector layers at constant radii
			//
			std::array<double, 6> genPar = allGenPar[iTrack];
			double curv(genPar[0]), phi0(genPar[1]), dca(genPar[2]), dzds(
					genPar[3]), z0(genPar[4]);
			// local constant (Bfield) helix
			GblSimpleHelix hlx = GblSimpleHelix(curv, phi0, dca, dzds, z0);

			std::vector<GblDetectorLayer> layers;

			// beam spot (required as first (and common) point for composed trajectory)
			double sArc = hlx.getArcLengthXY(beamPos[0], beamPos[1]);
			// add virtual layer
			layers.push_back(
					CreateImpactPar("impact", 0, beamPos[0], beamPos[1],
							beamPos[2], phi0 + sArc * curv, dzds, beamSize[0],
							beamSize[1], beamSize[2]));

			unsigned int cLayer = 0;
			for (unsigned int iSuper = 0; iSuper < nSuper; ++iSuper) {
				double radius = rInner[iSuper];
				double step = (rOuter[iSuper] - radius) / (nLayer[iSuper] - 1);
				bool abort = false; 
				for (unsigned int iLayer = 0; iLayer < nLayer[iSuper];
						++iLayer) {
					// check for |dca| < radius
					if (fabs(dca) < radius) {
						// arc length to layer
						double sArcAtLayer = hlx.getArcLengthR(radius);
						if (sArcAtLayer == 0.){
							abort = true;
							break; 
						}
						// phi of position at layer
						double phiPos = hlx.getPhi(radius);
						// (virtual) wire position
						double xPos(cos(phiPos) * radius), yPos(
								sin(phiPos) * radius), zPos(z0 + sArcAtLayer * dzds);
						// add virtual wire
						if (zPos > zStart[iSuper] and zPos < zEnd[iSuper])
							layers.push_back(
									CreateWireCdc("wire", ++cLayer, xPos, yPos,
											zPos, phi0 + sArcAtLayer * curv, dzds,
											stereo[iSuper], 0.015));
					}
					// next layer
					radius += step;
				}
				if (abort) break; 
			}

			//
			// create GBL trajectory (list of GBL points)
			//
			// seed (with true parameters)
			std::array<double, 6> seedPar = allGenPar[iTrack];
			// optionally distort seed
			//seedPar[0] += 0.01 * seedPar[0] * iTrack;
			//seedPar[1] += 0.001 * iTrack; // different distortions per track
			//seedPar[2] += 0.001 * iTrack; // different distortions per track
			double seedCurv(seedPar[0]), seedPhi0(seedPar[1]), seedDca(
					seedPar[2]), seedDzds(seedPar[3]), seedZ0(seedPar[4]);
			GblSimpleHelix seed = GblSimpleHelix(seedCurv, seedPhi0, seedDca,
					seedDzds, seedZ0);
			// (previous) arc-length
			double sOld = 0.;
			const double cosLambdaSeed = 1. / sqrt(1. + (seedDzds * seedDzds));
			// for multiple scattering (not (yet) implemented)
			const double qbyp = seedCurv / bfac * cosLambdaSeed;

			// transformations to common POSITION parameters (at vertex) of composed trajectory at first point
			Matrix<double, 2, 3> innerTrafo = Matrix<double, 2, 3>::Zero();
			// common parameters (offsets in curvilinear system)
			innerTrafo(0, 0) = -sin(seedPhi0); // u vs vtx
			innerTrafo(0, 1) = cos(seedPhi0); // u vs vty
			innerTrafo(1, 0) = -cos(seedPhi0) * cosLambdaSeed * seedDzds; // v vs vtx
			innerTrafo(1, 1) = -sin(seedPhi0) * cosLambdaSeed * seedDzds; // v vs vty
			innerTrafo(1, 2) = cosLambdaSeed; // v vs vtz

			listOfTrafos.push_back(innerTrafo);

			// list of points on trajectory
			std::vector<GblPoint> listOfPoints;
			for (unsigned int iLayer = 0; iLayer < layers.size(); ++iLayer) {
				GblDetectorLayer &layer = layers[iLayer];
				unsigned int lid = layer.getLayerID(); // layer ID (0 = vertex)
				// std::cout << " hit " << lid << " " << hits[iLayer].transpose() << std::endl;
				// prediction from seeding helix
				GblHelixPrediction pred = layer.intersectWithHelix(seed);
				double sArcAtLayer = pred.getArcLength();	// arc-length
				Vector2d measPrediction = pred.getMeasPred(); // measurement prediction
				Vector2d measPrecision = layer.getPrecision(); // measurement precision
				// residuals
				Vector2d res = measPrediction; // (virtual) wire positioned at measurement
				// smear u according to resolution
				Vector2d sigma = layer.getResolution();
				if (lid > 0) // wire layer?
					res[0] += sigma[0] * unrm();
				else
					std::cout << " impact par " << res[0] << " " << res[1]
							<< " " << seedCurv << " " << seedPhi0 << " "
							<< seedDca << " " << seedDzds << " " << seedZ0
							<< " " << sArcAtLayer << " " << layers.size() << std::endl;
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
						(sArcAtLayer - sOld) / cosLambdaSeed, cosLambdaSeed, bfac);
				sOld = sArcAtLayer;
				// point with (independent) measurements (in measurement system)
				GblPoint point(jacPointToPoint);
				if (lid > 0 or (iTrack == 0 and useBeamSpot)) // vertex only for first track (large correlations!)
					point.addMeasurement(proL2m, res, measPrecision);
				// global labels and parameters for rigid body alignment
				if (lid > 0) {

					std::vector<int> labGlobal(6);
					Vector3d pos = pred.getPosition();
					Vector3d dir = pred.getDirection();
					/* Layer alignment in local (measurement) system
					 for (int p = 0; p < 6; p++)
					 labGlobal[p] = (iLayer + 1) * 10 + p + 1;
					 Matrix<double, 2, 6> derGlobal = layer.getRigidBodyDerLocal(pos,
					 dir);
					 */
					// Chamber alignment in global system
					unsigned int layerID = layer.getLayerID();
					for (int p = 0; p < 6; p++)
						labGlobal[p] = layerID * 1000 + p + 1; // layer alignment
					Matrix<double, 2, 6> derGlobal =
							layer.getRigidBodyDerGlobal(pos, dir).block<2, 6>(0,
									0);
					point.addGlobals(labGlobal, derGlobal);
					// add scatterer to point
					double radlen = layer.getRadiationLength()
							/ fabs(pred.getCosIncidence());
					double errMs = gblMultipleScatteringError(qbyp, radlen); // simple model
					if (errMs > 0.) {
						Vector2d scat(0., 0.);
						Vector2d scatPrec(1. / (errMs * errMs),
								1. / (errMs * errMs)); // scattering precision matrix is diagonal in curvilinear system
						point.addScatterer(scat, scatPrec);
					}
				}
				// add point to trajectory
				listOfPoints.push_back(point);
			}
			listsOfPoints.push_back(listOfPoints);
		}

		// require 2 decent trajectories for composed trajectory
		unsigned int nAccepted = 0;
		for (unsigned int iTrack = 0; iTrack < nTrack; ++iTrack)
			// sufficient length?
			if (listsOfPoints[iTrack].size() >= 10)
				nAccepted++;
		if (nAccepted != 2)
			continue;

		// prepare composed trajectory
		std::vector<std::pair<std::vector<GblPoint>, MatrixXd> > aPointsAndTransList;
		for (unsigned int iTrack = 0; iTrack < nTrack; ++iTrack) {
			aPointsAndTransList.push_back(
					make_pair(listsOfPoints[iTrack], listOfTrafos[iTrack]));
		}
		//
		// fit composed GBL trajectory
		//

		// create trajectory
		GblTrajectory traj(aPointsAndTransList);
		//std::cout << " Trajectory contructed: " << traj.isValid() << std::endl;
		// fit trajectory
		double Chi2;
		int Ndf;
		double lostWeight;
		unsigned int ierr = traj.fit(Chi2, Ndf, lostWeight);
		std::cout << " Fit " << iTry << ": " << Chi2 << ", " << Ndf << ", "
				<< lostWeight << std::endl;
		traj.printTrajectory();
		//traj.printData();
		// successfully fitted?
		if (!ierr) {
			// write to MP binary file
			traj.milleOut(mille);
			// update statistics
			Chi2Sum += Chi2;
			NdfSum += Ndf;
			LostSum += lostWeight;
			numFit++;
			/* look at (track parameter) corrections
			 VectorXd aCorrection(5);
			 MatrixXd aCovariance(5, 5);
			 traj.getResults(1, aCorrection, aCovariance);
			 std::cout << " cor1 " << std::endl << aCorrection << std::endl;
			 std::cout << " cov1 " << std::endl << aCovariance << std::endl;
			 */
		}

	}

	clock_t endTime = clock();
	double diff = endTime - startTime;
	double cps = CLOCKS_PER_SEC;
	std::cout << " Time elapsed " << diff / cps << " s" << std::endl;
	std::cout << " Chi2/Ndf = " << Chi2Sum / NdfSum << std::endl;
	std::cout << " Tracks fitted " << numFit << std::endl;
	if (LostSum > 0.)
		std::cout << " Weight lost   " << LostSum << std::endl;
	return 0; 
}


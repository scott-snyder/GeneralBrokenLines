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
 *   - Multiple scattering in sensors (air inbetween ignored)
 *   - Curvilinear system (T,U,V) as local coordinate system and (q/p, slopes, offsets) as local track parameters
 *
 * \remark To exercise (mis)alignment different sets of layers (with different geometry)
 * for simulation and reconstruction can be used.
 */
void exampleSit() {

	// detector layers (ordered in X):
	// name, position (x,y,z), thickness (X/X_0), (1 or 2) measurements (direction in YZ, resolution)
	std::vector<GblSiliconLayer> layers;
	layers.push_back(
			GblSiliconLayer("PIX1", 2.0, 0., 0., 0.0033, 0., 0.0010, 90.,
					0.0020)); // pixel
	layers.push_back(
			GblSiliconLayer("PIX2", 3.0, 0., 0., 0.0033, 0., 0.0010, 90.,
					0.0020)); // pixel
	layers.push_back(
			GblSiliconLayer("PIX3", 4.0, 0., 0., 0.0033, 0., 0.0010, 90.,
					0.0020)); // pixel
	layers.push_back(
			GblSiliconLayer("S2D4", 6.0, 0., 0., 0.0033, 0., 0.0025, 5.0,
					0.0025)); // strip 2D, +5 deg stereo angle
	layers.push_back(
			GblSiliconLayer("S2D5", 8.0, 0., 0., 0.0033, 0., 0.0025, -5.,
					0.0025)); // strip 2D, -5 deg stereo angle
	layers.push_back(
			GblSiliconLayer("S2D6", 10., 0., 0., 0.0033, 0., 0.0025, 5.0,
					0.0025)); // strip 2D, +5 deg stereo angle
	layers.push_back(
			GblSiliconLayer("S2D7", 12., 0., 0., 0.0033, 0., 0.0025, -5.,
					0.0025)); // strip 2D, -5 deg stereo angle
	layers.push_back(GblSiliconLayer("S1D8", 15., 0., 0., 0.0033, 0., 0.0040)); // strip 1D

	/* print layers
	 for (unsigned int iLayer = 0; iLayer < layers.size(); ++iLayer) {
	 layers[iLayer].print();
	 } */

	unsigned int nTry = 1000; //: number of tries
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
			Vector3d measPos = pred.getMeasPos();
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
		const double sinLambdaSeed = seedDzds * cosLambdaSeed;
		// list of points on trajectory
		std::vector<GblPoint> listOfPoints;
		for (unsigned int iLayer = 0; iLayer < layers.size(); ++iLayer) {
			// std::cout << " hit " << iLayer << " " << hits[iLayer].transpose() << std::endl;
			// prediction from seeding helix
			GblHelixPrediction pred = layers[iLayer].intersectWithHelix(seed);
			double sArc = pred.getArcLength(); // arc-length
			Vector2d measPrediction = pred.getMeasPred(); // measurement prediction
			Vector2d measPrecision = layers[iLayer].getPrecision(); // measurement precision
			// residuals
			Vector2d res(hits[iLayer][0] - measPrediction[0],
					hits[iLayer][1] - measPrediction[1]);
			// Curvilinear system: track direction T, U = Z x T / |Z x T|, V = T x U
			// as local system
			const double phi = seedPhi0 + seedCurv * sArc;
			const double cosPhi = cos(phi);
			const double sinPhi = sin(phi);
			Vector3d uDir(-sinPhi, cosPhi, 0.); // U direction
			Vector3d vDir(-sinLambdaSeed * cosPhi, -sinLambdaSeed * sinPhi,
					cosLambdaSeed); // V direction
			// transformation global system to local (u,v) (matrix from row vectors)
			Matrix<double, 2, 3> transG2l;
			transG2l << uDir.transpose(), vDir.transpose();
			// transformation measurement system to global system
			Matrix3d transM2g = pred.getMeasSystemDirs().inverse();
			// projection matrix (measurement plane to local (u,v))
			Matrix2d proM2l = transG2l * transM2g.block<3, 2>(0, 0); // skip measurement normal
			// projection matrix (local (u,v) to measurement plane)
			Matrix2d proL2m = proM2l.inverse();
			// propagation
			Matrix5d jacPointToPoint = gblSimpleJacobian(sArc - sOld,
					cosLambdaSeed, bfac);
			sOld = sArc;
			// point with (independent) measurements (in measurement system)
			GblPoint point(jacPointToPoint);
			point.addMeasurement(proL2m, res, measPrecision);
			// global labels and parameters for rigid body alignment
			std::vector<int> labGlobal(6);
			for (int p = 0; p < 6; p++)
				labGlobal[p] = iLayer * 10 + p + 1;
			Vector3d rotCenter = layers[iLayer].getCenter();
			Matrix<double, 2, 6> derGlobal = pred.getRigidBodyDerGlobal(
					rotCenter).block<2, 6>(0, 0);
			point.addGlobals(labGlobal, derGlobal);
			// add scatterer to point
			double radlen = layers[iLayer].getRadiationLength()
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

/// Multiple scattering error
/**
 * Angular error in plane, simple model (Rossi, Greisen)
 * \param [in] qbyp    q/p [1/GeV]
 * \param [in] xbyx0   thickness / radiation length
 */
double gblMultipleScatteringError(double qbyp, double xbyx0) {
	return 0.015 * fabs(qbyp) * sqrt(xbyx0);
}

/// Simple jacobian: quadratic in arc length difference
/**
 * \param [in] ds    (3D) arc-length
 * \param [in] cosl  cos(lambda)
 * \param [in] bfac  Bz*c
 * \return jacobian
 */
Matrix5d gblSimpleJacobian(double ds, double cosl, double bfac) {
	Matrix5d jac;
	jac.setIdentity();
	jac(1, 0) = -bfac * ds * cosl;
	jac(3, 0) = -0.5 * bfac * ds * ds * cosl;
	jac(3, 1) = ds;
	jac(4, 2) = ds;
	return jac;
}

///  unit normal distribution, Box-Muller method, polar form
double unrm() {
	static double unrm2 = 0.0;
	static bool cached = false;
	if (!cached) {
		double x, y, r;
		do {
			x = 2.0 * static_cast<double>(rand())
					/ static_cast<double>(RAND_MAX) - 1;
			y = 2.0 * static_cast<double>(rand())
					/ static_cast<double>(RAND_MAX) - 1;
			r = x * x + y * y;
		} while (r == 0.0 || r > 1.0);
		// (x,y) in unit circle
		double d = sqrt(-2.0 * log(r) / r);
		double unrm1 = x * d;
		unrm2 = y * d;
		cached = true;
		return unrm1;
	} else {
		cached = false;
		return unrm2;
	}
}

///  uniform distribution [0..1]
double unif() {
	return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

/// Create helix prediction.
/**
 * \param [in] sArc     arc length
 * \param [in] aPred    prediction for measurement (u,v)
 * \param [in] tDir     track direction at prediction
 * \param [in] uDir     measurement direction for u
 * \param [in] vDir     measurement direction for v
 * \param [in] nDir     normal to measurement plane
 * \param [in] aPos     position at prediction
 */
GblHelixPrediction::GblHelixPrediction(double sArc, const Vector2d& aPred,
		const Vector3d& tDir, const Vector3d& uDir, const Vector3d& vDir,
		const Vector3d& nDir, const Vector3d& aPos) :
		sarc(sArc), pred(aPred), tdir(tDir), udir(uDir), vdir(vDir), ndir(nDir), pos(
				aPos) {
	global2meas << uDir.transpose(), vDir.transpose(), nDir.transpose();
}

GblHelixPrediction::~GblHelixPrediction() {
}

/// Get arc-length.
double GblHelixPrediction::getArcLength() const {
	return sarc;
}

/// Get prediction.
const Vector2d& GblHelixPrediction::getMeasPred() const {
	return pred;
}

/// Get position.
const Vector3d& GblHelixPrediction::getMeasPos() const {
	return pos;
}

/// Get directions of measurement system.
/**
 * Matrix from row vectors (transformation from global to measurement system)
 */
const Eigen::Matrix3d& GblHelixPrediction::getMeasSystemDirs() const {
	return global2meas;;
}

/// Get cosine of incidence
double GblHelixPrediction::getCosIncidence() const {
	return tdir.dot(ndir);
}

/// Get rigid body derivatives in global frame.
/**
 * \param [in] rotCenter  rotation center
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
 * \endcode
 */
const Matrix<double, 3, 6> GblHelixPrediction::getRigidBodyDerGlobal(
		Eigen::Vector3d& rotCenter) const {
// lever arms (for rotations)
	Vector3d dist = pos - rotCenter;
// dr/dm (residual vs measurement, 1-tdir*ndir^t/tdir*ndir)
	Matrix3d drdm = Matrix3d::Identity()
			- (tdir * ndir.transpose()) / (tdir.transpose() * ndir);
// dm/dg (measurement vs 6 rigid body parameters, global system)
	Matrix<double, 3, 6> dmdg = Matrix<double, 3, 6>::Zero();
	dmdg(0, 0) = 1.;
	dmdg(0, 4) = dist(2);
	dmdg(0, 5) = -dist(1);
	dmdg(1, 1) = 1.;
	dmdg(1, 3) = -dist(2);
	dmdg(1, 5) = dist(0);
	dmdg(2, 2) = 1.;
	dmdg(2, 3) = dist(1);
	dmdg(2, 4) = -dist(0);
// drl/dg (local residuals vs rigid body parameters)
	return global2meas * drdm * dmdg;
}

/// Create simple helix.
/**
 * \param [in] aRinv      curvature (1/R)
 * \param [in] aPhi0      azimuth at PCA
 * \param [in] aDca       XY distance at PCA
 * \param [in] aDzds      slope in ZS
 * \param [in] aZ0        offset in ZS
 */
GblSimpleHelix::GblSimpleHelix(double aRinv, double aPhi0, double aDca,
		double aDzds, double aZ0) :
		rinv(aRinv), phi0(aPhi0), dca(aDca), dzds(aDzds), z0(aZ0), cosPhi0(
				cos(phi0)), sinPhi0(sin(phi0)), xRelCenter(
				-(1. - dca * rinv) * sinPhi0), yRelCenter(
				(1. - dca * rinv) * cosPhi0) {
}

GblSimpleHelix::~GblSimpleHelix() {
}

/// Get (2D) arc length for given point.
/**
 * \param [in] xPos   X Position
 * \param [in] yPos   Y Position
 */
double GblSimpleHelix::getArcLengthXY(double xPos, double yPos) const {
// line
	if (rinv == 0)
		return cosPhi0 * xPos + sinPhi0 * yPos;
// helix
	double dx = xPos * rinv - xRelCenter;
	double dy = yPos * rinv - yRelCenter;
	double dphi = atan2(dx, -dy) - phi0;
	if (dphi > M_PI)
		dphi -= 2.0 * M_PI;
	else if (dphi < -M_PI)
		dphi += 2.0 * M_PI;
	return dphi / rinv;
}

/// Move to new reference point (X,y)
/**
 * \param [in] xPos      X Position
 * \param [in] yPos      Y Position
 * \param [out] newPhi0  new phi0
 * \param [out] newDca   new dca
 * \param [out] newZ0    new z0
 */
void GblSimpleHelix::moveToXY(double xPos, double yPos, double& newPhi0,
		double& newDca, double& newZ0) const {
// start values
	newPhi0 = phi0;
	newDca = dca;
	newZ0 = z0;
// Based on V. Karimaki, NIM A305 (1991) 187-191, eqn (19)
	const double u = 1. - rinv * dca;
	const double dp = -xPos * sinPhi0 + yPos * cosPhi0 + dca;
	const double dl = xPos * cosPhi0 + yPos * sinPhi0;
	const double sa = 2. * dp - rinv * (dp * dp + dl * dl);
	const double sb = rinv * xPos + u * sinPhi0;
	const double sc = -rinv * yPos + u * cosPhi0;
	const double sd = sqrt(1. - rinv * sa);
// transformed parameters
	double sArc;
	if (rinv == 0.) {
		newDca = dp;
		sArc = dl;
	} else {
		newPhi0 = atan2(sb, sc);
		newDca = sa / (1. + sd);
		double dphi = newPhi0 - phi0;
		if (dphi > M_PI)
			dphi -= 2.0 * M_PI;
		else if (dphi < -M_PI)
			dphi += 2.0 * M_PI;
		sArc = dphi / rinv;
	}
	newZ0 += sArc * dzds;
}

/// Get prediction
/*
 * \param [in] refPos  reference position on detector plane
 * \param [in] uDir    measurement direction 'u'
 * \param [in] vDir    measurement direction 'v'
 */
GblHelixPrediction GblSimpleHelix::getPrediction(const Eigen::Vector3d& refPos,
		const Eigen::Vector3d& uDir, const Eigen::Vector3d& vDir) const {
// normal to (u,v) measurement plane
	Vector3d nDir = uDir.cross(vDir).normalized();
// ZS direction
	const double cosLambda = 1. / sqrt(1. + dzds * dzds);
	const double sinLambda = dzds * cosLambda;
	double sArc2D;
	Vector3d dist, pos, tDir;
// line (or helix)
	if (rinv == 0.) {
		// track direction
		tDir << cosLambda * cosPhi0, cosLambda * sinPhi0, sinLambda;
		// distance (of point at dca to reference)
		Vector3d pca(dca * sinPhi0, -dca * cosPhi0, z0);
		dist = pca - refPos;
		// arc-length
		double sArc3D = -dist.dot(nDir) / tDir.dot(nDir);
		sArc2D = sArc3D * cosLambda;
		// position at prediction
		pos = pca + sArc3D * tDir;
		// distance (of point at sArc to reference)
		dist = pos - refPos;
	} else {
		// initial guess of 2D arc-length
		sArc2D = this->getArcLengthXY(refPos(0), refPos(1));
		unsigned int nIter = 0;
		while (nIter < 10) {
			nIter += 1;
			// track direction
			const double dPhi = sArc2D * rinv;
			const double cosPhi = cos(phi0 + dPhi);
			const double sinPhi = sin(phi0 + dPhi);
			tDir << cosLambda * cosPhi, cosLambda * sinPhi, sinLambda;
			// position at prediction
			pos << (xRelCenter + sinPhi) / rinv, (yRelCenter - cosPhi) / rinv, z0
					+ dzds * sArc2D;
			// distance (of point at sArc to reference)
			dist = pos - refPos;
			// arc-length correction (linearizing helix at sArc)
			const double sCorr3D = -dist.dot(nDir) / tDir.dot(nDir);
			if (fabs(sCorr3D) > 0.00001) {
				// iterate
				sArc2D += sCorr3D * cosLambda;
			} else {
				// converged
				break;
			}
		}
	}
// projections on measurement directions
	Vector2d pred(dist.dot(uDir), dist.dot(vDir));
	return GblHelixPrediction(sArc2D, pred, tDir, uDir, vDir, nDir, pos);
}

/// Create a silicon layer with 1D measurement.
/**
 * Create silicon layer with 1D measurement (u)
 * \param [in] aName      name
 * \param [in] xPos       X-position (of center)
 * \param [in] yPos       Y-position (of center)
 * \param [in] zPos       Z-position (of center)
 * \param [in] thickness  thickness / radiation_length
 * \param [in] uAngle     angle of u-direction in YZ plane
 * \param [in] uRes       resolution in u-direction
 */
GblSiliconLayer::GblSiliconLayer(const std::string aName, double xPos,
		double yPos, double zPos, double thickness, double uAngle, double uRes) :
		name(aName), measDim(1), xbyx0(thickness), center(xPos, yPos, zPos), resolution(
				uRes, 0.), precision(1. / (uRes * uRes), 0.), udir(0.,
				cos(uAngle / 180. * M_PI), sin(uAngle / 180. * M_PI)), vdir(0.,
				-udir(2), udir(1)) {
}

/// Create a silicon layer with 2D measurement.
/**
 * Create silicon layer with 2D measurement (u,v)
 * \param [in] aName      name
 * \param [in] xPos       X-position (of center)
 * \param [in] yPos       Y-position (of center)
 * \param [in] zPos       Z-position (of center)
 * \param [in] thickness  thickness / radiation_length
 * \param [in] uAngle     angle of u-direction in YZ plane
 * \param [in] uRes       resolution in u-direction
 * \param [in] vAngle     angle of v-direction in YZ plane
 * \param [in] vRes       resolution in v-direction
 */
GblSiliconLayer::GblSiliconLayer(const std::string aName, double xPos,
		double yPos, double zPos, double thickness, double uAngle, double uRes,
		double vAngle, double vRes) :
		name(aName), measDim(2), xbyx0(thickness), center(xPos, yPos, zPos), resolution(
				uRes, vRes), precision(1. / (uRes * uRes), 1. / (vRes * vRes)), udir(
				0., cos(uAngle / 180. * M_PI), sin(uAngle / 180. * M_PI)), vdir(
				0., cos(vAngle / 180. * M_PI), sin(vAngle / 180. * M_PI)) {
}

GblSiliconLayer::~GblSiliconLayer() {
}

/// Print GblSiliconlayer.
void GblSiliconLayer::print() const {
	IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::cout << " Layer " << name << " " << measDim << "D, " << xbyx0
			<< " X0, @ " << center.transpose().format(CleanFmt) << ", res "
			<< resolution.transpose().format(CleanFmt) << ", udir "
			<< udir.transpose().format(CleanFmt) << ", vdir "
			<< vdir.transpose().format(CleanFmt) << std::endl;
}

/// Get radiation length.
double GblSiliconLayer::getRadiationLength() const {
	return xbyx0;
}

/// Get resolution.
Eigen::Vector2d GblSiliconLayer::getResolution() const {
	return resolution;
}

/// Get precision.
Eigen::Vector2d GblSiliconLayer::getPrecision() const {
	return precision;
}

/// Get center.
Eigen::Vector3d GblSiliconLayer::getCenter() const {
	return center;
}

/// Intersect with helix.
/**
 * \param [in]  hlx  helix
 */
GblHelixPrediction GblSiliconLayer::intersectWithHelix(
		GblSimpleHelix hlx) const {
	return hlx.getPrediction(center, udir, vdir);
}

}

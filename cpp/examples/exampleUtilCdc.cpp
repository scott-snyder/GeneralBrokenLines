/*
 * exampleCdc.cpp
 *
 *  Created on: 23 Mar 2023
 *      Author: kleinwrt
 */

/** \file
 *  ExampleUtil(ities) for a Cylindrical Drift Chamber.
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

#include "exampleUtilCdc.h"

using namespace Eigen;

//! Namespace for the general broken lines package
namespace gbl {

/// Create a drift chamber wire with 1D measurement.
/**
 * Create drift chamber wire at given position with 1D measurement (u)
 * \param [in] aName       name
 * \param [in] layer       layer ID
 * \param [in] xPos        X-position (of center)
 * \param [in] yPos        Y-position (of center)
 * \param [in] zPos        Z-position (of center)
 * \param [in] phi         track direction in XY
 * \param [in] tanLambda   track direction in ZS
 * \param [in] stereoAngle stereo angle
 * \param [in] uRes        resolution in u-direction
 */
GblDetectorLayer CreateWireCdc(const std::string aName, unsigned int layer,
		double xPos, double yPos, double zPos, double phi, double tanLambda,
		double stereoAngle, double uRes) {
	Vector3d aCenter(xPos, yPos, zPos);
	Vector2d aResolution(uRes, 0.);
	Vector2d aPrecision(1. / (uRes * uRes), 0.);
// track direction (N)
	const double cosLambda = 1. / sqrt(1. + tanLambda * tanLambda);
	Vector3d tDir(cos(phi) * cosLambda, sin(phi) * cosLambda,
			tanLambda * cosLambda);
// wire direction (V)
	const double rad = sqrt(xPos * xPos + yPos * yPos);
	const double scale = sqrt(1. + stereoAngle * stereoAngle);
	Vector3d vDir(-yPos * stereoAngle / rad / scale,
			xPos * stereoAngle / rad / scale, 1. / scale);
// measurement direction perpendicular to track and wire direction (U)
	Vector3d uDir = vDir.cross(tDir).normalized();
// normal to measurement plane
	Vector3d nDir = uDir.cross(vDir).normalized();

	Matrix3d measTrafo;
	measTrafo << uDir[0], uDir[1], uDir[2], vDir[0], vDir[1], vDir[2], nDir[0], nDir[1], nDir[2]; // U,V,N
	Matrix3d alignTrafo = Matrix3d::Identity();			 // I,J,K
	return GblDetectorLayer(aName, layer, 1, 0.0, aCenter, aResolution,
			aPrecision, measTrafo, alignTrafo); // no mult. scat.
}

/// Create detector plane for impact parameters as 2D measurement.
/**
 * \param [in] aName       name
 * \param [in] layer       layer ID
 * \param [in] xPos        X-position (of center)
 * \param [in] yPos        Y-position (of center)
 * \param [in] zPos        Z-position (of center)
 * \param [in] phi         track direction in XY
 * \param [in] tanLambda   track direction in ZS
 * \param [in] xRes        resolution in x-direction
 * \param [in] yRes        resolution in y-direction
 * \param [in] zRes        resolution in z-direction
 */
GblDetectorLayer CreateImpactPar(const std::string aName, unsigned int layer,
		double xPos, double yPos, double zPos, double phi, double tanLambda,
		double xRes, double yRes, double zRes) {
// track direction (N)
	const double cosLambda = 1. / sqrt(1. + tanLambda * tanLambda);
	Vector3d tDir(cos(phi) * cosLambda, sin(phi) * cosLambda,
			tanLambda * cosLambda);
// beam direction (V) = Z
	Vector3d vDir(0., 0., 1.);
// measurement direction perpendicular to track and beam direction (U)
	Vector3d uDir = vDir.cross(tDir).normalized();
// normal to measurement plane
	Vector3d nDir = uDir.cross(vDir).normalized();
// beam variance in U direction
	const double uVar = uDir[0] * xRes * xRes * uDir[0]
			+ uDir[1] * yRes * yRes * uDir[1];
	Vector3d aCenter(xPos, yPos, zPos);
	Vector2d aResolution(sqrt(uVar), zRes);
	Vector2d aPrecision(1. / uVar, 1. / (zRes * zRes));

	Matrix3d measTrafo;
	measTrafo << uDir[0], uDir[1], uDir[2], vDir[0], vDir[1], vDir[2], nDir[0], nDir[1], nDir[2]; // U,V,N
	Matrix3d alignTrafo = Matrix3d::Identity();			 // I,J,K
	return GblDetectorLayer(aName, layer, 2, 0., aCenter, aResolution,
			aPrecision, measTrafo, alignTrafo);
}

}


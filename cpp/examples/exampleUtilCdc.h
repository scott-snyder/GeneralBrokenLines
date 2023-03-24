/*
 * exampleUtilCdc.h
 *
 *  Created on: 23 Mar 2023
 *      Author: kleinwrt
 */

/** \file
 *  Definitions for exampleUtil(ities) for a Cylindrical Drift Chamber.
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

#ifndef SRC_EXAMPLEUTILCDC_H_
#define SRC_EXAMPLEUTILCDC_H_

#include "GblUtilities.h"

//! Namespace for the general broken lines package
namespace gbl {

GblDetectorLayer CreateWireCdc(const std::string aName, unsigned int layer,
		double xPos, double yPos, double zPos, double phi, double tanLambda,
		double stereoAngle, double uRes);
GblDetectorLayer CreateImpactPar(const std::string aName, unsigned int layer,
		double xPos, double yPos, double zPos, double phi, double tanLambda,
		double xRes, double yRes, double zRes);

}

#endif /* SRC_EXAMPLEUTILCDC_H_ */

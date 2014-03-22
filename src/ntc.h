/**
 * @file ntc.h
 * @brief NTC thermistor library
 * @version 1.0
 * @copyright GNU Lesser General Public License version 3
 *            <http://www.gnu.org/licenses/lgpl.html>
 * Copyright (c) 2007, 2013 - SoftQuadrat GmbH, Germany
 * Contact: thermistor (at) softquadrat.de
 * Web site: thermistor.sourceforge.net
 *******************************************************************************
 * This program is free software: you can redistribute it and/or modif         *
 *    it under the terms of the GNU Lesser General Public License as published *
 *    by the Free Software Foundation, either version 3 of the License, or     *
 *    (at your option) any later version.                                      *
 *                                                                             *
 *    This program is distributed in the hope that it will be useful,          *
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *    GNU Lesser General Public License for more details.                      *
 *                                                                             *
 *    You should have received a copy of the GNU Lesser General Public License *
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.    *
 *******************************************************************************
 */
#ifndef _NTC_H_
#define _NTC_H_
#ifdef __cplusplus
extern "C" {
#endif
/* ========================================================================== */

/* internal public functions ================================================ */
/**
 * Conversion from temperature to resistance.
 * Calculates and returns resistance for given temperature.
 * @param t temperature (in degree Celsius).
 * @return corresponding resistance.
 */
double dNtcTempToRes (double dT, double dCoeff[]);

/**
 * Conversion from resistance to temperature.
 * Calculates and returns temperature for given resistance.
 * @param t resistance (in Ohm).
 * @return corresponding temperature.
 */
double dNtcResToTemp(double dR, double dCoeff[]);

/* ========================================================================== */
#ifdef __cplusplus
}
#endif
#endif /* _NTC_H_ defined */

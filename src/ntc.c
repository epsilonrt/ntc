/**
 * @file ntc.c
 * @brief NTC thermistor library (Implementation)
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
#include <math.h>
#include "ntc.h"

/* constants ================================================================ */
#define TABS (-273.15)

/* private functions ======================================================== */
/*
 * Evaluates p(x) for a polynom p.
 * Calculates the value of polynom p at x accordings to
 * Horners schema.
 * @param p polynom.
 * @param x value to be inserted into the polynom.
 * @return calculated polynom value.
 */
static double
poly(double x, int degree, double p[]) {
  double retval = 0.0;
  int i;

  for (i = degree; i >= 0; i--) {

    retval = retval * x + p[i];
  }
  return retval;
}

/* internal public functions ================================================ */
// -----------------------------------------------------------------------------
double
dNtcTempToRes (double dT, double dCoeff[]) {
  double r;
  double u, v, p, q, b, c, d;

  dT = dT - TABS;
  d = (dCoeff[0] - 1.0 / dT) / dCoeff[3];
  c = dCoeff[1] / dCoeff[3];
  b = dCoeff[2] / dCoeff[3];
  q = 2.0 / 27.0 * b * b * b - 1.0 / 3.0 * b * c + d;
  p = c - 1.0 / 3.0 * b * b;
  v = - pow(q / 2.0 + sqrt(q * q / 4.0 + p * p * p / 27.0), 1.0 / 3.0);
  u =   pow(-q / 2.0 + sqrt(q * q / 4.0 + p * p * p / 27.0), 1.0 / 3.0);
  r  = exp(u + v - b / 3.0);
  return r;
}

// -----------------------------------------------------------------------------
double
dNtcResToTemp(double dR, double dCoeff[])
{
  double ti;

  ti = poly(log(dR), 3, dCoeff);
  ti = 1.0 / ti + TABS;
  return ti;
}

/* ========================================================================== */

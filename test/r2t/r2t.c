/*
 * NTC thermistor library
 * Version 1.0
 * Copyright (C) 2007, 2013 - SoftQuadrat GmbH, Germany
 * Contact: thermistor (at) softquadrat.de
 * Web site: thermistor.sourceforge.net
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 * USA
 */

/** @file r2t.c
 * Program calculating the temperature from resistance value of an NTC thermistor.
 * The Steinhart-Hart polynom allows calculation of absolute temperature
 * from resistance of an NTC thermistor
 *
 * <center><i>
 * 1/t = 1/t<sub>0</sub>
 *       + c<sub>1</sub> &middot; ln(r/r<sub>0</sub>)
 *       + c<sub>2</sub> &middot; ln(r/r<sub>0</sub>)<sup>2</sup>
 *       + c<sub>3</sub> &middot; ln(r/r<sub>0</sub>)<sup>3</sup>
 * </i></center>
 *
 * where (<b>r<sub>0</sub></b>,<b>t<sub>0</sub></b>) is a fixed resistance temperature pair.
 *
 * By substitution
 *
 * <center><i>
 * ln(r/r<sub>0</sub>) = ln(r) - ln(r<sub>0</sub>)
 * </i></center>
 *
 * this leads to a polynom in <b>ln(r)</b>
 *
 * <center><i>
 * 1/t = a<sub>0</sub>
 *       + a<sub>1</sub> &middot; ln r
 *       + a<sub>2</sub> &middot; (ln r)<sup>2</sup>
 *       + a<sub>3</sub> &middot; (ln r)<sup>3</sup>
 * </i></center>
 *
 * The program calculates the calculates the temperature from a given resistance value of an NTC
 * according to that formula.
 */
#include <stdio.h>
#include <stdlib.h>
#include <ntc.h>

/**
 * Main function for conversion from resistance to temperature.
 * Calculates temperature for given resistance and prints result to console.
 * If called without parameters, the user is requested for resistance value,
 * otherwise all arguments are used as resistances, converted to temperature
 * and printed to console.
 * @param argc number of arguments.
 * @param argv argument list.
 * @return 0 indicating no error.
 */
int main(int argc, char *argv[])
{
  int i;
  double r;

  /* Coefficients of Steinhart-Hart polynom. */
  double a[] = {
    4.524024725919526e-004,
    3.934722516618191e-004,
    -7.642331765196044e-006,
    4.048572707661904e-007,
  };

  printf("Thermistor library version 1.0\n");
  printf("Copyright (C) 2007, 2013 - SoftQuadrat GmbH, Germany\n\n");
  if (argc > 1) {

    for (i = 1; i < argc; i++) {

      sscanf(argv[i], "%lf", &r);
      printf("tResistance : %f\tTemperature : %f\n", r, dNtcResToTemp(r, a));
    }
  }
  else {

    int e;

    do {
      printf("Resistance... : ");
      e = scanf("%lf", &r);
    } while (e == 0);
    printf("Temperature.. : %f\n", dNtcResToTemp(r, a));
  }
  return 0;
}

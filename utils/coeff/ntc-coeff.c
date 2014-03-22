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

/** @file ntc-coeff.c
 *
 * Program calculating the coefficients of an extended Steinhart-Hart polynom.
 *
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
 * @image html "rtdiagram.png"
 * @image latex rtdiagram.png
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
 * The program calculates the coefficients <b>a<sub>0</sub></b>, <b>a<sub>1</sub></b>, <b>a<sub>2</sub></b>
 * and <b>a<sub>3</sub></b> from a T-R table, minimizing the sum of squres
 *
 * <center><i>
 * Sum (1/t(r<sub>n</sub>) - 1/t<sub>n</sub>)<sup>2</sup>
 * </i></center>
 *
 */

/***********
* Includes *
***********/
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*********
* Macros *
*********/

/** Dimension of space U. */
#define M 4

/** Absolute Zero. */
#define TABS (-273.15)

/** Maximal line length. */
#define MAX_LENGTH 80

/** Debug flag. */
#define DEBUG 1

/***********
* Typedefs *
***********/
/** Type definition for a polynom. */
typedef double polynom[M];

/************
* Variables *
************/
//* Dimension of space V. */
int n;
static int verbose;

/** Canoncical base of U. */
polynom basis[] = {
  {1.0L, 0.0L, 0.0L, 0.0L},             /* basis[0] = 1                         */
  {0.0L, 1.0L, 0.0L, 0.0L},             /* basis[1] = x                         */
  {0.0L, 0.0L, 1.0L, 0.0L},             /* basis[2] = x^2                       */
  {0.0L, 0.0L, 0.0L, 1.0L}              /* basis[3] = x^3                       */
};

/** x-values calculated from r-values in T-R table. */
double *x;
/** y-values calculated from t-values in T-R table. */
double *y;

/**************
* Prototyping *
**************/
/* Evaluate p(x) */
double value(polynom p, double x);
/* Evaluate [p,q] */
double skalarpoly(polynom p, polynom q);
/* Evaluate p *= fact */
void mult(polynom p, double fact);
/* Evaluate p += fact*q */
void linear(polynom p, polynom q, double fact);
/* Build orthonormal base p from p */
void orthonormal(polynom p[]);
/* Evaluate [p, pf] */
double skalar(polynom p);
/* Evaluate approximating polynom. */
polynom *approx(void);
/* Reads all temperature- resistance pairs from an T-R table file. */
void readtable(const char *filename);
/* Tests the approximation polynom with all t-r pairs. */
void testresult(polynom *erg);
/* Exits with error message in case of errors. */
void errexit(char *format, ...);

void usage (const char * me);

/**
 * Main function for calculating the approximation polynom.
 * Calculation is done in several steps
 *   -# Read all t-r pairs and converting them to x-y values.
 *   -# Evaluate orthonormal base.
 *   -# Evaluate approximation polynom
 *   -# Test approximation polynom
 * @return 0 indicating no error.
 */
int main(int argc, char *argv[])
{
  polynom *erg;
  const char * f;

  printf("Thermistor library version 1.0\n");
  printf("Copyright (C) 2007, 2013 - SoftQuadrat GmbH, Germany\n\n");
  if (argc < 2) {

    usage (argv[0]);
  }
  if (argc > 2) {

    if (strcmp (argv[1], "-v") == 0) {

      f = argv[2];
      verbose = 1;
    }
    else if (strcmp (argv[2], "-v") == 0) {

      f = argv[1];
      verbose = 1;
    }
    else {

      usage (argv[0]);
    }
  }
  else {

    f = argv[1];
    verbose = 0;
  }

  readtable(f);
  orthonormal(basis);
  erg = approx();
  testresult(erg);
  return 0;
}


/************
* Functions *
************/

void
usage (const char * me) {

  fprintf(stderr, "usage : %s [ options ] file  [ options ]\n", me);
  fprintf(stderr,
  "Program calculating the coefficients of an extended Steinhart-Hart polynom.\n"
  " The Steinhart-Hart polynom allows calculation of absolute temperature\n"
  " from resistance of an NTC thermistor\n\n");

  fprintf(stderr,"valid options are :\n");
  fprintf(stderr,
  "  -v\tenables verbose output\n");
  exit(EXIT_FAILURE);
}

/**
 * Evaluates p(x) for a polynom p.
 * Calculates the value of polynom p at x accordings to
 * Horners schema.
 * @param p polynom.
 * @param x value to be inserted into the polynom.
 * @return calculated polynom value.
 */
double value(polynom p, double x)
{
  int i;
  double retval = 0.0L;

  for (i = M - 1; i >= 0; i--)
    retval = retval * x + p[i];
  return retval;
}

/**
 * Evaluates [p,q] for two polynoms p and q.
 * Calculates the scalar product of two polynoms p and q.
 * This is defined as sum
 * <center>[p, q] := Sum p(x<sub>i</sub>) * q(x<sub>i</sub>) &uuml;ber i = 0, .., N - 1</center>
 * @param p first polynom.
 * @param q second polynom.
 * @return calculated scalar product.
 */
double skalarpoly(polynom p, polynom q)
{
  int i;
  double retval = 0.0L;

  for (i = 0; i < n; i++)
    retval += value(p, x[i]) * value(q, x[i]);
  return retval;
}

/**
 * Evaluates p *= fact for a polynom p and a factor fact.
 * Multiplies p with a factor fact.
 * @param p polynom.
 * @param fact factor.
 */
void mult(polynom p, double fact)
{
  int i;

  for (i = 0; i < M; i++)
    p[i] *= fact;
}

/**
 * Evaluates p += fact*q for two polynoms p and q and a factor fact.
 * Adds a multiple of a polynom q to polynom p.
 * @param p polynom, to be added to.
 * @param q qolynom, added with a factor.
 * @param fact factor.
 */
void linear(polynom p, polynom q, double fact)
{
  int i;

  for (i = 0; i < M; i++)
    p[i] += q[i] * fact;
}

/**
 * Converts a base to an orthonormal base.
 * @param p base in form of an array of polynoms.
 */
void orthonormal(polynom p[])
{
  int i, j;
  double fact, norm;

  if (verbose)
  {
    printf("function orthonormal\n");
    printf("====================\n");
  }
  for (i = 0; i < M; i++) {
    if (verbose)
      printf("Evaluating polynom number %d\n", i);
    for (j = 0; j < i; j++) {
      fact = skalarpoly(p[i], p[j]);
      linear(p[i], p[j], -fact);
    }
    norm = skalarpoly(p[i], p[i]);
    mult(p[i], 1.0L / sqrt(norm));
    if (verbose)
    {
      printf("Polynom %d: ", i);
      for (j = 0; j < M; j++)
      {
        printf("%f ", p[i][j]);
      }
      printf("\n");
    }
  }
  if (verbose) {
    printf("Testing orthonormal base\n");
    for (i = 0; i < M; i++) {
      for (j = 0; j <= i; j++)
        printf("%.15f ", skalarpoly(basis[i], basis[j]));
      printf("\n");
    }
    printf("\n");
  }
}

/**
 * Evaluates [p, p<sub>f</sub>] for given polynom p and solving polynom p<sub>f</sub>.
 * @param p polynom.
 * @return calculated scalar product.
 */
double skalar(polynom p)
{
  int i;
  double retval = 0.0L;

  for (i = 0; i < n; i++)
    retval += y[i] * value(p, x[i]);
  return retval;
}

/**
 * Evaluate approximation polynom u<sub>f</sub>.
 * @return approximation polynom u<sub>f</sub>.
 */
polynom *approx(void)
{
  int i;
  double fact;
  polynom *erg = malloc(sizeof(polynom));

  if (verbose)
  {
    printf("function approx\n");
    printf("===============\n");
  }
  for (i = 0; i < M; i++)
    (*erg)[i] = 0.0L;
  for (i = 0; i < M; i++) {
    if (verbose)
      printf("Approximating with polynom number %d\n", i);
    fact = skalar(basis[i]);
    linear(*erg, basis[i], fact);
  }
  printf("Steinhart-Hart coefficients\n");
  for (i = 0; i < M; i++)
    printf("a[%d] = %.15e\n", i, (*erg)[i]);
  if (verbose)
    printf("\n");
  return erg;
}

/**
 * Reads all temperature- resistance pairs from an T-R table file.
 * The resulting values are converted from t-r pairs to x-y pairs
 * where x = ln(r) and y = 1 / (t - TABS).
 * @param filename name of file with all t-r pairs.
 */
void readtable(const char *filename)
{
  double temp;
  double res;
  FILE *fr;
  char *line;
  double *t, *r;
  int num;
  int i;

  if (verbose)
  {
    printf("function readtable\n");
    printf("==================\n");
  }
  num = 50;
  t = malloc(num * sizeof(double));
  r = malloc(num * sizeof(double));
  fr = fopen(filename, "r");
  if (fr == NULL)
    errexit("Cannot find file %s", filename);
  line = malloc(MAX_LENGTH);
  for (i = 0; fgets(line, MAX_LENGTH - 1, fr) != NULL; i++)
  {
    if (i >= num)
    {
      num *= 2;
      t = realloc(t, num * sizeof(double));
      r = realloc(r, num * sizeof(double));

    }
    sscanf(line, "%lf\t%lf", &temp, &res);
    t[i] = temp;
    r[i] = res;
    if (verbose)
      printf("t=%8.2f\tr=%8.2f\n", t[i], r[i]);
  }
  if (verbose)
    printf("\n");
  n = i;
  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  for (i = 0; i < n; i++) {
    x[i] = log(r[i]);
    y[i] = 1.0 / (t[i] - TABS);
    if (verbose)
      printf("x=%8.2f\ty=%9.4f\n", x[i], y[i]);
  }
  if (verbose)
    printf("\n");
}

/**
 * Tests the approximation polynom with all t-r pairs.
 * Prints out all calculated values and the maximal error.
 * The function will do nothing, if verbose mode is off
 * (Macro DEBUG = 0).
 * @param erg approximation polynom.
 */
void testresult(polynom *erg)
{
  double maxerr, val1, val2, err, temp = 0;
  int i;

  if (verbose)
  {
    printf("function testresult\n");
    printf("===================\n");
    maxerr = 0.0;
    for (i = 0; i < n; i++) {
      val1 = 1.0 / value(*erg, x[i]) + TABS;
      val2 = 1.0 / y[i] + TABS;
      err = val1 - val2;
      if (err < 0)
        err = -err;
      printf("%8.3f\t%8.1f\t%8.1f\n", val1, exp(x[i]), val2);
      if (err > maxerr)
      {
        temp = val2;
        maxerr = err;
      }
    }
    printf("\n");
    printf("Maximal error=%7.5f at temperature=%5.1f\n", maxerr, temp);
    printf("\n");
  }
}

/**
 * Exits with error message in case of errors.
 * @param format of error message.
 */
void errexit(char *format, ...)
{
  va_list ap;

  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  exit(EXIT_FAILURE);
}


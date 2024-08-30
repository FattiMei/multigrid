/******************************************************************************
 *                      Code generated with SymPy 1.12.1                      *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                       This file is part of 'project'                       *
 ******************************************************************************/
#include "symbolic.h"
#include <math.h>

int symbolic() {

   int symbolic_result;
   symbolic_result = 0;
   return symbolic_result;

}

double solution_1d(double x) {

   double solution_1d_result;
   solution_1d_result = exp(-pow(x, 2)*sin((1.0/10.0)*x));
   return solution_1d_result;

}

double forcing_term_1d(double x) {

   double forcing_term_1d_result;
   forcing_term_1d_result = -((1.0/100.0)*pow(x, 2)*pow(x*cos((1.0/10.0)*x) + 20*sin((1.0/10.0)*x), 2) + (1.0/100.0)*pow(x, 2)*sin((1.0/10.0)*x) - 2.0/5.0*x*cos((1.0/10.0)*x) - 2*sin((1.0/10.0)*x))*exp(-pow(x, 2)*sin((1.0/10.0)*x));
   return forcing_term_1d_result;

}

double solution_2d(double x, double y) {

   double solution_2d_result;
   solution_2d_result = exp(x)*exp(-2*y);
   return solution_2d_result;

}

double forcing_term_2d(double x, double y) {

   double forcing_term_2d_result;
   forcing_term_2d_result = -5*exp(x)*exp(-2*y);
   return forcing_term_2d_result;

}

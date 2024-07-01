#include "bessel.h"

#ifdef USE_BARE_C_BESSEL

#ifndef KMAX
#define KMAX 2000
#endif

static inline double bessel_i0( double x );
static inline double bessel_k0( double x );

/*
 * taken from https://doi.org/10.1103/PhysRevB.86.115447
 * =====================================================
 * Robert E. Throckmorton and Oskar Vafek
 * Phys. Rev. B 86, 115447 (Published 26 September 2012)
*/
double gate_screened_coulomb_shape( double dist_over_gate_dist ) {
    double val = 0.0;
    for (int k=0; k<KMAX; ++k)
        val += bessel_k0( (2.*k + 1.)*M_PI * dist_over_gate_dist );
    return val;
}

// Implementation based on [Numerical Recipes, Press et. al.].

static inline double bessel_i0( double x ) {
   double ax, res;
   double y;
   if ((ax=fabs(x)) < 3.75) {
      y = x/3.75;
      y = y*y;
      res = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 +
                    y*(0.360768e-1 + y*0.45813e-2)))));
   } else {
      y = 3.75/ax;
      res = (exp(ax)/sqrt(ax)) * (0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 +
                    y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 +
                    y*(0.2635537e-1 + y*(-0.1647633e-1 + y*0.392377e-2))))))));
   }
   return res;
}

static inline double bessel_k0( double x ) {
   double y, res;
   if (x <= 2.0) {
      y = x*x/4.0;
      res = (-log(x/2.0) * bessel_i0(x)) + (-0.57721566 + y*(0.42278420 +
                    y*(0.23069756 + y*(0.3488590e-1 + y*(0.262698e-2 +
                    y*(0.10750e-3 + y*0.74e-5))))));
   } else {
      y = 2.0/x;
      res = (exp(-x)/sqrt(x)) * (1.25331414 + y*(-0.7832358e-1 +
                    y*(0.2189568e-1 + y*(-0.1062446e-1 + y*(0.587872e-2 +
                    y*(-0.251540e-2 + y*0.53208e-3))))));
   }
   return res;
}

#endif // USE_BARE_C_BESSEL

#include "bessel.h"
#include <cmath>

#ifndef KMAX
#define KMAX 2000
#endif

/*
 * taken from https://doi.org/10.1103/PhysRevB.86.115447
 * =====================================================
 * Robert E. Throckmorton and Oskar Vafek
 * Phys. Rev. B 86, 115447 (Published 26 September 2012)
*/
double gate_screened_coulomb_shape( double dist_over_gate_dist ) {
    double val = 0.0;
    for (int k=0; k<KMAX; ++k)
        val += std::cyl_bessel_k( 0, (2.*k + 1.)*M_PI * dist_over_gate_dist );
    return val;
}


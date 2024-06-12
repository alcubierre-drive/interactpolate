#include "interactpolate.h"
#include "bessel.h"
#include "spline.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#ifndef POW2
#define POW2(x) ((x)*(x))
#endif

#ifndef INTERACTION_INTERP_NUM
#define INTERACTION_INTERP_NUM 2000
#endif

#ifndef INTERACTION_ZERO_SHIFT
#define INTERACTION_ZERO_SHIFT 1.e-9
#endif

#ifndef INTERACTION_GATE_DIST_ZERO
#define INTERACTION_GATE_DIST_ZERO 2.e-4
#endif

#ifndef wrn_printf
#define wrn_printf( ... ) fprintf( stderr, __VA_ARGS__ )
#endif

static int _nthr = -1;
void interactpolate_init_set_num_threads( int nthr ) {
    _nthr = nthr;
}

struct interactpolate_t {
    double ohno_dist;
    double gate_dist;
    double U;

    double x[INTERACTION_INTERP_NUM];
    double y[INTERACTION_INTERP_NUM];

    cspline_t* spline;

    double interp_maxdist;
    double interp_maxdist_lambda;
    double interp_maxdist_V0;
};

interactpolate_t* interactpolate_init( double U, double ohno_dist, double gate_dist ) {
    if (gate_dist < 0) gate_dist = fabs(gate_dist);
    if (ohno_dist < 0) ohno_dist = fabs(ohno_dist);

    if (gate_dist < ohno_dist) {
        gate_dist = 5.*ohno_dist;
        wrn_printf( "cannot do interactpolate where gate is closer than Ohno distance\n" );
        wrn_printf( "reset gate_dist=%.2e (ohno_dist=%.2e)\n", gate_dist, ohno_dist );
    }

    interactpolate_t* ii = calloc( 1, sizeof(interactpolate_t) );

    ii->U = U;
    ii->ohno_dist = ohno_dist;
    ii->gate_dist = gate_dist;
    ii->interp_maxdist = gate_dist * 5.;

    double R0 = ii->gate_dist * INTERACTION_GATE_DIST_ZERO; // here we should have the 1/R behavior
    double G0 = gate_screened_coulomb_shape( R0 / ii->gate_dist ) / ii->gate_dist; // this shuold be alpha/R0
    double alpha = G0 * R0;

    int nthr = _nthr <= 0 ? omp_get_max_threads() : _nthr;
    #pragma omp parallel for num_threads(nthr)
    for (int i=0; i<INTERACTION_INTERP_NUM; ++i) {
        double xrel = (double)i/(double)(INTERACTION_INTERP_NUM - 1);
        ii->x[i] = xrel * ii->interp_maxdist + INTERACTION_ZERO_SHIFT; // for safety
        double G_r = gate_screened_coulomb_shape( ii->x[i] / ii->gate_dist ) / ii->gate_dist / alpha;
        ii->y[i] = ii->U / sqrt( 1. + 1./(POW2(G_r * ii->ohno_dist)) );
    }

    ii->spline = cspline_init( INTERACTION_INTERP_NUM, ii->x, ii->y );

    ii->interp_maxdist = ii->x[INTERACTION_INTERP_NUM-1];

    double r1 = ii->x[INTERACTION_INTERP_NUM-2], r2 = ii->x[INTERACTION_INTERP_NUM-1],
           V1 = ii->y[INTERACTION_INTERP_NUM-2], V2 = ii->y[INTERACTION_INTERP_NUM-1];

    double lambda = log( V1*sqrt(r2) / (V2*sqrt(r1)) ) / (r2 - r1);
    double V01 = V1 / exp(-lambda * r1) * sqrt(r1);
    double V02 = V2 / exp(-lambda * r2) * sqrt(r2);

    ii->interp_maxdist_lambda = lambda;
    ii->interp_maxdist_V0 = 0.5*(V01+V02);

    return ii;
}

interactpolate_t* interactpolate_abinit( double U, double alpha, double eps, double gate_dist ) {
    if (gate_dist < 0) gate_dist = fabs(gate_dist);
    if (alpha < 0) alpha = fabs(alpha);
    if (eps < 0) eps = fabs(eps);

    interactpolate_t* ii = calloc( 1, sizeof(interactpolate_t) );

    double V0 = alpha / (eps * gate_dist),
           ohno_dist = alpha / (eps * U);

    ii->U = U;
    ii->ohno_dist = ohno_dist;
    ii->gate_dist = gate_dist;
    ii->interp_maxdist = gate_dist * 5.;

    int nthr = _nthr <= 0 ? omp_get_max_threads() : _nthr;
    #pragma omp parallel for num_threads(nthr)
    for (int i=0; i<INTERACTION_INTERP_NUM; ++i) {
        double xrel = (double)i/(double)(INTERACTION_INTERP_NUM - 1);
        ii->x[i] = xrel * ii->interp_maxdist - INTERACTION_ZERO_SHIFT; // yes, it's the other way 'round here
        double ii_kernel = gate_screened_coulomb_shape( sqrt(POW2(ii->x[i])+POW2(ii->ohno_dist))/ii->gate_dist );
        ii->y[i] = V0 * 4 * ii_kernel;
    }

    ii->spline = cspline_init( INTERACTION_INTERP_NUM, ii->x, ii->y );

    ii->interp_maxdist = ii->x[INTERACTION_INTERP_NUM-1];

    double r1 = ii->x[INTERACTION_INTERP_NUM-2], r2 = ii->x[INTERACTION_INTERP_NUM-1],
           V1 = ii->y[INTERACTION_INTERP_NUM-2], V2 = ii->y[INTERACTION_INTERP_NUM-1];

    double lambda = log( V1*sqrt(r2) / (V2*sqrt(r1)) ) / (r2 - r1);
    double V01 = V1 / exp(-lambda * r1) * sqrt(r1);
    double V02 = V2 / exp(-lambda * r2) * sqrt(r2);

    ii->interp_maxdist_lambda = lambda;
    ii->interp_maxdist_V0 = 0.5*(V01+V02);

    return ii;
}

double interactpolate_get_safe( const interactpolate_t* ii, double r ) {
    if (r < 0) r = fabs(r);
    if (r < 1.e-11)
        return ii->U;
    else if (r>=ii->interp_maxdist)
        return ii->interp_maxdist_V0 * exp( -ii->interp_maxdist_lambda * r ) / sqrt(r);
    else {
        return cspline_eval( ii->spline, r );
    }
}

double interactpolate_get( const interactpolate_t* ii, double r ) {
    if (r>=ii->interp_maxdist)
        return ii->interp_maxdist_V0 * exp( -ii->interp_maxdist_lambda * r ) / sqrt(r);
    else
        return cspline_eval( ii->spline, r );
}

void interactpolate_free( interactpolate_t* ii ) {
    cspline_free( ii->spline );
    free( ii );
}



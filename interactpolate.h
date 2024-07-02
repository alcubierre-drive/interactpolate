#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct interactpolate_t interactpolate_t;

// default: don't use any constraint
void interactpolate_init_set_num_threads( int nthr );
// default: use 2048 interpolation points
void interactpolate_init_set_num_points( int npts );

interactpolate_t* interactpolate_init( double U, double ohno_dist, double gate_dist );
interactpolate_t* interactpolate_abinit( double U, double alpha, double eps, double gate_dist );

double interactpolate_get_safe( const interactpolate_t* ii, double r );
double interactpolate_get( const interactpolate_t* ii, double r );
void interactpolate_free( interactpolate_t* ii );

#ifdef __cplusplus
}
#endif

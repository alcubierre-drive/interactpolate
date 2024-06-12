#pragma once

#ifdef __cplusplus
#include <cstdint>
extern "C" {
#else
#include <stdint.h>
#endif

typedef struct cspline_t cspline_t;

cspline_t* cspline_init( int64_t N, double* x, double* y );

double cspline_eval( const cspline_t* s, double x_new );
double cspline_eval_deriv( const cspline_t* s, double x_new );
double cspline_eval_integral( const cspline_t* s, double x_new );

void cspline_free( cspline_t* s );

#ifdef __cplusplus
}
#endif

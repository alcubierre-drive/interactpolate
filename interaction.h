#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct interaction_t interaction_t;

// default: don't use any constraint
void interaction_init_set_num_threads( int nthr );

interaction_t* interaction_init( double U, double ohno_dist, double gate_dist );
double interaction_get( const interaction_t* ii, double r );
void interaction_free( interaction_t* ii );

#ifdef __cplusplus
}
#endif

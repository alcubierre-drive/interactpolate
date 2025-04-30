#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** opaque handle for interpolation */
typedef struct interactpolate_t interactpolate_t;

/**
 * Set the number of threads in *generation* of the interpolation points
 * (default given by OpenMP). Not thread-safe, use single-threaded *before*
 * calling :c:func:`interactpolate_abinit` if non-default value desired.
 */
void interactpolate_init_set_num_threads( int nthr );

/**
 * Set the number of points for interpolation (default: 2048). Not thread-safe,
 * use single-threaded *before* calling :c:func:`interactpolate_abinit` if
 * non-default value desired.
 */
void interactpolate_init_set_num_points( int npts );

/**
 * Initialize interpolation handle given Hubbard interaction :math:`U`, fine
 * structure constant :math:`\alpha = 14.40\,\mathrm{eV}\text{Å}`, dielectric
 * constant :math:`\epsilon` (use :math:`\epsilon=1` for vacuum) and gate
 * distance :math:`\xi`. Returned handle must be freed with
 * :c:func:`interactpolate_free`.
 */
interactpolate_t* interactpolate_abinit( double U, double alpha, double eps, double gate_dist );

/**
 * Initialize interpolation handle, *deprecated* version that uses Ohno
 * distance, :math:`a`, Hubbard interaction :math:`U`, and gate distance
 * :math:`\xi`. Can fail in some instances of badly chosen parameters. Use
 * :c:func:`interactpolate_abinit` whenever possible!
 */
interactpolate_t* interactpolate_init( double U, double ohno_dist, double gate_dist );

/**
 * Get the value of the Coulomb interaction :math:`V(r)` for a distance
 * :math:`r`. Thread-safe (i.e. can be used inside a parallel region). Does not
 * check for domain of :math:`r > 0`.
 */
double interactpolate_get( const interactpolate_t* ii, double r );

/**
 * same as :c:func:`interactpolate_get` but with domain check (and correction)
 * for :math:`r`
 */
double interactpolate_get_safe( const interactpolate_t* ii, double r );

/** free all resources associated to a handle */
void interactpolate_free( interactpolate_t* ii );

#ifdef __cplusplus
}
#endif

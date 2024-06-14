#include <interactpolate.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

double wtime( void ) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1.e+9;
}

int main() {
    double U = 4.0,
           alpha = 14.40,
           epsilon = 12.0,
           xi = 200.0;

    interactpolate_t* ii = interactpolate_abinit( U, alpha, epsilon, xi );

    size_t N = 100000;
    double xmax = 1000.0;
    double* x = calloc( N, sizeof(double) );
    double* y = calloc( N, sizeof(double) );

    for (size_t i=0; i<N; ++i)
        x[i] = (double)i/(double)N * xmax;

    double tick, tock;

    printf( "starting interactpolate_get get serial\n" );
    tick = wtime();
    for (size_t i=0; i<N; ++i)
        y[i] = interactpolate_get( ii, x[i] );
    tock = wtime();
    printf( "%.4fs\n", tock-tick );

    printf( "starting interactpolate_get parallel (%i threads)\n", omp_get_max_threads() );
    tick = wtime();
    #pragma omp parallel for
    for (size_t i=0; i<N; ++i)
        y[i] = interactpolate_get( ii, x[i] );
    tock = wtime();
    printf( "%.4fs\n", tock-tick );

    printf( "timing done\n" );

    FILE* f = fopen("test.bin", "w");
    fwrite( x, sizeof(x[0]), N, f );
    fwrite( y, sizeof(y[0]), N, f );
    fclose( f );
    printf( "output to test.bin\n" );

    free( x );
    free( y );

    interactpolate_free( ii );
}

#include "interactpolate.h"
#include <stdio.h>

int main() {
    interactpolate_t* ii = interactpolate_abinit( 8.0, 14.0, 8.0, 20.0 );

    for (double x=0.0; x<=1000.0; x+=1.e-1)
        printf( "%.5e %.5e\n", x, interactpolate_get( ii, x ) );

    interactpolate_free( ii );
}

#include "spline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// adopted from
// https://medium.com/@martinmikkelsen/cubic-splines-in-the-c-language-4463fe9a3a83

struct cspline_t {
    int64_t N;
    double *x, *y;
    double *b, *c, *d;

    int reg;
    double x0, xN_x0;
};

static inline int64_t binsearch( int64_t N, double* x, double x_new ) {
    int64_t i = 0;
    int64_t j = N - 1;

    while (j - i > 1) {
        int64_t mid = (i + j) / 2;
        if (x_new > x[mid]) {
            i = mid;
        } else {
            j = mid;
        }
    }

    return i;
}

// use this 'search' in case spacing is regular. much faster.
static inline int64_t binsearch_reg( const cspline_t* s, double x_new ) {
    return (x_new - s->x0) / s->xN_x0 * s->N;
}

cspline_t* cspline_init( int64_t N, double* x, double* y ) {
    cspline_t* s = (cspline_t*)malloc(sizeof(cspline_t));
    s->N = N;
    s->x = calloc(N, sizeof(double));
    s->y = calloc(N, sizeof(double));
    s->b = calloc(N, sizeof(double));
    s->c = calloc((N - 1), sizeof(double));
    s->d = calloc((N - 1), sizeof(double));

    // check if spacing is regular and then set s->reg
    double dx_regular = x[1] - x[0];
    int not_regular = 0;
    for (int64_t i = 0; i < N; ++i) {
        s->x[i] = x[i];
        s->y[i] = y[i];
        if (i > 0)
            not_regular += (fabs(dx_regular - (x[i] - x[i-1])) > 1.e-9);
    }
    s->reg = !not_regular;
    s->x0 = x[0];
    s->xN_x0 = x[N-1] - x[0];

    int64_t i;
    double* dx = calloc(N - 1, sizeof(double));
    double* p = calloc(N - 1, sizeof(double));

    for (i = 0; i < N - 1; ++i) {
        dx[i] = x[i + 1] - x[i];
        p[i] = (y[i + 1] - y[i]) / dx[i];
    }
    double* D = calloc(N, sizeof(double));
    double* B = calloc(N, sizeof(double));
    double* Q = calloc(N - 1, sizeof(double));

    D[0] = 2;
    D[N - 1] = 2;
    B[0] = 3 * (p[0]);
    B[N - 1] = 3 * (p[N - 2]);
    Q[0] = 1;

    for (i = 0; i < N - 2; ++i) {
        D[i + 1] = 2 * dx[i] / dx[i + 1] + 2;
        B[i + 1] = 3 * (p[i] + p[i + 1] * dx[i] / dx[i + 1]);
        Q[i + 1] = dx[i] / dx[i + 1];
    }
    for (i = 1; i < N; ++i) {
        D[i] -= Q[i - 1] / D[i - 1];
        B[i] -= B[i - 1] / D[i - 1];
    }
    s->b[N - 1] = B[N - 1] / D[N - 1];

    for (i = N - 2; i >= 0; --i) {
        s->b[i] = (B[i] - Q[i] * s->b[i + 1]) / D[i];
    }

    for (i = 0; i < N - 1; ++i) {
        s->c[i] = (-2 * s->b[i] - s->b[i + 1] + 3 * p[i]) / dx[i];
        s->d[i] = (s->b[i] + s->b[i + 1] - 2 * p[i]) / dx[i] / dx[i];
    }

    free(dx);
    free(p);
    free(D);
    free(B);
    free(Q);

    return s;
}

#ifndef POW2
#define POW2( x ) ((x)*(x))
#endif
#ifndef POW3
#define POW3( x ) ((x)*(x)*(x))
#endif

double cspline_eval( const cspline_t* s, double x_new ) {
    int64_t i = s->reg ? binsearch_reg(s, x_new) : binsearch(s->N, s->x, x_new);
    double dx = x_new - s->x[i];
    double y_new = s->y[i] + s->b[i] * dx + s->c[i] * POW2(dx) + s->d[i] * POW3(dx);
    return y_new;
}

double cspline_eval_deriv( const cspline_t* s, double x_new ) {
    int64_t i = s->reg ? binsearch_reg(s, x_new) : binsearch(s->N, s->x, x_new);
    double dx = x_new - s->x[i];
    double deriv = s->b[i] + 2 * s->c[i] * dx + 3 * s->d[i] * pow(dx, 2);
    return deriv;
}

double cspline_eval_integral( const cspline_t* s, double x_new ) {
    int64_t i = s->reg ? binsearch_reg(s, x_new) : binsearch(s->N, s->x, x_new);
    double integral = 0;
    for (int64_t j = 0; j < i; ++j) {
        double dx = s->x[j + 1] - s->x[j];
        integral += s->y[j] * dx + s->b[j] * pow(dx, 2) / 2 +
                    s->c[j] * pow(dx, 3) / 3 + s->d[j] * pow(dx, 4) / 4;
    }
    double dx = x_new - s->x[i];
    integral += s->y[i] * dx + s->b[i] * pow(dx, 2) / 2 +
                s->c[i] * pow(dx, 3) / 3 + s->d[i] * pow(dx, 4) / 4;
    return integral;
}

void cspline_free( cspline_t* s ) {
    free(s->x);
    free(s->y);
    free(s->b);
    free(s->c);
    free(s->d);
    free(s);
}

/*
int main() {
    double x[100];
    double y[100];
    for (int64_t i=0; i<100; ++i) {
        x[i] = (double)i/100.;
        y[i] = sin(3*x[i]);
    }
    cspline_t* c = cspline_init(100, x, y);

    for (double xx=0.1; xx<0.9; xx+=1.e-6)
        printf( "%.5e %.5e %.5e\n", xx, sin(3*xx), cspline_eval(c, xx) );
    cspline_free(c);
    return 0;
}
*/

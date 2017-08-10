#ifndef POLY_H
#define POLY_H

#if !defined(NAME_H)
#warning "poly.h" requires "name.h"
#endif

#define lagrange_eval  PREFIXED_NAME(lagrange_eval )
#define lagrange_setup PREFIXED_NAME(lagrange_setup)
#define gll_lag_setup  PREFIXED_NAME(gll_lag_setup )
#define gauss_nodes    PREFIXED_NAME(gauss_nodes   )
#define gauss_quad     PREFIXED_NAME(gauss_quad    )
#define lobatto_nodes  PREFIXED_NAME(lobatto_nodes )
#define lobatto_quad   PREFIXED_NAME(lobatto_quad  )

/*--------------------------------------------------------------------------
   Quadrature Nodes and Weights Calculation

    Gauss   -> Gauss-Legendre quadrature (open)
    Lobatto -> Gauss-Lobatto-Legendre quadrature (closed at both ends)
   
   the _quad functions compute both nodes and weights
  --------------------------------------------------------------------------*/

void   gauss_nodes(double *restrict z, int n); /* n nodes (order = 2n-1) */
void lobatto_nodes(double *restrict z, int n); /* n nodes (order = 2n-3) */

void   gauss_quad(double *restrict z, double *restrict w, int n);
void lobatto_quad(double *restrict z, double *restrict w, int n);

/*--------------------------------------------------------------------------
   Lagrangian basis function evaluation
   
   Usage:
   
   double z[N] = ..., x = ...; // nodes and evaluation point
   double w[N]; // Coefficients for Lagrange polynomial
   double p[3*N];
   
   lagrange_setup(w, z,n);  // Compute w
   
   int d = ...; // 0, 1, or 2  --- the highest derivative to compute
   lagrange_eval(p, z,w,N,d, x);
   // now p[i] = h_i(x), 0 <= i < N 
   // if d>=1, p[N+i] = h_i'(x)
   // if d>=2, p[2*N+i] = h_i''(x)
  --------------------------------------------------------------------------*/

void lagrange_setup(double *restrict w, const double *restrict z, unsigned n);

void lagrange_eval(
  double *restrict p,
  const double *restrict z, const double *restrict w, unsigned n,
  int der, double x);

/*--------------------------------------------------------------------------
   Lagrangian basis function evaluation
   for Gauss-Lobatto-Legendre quadrature nodes
   
   Usage:
   double x;  // evaluation point
   double p[3*N];

   gll_lag_fun *gll_lag = gll_lag_setup(N);
   
   int d = ...; // 0, 1, or 2  --- the highest derivative to compute
   gll_lag(p, N,d, x);
   // now p[i] = h_i(x), 0 <= i < N 
   // if d>=1, p[N+i] = h_i'(x)
   // if d>=2, p[2*N+i] = h_i''(x)
  --------------------------------------------------------------------------*/

typedef void gll_lag_fun(double *restrict p, unsigned n, int der, double x);

gll_lag_fun *gll_lag_setup(unsigned n);

#endif


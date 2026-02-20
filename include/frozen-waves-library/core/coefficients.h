/*
  Frozen Waves Library: A C library providing routines for computing quantities
  related to Frozen Waves

  File: include/core/coefficients.h
  Language standards: C99
  References: include/frozen-waves-library/references.txt
  License: include/frozen-waves-library/license.txt
  Repository: <https://github.com/jodesarro/frozen-waves-library>

  Description: Callable functions related to Frozen Wave coefficients.
*/

#ifndef FROZEN_WAVES_LIBRARY_COEFFICIENTS_H
#define FROZEN_WAVES_LIBRARY_COEFFICIENTS_H

#include "../impl/api_impl_.h"
#include <complex.h> /* For 'double complex' type */
#include <math.h>    /* Maths and constants */
#include <stdbool.h> /* For 'bool' type */

#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

/*
  Evaluates the A_q coefficients of a Frozen Wave (FW), resistant or not to a
  possible absorption due to a medium with complex refractive index. See Refs
  [1], [2], [3] and [4].

  Parameters:
  - N, parameter N of the FW.
  - L, parameter L of the FW.
  - beta, array of size iqmax = 2 * N + 1 containing the longitudinal wavenumber
  beta_q, where beta_q -> beta[iq], q -> iq - N, and 0<=iq<iqmax.
  - F, array of size izmax containing the morphological function F(z), where
  F(z) -> F[iz] and z -> iz * L / (izmax - 1), and 0<=iz<izmax.
  - izmax, size for the F[iz] array, where 0<=iz<izmax.
  - A, array of size iqmax = 2 * N + 1 to output the A_q coefficients, where A_q
  -> A[iq], q -> iq - N, and 0<=iq<iqmax.
  - absorption_resistant, if false, the FW is not absorption-resistant.

  Implementation: The A_q = int_0^L F(z) exp(-betabar_0) exp(j*2*pi*q*z/L) dz,
  where betabar_0 = imag(beta_0) if absorption_resistant = true or betabar_0 = 0
  if absorption_resistant = false, is computed by means of an approximation of
  the integral by the trapezoidal method for equally spaced z for a given
  morphologial function F(z). Notice that for purely real refractive indices,
  the absorption-resistant concerns may be ignored.
*/
FROZEN_WAVES_LIBRARY_API_IMPL_
void fw_A_coefficient(int N, double L, const double complex *beta,
                      const double *F, int izmax, double complex *A,
                      bool absorption_resistant)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Constants */
  double betabar_0 = absorption_resistant ? cimag(beta[N]) : 0.0;
  double c_2pi_L = 2.0 * M_PI / L;
  int iqmax = 2 * N + 1;

  for (int iq = 0; iq < iqmax; iq++) {

    double q = (double)(iq - N);

    /* Initial trapezoidal contribution (endpoints) */
    double complex aux = 0.5 * (F[izmax - 1] * exp(-betabar_0 * L) + F[0]);

    /* Sum over intermediate points */
    for (int iz = 1; iz < izmax - 1; iz++) {
      double z = L * (double)iz / (double)(izmax - 1);
      aux += F[iz] * cexp(I * c_2pi_L * q * z) * exp(-betabar_0 * z);
    }

    /* Normalize */
    A[iq] = aux / (double)(izmax - 1);
  }
}
#else
    ;
#endif

/*
  Evaluates the A_q,p coefficients of a Frozen Wave (FW) 2D made by P FWs,
  restricted to the case where all FWs share the same N, Q and P (i.e., for 1 <=
  p <= P, N_p = N, Q_p = Q, and L_p = L), resistant or not to a possible
  absorption due to a medium with complex refractive index. See Refs. [5], [6]
  and [7].

  Parameters:
  - P, number P of FWs in the transverse direction y.
  - N, parameter N of the FW.
  - L, parameter L of the FW.
  - beta, array of size iqmax = 2 * N + 1 containing the longitudinal wavenumber
  beta_q, where beta_q -> beta[iq], q -> iq - N, and 0<=iq<iqmax.
  - F, array of size iymax*izmax containing the morphological function F(y,z),
  where F(y,z) -> F[iz + izmax*iy], y -> iy, z -> iz * L / (izmax - 1),
  0<=iy<iymax, and 0<=iz<izmax.
  - iymax, size contribution in iy for the F[iz + izmax*iy] array, where
  0<=iy<iymax.
  - izmax, size contribution in iz for the F[iz + izmax*iy] array, where
  0<=iz<izmax.
  - A, array of size ipmax*iqmax = (P) * (2 * N + 1) to output the A_q,p
  coefficients, where A_q,p -> A[iq + iqmax*ip], p -> ip + 1, q -> iq - N,
  0<=ip<ipmax, and 0<=iq<iqmax.
  - absorption_resistant, if false, the FW is not absorption-resistant.

  Implementation: The A_q = int_0^L F(y=y_0p, z) exp(-betabar_0)
  exp(j*2*pi*q*z/L) dz, where y_0p -> iy * (ipmax - 1) / (iymax - 1), betabar_0
  = imag(beta_0) if absorption_resistant = true or betabar_0 = 0 if
  absorption_resistant = false, is computed by means of an approximation of the
  integral by the trapezoidal method for equally spaced z for a given
  morphologial function F(y,z). Notice that for purely real refractive indices,
  the absorption-resistant concerns may be ignored.
*/
FROZEN_WAVES_LIBRARY_API_IMPL_
void fw2d_A_coefficient_restricted(int P, int N, double L,
                                   const double complex *beta, const double *F,
                                   int iymax, int izmax, double complex *A,
                                   bool absorption_resistant)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Constants */
  double betabar_0 = absorption_resistant ? cimag(beta[N]) : 0.0;
  double c_2pi_L = 2.0 * M_PI / L;
  int iqmax = 2 * N + 1;
  int ipmax = P;

  for (int ip = 0; ip < ipmax; ip++) {

    for (int iq = 0; iq < iqmax; iq++) {

      double q = (double)(iq - N);

      /* Map ip-index to iy-index (linear scaling) */
      int iy =
          (int)floor((double)ip * (double)(iymax - 1) / (double)(ipmax - 1));

      /* Initial trapezoidal contribution (endpoints) */
      double complex aux =
          0.5 * (F[(izmax - 1) + izmax * iy] * exp(-betabar_0 * L) +
                 F[0 + izmax * iy]);

      /* Sum over intermediate z points */
      for (int iz = 1; iz < izmax - 1; iz++) {
        double z = L * (double)iz / (double)(izmax - 1);
        aux += F[iz + izmax * iy] * cexp(I * c_2pi_L * q * z) *
               exp(-betabar_0 * z);
      }

      /* Normalize */
      A[iq + iqmax * ip] = aux / (double)(izmax - 1);
    }
  }
}
#else
    ;
#endif

/*
  Evaluates the A_q,sp coefficients of a Frozen Wave (FW) 3D made by SxP FWs,
  restricted to the case where all FWs share the same N, Q and P (i.e., for 1 <=
  p <= P and 1 <= s <= S, N_sp = N, Q_sp = Q, and L_sp = L), resistant or not to
  a possible absorption due to a medium with complex refractive index. See Refs.
  [7] and [8].

  Parameters:
  - S, number S of FWs in the transverse direction x.
  - P, number P of FWs in the transverse direction y.
  - N, parameter N of the FW.
  - L, parameter L of the FW.
  - beta, array of size iqmax = 2 * N + 1 containing the longitudinal wavenumber
  beta_q, where beta_q -> beta[iq], q -> iq - N, and 0<=iq<iqmax.
  - F, array of size ixmax*iymax*izmax containing the morphological function
  F(x,y,z), where F(x,y,z) -> F[iz + izmax*iy + izmax*iymax*ix], x -> ix, y ->
  iy, z -> iz * L / (izmax - 1), 0<=ix<ixmax, 0<=iy<iymax, and 0<=iz<izmax.
  - ixmax, size contribution in ix for the F[iz + izmax*iy + izmax*iymax*ix]
  array, where 0<=ix<ixmax.
  - iymax, size contribution in iy for the F[iz + izmax*iy + izmax*iymax*ix]
  array, where 0<=iy<iymax.
  - izmax, size contribution in iz for the F[iz + izmax*iy + izmax*iymax*ix]
  array, where 0<=iz<izmax.
  - A, array of size ismax*ipmax*iqmax = (S) * (P) * (2 * N + 1) to output the
  A_q,sp coefficients, where A_q,sp -> A[iq + iqmax*ip + iqmax*ipmax*is], s ->
  is + 1, p -> ip + 1, q -> iq - N, 0<=is<ismax, 0<=ip<ipmax, and 0<=iq<iqmax.
  - absorption_resistant, if false, the FW is not absorption-resistant.

  Implementation: The A_q = int_0^L F(x=x_0s, y=y_0p, z) exp(-betabar_0)
  exp(j*2*pi*q*z/L) dz, where x_0s -> ix * (ismax - 1) / (ixmax - 1), y_0p -> iy
  * (ipmax - 1) / (iymax - 1), betabar_0 = imag(beta_0) if absorption_resistant
  = true or betabar_0 = 0 if absorption_resistant = false, is computed by means
  of an approximation of the integral by the trapezoidal method for equally
  spaced z for a given morphologial function F(x,y,z). Notice that for purely
  real refractive indices, the absorption-resistant concerns may be ignored.
*/
FROZEN_WAVES_LIBRARY_API_IMPL_
void fw3d_A_coefficient_restricted(int S, int P, int N, double L,
                                   const double complex *beta, const double *F,
                                   int ixmax, int iymax, int izmax,
                                   double complex *A, bool absorption_resistant)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Constants */
  double betabar_0 = absorption_resistant ? cimag(beta[N]) : 0.0;
  double c_2pi_L = 2.0 * M_PI / L;
  int iqmax = 2 * N + 1;
  int ismax = S;
  int ipmax = P;

  for (int is = 0; is < ismax; is++) {
    for (int ip = 0; ip < ipmax; ip++) {
      for (int iq = 0; iq < iqmax; iq++) {

        double q = (double)(iq - N);

        /* Map is-index and ip-index to ix and iy indices */
        int ix =
            (int)floor((double)is * (double)(ixmax - 1) / (double)(ismax - 1));
        int iy =
            (int)floor((double)ip * (double)(iymax - 1) / (double)(ipmax - 1));

        /* Initial trapezoidal contribution (endpoints in z) */
        double complex aux =
            0.5 *
            (F[(izmax - 1) + izmax * iy + izmax * iymax * ix] *
             exp(-betabar_0 * L) * +F[0 + izmax * iy + izmax * iymax * ix]);

        /* Sum over intermediate z points */
        for (int iz = 1; iz < izmax - 1; iz++) {
          double z = L * (double)iz / (double)(izmax - 1);
          aux += F[iz + izmax * iy + izmax * iymax * ix] *
                 cexp(I * c_2pi_L * q * z) * exp(-betabar_0 * z);
        }

        /* Normalize */
        A[iq + iqmax * ip + iqmax * ipmax * is] = aux / (double)(izmax - 1);
      }
    }
  }
}
#else
    ;
#endif

#endif /* FROZEN_WAVES_LIBRARY_COEFFICIENTS_H */
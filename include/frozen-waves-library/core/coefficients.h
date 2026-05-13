/*
  Frozen Waves Library: A C library providing routines for computing quantities
  related to Frozen Waves

  Repository: <https://github.com/jodesarro/frozen-waves-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Callable functions related to Frozen Wave coefficients.
*/

#ifndef FROZEN_WAVES_LIBRARY_COEFFICIENTS_H
#define FROZEN_WAVES_LIBRARY_COEFFICIENTS_H

#include "../impl/api_impl_.h"
#include <complex.h> /* For 'double complex' type */
#include <math.h>    /* Maths and constants */
#include <omp.h>     /* For parallelization */
#include <stdbool.h> /* For 'bool' type */

#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Evaluates the A_q,ps coefficients of a Frozen Wave (FW) 3D made by P times S
  FWs, restricted to the case where all FWs share the same N, Q and L (i.e.,
  for 1 <= p <= P and 1 <= s <= S, N_ps = N, Q_ps = Q, and L_ps = L),
  resistant or not to a possible absorption due to a medium with complex
  refractive index. See Refs. [7] and [8].

  Parameters:
  - P, number P of FWs in the transverse direction x.
  - S, number S of FWs in the transverse direction y.
  - N, parameter N of the FW.
  - L, parameter L of the FW.
  - beta, array of size iqmax = 2 * N + 1 containing the longitudinal
  wavenumber beta_q, where beta_q -> beta[iq], q -> iq - N, and 0 <= iq < iqmax.
  - F, array of size ixmax * iymax * izmax containing the morphological function
  F(x,y,z), where F(x,y,z) -> F[iz + izmax*iy + izmax*iymax*ix], x -> ix, y
  -> iy, z -> iz * L / (izmax - 1), 0 <= ix < ixmax, 0 <= iy < iymax, and
  0 <= iz < izmax.
  - ixmax, size contribution in ix for the F[iz + izmax*iy + izmax*iymax*ix]
  array, where 0 <= ix < ixmax.
  - iymax, size contribution in iy for the F[iz + izmax*iy + izmax*iymax*ix]
  array, where 0 <= iy < iymax.
  - izmax, size contribution in iz for the F[iz + izmax*iy + izmax*iymax*ix]
  array, where 0 <= iz < izmax.
  - A, array of size ipmax * ismax * iqmax = P * S * (2 * N + 1) to output
  the A_q,ps coefficients, where A_q,ps -> A[iq + iqmax*is +
  iqmax*ismax*ip], p -> ip + 1, s -> is + 1, q -> iq - N, 0 <= ip < ipmax = P,
  0 <= is < ismax = S, and 0 <= iq < iqmax = 2 * N + 1.
  - absorption_resistant, if false, the FW is not absorption-resistant.

  Implementation: The A_q,ps = int_0^L F(x=x_0p, y=y_0s, z) exp(-betabar_0)
  exp(i 2 pi q z/L) dz, where x_0p -> ix * (ipmax - 1) / (ixmax - 1), y_0s
  -> iy * (ismax - 1) / (iymax - 1), betabar_0 = imag(beta_0) if
  absorption_resistant = true or betabar_0 = 0 if absorption_resistant =
  false, is computed by means of an approximation of the integral by the
  trapezoidal method for equally spaced z for a given morphologial function
  F(x,y,z). Notice that for purely real refractive indices, the
  absorption-resistant concerns may be ignored.
*/
void fw3d_A_restricted(int P, int S, int N, double L,
                       const double complex *beta, const double *F, int ixmax,
                       int iymax, int izmax, double complex *A,
                       bool absorption_resistant)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Constants */
  double betabar_0 = absorption_resistant ? cimag(beta[N]) : 0.0;
  double c_2pi_L = 2.0 * M_PI / L;
  int ipmax = P;
  int ismax = S;
  int iqmax = 2 * N + 1;

#pragma omp parallel for collapse(2)

  for (int ip = 0; ip < ipmax; ip++) {

    /* Map ip-index to ix-index */
    int ix = (ipmax == 1) ? 0
                          : (int)floor((double)ip * (double)(ixmax - 1) /
                                       (double)(ipmax - 1));

    for (int is = 0; is < ismax; is++) {

      /* Map is-index to iy-index */
      int iy = (ismax == 1) ? 0
                            : (int)floor((double)is * (double)(iymax - 1) /
                                         (double)(ismax - 1));

      for (int iq = 0; iq < iqmax; iq++) {

        double q = (double)(iq - N);

        /* Initial trapezoidal contribution (endpoints in z) */
        double complex aux =
            0.5 * (F[(izmax - 1) + izmax * iy + izmax * iymax * ix] *
                       exp(-betabar_0 * L) +
                   F[0 + izmax * iy + izmax * iymax * ix]);

        /* Sum over intermediate z points */
        for (int iz = 1; iz < izmax - 1; iz++) {
          double z = L * (double)iz / ((double)(izmax - 1));
          aux += F[iz + izmax * iy + izmax * iymax * ix] * exp(-betabar_0 * z) *
                 cexp(I * c_2pi_L * q * z);
        }

        /* Normalize and store */
        A[iq + iqmax * is + iqmax * ismax * ip] = aux / ((double)(izmax - 1));
      }
    }
  }
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Evaluates the A_q coefficients of a single Frozen Wave (FW),
  resistant or not to a possible absorption due to a medium with complex
  refractive index. See Refs. [1], [2] and [3].

  Parameters:
  - N, parameter N of the FW.
  - L, parameter L of the FW.
  - beta, array of size iqmax = 2 * N + 1 containing the longitudinal
  wavenumber beta_q, where beta_q -> beta[iq], q -> iq - N, and 0 <= iq < iqmax.
  - F, array of size izmax containing the morphological function
  F(z), where F(z) -> F[iz], z -> iz * L / (izmax - 1), and 0 <= iz < izmax.
  - izmax, size contribution in iz for the F[iz] array, where 0 <= iz < izmax.
  - A, array of size iqmax = 2 * N + 1 to output the A_q coefficients, where
  A_q -> A[iq], q -> iq - N, and 0 <= iq < iqmax.
  - absorption_resistant, if false, the FW is not absorption-resistant.

  Implementation: The A_q = int_0^L F(z) exp(-betabar_0)
  exp(i 2 pi q z/L) dz, where betabar_0 = imag(beta_0) if
  absorption_resistant = true or betabar_0 = 0 if absorption_resistant =
  false, is computed through the function fw3d_A_restricted()
  with P = S = ixmax = iymax = 1.
*/
void fw_A(int N, double L, const double complex *beta, const double *F,
          int izmax, double complex *A, bool absorption_resistant)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  fw3d_A_restricted(1, 1, N, L, beta, F, 1, 1, izmax, A, absorption_resistant);
}
#else
    ;
#endif

#endif /* FROZEN_WAVES_LIBRARY_COEFFICIENTS_H */
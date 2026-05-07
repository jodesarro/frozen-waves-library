/*
  Frozen Waves Library: A C library providing routines for computing quantities
  related to Frozen Waves

  Repository: <https://github.com/jodesarro/frozen-waves-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Callable functions related to Frozen Wave wavenumbers.
*/

#ifndef FROZEN_WAVES_LIBRARY_WAVENUMBERS_H
#define FROZEN_WAVES_LIBRARY_WAVENUMBERS_H

#include "../impl/api_impl_.h"
#include <complex.h> /* For 'double complex' type */
#include <stdbool.h> /* For 'bool' type */

#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Fallback for M_SQRT2 if not defined by <math.h> */
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Evaluates the wavenumbers k, beta_q and h_q of a Frozen Wave (FW) using the
  traditional method. See Refs. [1] and [2].

  Parameters:
  - N, parameter N of the FW.
  - Q, parameter Q of the FW.
  - L, parameter L of the FW.
  - k0, wavenumber k_0 with respect to the vacuum (k_0 = omega_0/c, where
  omega_0 is the angular frequency and c the speed of light).
  - nref, refractive index.
  - k, to output the wavenumber k.
  - beta, array of size iqmax = 2 * N + 1 to output the longitudinal wavenumber
  beta_q, where beta_q -> beta[iq], q -> iq - N, and 0 <= iq < iqmax.
  - h, array of size iqmax = 2 * N + 1 to output the transverse wavenumber h_q,
  where h_q -> h[iq], q -> iq - N, and 0 <= iq < iqmax.

  Implementation: It is computed using the dispersion relationship
  k^2=h_q^2+beta_q^2. In the traditional method, it is assumed real axicon
  angles for the Bessel beams, [i.e, 0 <= Re(beta_q) <= Re(k) and 0 <= Re(h_q)
  <= Re(k) must be satisfied], to ensure Re(beta_q) >= 0 and to finally find the
  relations Im(k) / Re(k) = Im(beta_q) / Re(beta_q) = Im(h_q) / Re(h_q) =
  Im(n_ref) / Re(n_ref). Notice that the condition Im(n_ref) << Re(n_ref), or at
  least Im(n_ref) <= 2Re(n_ref)/3pi must be satisfied in order to avoid the
  infinite behavior of the Bessel functions.
*/
void fw_wavenumbers_traditional(int N, double Q, double L, double k0,
                                double complex nref, double complex *k,
                                double complex *beta, double complex *h)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Angular wavenumber */
  *k = nref * k0;

  /* Constants */
  double c_2pi_L = 2.0 * M_PI / L;
  double re_nref = creal(nref);
  double im_nref = cimag(nref);
  double complex k2 = (*k) * (*k);
  int iqmax = 2 * N + 1;

  /* Other wavenumbers */
  for (int iq = 0; iq < iqmax; iq++) {

    double q = (double)(iq - N);

    /* Longitudinal wavenumber */
    double re_b_q = Q + c_2pi_L * q;
    double im_b_q = re_b_q * im_nref / re_nref;
    beta[iq] = re_b_q + I * im_b_q;

    /* Transverse wavenumber */
    double complex h2 = k2 - beta[iq] * beta[iq];
    h[iq] = csqrt(h2);
  }
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Evaluates the wavenumbers k, beta_q and h_q of a Frozen Wave (FW)
  considering the purely real transverse wavenumber h_q method. See Ref.
  [3].

  Parameters:
  - N, parameter N of the FW.
  - Q, parameter Q of the FW.
  - L, parameter L of the FW.
  - k0, wavenumber k_0 with respect to the vacuum (k_0 = omega_0/c, where
  omega_0 is the angular frequency and c the speed of light).
  - nref, refractive index.
  - k, to output the wavenumber k.
  - beta, array of size iqmax = 2 * N + 1 to output the longitudinal
  wavenumber beta_q, where beta_q -> beta[iq], q -> iq - N, and 0 <= iq < iqmax.
  - h, array of size iqmax = 2 * N + 1 to output the transverse wavenumber
  h_q, where h_q -> h[iq], q -> iq - N, and 0 <= iq < iqmax.

  Implementation: It is computed using the dispersion relationship
  k^2=h_q^2+beta_q^2 for purely real h_q. In the purely real h_q method, the
  condition 0 <= Re(beta_q) <= Re(k) must be also satisfied in order to
  ensure Re(beta_q) >= 0 and Im(h_q) = 0.
*/
void fw_wavenumbers_purely_real_h(int N, double Q, double L, double k0,
                                  double complex nref, double complex *k,
                                  double complex *beta, double complex *h)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Angular wavenumber */
  *k = nref * k0;

  /* Constants */
  double re_nref = creal(nref);
  double re_nref2 = re_nref * re_nref;
  double im_nref = cimag(nref);
  double im_nref2 = im_nref * im_nref;
  double c_2pi_L = 2.0 * M_PI / L;
  double k02 = k0 * k0;
  int iqmax = 2 * N + 1;

  /* Other wavenumbers */
  for (int iq = 0; iq < iqmax; iq++) {

    double q = (double)(iq - N);

    /* Longitudinal wavenumber */
    double re_b_q = Q + c_2pi_L * q;
    double im_b_q = k02 * re_nref * im_nref;
    beta[iq] = re_b_q + I * im_b_q;

    /* Transverse wavenumber */
    double h2 = (re_nref2 - im_nref2) * k02 - re_b_q * re_b_q + im_b_q * im_b_q;
    h[iq] = csqrt(h2);
  }
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Evaluates the wavenumbers k, beta_q and h_q of a Frozen Wave (FW) using
  the paraxial approximation for the transverse wavenumber h_q. See Ref.
  [3].

  Parameters:
  - N, parameter N of the FW.
  - Q, parameter Q of the FW.
  - L, parameter L of the FW.
  - k0, wavenumber k_0 with respect to the vacuum (k_0 = omega_0/c, where
  omega_0 is the angular frequency and c the speed of light).
  - nref, refractive index.
  - k, to output the wavenumber k.
  - beta, array of size iqmax = 2 * N + 1 to output the longitudinal
  wavenumber beta_q, where beta_q -> beta[iq], q -> iq - N, and 0 <= iq < iqmax.
  - h, array of size iqmax = 2 * N + 1 to output the transverse wavenumber
  h_q, where h_q -> h[iq], q -> iq - N, and 0 <= iq < iqmax.

  Implementation: It is computed using the dispersion relationship
  h_q=sqrt(2)*k*sqrt(1-beta_q/k) (paraxial approximation). In the paraxial
  h_q method, the condition 0 <= Re(beta_q) <= Re(k) must be also satisfied
  in order to ensure Re(beta_q) >= 0 and Im(h_q) = 0 (notice this also
  implies in purely real h_q).
*/
void fw_wavenumbers_paraxial_h(int N, double Q, double L, double k0,
                               double complex nref, double complex *k,
                               double complex *beta, double complex *h)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Angular wavenumber */
  *k = nref * k0;

  /* Constants */
  double re_k = creal(*k);
  double im_k = cimag(*k);
  double c_2pi_L = 2.0 * M_PI / L;
  int iqmax = 2 * N + 1;

  /* Other wavenumbers */
  for (int iq = 0; iq < iqmax; iq++) {

    double q = (double)(iq - N);

    /* Longitudinal wavenumber */
    double re_b_q = Q + c_2pi_L * q;
    double im_b_q = im_k * (2.0 - re_b_q / re_k);
    beta[iq] = re_b_q + I * im_b_q;

    /* Transverse wavenumber */
    h[iq] = M_SQRT2 * (*k) * csqrt(1.0 - beta[iq] / (*k));
  }
}
#else
    ;
#endif

#endif /* FROZEN_WAVES_LIBRARY_WAVENUMBERS_H */
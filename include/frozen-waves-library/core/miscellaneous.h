/*
  Frozen Waves Library: A C library providing routines for computing quantities
  related to Frozen Waves

  Repository: <https://github.com/jodesarro/frozen-waves-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Callable functions related to Frozen Wave miscellaneous.
*/

#ifndef FROZEN_WAVES_LIBRARY_MISCELLANEOUS_H
#define FROZEN_WAVES_LIBRARY_MISCELLANEOUS_H

#include "../impl/api_impl_.h"
#include <complex.h> /* For 'double complex' type */
#include <float.h>   /* For DBL_EPSILON */
#include <math.h>    /* For NAN and maths */
#include <stdbool.h> /* For 'bool' type */

#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Fallback for M_1_PI if not defined by <math.h> */
#ifndef M_1_PI
#define M_1_PI 0.3183098861837906715377675
#endif

/* First zero of cylindrical Bessel function of the 1st kind and order 0 */
#define C_CYLJ_01 2.4048255576957727686216
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the spot radius of a Bessel beam (BB) of transverse wavenumber h.

  Parameters:
  - h, transverse wavenumber of the BB.
  - asymptotic, if true, uses the asymptotic approximation.

  Implementation: It is computed through c_0 / Re(h), where c_0 = 2.4048... in
  general or c_0 = 3pi/4 in the asymptotic expansion approximation of the Bessel
  functions. Notice that for complex transverse wavenumber h, Im(h) << Re(h), or
  at least Im(h) <= 2Re(h)/3pi must be satisfied to guarantee finite behavior
  for the Bessel functions.
*/
double bb_spot_radius(double complex h, bool asymptotic)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  double re_h = creal(h);
  double c0 = (asymptotic) ? (0.75 * M_PI) : C_CYLJ_01;
  return c0 / re_h;
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the penetration depth of a Bessel beam (BB) for a given longitudinal
  wavenumber beta.

  Parameter:
  - beta, longitudinal wavenumber of the BB.

  Implementation: It is computed using the expression 0.5 / alpha, where alpha =
  -Im(beta).
*/
double bb_penetration_depth(double complex beta)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  double alpha = -cimag(beta);
  return 0.5 / alpha;
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the axicon angle theta of a Bessel beam (BB) for a given wavenumber k
  and transverse wavenumber h.

  Parameters:
  - k, angular wavenumber of the BB.
  - h, transverse wavenumber of the BB.
  - in_degree, returns in radians if false, and in degrees if true.

  Implementation: It is computed using the expression theta = arcsin(h/k).
*/
double complex bb_axicon_angle(double complex k, double complex h,
                               bool in_degree)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  double complex theta = casin(h / k);
  if (in_degree) {
    theta *= (180.0 * M_1_PI); /* 180/pi */
  }
  return theta;
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the aperture radius R for generating a Bessel beam (BB).

  Parameters:
  - beta, longitudinal wavenumber of the BB.
  - h, transverse wavenumber of the BB.
  - L, longitudinal range.

  Implementation: It is computed using R = (h/beta)L, where beta is the
  longitudinal wavenumber, h the transverse wavenumber, and L the longitudinal
  range. It returns NAN if the result is not purely real.
*/
double bb_aperture_radius(double complex beta, double complex h, double L)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  double complex c_h_beta = h / beta;
  if (fabs(cimag(c_h_beta)) < DBL_EPSILON) {
    return L * creal(c_h_beta);
  } else {
    return NAN;
  }
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns R_max, the maximum aperture radius possible for an experimental
  generation of a Bessel Beam (BB) of complex transverse wavenumber h.

  Parameter:
  - h, transverse wavenumber of the BB.

  Implementation: It is computed using R_max = -0.5/Im(h).
*/
double bb_aperture_radius_max(double complex h)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  double im_h = cimag(h);
  if (fabs(im_h) < DBL_EPSILON) {
    return INFINITY;
  } else {
    return -0.5 / im_h;
  }
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns R_min, the minimum aperture radius possible for an experimental
  generation of a Bessel Beam (BB) of complex transverse wavenumber h.

  Parameter:
  - h, transverse wavenumber of the BB.
  - asymptotic, if true, uses the asymptotic approximation.

  Implementation: Returns bb_spot_radius().
*/
double bb_aperture_radius_min(double complex h, bool asymptotic)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  return bb_spot_radius(h, asymptotic);
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the maximum value possible for the integer parameter N of a Frozen
  Wave (FW). See Refs. [1], [2] and [3].

  Parameters:
  - Q, parameter Q of the FW.
  - L, parameter L of the FW.
  - k, wavenumber.

  Implementation: It is computed by means of the condition 0 <= Q +/- 2 pi N/L
  <= Re(k).
*/
int fw_N_max(double Q, double L, double complex k)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  double re_k = creal(k);

  /* Minimum between Q and re_k - Q */
  double min_Q_re_k = (re_k - Q < Q) ? re_k - Q : Q;

  /* Maximum value for N */
  double N_max = L * 0.5 * M_1_PI * min_Q_re_k;
  return (int)floor(N_max);
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the parameter Q of the Frozen Wave (FW) technique for a given spot
  radius Delta_rho for the FW beam. See Refs. [1] and [2].

  Parameters:
  - k, wavenumber.
  - spot_radius, FW spot radius.
  - asymptotic, if true, uses the asymptotic approximation.

  Implementation: It is computed using the expression Q = sqrt(k^2 - c_0^2 /
  Delta_rho^2) of the traditional method for wavenumbers calculation, where c_0
  = 2.4048... in general or c_0 = 3pi/4 for asymptotic expansion approximation
  of the spot radius. Notice that for complex transverse wavenumber h_q, Im(h_q)
  << Re(h_q), or at least Im(h_q) <= 2Re(h_q)/3pi must be satisfied to avoid the
  infinite behavior of the Bessel functions.
*/
double fw_Q_from_spot_radius_traditional(double complex k, double spot_radius,
                                         bool asymptotic)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Constants */
  double c0 = (asymptotic) ? (0.75 * M_PI) : C_CYLJ_01;
  double re_k = creal(k);
  double spot_radius2 = spot_radius * spot_radius;

  /* Return Q from spot radius */
  return sqrt(re_k * re_k - c0 * c0 / spot_radius2);
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the parameter Q of the Frozen Wave (FW) technique for a given spot
  radius Delta_rho for the FW beam, considering purely real transverse
  wavenumbers. See Ref. [3].

  Parameters:
  - k, wavenumber.
  - spot_radius, FW spot radius.
  - asymptotic, if true, uses the asymptotic approximation.

  Implementation: Returns fw_Q_from_spot_radius_traditional().
*/
double fw_Q_from_spot_radius_purely_real_h(double complex k, double spot_radius,
                                           bool asymptotic)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Same case as of the traditional method */
  return fw_Q_from_spot_radius_traditional(k, spot_radius, asymptotic);
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the parameter Q of the Frozen Wave (FW) technique for a given spot
  radius Delta_rho for the FW beam, considering the paraxial approximation. See
  Ref. [3].

  Parameters:
  - k, wavenumber.
  - spot_radius, FW spot radius.
  - asymptotic, if true, uses the asymptotic approximation.

  Implementation: It is computed using the expression Q = Re(k) -
  0.5Re(k)c_0^2/(|k|Delta_rho)^2 of the paraxial h_q method for wavenumbers
  calculation, where c_0 = 2.4048 in general or c_0 = 3pi / 4 for asymptotic
  expansion approximation of the spot radius.
*/
double fw_Q_from_spot_radius_paraxial_h(double complex k, double spot_radius,
                                        bool asymptotic)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  /* Constants */
  double c0 = (asymptotic) ? (0.75 * M_PI) : C_CYLJ_01;
  double re_k = creal(k);
  double spot_radius2 = spot_radius * spot_radius;
  double abs_k2 = cabs(k) * cabs(k);

  /* Return Q from spot radius */
  return re_k * (1.0 - 0.5 * c0 * c0 / (spot_radius2 * abs_k2));
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Returns the Frozen Wave (FW) absorption-resistant condition. See Ref. [4].

  Parameters:
  - N, parameter N of the FW.
  - beta, array of size iqmax = 2 * N + 1, containing the longitudinal
  wavenumber beta_q, where beta_q -> beta[iq] with q -> iq - N, and 0 <= iq <
  iqmax.

  Implementation: It is computed using (Im(beta_N) - Im(beta_-N)) / Im(beta_0),
  where beta_q is the longitudinal wavenumber of the FW. The result must be <<1
  in order to validate the approximation Im(beta_q) ~ Im(beta_0) for all q, and
  so to guarantee the absorption-resistant property of the FW.
*/
double fw_absorption_resistant_condition(int N, const double complex *beta)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  int iqmax = 2 * N + 1;
  return (beta[iqmax - 1] - beta[0]) / beta[N];
}
#else
    ;
#endif

#endif /* FROZEN_WAVES_LIBRARY_MISCELLANEOUS_H */
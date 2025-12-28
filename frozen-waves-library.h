/*
    Frozen Waves Library: A C library providing routines for computing
    quantities related to Frozen Waves

    MIT License Copyright (c) 2025 Jhonas Olivati de Sarro

    https://github.com/jodesarro/frozen-waves-library
*/

/* CHANGELOG
    0.0.0 until Dec 27, 2025 (current version)
        - Inclusion of the first functions: bb_spot_radius,
        bb_penetration_depth, bb_axicon_angle, bb_aperture_radius,
        bb_aperture_radius_max, bb_aperture_radius_min,
        fw_wavenumbers_traditional, fw_wavenumbers_purely_real_h,
        fw_wavenumbers_paraxial_h, fw_N_max, fw_Q_from_spot_radius_traditional,
        fw_Q_from_spot_radius_purely_real_h, fw_Q_from_spot_radius_paraxial_h.
*/

/* SCIENTIFIC NOTES:
    - The harmonic convention adopted is exp(+j omega_0 t) where omega_0 is
    the angular frequency, which means for the refractive index n_ref,
    Re(n_ref) >= 0 and Im(n_ref) = -kappa <= 0, where kappa >= 0, for the
    angular wavenumber k, Re(k) >= 0 and Im(k) <= 0, for the
    transverse wavenumber h, Re(h) >= 0 and Im(h) <= 0, and for the
    longitudinal wavenumber beta, Re(beta) >= 0 and Im(beta) = -alpha <= 0,
    where alpha >= 0.
*/

/* REFERENCES
    [1] M. Zamboni-Rached, "Stationary optical wave fields with arbitrary
    longitudinal shape by superposing equal frequency Bessel beams:
    Frozen Waves," Optics Express, vol. 12, no. 17, pp. 4001--4006,
    Aug. 2004, doi: 10.1364/OPEX.12.004001.

    [2] M. Zamboni-Rached, E. Recami, and H. E. Hern ndez-Figueroa, "Theory of
    'frozen waves': modeling the shape of stationary wave fields," Journal of
    the Optical Society of America A, vol. 22, no. 11, pp. 2465--2475,
    Nov. 2005, doi: 10.1364/JOSAA.22.002465.

    [3] M. Zamboni-Rached and M. Mojahedi, "Shaping finite-energy diffraction-
    and attenuation-resistant beams through Bessel-Gauss beam superposition,"
    Physical Review A, vol. 92, no. 4, p. 043839, Oct. 2015,
    doi: 10.1103/PhysRevA.92.043839.
*/

#ifndef FROZEN_WAVES_LIBRARY_H
#define FROZEN_WAVES_LIBRARY_H

// C++ compilers: Treat the whole code as a C code
#ifdef __cplusplus
extern "C" {
#endif


// -v--------------------------------------------------------------------------
// EXTERNAL LIBRARIES
#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdbool.h>
// EXTERNAL LIBRARIES
// -^--------------------------------------------------------------------------


// -v--------------------------------------------------------------------------
// MISC

// Explicit use of the imaginary unit J instead of I
#ifndef J
#define J I
#endif

// Fallback for M_PI if not defined by <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Fallback for M_1_PI if not defined by <math.h>
#ifndef M_1_PI
#define M_1_PI 0.3183098861837906715377675
#endif

// First zero of cylindrical Bessel function of the 1st kind and order 0
#define C_CYLJ_01 2.4048255576957727686216

// MISC
// -^--------------------------------------------------------------------------


// -v--------------------------------------------------------------------------
// BESSEL BEAMS

/* BB SPOT RADIUS  
    Returns the spot radius of a Bessel beam (BB) of transverse wavenumber
    h. It is computed through c_0 / Re(h), where c_0 = j_01 = 2.4048... in
    general or c_0 = 3pi/4 in the asymptotic expansion approximation of the
    Bessel functions. Notice that for complex transverse wavenumber h,
    Im(h) << Re(h), or at least Im(h) <= 2Re(h)/3pi must be satisfied to
    guarantee finite behavior of the Bessel functions.

    Parameter:
    - h, transverse wavenumber of the BB.
    - asymptotic, if true, uses the asymptotic approximation.
*/
double bb_spot_radius(double complex h, bool asymptotic);

/* BB PENETRATION DEPTH
    Returns the penetration depth of a Bessel beam (BB) for a given
    longitudinal wavenumber beta, using the expression 0.5 / alpha,
    where alpha = -Im(beta).
    
    Parameter:
    - beta, longitudinal wavenumber of the BB.
*/
double bb_penetration_depth(double complex beta);

/* BB AXICON ANGLE
    Returns the axicon angle theta of a Bessel beam (BB) for a given wavenumber
    k and transverse wavenumber h, using the expression theta = asin(h/k).

    Parameters:
    - k, angular wavenumber of the BB;
    - h, transverse wavenumber of the BB;
    - in_degree, returns in radians if false, and in degrees if true.
*/
double complex bb_axicon_angle(double complex k, double complex h,
bool in_degree);

/* BB APERTURE RADIUS
    Returns the aperture radius R for generating a Bessel beam (BB)
    using R = (h/beta)L, where beta is the longitudinal wavenumber,
    h the transverse wavenumber, and L the longitudinal range.

    Parameters:
    - beta, longitudinal wavenumber of the BB;
    - h, transverse wavenumber of the BB;
    - L, longitudinal range.
*/
double bb_aperture_radius(double complex beta, double complex h, double L);

/* BB APERTURE RADIUS MAX
    Returns the maximum aperture radius R_max = -0.5Im(h) possible for an
    experimental generation of a Bessel Beam (BB) of complex transverse
    wavenumber h.

    Parameter:
    - h, transverse wavenumber of the BB;
*/
double bb_aperture_radius_max(double complex h);

/* BB APERTURE RADIUS MIN
    Returns the minimum aperture radius R_min = 3pi/4Re(h) possible
    for an experimental generation of a Bessel Beam (BB).

    Parameter:
    - h, transverse wavenumber of the BB;
*/
double bb_aperture_radius_min(double complex h);

// BESSEL BEAMS
// -^--------------------------------------------------------------------------


// -v--------------------------------------------------------------------------
// FROZEN WAVES

/* FW WAVENUMBERS TRADITIONAL
    Evaluates the wavenumbers k, beta_q and h_q of a Frozen Wave (FW) using the
    traditional method together with the dispersion relationship
    k^2=h_q^2+beta_q^2. See Refs. [1,2]. In such method, it is assumed real
    axicon angles for the Bessel beams (BBs), [i.e, 0 <= Re(beta_q) <= Re(k)
    and 0 <= Re(h_q) <= Re(k) must be satisfied], to ensure Re(beta_q) >= 0 and
    to finally find the relations Im(k) / Re(k) = Im(beta_q) / Re(beta_q)
    = Im(h_q) / Re(h_q) = Im(n_ref) / Re(n_ref). Moreover, the condition
    Im(n_ref) << Re(n_ref), or at least Im(n_ref) <= 2Re(n_ref)/3pi must be
    satisfied in order to avoid the infinite behavior of the Bessel functions.

    Parameters:
    - N, parameter N related to the total of 2N+1 BBs;
    - Q, parameter Q of the FW;
    - L, parameter L of the FW;
    - k0, wavenumber k_0 with respect to the vacuum (k_0 = omega_0/c, where
    omega_0 is the angular frequency and c the speed of light);
    - nref, refractive index;
    - IQMAX, size of the h_q and b_q arrays;
    - *k, wavenumber;
    - b[], array of the longitudinal wavenumber beta_q;
    - h[], array of the transverse wavenumber h_q.
*/
void fw_wavenumbers_traditional(int N, double Q, double L, double k0,
double complex nref, int IQMAX, double complex *k, double complex b[],
double complex h[]);

/* FW WAVENUMBERS PURELY REAL H
    Evaluates the wavenumbers k, beta_q and h_q of a Frozen Wave (FW)
    considering purely real transverse wavenumber h_q together with the
    dispersion relationship k^2=h_q^2+beta_q^2. See Ref. [3]. In such method,
    the condition 0 <= Re(beta_q) <= Re(k) must be also satisfied in order to
    ensure Re(beta_q) >= 0 and Im(h_q) = 0.

    Parameters:
    - N, parameter N related to the total of 2N+1 BBs;
    - Q, parameter Q of the FW;
    - L, parameter L of the FW;
    - k0, wavenumber k_0 with respect to the vacuum (k_0 = omega_0/c, where
    omega_0 is the angular frequency and c the speed of light);
    - nref, refractive index;
    - IQMAX, size of the h_q and b_q arrays;
    - *k, wavenumber;
    - b[], array of the longitudinal wavenumber beta_q;
    - h[], array of the transverse wavenumber h_q.
*/
void fw_wavenumbers_purely_real_h(int N, double Q, double L, double k0,
double complex nref, int IQMAX, double complex *k, double complex b[],
double complex h[]);

/* FW WAVENUMBERS PARAXIAL H
    Evaluates the wavenumbers k, beta_q and h_q of a Frozen Wave (FW) using the
    paraxial approximation for the transverse wavenumber h_q. See Ref. [3]. In
    such method, the condition 0 <= Re(beta_q) <= Re(k) must be also satisfied
    in order to ensure Re(beta_q) >= 0 and Im(h_q) = 0 (notice this also
    implies in purely real h_q).

    Parameters:
    - N, parameter N related to the total of 2N+1 BBs;
    - Q, parameter Q of the FW;
    - L, parameter L of the FW;
    - k0, wavenumber k_0 with respect to the vacuum (k_0 = omega_0/c, where
    omega_0 is the angular frequency and c the speed of light);
    - nref, refractive index;
    - IQMAX, size of the h_q and b_q arrays;
    - *k, wavenumber;
    - b[], array of the longitudinal wavenumber beta_q;
    - h[], array of the transverse wavenumber h_q.
*/
void fw_wavenumbers_paraxial_h(int N, double Q, double L, double k0,
double complex nref, int IQMAX, double complex *k, double complex b[],
double complex h[]);

/* FW N MAX
    Returns the maximum value possible for the parameter N of a Frozen Wave
    (FW) by means of the condition 0 <= Re(beta_q) <= Re(k). See Refs. [1,2,3].

    Parameters:
    - Q, parameter Q of the FW;
    - L, parameter L of the FW;
    - k, wavenumber.
*/
int fw_N_max(double Q, double L, double complex k);

/* FW Q FROM SPOT RADIUS TRADITIONAL
    Returns the parameter Q of the Frozen Wave (FW) technique for a given spot
    radius for the FW beam, using the expression
    Q = sqrt(k^2 - c_0^2 / spot_radius^2) of the traditional method of
    wavenumbers calculation, where c0 = j_01 = 2.4048... in general or
    c_0 = 3pi / 4 for asymptotic expansion approximation of the spot radius.
    See Refs. [1,2]. Notice that for complex transverse wavenumber h_q,
    Im(h_q) << Re(h_q), or at least Im(h_q) <= 2Re(h_q)/3pi must be satisfied
    to avoid the infinite behavior of the Bessel functions.

    Parameters:
    - k, wavenumber;
    - spot_radius, FW spot radius;
    - asymptotic, if true, uses the asymptotic approximation.
*/
double fw_Q_from_spot_radius_traditional(double complex k, double spot_radius,
bool asymptotic);

/* FW Q FROM SPOT RADIUS PURELY REAL H
    Returns the parameter Q of the Frozen Wave (FW) technique for a given spot
    radius for the FW beam, using the expression
    Q = sqrt(k^2 - c_0^2 / spot_radius^2) of the purely real h_q method of
    wavenumbers calculation, where c0 = j_01 = 2.4048... in general or
    c_0 = 3pi / 4 for asymptotic expansion approximation of the spot radius.
    See Ref. [3].

    Parameters:
    - k, wavenumber;
    - spot_radius, FW spot radius;
    - asymptotic, if true, uses the asymptotic approximation.
*/
double fw_Q_from_spot_radius_purely_real_h(double complex k,
double spot_radius, bool asymptotic);

/* FW Q FROM SPOT RADIUS PARAXIAL H
    Returns the parameter Q of the Frozen Wave (FW) technique for a given spot
    radius for the FW beam, using the expression
    Q = Re(k) - 0.5Re(k)c_0^2/(|k|spot_radius)^2 of the paraxial h_q method of
    wavenumbers calculation, where c_0 = j_01 = 2.4048... in general or
    c_0 = 3pi / 4 for asymptotic expansion approximation of the spot radius.
    See Ref. [3].

    Parameters:
    - k, wavenumber;
    - spot_radius, FW spot radius;
    - asymptotic, if true, uses the asymptotic approximation.
*/
double fw_Q_from_spot_radius_paraxial_h(double complex k,
double spot_radius, bool asymptotic);

// FROZEN WAVES
// -^--------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif // C code

#endif // FROZEN_WAVES_LIBRARY_H

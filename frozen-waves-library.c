/*
    Frozen Waves Library: A C library providing routines for computing
    quantities related to Frozen Waves

    MIT License Copyright (c) 2025 Jhonas Olivati de Sarro

    https://github.com/jodesarro/frozen-waves-library
*/

/* CURRENT VERSION
    0.0.0 Dec 27, 2025
*/

/* DOCUMENTATION
    Refer to the header file for CHANGELOG, SCIENTIFIC NOTES, REFERENCES,
    function parameters and descriptions, and more.
*/

#include "frozen-waves-library.h"


// -v--------------------------------------------------------------------------
// BESSEL BEAMS

double bb_spot_radius(double complex h, bool asymptotic) {
    double re_h = creal(h);
    double c0 = (asymptotic) ? (0.75 * M_PI) : C_CYLJ_01;
    return c0 / re_h;
}

double bb_penetration_depth(double complex beta) {
    double alpha = -cimag(beta);
    return 0.5 / alpha;
}

double complex bb_axicon_angle(double complex k, double complex h,
    bool in_degree)
{
    double complex theta = casin(h / k);
    if (in_degree) {
        theta *= (180.0 * M_1_PI); // 180/pi
    }
    return theta;
}

double bb_aperture_radius(double complex beta, double complex h, double L) {
    double complex c0 = h/beta;
    if (fabs(cimag(c0)) < DBL_EPSILON) {
        return L * creal(c0);
    } else {
        return NAN;
    }
}

double bb_aperture_radius_max(double complex h) {
    double im_h = cimag(h);
    return -0.5 / im_h;
}

double bb_aperture_radius_min(double complex h) {
    double re_h = creal(h);
    double c0 = 0.75 * M_PI; // 3pi/4
    return c0 / re_h;
}

// BESSEL BEAMS
// -^--------------------------------------------------------------------------


// -v--------------------------------------------------------------------------
// FROZEN WAVES

void fw_wavenumbers_traditional(int N, double Q, double L, double k0,
    double complex nref, int IQMAX, double complex *k, double complex b[],
    double complex h[])
{
    *k = nref * k0;
    double re_nref = creal(nref);
    double im_nref = cimag(nref);
    for (int iq = 0; iq < IQMAX; iq++) {
        double q = (double)(iq - N);
        double re_b_q = Q + 2.0 * M_PI * q / L;
        double im_b_q = re_b_q * im_nref / re_nref;
        b[iq] = CMPLX(re_b_q, im_b_q);
        double complex h2 = (*k) * (*k) - b[iq] * b[iq];
        h[iq] = csqrt(h2);
    }
}

void fw_wavenumbers_purely_real_h(int N, double Q, double L, double k0,
    double complex nref, int IQMAX, double complex *k, double complex b[],
    double complex h[])
{
    *k = nref * k0;
    double re_nref = creal(nref);
    double im_nref = cimag(nref);
    for (int iq = 0; iq < IQMAX; iq++) {
        double q = (double)(iq - N);
        double re_b_q = Q + 2.0 * M_PI * q / L;
        double im_b_q = k0 * k0 * re_nref * im_nref;
        b[iq] = CMPLX(re_b_q, im_b_q);
        double h2 = (re_nref * re_nref - im_nref * im_nref) * k0 * k0
                      - re_b_q * re_b_q + im_b_q * im_b_q;
        h[iq] = csqrt(h_2);
    }
}

void fw_wavenumbers_paraxial_h(int N, double Q, double L, double k0,
    double complex nref, int IQMAX, double complex *k, double complex b[],
    double complex h[])
{
    *k = nref * k0;
    double re_k = creal(*k);
    double im_k = cimag(*k);
    for (int iq = 0; iq < IQMAX; iq++) {
        double q = (double)(iq - N);
        double re_b_q = Q + 2.0 * M_PI * q / L;
        double im_b_q = im_k * (2.0 - re_b_q / re_k);
        b[iq] = CMPLX(re_b_q, im_b_q);
        h[iq] = M_SQRT2 * (*k) * csqrt(1.0 - b[iq] / (*k));
    }
}

int fw_N_max(double Q, double L, double complex k){
    double re_k = creal(k);
    double c0 = (re_k - Q < Q) ? re_k - Q : Q;
    double N_max = L * 0.5 * M_1_PI * c0;
    return (int)floor(N_max);
}

double fw_Q_from_spot_radius_traditional(double complex k, double spot_radius,
    bool asymptotic)
{
    double c0 = (asymptotic) ? (0.75 * M_PI) : C_CYLJ_01;
    double re_k = creal(k);
    double spot_radius2 = spot_radius * spot_radius;
    double Q = sqrt(re_k * re_k - c0 * c0 / spot_radius2);
    if (fabs(Q) > fabs(re_k)) return NAN; // Must obey Q <= re_k
    return Q;
}

double fw_Q_from_spot_radius_purely_real_h(double complex k,
    double spot_radius, bool asymptotic)
{
    return fw_Q_from_spot_radius_traditional(k, spot_radius, asymptotic);
}

double fw_Q_from_spot_radius_paraxial_h(double complex k, double spot_radius,
     bool asymptotic)
{
    double c0 = (asymptotic) ? (0.75 * M_PI) : C_CYLJ_01;
    double re_k = creal(k);
    double spot_radius2 = spot_radius * spot_radius;
    double abs_k = cabs(k);
    double Q = re_k * ( 1.0 - 0.5 * c0 * c0 / (spot_radius2 * abs_k * abs_k) );
    if (fabs(Q) > fabs(re_k)) return NAN; // Must obey Q <= re_k
    return Q;
}

// FROZEN WAVES
// -^--------------------------------------------------------------------------

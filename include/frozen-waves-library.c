/*
    Frozen Waves Library: A C library providing routines for computing
    quantities related to Frozen Waves

    MIT License Copyright (c) 2025 Jhonas Olivati de Sarro

    <https://github.com/jodesarro/frozen-waves-library>
*/

/* CURRENT VERSION
    0.0.1 Dec 28, 2025
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
bool in_degree) {
    double complex theta = casin(h / k);
    if (in_degree) {
        theta *= (180.0 * M_1_PI); // 180/pi
    }
    return theta;
}

double bb_aperture_radius(double complex beta, double complex h, double L) {
    double complex c_h_beta = h/beta;
    if (fabs(cimag(c_h_beta)) < DBL_EPSILON) {
        return L * creal(c_h_beta);
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
double complex nref, double complex *k, double complex beta[],
double complex h[], int IQMAX) {

    // Angular wavenumber
    *k = nref * k0;

    // Constants
    double c_2pi_L = 2.0 * M_PI / L;
    double re_nref = creal(nref);
    double im_nref = cimag(nref);
    double complex k2 = (*k) * (*k)
    
    // Other wavenumbers
    for (int iq = 0; iq < IQMAX; iq++) {
        
        double q = (double)(iq - N);
        
        // Longitudinal wavenumber
        double re_b_q = Q + c_2pi_L * q;
        double im_b_q = re_b_q * im_nref / re_nref;        
        beta[iq] = CMPLX(re_b_q, im_b_q);
        
        // Transverse wavenumber
        double complex h2 = k2 - beta[iq] * beta[iq];
        h[iq] = csqrt(h2);
    }
}

void fw_wavenumbers_purely_real_h(int N, double Q, double L, double k0,
double complex nref, double complex *k, double complex beta[],
double complex h[], int IQMAX) {
    
    // Angular wavenumber
    *k = nref * k0;
    
    // Constants
    double re_nref = creal(nref);
    double re_nref2 = re_nref * re_nref;
    double im_nref = cimag(nref);
    double im_nref2 = im_nref * im_nref;
    double c_2pi_L = 2.0 * M_PI / L;
    double k02 = k0 * k0;
    
    // Other wavenumbers
    for (int iq = 0; iq < IQMAX; iq++) {
        
        double q = (double)(iq - N);

        // Longitudinal wavenumber
        double re_b_q = Q + c_2pi_L * q;
        double im_b_q = k02 * re_nref * im_nref;
        beta[iq] = CMPLX(re_b_q, im_b_q);

        // Transverse wavenumber
        double h2 = (re_nref2 - im_nref2) * k02
                    - re_b_q * re_b_q + im_b_q * im_b_q;
        h[iq] = csqrt(h2);
    }
}

void fw_wavenumbers_paraxial_h(int N, double Q, double L, double k0,
double complex nref, double complex *k, double complex beta[],
double complex h[], int IQMAX) {

    // Angular wavenumber
    *k = nref * k0;

    // Constants
    double re_k = creal(*k);
    double im_k = cimag(*k);
    double c_2pi_L = 2.0 * M_PI / L;

    // Other wavenumbers
    for (int iq = 0; iq < IQMAX; iq++) {

        double q = (double)(iq - N);

        // Longitudinal wavenumber
        double re_b_q = Q + c_2pi_L * q;
        double im_b_q = im_k * (2.0 - re_b_q / re_k);
        beta[iq] = CMPLX(re_b_q, im_b_q);
        
        // Transverse wavenumber
        h[iq] = M_SQRT2 * (*k) * csqrt(1.0 - beta[iq] / (*k));
    }
}

int fw_N_max(double Q, double L, double complex k) {
    
    double re_k = creal(k);

    // Minimum between Q and re_k - Q
    double min_Q_re_k = (re_k - Q < Q) ? re_k - Q : Q;

    // Maximum value for N
    double N_max = L * 0.5 * M_1_PI * min_Q_re_k;
    return (int)floor(N_max);
}

double fw_Q_from_spot_radius_traditional(double complex k, double spot_radius,
bool asymptotic) {
    
    // Constants
    double c0 = (asymptotic) ? (0.75 * M_PI) : C_CYLJ_01;
    double re_k = creal(k);
    double spot_radius2 = spot_radius * spot_radius;

    // Return Q from spot radius
    return sqrt(re_k * re_k - c0 * c0 / spot_radius2);
}

double fw_Q_from_spot_radius_purely_real_h(double complex k,
double spot_radius, bool asymptotic) {
    // Same case as of traditional method
    return fw_Q_from_spot_radius_traditional(k, spot_radius, asymptotic);
}

double fw_Q_from_spot_radius_paraxial_h(double complex k, double spot_radius,
bool asymptotic) {

    // Constants
    double c0 = (asymptotic) ? (0.75 * M_PI) : C_CYLJ_01;
    double re_k = creal(k);
    double spot_radius2 = spot_radius * spot_radius;
    double abs_k2 = cabs(k) * cabs(k);
        
    // Return Q from spot radius
    return re_k * ( 1.0 - 0.5 * c0 * c0 / (spot_radius2 * abs_k2) );
}

double fw_absorption_resistant_condition(int N, const double complex beta[],
int IQMAX) {
    return (beta[IQMAX-1] - beta[0]) / beta[N];
}

void fw_A_coefficient(int N, int L, const double complex beta[], int IQMAX,
const double F[], int IZMAX, double complex A[], bool absorption_resistant) {
    
    // Constants
    double betabar_0 = absorption_resistant ? cimag(beta[N]) : 0.0;
    double c_2pi_L = 2.0 * M_PI / L;

    for (int iq = 0; iq < IQMAX; iq++) {
        
        double q = (double)(iq - N);

        // Initial trapezoidal contribution (endpoints)
        double complex aux = 0.5 * (F[IZMAX - 1] * exp(-betabar_0 * L) + F[0]);

        // Sum over intermediate points
        for (int iz = 1; iz < IZMAX - 1; iz++) {
            double z = L * (double)iz / (double)(IZMAX - 1);
            aux += F[iz] * cexp(J * c_2pi_L * q * z )
                   * exp(-betabar_0 * z);
        }

        // Normalize
        A[iq] = aux / (double)(IZMAX - 1);
    }
}

void fw2d_A_coefficient_restricted(int N, int L, const double complex beta[],
int IQMAX, const double F[], int IYMAX, int IZMAX, double complex A[],
int IPMAX, bool absorption_resistant) {

    // Constants
    double betabar_0 = absorption_resistant ? cimag(beta[N]) : 0.0;
    double c_2pi_L = 2.0 * M_PI / L;

    for (int ip = 0; ip < IPMAX; ip++) {
        
        for (int iq = 0; iq < IQMAX; iq++) {
            
            double q = (double)(iq - N);

            // Map ip-index to iy-index (linear scaling)
            int iy = (int)floor((double)ip * (double)(IYMAX - 1)
                     / (double)(IPMAX - 1));

            // Initial trapezoidal contribution (endpoints)
            double complex aux = 0.5 * (F[iy + IYMAX * (IZMAX - 1)]
                                 * exp(-betabar_0 * L) + F[iy + 0]);

            // Sum over intermediate z points
            for (int iz = 1; iz < IZMAX - 1; iz++) {
                double z = L * (double)iz / (double)(IZMAX - 1);
                aux += F[iy + IYMAX * iz] * cexp(J * c_2pi_L * q * z )
                       * exp(-betabar_0 * z);
            }

            // Normalize
            A[iq + IQMAX * ip] = aux / (double)(IZMAX - 1);
        }
    }
}

void fw3d_A_coefficient_restricted(int N, int L, const double complex beta[],
int IQMAX, const double F[], int IXMAX, int IYMAX, int IZMAX,
double complex A[], int ISMAX, int IPMAX, bool absorption_resistant) {
    
    // Constants
    double betabar_0 = absorption_resistant ? cimag(beta[N]) : 0.0;
    double c_2pi_L = 2.0 * M_PI / L;

    for (int is = 0; is < ISMAX; is++) {
        for (int ip = 0; ip < IPMAX; ip++) {
            for (int iq = 0; iq < IQMAX; iq++) {
                
                double q = (double)(iq - N);

                // Map is-index and ip-index to ix and iy indices
                int ix = (int)floor((double)is * (double)(IXMAX - 1)
                         / (double)(ISMAX - 1));
                int iy = (int)floor((double)ip * (double)(IYMAX - 1)
                         / (double)(IPMAX - 1));

                // Initial trapezoidal contribution (endpoints in z)
                double complex aux = 0.5
                           * ( F[ix + IXMAX * iy + IXMAX * IYMAX * (IZMAX - 1)]
                           * exp(-betabar_0 * L) + F[ix + IXMAX * iy + 0]);

                // Sum over intermediate z points
                for (int iz = 1; iz < IZMAX - 1; iz++) {
                    double z = L * (double)iz / (double)(IZMAX - 1);
                    aux += F[ix + IXMAX * iy + IXMAX * IYMAX * iz]
                           * cexp(J * c_2pi_L * q * z)
                           * exp(-betabar_0 * z);
                }

                // Normalize
                A[iq + IQMAX * is + IQMAX * ISMAX * ip]
                                                   = aux / (double)(IZMAX - 1);
            }
        }
    }
}

// FROZEN WAVES
// -^--------------------------------------------------------------------------

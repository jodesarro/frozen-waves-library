/* 
    Frozen Waves Library: A C library providing routines for computing
    quantities related to Frozen Waves

    File: include/frozen-waves-library/impl/cplx_c_cpp_impl_.h
    Version: include/frozen-waves-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 and C++98 guards, macros and typedefs
    License: include/frozen-waves-library/license.txt

    Description:
        Defines macros and typedefs for ensuring C++ compatibility
*/

#ifndef FROZEN_WAVES_LIBRARY_CPLX_C_CPP_IMPL_H
#define FROZEN_WAVES_LIBRARY_CPLX_C_CPP_IMPL_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */

#include <complex> /* For complex numbers */
#define I_IMPL_ std::complex<double>(0.0, 1.0)
#define creal(z) std::real(z)
#define cimag(z) std::imag(z)
#define cabs(z) std::abs(z)
#define cexp(z) std::exp(z)
#define sin(z) std::sin(z)
#define cos(z) std::cos(z)
#define casin(z) std::asin(z)
typedef std::complex<double> tpdfcplx_impl_;

#else

/* C99 */

#include <complex.h> /* For complex numbers */
#define I_IMPL_ I
typedef double complex tpdfcplx_impl_;

#endif /* __cplusplus */

#endif /* FROZEN_WAVES_LIBRARY_CPLX_C_CPP_IMPL_H */
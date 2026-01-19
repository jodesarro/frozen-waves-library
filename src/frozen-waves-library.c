/* 
    Frozen Waves Library: A C library providing routines for computing
    quantities related to Frozen Waves

    File: src/frozen-waves-library.c
    Version: include/frozen-waves-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with C++ guards
    References: include/frozen-waves-library/references.txt
    License: include/frozen-waves-library/license.txt

    Description:
        Wrapper for compiling the include/frozen-waves-library.h
*/

/* Overwrite 'static inline' */
/* Overwrite 'static inline' */
#if defined(_WIN32) || defined(_WIN64)
    #ifdef __cplusplus
        #define FROZEN_WAVES_LIBRARY_STATIC_INLINE_IMPL_ __declspec(dllexport) extern "C"
    #else
        #define FROZEN_WAVES_LIBRARY_STATIC_INLINE_IMPL_ __declspec(dllexport)
    #endif
#else
    #ifdef __cplusplus
        #define FROZEN_WAVES_LIBRARY_STATIC_INLINE_IMPL_ extern "C"
    #else
        #define FROZEN_WAVES_LIBRARY_STATIC_INLINE_IMPL_
    #endif
#endif

#include "../include/frozen-waves-library.h"
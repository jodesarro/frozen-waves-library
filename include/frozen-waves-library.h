/* 
    Frozen Waves Library: A C library providing routines for computing
    quantities related to Frozen Waves

    File: include/frozen-waves-library.h
    Version: include/frozen-waves-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/frozen-waves-library/references.txt
    License: include/frozen-waves-library/license.txt

    Description:
        Include headers of all core functions of the
        include/frozen-waves-library/core folder.

    Scientific notes:
        The harmonic convention adopted is exp(+j omega_0 t) where omega_0 is
        the angular frequency, which means for the refractive index n_ref,
        Re(n_ref) >= 0 and Im(n_ref) = -kappa <= 0, where kappa >= 0, for the
        angular wavenumber k, Re(k) >= 0 and Im(k) <= 0, for the
        transverse wavenumber h, Re(h) >= 0 and Im(h) <= 0, and for the
        longitudinal wavenumber beta, Re(beta) >= 0 and Im(beta) = -alpha <= 0
        where alpha >= 0.
*/

#ifndef FROZEN_WAVES_LIBRARY_H
#define FROZEN_WAVES_LIBRARY_H

#include "frozen-waves-library/version.h"
#include "frozen-waves-library/core/coefficients.h"
#include "frozen-waves-library/core/fields.h"
#include "frozen-waves-library/core/miscellaneous.h"
#include "frozen-waves-library/core/wavenumbers.h"

#endif /* FROZEN_WAVES_LIBRARY_H */
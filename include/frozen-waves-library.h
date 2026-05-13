/*
  Frozen Waves Library: A C library providing routines for computing quantities
  related to Frozen Waves

  Repository: <https://github.com/jodesarro/frozen-waves-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99
  Last updated: 2026-05-13

  Description: Include headers of all core functions of the
  include/frozen-waves-library/core folder.

  Scientific notes: The harmonic convention adopted is exp(+i omega_0 t) where
  omega_0 is the angular frequency, which means for the refractive index n_ref,
  Re(n_ref) >= 0 and Im(n_ref) = -kappa <= 0, where kappa >= 0, for the angular
  wavenumber k, Re(k) >= 0 and Im(k) <= 0, for the transverse wavenumber h,
  Re(h) >= 0 and Im(h) <= 0, and for thelongitudinal wavenumber beta, Re(beta)
  >= 0 and Im(beta) = -alpha <= 0 where alpha >= 0.
*/

#ifndef FROZEN_WAVES_LIBRARY_H
#define FROZEN_WAVES_LIBRARY_H

#include "frozen-waves-library/core/coefficients.h"
#include "frozen-waves-library/core/miscellaneous.h"
#include "frozen-waves-library/core/psi.h"
#include "frozen-waves-library/core/wavenumbers.h"

#endif /* FROZEN_WAVES_LIBRARY_H */
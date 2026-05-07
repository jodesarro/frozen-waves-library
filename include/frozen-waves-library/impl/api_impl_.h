/*
  Frozen Waves Library: A C library providing routines for computing quantities
  related to Frozen Waves

  Repository: <https://github.com/jodesarro/frozen-waves-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Define the API with macros for compilation, and for
  header-only or compiled library usage.
*/

#ifndef FROZEN_WAVES_LIBRARY_API_IMPL_H
#define FROZEN_WAVES_LIBRARY_API_IMPL_H

#ifdef FROZEN_WAVES_LIBRARY_EXPORTS_IMPL_
#if defined(_WIN32) || defined(_WIN64)
#define FROZEN_WAVES_LIBRARY_VISIBILITY_IMPL_ __declspec(dllexport)
#else
#define FROZEN_WAVES_LIBRARY_VISIBILITY_IMPL_                                  \
  __attribute__((visibility("default")))
#endif
#elif defined(FROZEN_WAVES_LIBRARY_IMPORTS)
#if defined(_WIN32) || defined(_WIN64)
#define FROZEN_WAVES_LIBRARY_VISIBILITY_IMPL_ __declspec(dllimport)
#else
#define FROZEN_WAVES_LIBRARY_VISIBILITY_IMPL_
#endif
#else
#define FROZEN_WAVES_LIBRARY_VISIBILITY_IMPL_ static inline
#endif

#define FROZEN_WAVES_LIBRARY_API_IMPL_ FROZEN_WAVES_LIBRARY_VISIBILITY_IMPL_

#endif /* FROZEN_WAVES_LIBRARY_API_IMPL_H */
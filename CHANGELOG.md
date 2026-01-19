# CHANGELOG

Refer to the [References](README.md#references) section of the
[README.md](README.md) file for a complete list of references.

## 0.1.1 Jan 19, 2026
- Library reorganized using specific folders: src, include, core, extern,
and impl folders.
- The new file src/frozen-waves-library.c is for compiling purposes.
- Added files src/frozen-waves-library-declarations.c and
src/frozen-waves-library-declarations.cpp, containing a list of all the
functions declarations.
- The main header include/frozen-waves-library.h now is responsible for
including the headers of the core folder.
- Creation of CHANGELOG.md and version.h for versioning control.
- README.md file created.

## 0.1.0 Dec 28, 2025
- Improvement on the comments in the all files.
- Size of arrays placed after the arrays in the parameters of all
functions for standardization purposes.
- `b[]` arrays changed to `beta[]` arrays for standardization purposes.
- Inclusion of the functions: `fw_absorption_resistant_condition`,
`fw_A_coefficient`, `fw2d_A_coefficient`, `fw3d_A_coefficient_restricted`.

## 0.0.0 until Dec 27, 2025
- Inclusion of the first functions: `bb_spot_radius`,
`bb_penetration_depth`, `bb_axicon_angle`, `bb_aperture_radius`,
`bb_aperture_radius_max`, `bb_aperture_radius_min`,
`fw_wavenumbers_traditional`, `fw_wavenumbers_purely_real_h`,
`fw_wavenumbers_paraxial_h`, `fw_N_max`, `fw_Q_from_spot_radius_traditional`,
`fw_Q_from_spot_radius_purely_real_h`, `fw_Q_from_spot_radius_paraxial_h`.
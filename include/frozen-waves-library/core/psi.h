/*
  Frozen Waves Library: A C library providing routines for computing quantities
  related to Frozen Waves

  Repository: <https://github.com/jodesarro/frozen-waves-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Callable functions related to Frozen Wave scalar fields, i.e.,
  to the wave function psi.
*/

#ifndef FROZEN_WAVES_LIBRARY_PSI_H
#define FROZEN_WAVES_LIBRARY_PSI_H

#include "../../bessel-library.h" /* A library to evaluate Bessel functions */
#include "../../console-progress-bar.h" /* Routines to display progress bars */
#include "../impl/api_impl_.h"
#include <complex.h> /* For 'double complex' type */
#include <math.h>    /* Maths and constants */
#include <omp.h>     /* For parallelization */

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Evaluates the scalar wave function psi(x,y,z) of a Frozen Wave (FW) 3D made by
  P times S FWs, restricted to the case where all FWs share the same N, Q and L
  (i.e., for 1 <= p <= P and 1 <= s <= S, N_ps = N, Q_ps = Q, and L_ps = L).
  See Refs. [5], [6] and [7].

  Parameters:
  - P, number S of FWs in the transverse direction x.
  - S, number P of FWs in the transverse direction y.
  - N, parameter N of the FW.
  - beta, array of size iqmax = 2 * N + 1 containing the longitudinal
  wavenumber beta_q, where beta_q -> beta[iq], q -> iq - N, and 0 <= iq < iqmax.
  - h, array of size iqmax = 2 * N + 1 containing the transverse
  wavenumber h_q, where h_q -> h[iq], q -> iq - N, and 0 <= iq < iqmax.
  - A, array of size ipmax * ismax * iqmax = P * S * (2 * N + 1) containing the
  A_q,ps coefficients, where A_q,ps -> A[iq + iqmax*is + iqmax*ismax*ip], p ->
  ip + 1, s -> is + 1, q -> iq - N, 0 <= ip < ipmax = P, 0 <= is < ismax = S,
  and 0 <= iq < iqmax = 2 * N + 1.
  - x0, array of size ipmax = P containing the x-coordinates of the origins of
  the FWs, where x_0,p -> x0[ip], p -> ip + 1, 0 <= ip < ipmax = P.
  - y0, array of size ismax = S containing the y-coordinates of the origins of
  the FWs, where y_0,s -> y0[is], s -> is + 1, 0 <= is < ismax = S.
  - z0, array of size ipmax * ismax = P * S containing the z-coordinates of the
  origins of the FWs, where z_0,ps -> z0[is + ismax*ip], p -> ip + 1, s -> is +
  1, 0 <= ip < ipmax = P, 0 <= is < ismax = S.
  - xmin, minimum value for the x coordinate, where x >= xmin.
  - xmax, maximum value for the x coordinate, where x <= xmax.
  - xpoints, number of points in the x direction, where 0 <= ix < xpoints.
  - ymin, minimum value for the y coordinate, where y >= ymin.
  - ymax, maximum value for the y coordinate, where y <= ymax.
  - ypoints, number of points in the y direction, where 0 <= iy < ypoints.
  - zmin, minimum value for the z coordinate, where z >= zmin.
  - zmax, maximum value for the z coordinate, where z <= zmax.
  - zpoints, number of points in the z direction, where 0 <= iz < zpoints.
  - psi, array of size xpoints * ypoints * zpoints to output the scalar wave
  function psi(x,y,z), where psi(x,y,z) -> psi[iz + zpoints*iy +
  zpoints*ypoints*ix], x -> xmin + ix * (xmax - xmin) / (xpoints - 1), y -> ymin
  + iy * (ymax - ymin) / (ypoints - 1), z -> zmin + iz * (zmax - zmin) /
  (zpoints - 1).

  Implementation: The scalar wave function is computed through the formula
  psi(x,y,z) = sum_{p=1}^P sum_{s=1}^S sum_{q=-N}^N A_q,ps cyl_j(0, h_q
  sqrt{(x-x_0,p)^2 + (y-y_0,s)^2}) exp(-i beta_q (z - z_0,ps)).

  Details: The computation is parallelized with OpenMP. If the macro
  FROZEN_WAVES_LIBRARY_PROGRESS_BAR is defined, a progress bar is displayed on
  the console.
*/
void fw3d_psi_restricted(int P, int S, int N, const double complex *beta,
                         const double complex *h, const double complex *A,
                         const double *x0, const double *y0, const double *z0,
                         double xmin, double xmax, int xpoints, double ymin,
                         double ymax, int ypoints, double zmin, double zmax,
                         int zpoints, double complex *psi)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  int ipmax = P;
  int ismax = S;
  int iqmax = N * 2 + 1;

  double dx, dy, dz;
  dx = (xpoints == 1) ? (0.) : (xmax - xmin) / (double)(xpoints - 1);
  dy = (ypoints == 1) ? (0.) : (ymax - ymin) / (double)(ypoints - 1);
  dz = (zpoints == 1) ? (0.) : (zmax - zmin) / (double)(zpoints - 1);

#ifdef FROZEN_WAVES_LIBRARY_PROGRESS_BAR
  long int completed = 0;
  long int total = (long int)xpoints * (long int)ypoints * (long int)zpoints;
  print_progress_bar_empty();
#pragma omp parallel for collapse(3) shared(completed)
#else
#pragma omp parallel for collapse(3)
#endif

  for (int ix = 0; ix < xpoints; ix++) {
    double x = xmin + ((double)ix) * dx;

    for (int iy = 0; iy < ypoints; iy++) {
      double y = ymin + ((double)iy) * dy;

      for (int iz = 0; iz < zpoints; iz++) {
        double z = zmin + ((double)iz) * dz;
        double complex aux = 0.;

        for (int ip = 0; ip < ipmax; ip++) {
          double square_of_diff_x_x0 = (x - x0[ip]) * (x - x0[ip]);

          for (int is = 0; is < ismax; is++) {
            double square_of_diff_y_y0 = (y - y0[is]) * (y - y0[is]);
            double varrho = sqrt(square_of_diff_x_x0 + square_of_diff_y_y0);

            for (int iq = 0; iq < iqmax; iq++) {
              aux += A[iq + iqmax * is + iqmax * ismax * ip] *
                     cyl_j(0., h[iq] * varrho) *
                     cexp(-I * beta[iq] * (z - z0[is + ismax * ip]));
            }
          }
        }
        psi[iz + zpoints * iy + zpoints * ypoints * ix] = aux;
#ifdef FROZEN_WAVES_LIBRARY_PROGRESS_BAR
#pragma omp atomic
        completed++;
#pragma omp critical
        print_progress_bar_every_percent(completed, total, 10);
#endif
      }
    }
  }
#ifdef FROZEN_WAVES_LIBRARY_PROGRESS_BAR
  print_progress_bar_full();
#endif
}
#else
    ;
#endif

FROZEN_WAVES_LIBRARY_API_IMPL_
/*
  Evaluates the scalar wave function psi(x,y,z) of a single Frozen Wave (FW).
  See Refs. [1], [2], [3] and [4].

  Parameters:
  - N, parameter N of the FW.
  - beta, array of size iqmax = 2 * N + 1 containing the longitudinal
  wavenumber beta_q, where beta_q -> beta[iq], q -> iq - N, and 0 <= iq < iqmax.
  - h, array of size iqmax = 2 * N + 1 containing the transverse
  wavenumber h_q, where h_q -> h[iq], q -> iq - N, and 0 <= iq < iqmax.
  - A, array of size iqmax = 2 * N + 1 containing the
  A_q coefficients, where A_q -> A[iq], q -> iq - N, and 0 <= iq < iqmax.
  - xmin, minimum value for the x coordinate, where x >= xmin.
  - xmax, maximum value for the x coordinate, where x <= xmax.
  - xpoints, number of points in the x direction, where 0 <= ix < xpoints.
  - ymin, minimum value for the y coordinate, where y >= ymin.
  - ymax, maximum value for the y coordinate, where y <= ymax.
  - ypoints, number of points in the y direction, where 0 <= iy < ypoints.
  - zmin, minimum value for the z coordinate, where z >= zmin.
  - zmax, maximum value for the z coordinate, where z <= zmax.
  - zpoints, number of points in the z direction, where 0 <= iz < zpoints.
  - psi, array of size xpoints * ypoints * zpoints to output the scalar wave
  function psi(x,y,z), where psi(x,y,z) -> psi[iz + zpoints*iy +
  zpoints*ypoints*ix], x -> xmin + ix * (xmax - xmin) / (xpoints - 1), y -> ymin
  + iy * (ymax - ymin) / (ypoints - 1), z -> zmin + iz * (zmax - zmin) /
  (zpoints - 1).

  Implementation: The psi(x,y,z) = sum_{q=-N}^N A_q cyl_j(0, h_q sqrt{x^2 +
  y^2}) exp(-i beta_q z) is computed through the function fw3d_psi_restricted()
  with S = P = 1 and x0[0] = y0[0] = z0[0] = 0.

  Details: The computation is parallelized with OpenMP. If the macro
  FROZEN_WAVES_LIBRARY_PROGRESS_BAR is defined, a progress bar is displayed on
  the console.
*/
void fw_psi(int N, const double complex *beta, const double complex *h,
            const double complex *A, double xmin, double xmax, int xpoints,
            double ymin, double ymax, int ypoints, double zmin, double zmax,
            int zpoints, double complex *psi)
#ifndef FROZEN_WAVES_LIBRARY_IMPORTS
{
  fw3d_psi_restricted(1, 1, N, beta, h, A, (double[]){0.}, (double[]){0.},
                      (double[]){0.}, xmin, xmax, xpoints, ymin, ymax, ypoints,
                      zmin, zmax, zpoints, psi);
}
#else
    ;
#endif

#endif /* FROZEN_WAVES_LIBRARY_PSI_H */
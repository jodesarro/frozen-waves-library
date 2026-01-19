# Frozen Waves Library: A C library providing routines for computing quantities related to Frozen Waves

This project provides a C library
([compatible with C++](#compatibility-with-c)) designed for the
computation of several quantities related to Frozen Waves (FWs).
If you are unfamiliar with FWs, see the works listed in the
[References](#references) for an overview.

## Available features

### Miscellaneous

<details>
  <summary>
    <code><b>bb_spot_radius(h, asymptotic)</b></code>
  </summary>

  - **Description:** Returns the spot radius of a Bessel beam (BB) of
  transverse wavenumber $h$.
  - **Parameters:**
    - `h`, transverse wavenumber of the BB.
    - `asymptotic`, if true, uses the asymptotic approximation.
  - **Implementation:** It is computed through $c_0/\mathrm{Re}(h)$, where
    $c_0 = 2.4048$ in general or $c_0 = 3\pi/4$ in the asymptotic
    expansion approximation of the Bessel functions.
    Notice that for complex transverse wavenumber $h$,
    $\mathrm{Im}(h) \ll \mathrm{Re}(h)$,
    or at least $\mathrm{Im}(h) \leq 2\mathrm{Re}(h)/3\pi$ must be satisfied
    to guarantee finite behavior for the Bessel functions.
</details>

<details>
  <summary>
    <code><b>bb_penetration_depth(beta)</b></code>
  </summary>

  - **Description:** Returns the penetration depth of a Bessel beam (BB) for a
  given longitudinal wavenumber $\beta$.
  - **Parameter:**
    - `beta`, longitudinal wavenumber of the BB.
  - **Implementation:** It is computed using the expression $0.5 / \alpha$,
  where $\alpha = -\mathrm{Im}(\beta)$.
</details>

<details>
  <summary>
    <code><b>bb_axicon_angle(k, h, in_degree)</b></code>
  </summary>

  - **Description:** Returns the axicon angle $\theta$ of a Bessel beam (BB)
  for a given wavenumber $k$ and transverse wavenumber $h$.
  - **Parameters:**
    - `k`, angular wavenumber of the BB.
    - `h`, transverse wavenumber of the BB.
    - `in_degree`, returns in radians if false, and in degrees if true.
  - **Implementation:** It is computed using the expression
  $\theta = \arcsin(h/k)$.
</details>

<details>
  <summary>
    <code><b>bb_aperture_radius(beta, h, L)</b></code>
  </summary>

  - **Description:** Returns the aperture radius $R$ for generating a Bessel
  beam (BB).
  - **Parameters:**
    - `beta`, longitudinal wavenumber of the BB.
    - `h`, transverse wavenumber of the BB.
    - `L`, longitudinal range.
  - **Implementation:** It is computed using $R = (h/\beta)L$, where $\beta$
  is the longitudinal wavenumber, $h$ the transverse wavenumber, and $L$ the
  longitudinal range.
</details>

<details>
  <summary>
    <code><b>bb_aperture_radius_max(h)</b></code>
  </summary>

  - **Description:** Returns $R_{max}$, the maximum aperture radius possible for
  an experimental generation of a Bessel Beam (BB) of complex transverse
  wavenumber $h$.
  - **Parameter:**
    - `h`, transverse wavenumber of the BB.
  - **Implementation:** It is computed using $R_{max} = -0.5/\mathrm{Im}(h)$.
</details>

<details>
  <summary>
    <code><b>bb_aperture_radius_min(h)</b></code>
  </summary>

  - **Description:** Returns $R_{min}$, the minimum aperture radius possible
  for an experimental generation of a Bessel Beam (BB) of complex transverse
  wavenumber $h$.
  - **Parameter:**
    - `h`, transverse wavenumber of the BB.
  - **Implementation:** It is computed using $R_{min} = 3\pi/4\mathrm{Re}(h)$.
</details>

<details>
  <summary>
    <code><b>fw_N_max(Q, L, k)</b></code>
  </summary>

  - **Description:** Returns the maximum value possible for the integer
  parameter $N$ of a Frozen Wave (FW). See Refs. [[1]](#ref1), [[2]](#ref2)
  and [[3]](#ref3).
  - **Parameters:**
    - `Q`, parameter $Q$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `k`, wavenumber $k$.
  - **Implementation:** It is computed by means of the condition
  $0 \leq \mathrm{Re}(\beta_q) \leq \mathrm{Re}(k)$, where $\beta_q$ is the
  longitudinal wavenumber.
</details>

<details>
  <summary>
    <code><b>fw_Q_from_spot_radius_traditional(k, spot_radius, asymptotic)</b></code>
  </summary>

  - **Description:** Returns the parameter $Q$ of the Frozen Wave (FW)
  technique for a given spot radius $\Delta_\rho$ for the FW beam. See
  Refs. [[1]](#ref1) and [[2]](#ref2).
  - **Parameters:**
    - `k`, wavenumber $k$.
    - `spot_radius`, FW spot radius.
    - `asymptotic`, if true, uses the asymptotic approximation.
  - **Implementation:** It is computed using the expression
    $Q = \sqrt{(k^2 - c_0^2 / \Delta_\rho^2)}$ of the traditional method for
    wavenumbers calculation, where $c_0 = 2.4048$ in general or
    $c_0 = 3\pi / 4$ for asymptotic expansion approximation of the spot radius.
    Notice that for complex transverse wavenumber $h_q$,
    $\mathrm{Im}(h_q) \ll \mathrm{Re}(h_q)$, or at least
    $\mathrm{Im}(h_q) \leq 2\mathrm{Re}(h_q)/3\pi$ must be satisfied to avoid
    the infinite behavior of the Bessel functions.
</details>

<details>
  <summary>
    <code><b>fw_Q_from_spot_radius_purely_real_h(k, spot_radius, asymptotic)</b></code>
  </summary>

  - **Description:** Returns the parameter $Q$ of the Frozen Wave (FW)
  technique for a given spot radius $\Delta_\rho$ for the FW beam, considering
  purely real transverse wavenumbers. See Ref. [[3]](#ref3).
  - **Parameters:**
    - `k`, wavenumber $k$.
    - `spot_radius`, FW spot radius.
    - `asymptotic`, if true, uses the asymptotic approximation.
  - **Implementation:** It is computed using the expression
    $Q = \sqrt{(k^2 - c_0^2 / \Delta_\rho^2)}$ of the purely real $h_q$ method
    for wavenumbers calculation, where $c_0 = 2.4048$ in general or
    $c_0 = 3\pi /4$ for asymptotic expansion approximation of the spot radius.
</details>

<details>
  <summary>
    <code><b>fw_Q_from_spot_radius_paraxial_h(k, spot_radius, asymptotic)</b></code>
  </summary>

  - **Description:** Returns the parameter $Q$ of the Frozen Wave (FW)
  technique for a given spot radius $\Delta_\rho$ for the FW beam, considering
  the paraxial approximation. See Ref. [[3]](#ref3).
  - **Parameters:**
    - `k`, wavenumber $k$.
    - `spot_radius`, FW spot radius.
    - `asymptotic`, if true, uses the asymptotic approximation.
  - **Implementation:** It is computed using the expression
    $Q = \mathrm{Re}(k) - 0.5\mathrm{Re}(k)c_0^2/(|k|\Delta_\rho)^2$ of the
    paraxial $h_q$ method for wavenumbers calculation, where $c_0 = 2.4048$
    in general or $c_0 = 3\pi / 4$ for asymptotic expansion approximation of
    the spot radius.
</details>

<details>
  <summary>
    <code><b>fw_absorption_resistant_condition(N, beta)</b></code>
  </summary>

  - **Description:** Returns the Frozen Wave (FW) absorption-resistant
  condition. See Ref. [[4]](#ref4).
  - **Parameters:**
    - `N`, parameter $N$ of the FW.
    - `beta`, array of size `iqmax = 2 * N + 1`, containing the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
  - **Implementation:** It is computed using
    $[\mathrm{Im}(\beta_N) - \mathrm{Im}(\beta_{-N})] / \mathrm{Im}(\beta_0)$,
    where $\beta_q$ is the longitudinal wavenumber of the FW.
    The result must be $\ll 1$ in order to validate the approximation
    $\mathrm{Im}(\beta_q) \approx \mathrm{Im}(\beta_0)$ for all $q$, and so to
    guarantee the absorption-resistant property of the FW.
</details>

### Wavenumbers

<details>
  <summary>
    <code><b>fw_wavenumbers_traditional(N, Q, L, k0, nref, &k, beta, h)</b></code>
  </summary>

  - **Description:** Evaluates the wavenumbers $k$, $\beta_q$ and $h_q$ of a
  Frozen Wave (FW) using the traditional method. See Refs. [[1]](#ref1) and
  [[2]](#ref2).
  - **Parameters:**
    - `N`, parameter $N$ of the FW.
    - `Q`, parameter $Q$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `k0`, wavenumber $k_0$ with respect to the vacuum ($k_0 = \omega_0/c$,
    where $\omega_0$ is the angular frequency and $c$ the speed of light).
    - `nref`, refractive index.
    - `k`, to output the wavenumber $k$.
    - `beta`, array of size `iqmax = 2 * N + 1` to output the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
    - `h`, array of size `iqmax = 2 * N + 1` to output the transverse
    wavenumber $h_q$, where $h_q$ ↦ `h[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
  - **Implementation:** It is computed using the dispersion relationship
    $k^2=h_q^2+\beta_q^2$.
    In the traditional method, it is assumed real axicon angles for the
    Bessel beams, [i.e, $0 \leq \mathrm{Re}(\beta_q) \leq \mathrm{Re}(k)$ and
    $0 \leq \mathrm{Re}(h_q) \leq \mathrm{Re}(k)$ must be satisfied], to
    ensure $\mathrm{Re}(\beta_q) \geq 0$ and to finally find the
    relations $\mathrm{Im}(k) / \mathrm{Re}(k)$
    $= \mathrm{Im}(\beta_q) / \mathrm{Re}(\beta_q)$
    $= \mathrm{Im}(h_q) / \mathrm{Re}(h_q)$
    $= \mathrm{Im}(n_{ref}) / \mathrm{Re}(n_{ref})$.
    Notice that the condition $\mathrm{Im}(n_{ref}) \ll \mathrm{Re}(n_{ref})$,
    or at least $\mathrm{Im}(n_{ref}) \leq 2\mathrm{Re}(n_{ref})/3\pi$ must be
    satisfied in order to avoid the infinite behavior of the Bessel functions.
</details>

<details>
  <summary>
    <code><b>fw_wavenumbers_purely_real_h(N, Q, L, k0, nref, &k, beta, h)</b></code>
  </summary>

  - **Description:** Evaluates the wavenumbers $k$, $\beta_q$ and $h_q$ of a
  Frozen Wave (FW) considering the purely real transverse wavenumber $h_q$ method.
  See Ref. [[3]](#ref3).
  - **Parameters:**
    - `N`, parameter $N$ of the FW.
    - `Q`, parameter $Q$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `k0`, wavenumber $k_0$ with respect to the vacuum ($k_0 = \omega_0/c$,
    where $\omega_0$ is the angular frequency and $c$ the speed of light).
    - `nref`, refractive index.
    - `k`, to output the wavenumber $k$.
    - `beta`, array of size `iqmax = 2 * N + 1` to output the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
    - `h`, array of size `iqmax = 2 * N + 1` to output the transverse
    wavenumber $h_q$, where $h_q$ ↦ `h[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
  - **Implementation:** It is computed using the dispersion relationship
    $k^2=h_q^2+\beta_q^2$ for purely real $h_q$.
    In the purely real $h_q$ method, the condition
    $0 \leq \textrm{Re}(\beta_q) \leq \textrm{Re}(k)$ must be also satisfied
    in order to ensure $\textrm{Re}(\beta_q) \geq 0$ and
    $\textrm{Im}(h_q) = 0$.
</details>

<details>
  <summary>
    <code><b>fw_wavenumbers_paraxial_h(N, Q, L, k0, nref, &k, beta, h)</b></code>
  </summary>

  - **Description:** Evaluates the wavenumbers $k$, $\beta_q$ and $h_q$ of a
  Frozen Wave (FW) using the paraxial approximation for the transverse
  wavenumber $h_q$. See Ref. [[3]](#ref3).
  - **Parameters:**
    - `N`, parameter $N$ of the FW.
    - `Q`, parameter $Q$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `k0`, wavenumber $k_0$ with respect to the vacuum ($k_0 = \omega_0/c$,
    where $\omega_0$ is the angular frequency and $c$ the speed of light).
    - `nref`, refractive index.
    - `k`, to output the wavenumber $k$.
    - `beta`, array of size `iqmax = 2 * N + 1` to output the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
    - `h`, array of size `iqmax = 2 * N + 1` to output the transverse
    wavenumber $h_q$, where $h_q$ ↦ `h[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
  - **Implementation:** In the paraxial $h_q$ method, the condition
    $0 \leq \mathrm{Re}(\beta_q) \leq \mathrm{Re}(k)$ must be also satisfied
    in order to ensure $\mathrm{Re}(\beta_q) \geq 0$ and
    $\mathrm{Im}(h_q) = 0$ (notice this also implies in purely real $h_q$).
</details>

### Coefficients

<details>
  <summary>
    <code><b>fw_A_coefficient(N, L, beta, F, izmax, A, absorption_resistant)</b></code>
  </summary>

  - **Description:** Evaluates the $A_q$ coefficients of a Frozen Wave (FW),
  resistant or not to a possible absorption due to a medium with complex
  refractive index. See Refs [[1]](#ref1), [[2]](#ref2), [[3]](#ref3) and
  [[4]](#ref1).
  - **Parameters:**
    - `N`, parameter $N$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `beta`, array of size `iqmax = 2 * N + 1` containing the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
    - `F`, array of size `izmax` containing the morphological function $F(z)$,
    where $F(z)$ ↦ `F[iz]`, $z$ ↦ `iz * L / (izmax - 1)`, and
    `0<=iz<izmax`.
    - `izmax`, size for the `F[iz]` array, where `0<=iz<izmax`.
    - `A`, array of size `iqmax = 2 * N + 1` to output the $A_q$ coefficients,
    where $A_q$ ↦ `A[iq]`, $q$ ↦ `iq - N`, and `0<=iq<iqmax`.
    - `absorption_resistant`, if false, the FW is not absorption-resistant.
  - **Implementation:** The
  $` A_q = \int_0^L F(z) e^{-\bar{\beta_0}} e^{j 2 \pi q z / L} \mathrm{d}z `$,
  where $\bar{\beta}_0 = \mathrm{Im}(\beta_0)$ if
  `absorption_resistant = true` or
  $\bar{\beta}_0 = 0$ if `absorption_resistant = false`,
  is computed by means of an approximation of the
  integral by the trapezoidal method for equally spaced $z$ for a given
  morphologial function $F(z)$.
  Notice that for purely real refractive indices, the absorption-resistant
  concerns may be ignored.
</details>

<details>
  <summary>
    <code><b>fw2d_A_coefficient_restricted(P, N, L, beta, F, ixmax, iymax, izmax, A, absorption_resistant)</b></code>
  </summary>

  - **Description:** Evaluates the $A_{q,p}$ coefficients of a Frozen Wave (FW)
  2D made by $P$ FWs, restricted to the case where all FWs share the same $N$,
  $Q$ and $L$ (i.e., for $1 \leq p \leq P$, $N_p = N$, $Q_p = Q$, and
  $L_p = L$), resistant or not to a possible absorption due to a medium with
  complex refractive index. See Refs. [[5]](#ref5), [[6]](#ref6) and
  [[7]](#ref7).
  - **Parameters:**
    - `P`, number $P$ of FWs in the transverse direction $y$.
    - `N`, parameter $N$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `beta`, array of size `iqmax = 2 * N + 1` containing the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
    - `F`, array of size `iymax*izmax` containing the morphological function
    $F(y,z)$, where $F(y,z)$ ↦ `F[iz + izmax*iy]`,
    $y$ ↦ `iy`, $z$ ↦ `iz * L / (izmax - 1)`, `0<=iy<iymax`,
    and `0<=iz<izmax`.
    - `iymax`, size contribution in `iy` for the `F[iz + izmax*iy]` array,
    where `0<=iy<iymax`.
    - `izmax`, size contribution in `iz` for the `F[iz + izmax*iy]` array,
    where `0<=iz<izmax`.
    - `A`, array of size `ipmax*iqmax = (P) * (2 * N + 1)` to output
    the $A_{q,p}$ coefficients, where
    $A_{q,p}$ ↦ `A[iq + iqmax*ip]`, $p$ ↦ `ip + 1`, $q$ ↦ `iq - N`,
    `0<=ip<ipmax`, and `0<=iq<iqmax`.
    - `absorption_resistant`, if false, the FW is not absorption-resistant.
  - **Implementation:**  The
  $` A_{q,p} = \int_0^L F(y=y_{0p},z) e^{-\bar{\beta_0}} e^{j 2 \pi q z / L} \mathrm{d}z `$,
  where $y_{0p}$ ↦ `iy * (ipmax - 1) / (iymax - 1)`,
  $\bar{\beta}_0 = \mathrm{Im}(\beta_0)$ if
  `absorption_resistant = true` or
  $\bar{\beta}_0 = 0$ if `absorption_resistant = false`,
  is computed by means of an approximation of the
  integral by the trapezoidal method for equally spaced $z$ for a given
  morphologial function $F(y,z)$.
  Notice that for purely real refractive indices, the absorption-resistant
  concerns may be ignored.
</details>

<details>
  <summary>
    <code><b>fw3d_A_coefficient_restricted(S, P, N, L, beta, F, ixmax, iymax, izmax, A, absorption_resistant)</b></code>
  </summary>

  - **Description:** Evaluates the $A_{q,ps}$ coefficients of a Frozen Wave
  (FW) 3D made by $S \times P$ FWs, restricted to the case where all FWs
  share the same $N$, $Q$ and $L$ (i.e., for $1 \leq p \leq P$ and
  $1 \leq s \leq S$, $N_{ps} = N$, $Q_{ps} = Q$, and $L_{ps} = L$), resistant
  or not to a possible absorption due to a medium with complex refractive
  index. See Refs. [[7]](#ref7) and [[8]](#ref8).
  - **Parameters:**
    - `S`, number $S$ of FWs in the transverse direction $x$.
    - `P`, number $P$ of FWs in the transverse direction $y$.
    - `N`, parameter $N$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `beta`, array of size `iqmax = 2 * N + 1` containing the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0<=iq<iqmax`.
    - `F`, array of size `ixmax*iymax*izmax` containing the morphological
    function $F(x,y,z)$,
    where $F(x,y,z)$ ↦ `F[iz + izmax*iy + izmax*iymax*ix]`,
    $x$ ↦ `ix`, $y$ ↦ `iy`, $z$ ↦ `iz * L / (izmax - 1)`,
    `0<=ix<ixmax`, `0<=iy<iymax`, and `0<=iz<izmax`.
    - `ixmax`, size contribution in `ix` for the
    `F[iz + izmax*iy + izmax*iymax*ix]` array, where `0<=ix<ixmax`.
    - `iymax`, size contribution in `iy` for the
    `F[iz + izmax*iy + izmax*iymax*ix]` array, where `0<=iy<iymax`.
    - `izmax`, size contribution in `iz` for the
    `F[iz + izmax*iy + izmax*iymax*ix]` array, where `0<=iz<izmax`.
    - `A`, array of size `ismax*ipmax*iqmax = (S) * (P) * (2 * N + 1)`
    to output the $A_{q,sp}$ coefficients,
    where $A_{q,sp}$ ↦ `A[iq + iqmax*ip + iqmax*ipmax*is]`, $s$ ↦ `is + 1`,
    $p$ ↦ `ip + 1`, $q$ ↦ `iq - N`, `0<=is<ismax`, `0<=ip<ipmax`,
    and `0<=iq<iqmax`.
    - `absorption_resistant`, if false, the FW is not absorption-resistant.
  - **Implementation:**  The
  $` A_{q,ps} = \int_0^L F(x=x_{0s}, y=y_{0p}, z) e^{-\bar{\beta_0}} e^{j 2 \pi q z / L} \mathrm{d}z `$,
  where $x_{0s}$ ↦ `ix * (ismax - 1) / (ixmax - 1)`,
  $y_{0p}$ ↦ `iy * (ipmax - 1) / (iymax - 1)`,
  $\bar{\beta_0} = \mathrm{Im}(\beta_0)$ if absorption_resistant = true or
  $\bar{\beta_0} = 0$ if absorption_resistant = false,
  is computed by means of an approximation of the
  integral by the trapezoidal method for equally spaced $z$ for a given
  morphologial function $F(x,y,z)$.
  Notice that for purely real refractive indices, the absorption-resistant
  concerns may be ignored.
</details>

## How to use

This library is in a header-only style, i.e., there is nothing to build
(see the section [Compiling the library](#compiling-the-library) if you still
want to compile it).
Therefore, you only need to paste all the content of the
[include](include/) folder inside the include folder of your project (if you
do not have an include
folder in your project, paste the content inside the root folder of your
project). Finally, just write `#include "frozen-waves-library.h"` at the very
beginning of your code and you shall be ready to use the functions.

## Some C details

In this library, the implementation is carried out in terms of the C99
standards. Therefore, all the complex variables are handled using the
`double complex` type of the C `<complex.h>` library.

Notice that all functions, macros, constants and files whose names contain
the suffix `_impl_` are internal and are not intended to be used by users.

## Compatibility with C++

This library uses `__cplusplus` compiler guards with `extern "C"` and
macros to ensure C++ compatibility (C++98 standard at least).

In this sense, when using C++ compilers, the following C functions are
automatically mapped to their C++ equivalent: `creal(z)`↦`std::real(z)`,
`cimag(z)`↦`std::imag(z)`, `cabs(z)`↦`std::abs(z)`,
`cexp(z)`↦`std::exp(z)`, `sin(z)`↦`std::sin(z)`, `cos(z)`↦`std::cos(z)`, and
`casin(z)`↦`std::asin(z)`; and all the complex values are handled by means of
the `std::complex<double>` type of the C++ `<complex>` library.

## Compiling the library

As aforementioned, usually it is not necessary to compile the library.
However, in any case, the [src](src/) folder contains the files
[frozen-waves-library-declarations.c](src/frozen-waves-library-declarations.c) and
[frozen-waves-library-declarations.cpp](src/frozen-waves-library-declarations.cpp), with
declarations of all functions in respectively C and C++, and the file
[frozen-waves-library.c](src/frozen-waves-library.c), which is a C wrapper
([compatible with C++](#compatibility-with-c)) that may be used for
compilation.

The following are examples of how to compile this library using C and C++
compilers.

<details>
  <summary>
    <b>Compiling on Windows with MinGW gcc</b>
  </summary>

  ```bash
  gcc -shared -o src/frozen-waves-library.dll src/frozen-waves-library.c -Iinclude
  ```
</details>

<details>
  <summary>
    <b>Compiling on Windows with MinGW g++</b>
  </summary>

  ```bash
  g++ -shared -o src/frozen-waves-library.dll src/frozen-waves-library.c -Iinclude
  ```
</details>

<details>
  <summary>
    <b>Compiling on Linux/macOS with gcc</b>
  </summary>

  ```bash
  gcc -shared -fPIC -o src/frozen-waves-library.so src/frozen-waves-library.c -Iinclude
  ```
</details>

<details>
  <summary>
    <b>Compiling on Linux/macOS with g++</b>
  </summary>

  ```bash
  g++ -shared -fPIC -o src/frozen-waves-library.so src/frozen-waves-library.c -Iinclude
  ```
</details>

<details>
  <summary>
    <b>Compiling on Windows with MSVC cl</b>
  </summary> 

  Compiling this library with MSVC targeting the C language is discouraged
  because MSVC does not support the `double complex` type from the C
  `<complex.h>` library. However, you may build it in MSVC
  targeting the C++ language (e.g., using the `/TP` flag).

  ```bash
  cl /TP src/frozen-waves-library.c
  ```
</details>

## Other programming languages

Once compiled, it is also possible to use this library together with other
programming languages.

The following is an example on how to load the C compiled library in Python
using `numpy` and `cffi`.

<details>
  <summary>
    <b>Example of usage in Python</b>
  </summary>

```python
import numpy as np
from cffi import FFI

ffi = FFI()

# Read the C functions declarations
with open("src/frozen-waves-library-declarations.c", "r") as f:
    ffi.cdef(f.read())

# Import the compiled file
fwl = ffi.dlopen("src/frozen-waves-library.so") # for Linux/macOS
# fwl = ffi.dlopen("src/frozen-waves-library.dll") # for Windows

# Use NumPy complex128 to declare k
k = np.complex128(9.94e6 - 3.18j)

# Call a function of the library and print the result
Q = fwl.fw_Q_from_spot_radius_traditional(k, 10e-6, 0)
print("Result: ", Q)
```
</details>

## Change log

Refer to the [CHANGELOG.md](CHANGELOG.md) file for the latest updates.

## Authorship

The codes and routines were developed and are updated by
<a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">
Jhonas O. de Sarro</a> ([@jodesarro](https://github.com/jodesarro)).

## Acknowledgements

The author is very grateful for the collaborations with Professor
<a href="https://www.researchgate.net/profile/Leonardo-Ambrosio">Leonardo A.
Ambrosio</a> of <a href="http://www.sel.eesc.usp.br/leonardo">Applied
Electromagnetics Group (AEG)</a> from University of São Paulo (USP).

## Licensing

This project is protected under [MIT License](LICENSE). A copy of the license
is available at
[include/frozen-waves-library/license.txt](include/frozen-waves-library/license.txt).

## References

<a id="ref1"></a>
[1] M. Zamboni-Rached, "Stationary optical wave fields with arbitrary
longitudinal shape by superposing equal frequency Bessel beams:
Frozen Waves," Optics Express, vol. 12, no. 17, pp. 4001--4006, Aug. 2004, 
[doi: 10.1364/OPEX.12.004001](https://doi.org/10.1364/OPEX.12.004001).

<a id="ref2"></a>
[2] M. Zamboni-Rached, E. Recami, and H. E. Hernández-Figueroa, "Theory of
'frozen waves': modeling the shape of stationary wave fields," Journal of the
Optical Society of America A, vol. 22, no. 11, pp. 2465–2475, Nov. 2005,
[doi: 10.1364/JOSAA.22.002465](https://doi.org/10.1364/JOSAA.22.002465).

<a id="ref3"></a>
[3] M. Zamboni-Rached and M. Mojahedi, "Shaping finite-energy diffraction-
and attenuation-resistant beams through Bessel-Gauss beam superposition,"
Physical Review A, vol. 92, no. 4, p. 043839, Oct. 2015,
[doi: 10.1103/PhysRevA.92.043839](https://doi.org/10.1103/PhysRevA.92.043839).

<a id="ref4"></a>
[4] M. Zamboni-Rached, "Diffraction-Attenuation resistant beams in
absorbing media," Optics Express, vol. 14, no. 5, pp. 1804–1809, Mar. 2006,
[doi: 10.1364/OE.14.001804](https://doi.org/10.1364/OE.14.001804).

<a id="ref5"></a>
[5] L. A. Ambrosio, "Millimeter-structured nondiffracting surface beams,"
Journal of the Optical Society of America B, vol. 36, no. 3, pp. 638–645,
Feb. 2019,
[doi: 10.1364/JOSAB.36.000638](https://doi.org/10.1364/JOSAB.36.000638).

<a id="ref6"></a>
[6] J. O. de Sarro and L. A. Ambrosio, "Surface beams resistant to
diffraction and attenuation and structured at the millimeter scale,"
Journal of the Optical Society of America B, vol. 38, no. 3, pp. 677–684,
Mar. 2021, [doi: 10.1364/JOSAB.412756](https://doi.org/10.1364/JOSAB.412756).

<a id="ref7"></a>
[7] A. H. Dorrah et al., "Light sheets for continuous-depth holography and
three-dimensional volumetric displays," Nature Photonics, vol. 17, pp.
427–434, Apr. 2023,
[doi: 10.1038/s41566-023-01188-y](https://doi.org/10.1038/s41566-023-01188-y).

<a id="ref8"></a>
[8] Still needs an official reference, but FW 3D is a superposition of FWs
with origins displaced in the plane transverse to the propagation
direction.

A copy of this list is available at
[include/frozen-waves-library/references.txt](include/frozen-waves-library/references.txt).

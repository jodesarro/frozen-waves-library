# Frozen Waves Library: A C library providing routines for computing quantities related to Frozen Waves

This project provides a C library designed for the
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
    $c_0 = 2.4048...$ in general or $c_0 = 3\pi/4$ in the asymptotic
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
  longitudinal range. It returns `NAN` if the result is not purely real.
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
    <code><b>bb_aperture_radius_min(h, asymptotic)</b></code>
  </summary>

  - **Description:** Returns $R_{min}$, the minimum aperture radius possible
  for an experimental generation of a Bessel Beam (BB) of complex transverse
  wavenumber $h$.
  - **Parameter:**
    - `h`, transverse wavenumber of the BB.
    - `asymptotic`, if true, uses the asymptotic approximation.
  - **Implementation:** Returns `bb_spot_radius()`.
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
  $0 \leq Q \pm 2\pi N/L \leq \mathrm{Re}(k)$.
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
    wavenumbers calculation, where $c_0 = 2.4048...$ in general or
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
  - **Implementation:** Returns `fw_Q_from_spot_radius_traditional()`.
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
    `0 <= iq < iqmax`.
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
    `0 <= iq < iqmax`.
    - `h`, array of size `iqmax = 2 * N + 1` to output the transverse
    wavenumber $h_q$, where $h_q$ ↦ `h[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
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
    `0 <= iq < iqmax`.
    - `h`, array of size `iqmax = 2 * N + 1` to output the transverse
    wavenumber $h_q$, where $h_q$ ↦ `h[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
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
    `0 <= iq < iqmax`.
    - `h`, array of size `iqmax = 2 * N + 1` to output the transverse
    wavenumber $h_q$, where $h_q$ ↦ `h[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
  - **Implementation:** It is computed using the dispersion relationship
    $h_q=\sqrt{2} k \sqrt{1-\beta_q/k}$ (paraxial approximation).
    In the paraxial $h_q$ method, the condition
    $0 \leq \mathrm{Re}(\beta_q) \leq \mathrm{Re}(k)$ must be also satisfied
    in order to ensure $\mathrm{Re}(\beta_q) \geq 0$ and
    $\mathrm{Im}(h_q) = 0$ (notice this also implies in purely real $h_q$).
</details>

### Coefficients

<details>
  <summary>
    <code><b>fw_A(N, L, beta, F, izmax, A, absorption_resistant)</b></code>
  </summary>

  - **Description:** Evaluates the $A_q$ coefficients of a single Frozen Wave
  (FW), resistant or not to a possible absorption due to a medium with complex refractive index. See Refs. [[1]](#ref1), [[2]](#ref2) and [[3]](#ref3).
  - **Parameters:**
    - `N`, parameter $N$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `beta`, array of size `iqmax = 2 * N + 1` containing the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
    - `F`, array of size `izmax` containing the morphological
    function $F(z)$, where $F(z)$ ↦ `F[iz]`, and `0 <= iz < izmax`.
    - `izmax`, size contribution in `iz` for the `F[iz]` array, where `0 <= iz < izmax`.
    - `A`, array of size `iqmax = 2 * N + 1` to output the $A_q$ coefficients,
    where $A_q$ ↦ `A[iq]`, $q$ ↦ `iq - N`, and `0 <= iq < iqmax`.
    - `absorption_resistant`, if false, the FW is not absorption-resistant.
  - **Implementation:**  The
  $` A_q = \int_0^L F(z) e^{-\bar{\beta_0}} e^{i 2 \pi q z / L} \mathrm{d}z `$,
  where $\bar{\beta_0} = \mathrm{Im}(\beta_0)$ if `absorption_resistant = true` or $\bar{\beta_0} = 0$ if `absorption_resistant = false`,
  is computed through the function `fw3d_A_restricted()` with `P = S = ixmax = iymax = 1`.
</details>

<details>
  <summary>
    <code><b>fw3d_A_restricted(P, S, N, L, beta, F, ixmax, iymax, izmax, A, absorption_resistant)</b></code>
  </summary>

  - **Description:** Evaluates the $A_{q,ps}$ coefficients of a Frozen Wave
  (FW) 3D made by $P \times S$ FWs, restricted to the case where all FWs
  share the same $N$, $Q$ and $L$ (i.e., for $1 \leq p \leq P$ and
  $1 \leq s \leq S$, $N_{ps} = N$, $Q_{ps} = Q$, and $L_{ps} = L$), resistant
  or not to a possible absorption due to a medium with complex refractive
  index. See Refs. [[7]](#ref7) and [[8]](#ref8).
  - **Parameters:**
    - `P`, number $P$ of FWs in the transverse direction $x$.
    - `S`, number $S$ of FWs in the transverse direction $y$.
    - `N`, parameter $N$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `beta`, array of size `iqmax = 2 * N + 1` containing the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
    - `F`, array of size `ixmax * iymax * izmax` containing the morphological
    function $F(x,y,z)$,
    where $F(x,y,z)$ ↦ `F[iz + izmax*iy + izmax*iymax*ix]`,
    $x$ ↦ `ix`, $y$ ↦ `iy`, $z$ ↦ `iz * L / (izmax - 1)`,
    `0 <= ix < ixmax`, `0 <= iy < iymax`, and `0 <= iz < izmax`.
    - `ixmax`, size contribution in `ix` for the
    `F[iz + izmax*iy + izmax*iymax*ix]` array, where `0 <= ix < ixmax`.
    - `iymax`, size contribution in `iy` for the
    `F[iz + izmax*iy + izmax*iymax*ix]` array, where `0 <= iy < iymax`.
    - `izmax`, size contribution in `iz` for the
    `F[iz + izmax*iy + izmax*iymax*ix]` array, where `0 <= iz < izmax`.
    - `A`, array of size `ipmax * ismax * iqmax = P * S * (2 * N + 1)`
    to output the $A_{q,ps}$ coefficients,
    where $A_{q,ps}$ ↦ `A[iq + iqmax*is + iqmax*ismax*ip]`, $p$ ↦ `ip + 1`,
    $s$ ↦ `is + 1`, $q$ ↦ `iq - N`, `0 <= ip < ipmax = P`, `0 <= is < ismax = S`,
    and `0 <= iq < iqmax = 2 * N + 1`.
    - `absorption_resistant`, if false, the FW is not absorption-resistant.
  - **Implementation:**  The
  $` A_{q,ps} = \int_0^L F(x=x_{0p}, y=y_{0s}, z) e^{-\bar{\beta_0}} e^{i 2 \pi q z / L} \mathrm{d}z `$,
  where $x_{0p}$ ↦ `ix * (ipmax - 1) / (ixmax - 1)`,
  $y_{0s}$ ↦ `iy * (ismax - 1) / (iymax - 1)`,
  $\bar{\beta_0} = \mathrm{Im}(\beta_0)$ if `absorption_resistant = true` or
  $\bar{\beta_0} = 0$ if `absorption_resistant = false`,
  is computed by means of an approximation of the
  integral by the trapezoidal method for equally spaced $z$ for a given
  morphologial function $F(x,y,z)$.
  Notice that for purely real refractive indices, the absorption-resistant
  concerns may be ignored.
</details>

### Wave functions

<details>
  <summary>
    <code><b>fw_psi(N, beta, h, A, xmin, xmax, xpoints, ymin, ymax, ypoints, zmin, zmax, zpoints, psi)</b></code>
  </summary>

  - **Description:** Evaluates the scalar wave function $\psi(x,y,z)$ of a single Frozen Wave (FW).
  See Refs. [[1]](#ref1), [[2]](#ref2) and [[3]](#ref3).
  - **Parameters:**
    - `N`, parameter $N$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `beta`, array of size `iqmax = 2 * N + 1` containing the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
    - `h`, array of size `iqmax = 2 * N + 1` containing the transverse
    wavenumber $h_q$, where $h_q$ ↦ `h[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
    - `A`, array of size `iqmax = 2 * N + 1` containing the
    $A_q$ coefficients, where $A_q$ ↦ `A[iq]`, $q$ ↦ `iq - N`, and `0 <= iq < iqmax`.
    - `xmin`, minimum value for the $x$ coordinate, where $x \geq$ `xmin`.
    - `xmax`, maximum value for the $x$ coordinate, where $x \leq$ `xmax`.
    - `xpoints`, number of points in the $x$ direction, where `0 <= ix < xpoints`.
    - `ymin`, minimum value for the $y$ coordinate, where $y \geq$ `ymin`.
    - `ymax`, maximum value for the $y$ coordinate, where $y \leq$ `ymax`.
    - `ypoints`, number of points in the $y$ direction, where `0 <= iy < ypoints`.
    - `zmin`, minimum value for the $z$ coordinate, where $z \geq$ `zmin`.
    - `zmax`, maximum value for the $z$ coordinate, where $z \leq$ `zmax`.
    - `zpoints`, number of points in the $z$ direction, where `0 <= iz < zpoints`.
    - `psi`, array of size `xpoints * ypoints * zpoints` to output the scalar wave
    function $\psi(x,y,z)$, where $\psi(x,y,z)$ ↦ `psi[iz + zpoints*iy + zpoints*ypoints*ix]`, $x$ ↦ `xmin + ix * (xmax - xmin) / (xpoints - 1)`, $y$ ↦ `ymin + iy * (ymax - ymin) / (ypoints - 1)`, $z$ ↦ `zmin + iz * (zmax - zmin) / (zpoints - 1)`.
  - **Implementation:** The $` \psi(x,y,z) = \sum_{q=-N}^N A_q \ J_0(h_q\sqrt{x^2 + y^2}) \ \exp(-i \beta_q z) `$ is computed through the function `fw3d_psi_restricted()` with `S = P = 1` and
  `x0[0] = y0[0] = z0[0] = 0`.
  - **Details:** The computation is parallelized with OpenMP. If the macro `FROZEN_WAVES_LIBRARY_PROGRESS_BAR` is defined, a progress bar is displayed on
  the console.
</details>

<details>
  <summary>
    <code><b>fw3d_psi_restricted(P, S, N, beta, h, A, x0, y0, z0, xmin, xmax, xpoints, ymin, ymax, ypoints, zmin, zmax, zpoints, psi)</b></code>
  </summary>

  - **Description:** Evaluates the scalar wave function $\psi(x,y,z)$ of a Frozen Wave (FW) 3D made by
  $P \times S$ FWs, restricted to the case where all FWs share the same $N$, $Q$ and $L$
  (i.e., for $1 \leq p \leq P$ and $1 \leq s \leq S$, $N_{ps} = N$, $Q_{ps} = Q$, and $L_{ps} = L$).
  See Refs. [[7]](#ref7) and [[8]](#ref8).
  - **Parameters:**
    - `P`, number $P$ of FWs in the transverse direction $x$.
    - `S`, number $S$ of FWs in the transverse direction $y$.
    - `N`, parameter $N$ of the FW.
    - `L`, parameter $L$ of the FW.
    - `beta`, array of size `iqmax = 2 * N + 1` containing the longitudinal
    wavenumber $\beta_q$, where $\beta_q$ ↦ `beta[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
    - `h`, array of size `iqmax = 2 * N + 1` containing the transverse
    wavenumber $h_q$, where $h_q$ ↦ `h[iq]`, $q$ ↦ `iq - N`, and
    `0 <= iq < iqmax`.
    - `A`, array of size `ipmax * ismax * iqmax = P * S * (2 * N + 1)` containing the
    $A_{q,ps}$ coefficients, where $A_{q,ps}$ ↦ `A[iq + iqmax*is + iqmax*ismax*ip]`, $p$ ↦
    `ip + 1`, $s$ ↦ `is + 1`, $q$ ↦ `iq - N`, `0 <= ip < ipmax = P`, `0 <= is < ismax = S`,
    and `0 <= iq < iqmax = 2 * N + 1`.
    - `x0`, array of size `ipmax = P` containing the $x$-coordinates of the origins of
    the FWs, where $x_{0p}$ ↦ `x0[ip]`, $p$ ↦ `ip + 1`, `0 <= ip < ipmax = P`.
    - `y0`, array of size `ismax = S` containing the $y$-coordinates of the origins of
    the FWs, where $y_{0s}$ ↦ `y0[is]`, $s$ ↦ `is + 1`, `0 <= is < ismax = S`.
    - `z0`, array of size `ipmax * ismax = P * S` containing the $z$-coordinates of the
    origins of the FWs, where $z_{0ps}$ ↦ `z0[is + ismax*ip]`, $p$ ↦ `ip + 1`, $s$ ↦ `is + 1`, `0 <= ip < ipmax = P`, `0 <= is < ismax = S`.
    - `xmin`, minimum value for the $x$ coordinate, where $x \geq$ `xmin`.
    - `xmax`, maximum value for the $x$ coordinate, where $x \leq$ `xmax`.
    - `xpoints`, number of points in the $x$ direction, where `0 <= ix < xpoints`.
    - `ymin`, minimum value for the $y$ coordinate, where $y \geq$ `ymin`.
    - `ymax`, maximum value for the $y$ coordinate, where $y \leq$ `ymax`.
    - `ypoints`, number of points in the $y$ direction, where `0 <= iy < ypoints`.
    - `zmin`, minimum value for the $z$ coordinate, where $z \geq$ `zmin`.
    - `zmax`, maximum value for the $z$ coordinate, where $z \leq$ `zmax`.
    - `zpoints`, number of points in the $z$ direction, where `0 <= iz < zpoints`.
    - `psi`, array of size `xpoints * ypoints * zpoints` to output the scalar wave
    function $\psi(x,y,z)$, where $\psi(x,y,z)$ ↦ `psi[iz + zpoints*iy + zpoints*ypoints*ix]`, $x$ ↦ `xmin + ix * (xmax - xmin) / (xpoints - 1)`, $y$ ↦ `ymin + iy * (ymax - ymin) / (ypoints - 1)`, $z$ ↦ `zmin + iz * (zmax - zmin) / (zpoints - 1)`.
  - **Implementation:** The scalar wave function is computed through the formula
  $` \psi(x,y,z) = \sum_{p=1}^P \sum_{s=1}^S \sum_{q=-N}^N A_{q,ps} \ J_0\left[h_q\sqrt{(x-x_{0p})^2 + (y-y_{0s})^2}\right] \ \exp\left[-i \beta_q (z - z_{0ps})\right] `$.
  - **Details:** The computation is parallelized with OpenMP. If the macro `FROZEN_WAVES_LIBRARY_PROGRESS_BAR` is defined, a progress bar is displayed on
  the console.
</details>

## How to use

This library is in a header-only style, i.e., there is nothing to build.
If you need a binary, you may follow the instructions in the
[Compiling the library](#compiling-the-library).

Otherwise, you only need to paste all the content of the
[include](include/) folder
inside the include folder of your project (if you do not have an include
folder in your project, paste the content inside the root folder of your
project).

Finally, just write `#include "frozen-waves-library.h"` at the very
beginning of your code and you shall be ready to use the functions.

## Some C details

In this library, the implementation is carried out in terms of the C99
standards. Therefore, all the complex variables are handled using the
`double complex` type of the C `<complex.h>` library.

Notice that all functions, macros, constants and files whose names contain
the suffix `_impl_` are internal and are not intended to be used by users.

The computation of some functions is parallelized with OpenMP. In this sense,
you need to use the flag `-fopenmp` when compiling your code.

## Compiling the library

As aforementioned, usually it is not necessary to compile the library.

However, in any case, the [src](src/) folder contains the file
[frozen-waves-library.c](src/frozen-waves-library.c), which is a C wrapper
that may be used for compilation.

The following are examples on how to compile this library using C compilers.

<details>
  <summary>
    <b>Compiling on Windows with MinGW gcc</b>
  </summary>

  ```bash
  gcc -fopenmp -shared -Iinclude src/frozen-waves-library.c -o frozen-waves-library.dll -Wl,--out-implib,libfrozen-waves-library.dll.a
  ```
</details>

<details>
  <summary>
    <b>Compiling on Linux with gcc</b>
  </summary>

  ```bash
  gcc -fopenmp -fPIC -shared -Iinclude src/frozen-waves-library.c -o libfrozen-waves-library.so
  ```
</details>

<details>
  <summary>
    <b>Compiling on Windows with MSVC cl</b>
  </summary> 

  Compiling this library with MSVC targeting the C language is discouraged
  because MSVC does not support the `double complex` type from the C
  `<complex.h>` library.
</details>

When using a compiled file of the library into C projects,
you must paste all the content of the
[include](include/) folder inside the include folder of your project,
and then write `#define FROZEN_WAVES_LIBRARY_IMPORTS`
before the `#include "frozen-waves-library.h"`.

## Macros

There are some macros that may be used in this library.

- `FROZEN_WAVES_LIBRARY_PROGRESS_BAR`: define this macro to print a progress bar on the console while running the function `fw3d_psi_restricted()`.
- `FROZEN_WAVES_LIBRARY_IMPORTS`: always and only define this macro when using this library through a compiled file.

All macros must be defined before the inclusion of the header of this library, i.e., you must write `#define FROZEN_WAVES_LIBRARY_PROGRESS_BAR` and/or `#define FROZEN_WAVES_LIBRARY_IMPORTS` before `#include "frozen-waves-library.h"`.

## Third Parties

This library makes use of the following third-party libraries, codes, or routines, already incorporated:

- Bessel Library: A C library with routines for computing Bessel functions. More information available at [<https://github.com/jodesarro/bessel-library>](https://github.com/jodesarro/bessel-library).

- Console Progress Bar: A C code for printing the progress of numeric iterations on the console. More information available at [<https://github.com/jodesarro/console-progress-bar>](https://github.com/jodesarro/console-progress-bar).

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

This project is protected under [MIT License](LICENSE).

[Third parties](#third-parties) may have their own [license](LICENSE).

## References

<a id="ref1"></a>
[1] M. Zamboni-Rached, "Stationary optical wave fields with arbitrary
longitudinal shape by superposing equal frequency Bessel beams:
Frozen Waves," Optics Express, vol. 12, no. 17, pp. 4001–4006, Aug. 2004, 
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
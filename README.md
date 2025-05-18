# Outline

<!-- toc -->

- [Outline](#outline)
- [DalitzPlot.jl](#dalitzplotjl)
  - [Installation](#installation)
  - [Usage](#usage)
- [Xs Package: for cross section and Dalitz plot](#xs-package-for-cross-section-and-dalitz-plot)
  - [Define amplitudes with factors for the calculation](#define-amplitudes-with-factors-for-the-calculation)
  - [Define the masses of initial and final particles](#define-the-masses-of-initial-and-final-particles)
  - [Define the momentum or total energy](#define-the-momentum-or-total-energy)
  - [Calculate](#calculate)
  - [Plot Dalitz Plot](#plot-dalitz-plot)
- [GEN Package: for Generating Events](#gen-package-for-generating-events)
- [FR Package: for Numerical Calculation of Feynman Rules](#fr-package-for-numerical-calculation-of-feynman-rules)
- [Basic conventions.](#basic-conventions)
  - [Dirac Gamma matrices.](#dirac-gamma-matrices)
  - [Functions for particles with different spins](#functions-for-particles-with-different-spins)
    - [Spinor for $S=1/2$](#spinor-for-s12)
    - [Polarized vector for $S=1$](#polarized-vector-for-s1)
    - [Rarita-Schwinger Spinor for $S=3/2$](#rarita-schwinger-spinor-for-s32)
    - [Levi-Civita tensor](#levi-civita-tensor)
  - [Additional defintions to \*](#additional-defintions-to-)

<!-- tocstop -->

# DalitzPlot.jl

This Julia package is designed for high-energy physics applications. It was originally developed for visualizing and analyzing particle decays using Dalitz plots, but its functionality has since expanded beyond Dalitz plots to support a wide range of features. It consists of the following subpackages (after v0.1.8, the package was divided into several subpackages):

- `Xs`: A module that provides utilities for cross section calculations and the generation of Dalitz or Dalitz-like plots, which are crucial for visualizing three-body and multi-body decays (up to 10 particles). Users can specify amplitudes to generate these plots.
- `GEN`: Stands for Event Generation. This subpackage is responsible for generating events, which users can then utilize for tasks such as cross section calculations or decay analysis. It is used internally by `Xs`.
- `FR`: Focuses on Feynman rules. This subpackage includes functions for calculating spinor or polarized vectors with momentum and spin, gamma matrices, and other related calculations.
- `qBSE`: Refers to the Quasipotential Bethe-Salpeter Equation. This subpackage provides tools for solving the Bethe-Salpeter equation, which is used in the study of bound states in quantum field theory. This is a specific theoretical model, and those not interested can disregard it. For those interested, please refer to `doc/qBSE.md`.
- `AUXs`: A collection of auxiliary functions designed to support and enhance the functionality of the main subpackages.

> **Note**: This package is primarily intended for theoretical and phenomenological studies. For experimental data analysis, please consider using specialized tools designed for that purpose.

## Installation

To install the "DalitzPlot" package, you can follow the standard Julia package manager procedure. Open Julia and use the following commands:

Using the Julia REPL

```julia
Pkg.add("DalitzPlot")
```

Alternatively, if you want to install it directly from the GitHub repository:

Using the Julia REPL:

```julia
import Pkg
Pkg.add(url="https://github.com/junhe1979/DalitzPlot.jl")
```

These commands will install the "DalitzPlot" package and allow you to use it in your Julia environment.

## Usage

After installation, the package can be used as:

```julia
using DalitzPlot
```

To use the subpackages:

```julia
using DalitzPlot.Xs
using DalitzPlot.GEN
using DalitzPlot.FR
using DalitzPlot.qBSE
using DalitzPlot.PLOT
using DalitzPlot.AUX
```

# Xs Package: for cross section and Dalitz plot

The cross section, denoted by $d\sigma$, can be expressed in terms of amplitudes, ${\mathcal M}$, as follows:

$d\sigma=F\frac{1}{S}|{\mathcal M}|^2d\Phi=(2\pi)^{4-3n}F\frac{1}{S}|{\mathcal M}|^2dR$

The decay width in CMS is given as

$\Gamma = \frac{1}{2M}\int d\Phi |{\mathcal M}|^2=\frac{ (2\pi)^{4-3n}}{2M}\int dR |{\mathcal M}|^2$.

Here, $dR = (2\pi)^{3n-4} d\Phi = \prod_{i}\frac{d^3k_i}{2E_i}\delta^4(\sum_{i}k_i-P)$ represents the Lorentz-invariant phase space for $n$ particles, and it is generated using the Monte-Carlo method described in Ref. [F. James, CERN 68-12].

The flux factor $F$ for the cross section is given by: $F=\frac{1}{2E2E'v_{12}}=\frac{1}{4[(p_1\cdot p_2)^2-m_1^2m_2^2]^{1/2}}    \frac{|p_1\cdot p_2|}{p_1^0p_2^0}$. In the laboratory or center of mass frame, the relation $\vec{p}_1^2 \vec{p}_2^2 = (\vec{p}_1 \cdot \vec{p}_2)^2$ is utilized. In the laboratory frame, the term $\frac{|p_1\cdot p_2|}{p_1^0 p_2^0}$ simplifies to 1.
Additionally, if a boson or zero-mass spinor particle is replaced with a non-zero mass spinor particle, the factor $1/2$ is replaced with the mass of the particle, $m$.
The total symmetry factor $S$ is given by $\prod_i m_i!$ if there are $m_i$ identical particles.

For decay width, the flux factor is modified to $F = \frac{1}{2E}$.

## Define amplitudes with factors for the calculation

Users are required to supply amplitudes with factors $(2\pi)^{4-3n}F|{\mathcal M}|^2\frac{1}{S}$ within the function, named `amp` here.

The simplest case is to assume it equals 1.

```julia
amp(tecm, kf, ch, para, p0)=1.
```

Define more intricate amplitudes for a 2->3 process.

This function, named `amp`, calculates amplitudes with factors for a 2->3 process. The input parameters are:

- `tecm`: Total energy in the center-of-mass frame.
- `kf`: Final momenta generated.
- `ch`: Information about the process (to be defined below).
- `para`: Additional parameters.
- `p0`: Additional parameters for possible fitting.

Users are expected to customize the amplitudes within this function according to their specific requirements.

```julia

function amps(tecm, kf, ch, para, p0)

    # get kf as momenta in the center-of-mass ,
    #k1,k2,k3=kf   
    #get kf as momenta in laboratory frame
    k1, k2, k3 = Xs.getkf(para.p, kf, ch)

    # Incoming particle momentum
    # Center-of-mass frame: p1 = [p 0.0 0.0 E1]
    #p1, p2 = pcm(tecm, ch.mi)
    # Laboratory frame
    p1, p2 = Xs.plab(para.p, ch.mi)

    #flux
    #flux factor for cross section
    fac = 1e9 / (4 * para.p * ch.mi[2] * (2 * pi)^5)

    k12 = k1 + k2
    s12 = k12 * k12
    m = 6.0
    A = 1 / (s12 - m^2 + im * m * 0.1)

    total = abs2(A) * fac * 0.389379e-3

    return total
end
```

## Define the masses of initial and final particles

The mass of initial and final particles is specified in a NamedTuple (named `ch` here) with fields `mi` and `mf`.
Particle names can also be provided for PlotD as `namei` and `namef`.

The function for amplitudes with factors is saved as `amp`.

Example usage:

```julia
ch = (pf=["p1","p2","p3"],mi=[mass_i_1, mass_i_2], mf=[mass_f_1, mass_f_2, mass_f_3], namei=["p^i_{1}", "p^i_{2}"], namef=["p^f_{1}", "p^f_{2}", "p^f_{3}"], amp=amp)
```

Make sure to replace `mass_i_1`, `mass_i_2`, `mass_f_1`, `mass_f_2`, and `mass_f_3` with the actual masses of the particles (1.0, 1.0, 1.0, 2.0, 3.0 here).

## Define the momentum or total energy

Momentum in the Laboratory frame and transfer it to the total energy in the center-of-mass frame.

Example usage:

```julia
p_lab = 20.0
tecm = Xs.pcm(p_lab, ch.mi)
tecm=10.
```

## Calculate

Calculate the cross section and related spectra using the GENEV function.

The function `Xsection` takes the momentum of the incoming particle in the Laboratory frame (`p_lab`), the information about the particles (`ch`), the axes representing invariant masses (`axes`), the total number of events (`nevtot`), the number of bins (`Nbin`), and additional parameters (`para`). The function uses the plab2pcm function to transform the momentum from the Laboratory frame to the center-of-mass frame.

Example usage:

```julia
using ProgressBars
function progress_callback(pb)
   ProgressBars.update(pb)  
end
nevtot=Int64(1e7)
pb = ProgressBar(1:nevtot)  
callback = i -> progress_callback(pb)  
res = Xs.Xsection(tecm, ch, callback,axes=[["p2","p3"], ["p1","p2"]], nevtot=Int64(1e7), Nbin=500, para=(p=p_lab, l=1.0),stype=2)
```

The calculation results are stored in the variable `res` as a `NamedTuple` with the following fields:

- `res.cs0`: The total cross section or decay width, calculated as $\int \text{amp} \, dR$. The result is given in units of GeV$^{-2}$. To convert to barns, use the relation $1\,\text{GeV}^{-2} = 0.3894\,\text{mb}$.
- `res.cs1`: The invariant mass spectrum, given by $d\int \text{amp} \, dR/dx$, where $x$ is the invariant mass $m_{ij}$ of the first two particles (`stype=1`) or the squared invariant mass $s_{ij} = m_{ij}^2$ (`stype=2`). The spectrum is binned according to the specified range and number of bins (`Nbin`).
- `res.cs2`: The Dalitz plot as $d\int \text{amp} \, dR/dxdy$.

These results allow you to analyze the total cross section, invariant mass distributions, and Dalitz plot for your process.

## Plot Dalitz Plot

```julia
DalitzPlot.PLOT.plotD(res)
```

<img src="test/DP.png" alt="描述文字" width="500" height="500">

# GEN Package: for Generating Events

The GEN package is used for generating events for cross-section calculations and Dalitz plots. The Lorentz-invariant phase space used here is defined as:

$dR = (2\pi)^{3n-4} d\Phi = \prod_{i}\frac{d^3k_i}{2E_i}\delta^4(\sum_{i}k_i-P)$

for $n$ particles. The events are generated using the Monte-Carlo method described in Ref. [F. James, CERN 68-15].

The primary function provided by this package is `GENEV`, which can be used as follows:

```julia
PCM, WT=GEN.GENEV(tecm,EM)
```

- Input:  total momentum in center of mass frame `tecm`, and the mass of particles `EM`.
  * `tecm`: a `Float64` value representing the total momentum in the center of mass frame.
  * `EM`: a `Vector{Float64}` containing the masses of the particles.
- Output: the momenta of the particles `PCM`, and a weight `WT`.
  * `PCM`: a StaticArrays `Vector{SVector{5,Float64}}` storing the momenta of the particles. Note that at most 18 particles can be considered.
  * `WT`:  `Float64` value representing the weight for this event.

# FR Package: for Numerical Calculation of Feynman Rules

# Basic conventions.

Since arrays in Julia are 1-indexed, a covariant 4-vector is represented as an `SVector{5, Type}(v1, v2, v3, v0, v5)`. The first three elements correspond to a 3-vector, the fourth element represents the time component (or the 0th component), and the fifth element represents mass in the case of momentum, but it is typically meaningless in most other contexts. The addition of the fifth element serves to distinguish it from the four-dimensional Dirac gamma matrices.

For example, a momentum is `SVector{5, Type}(kx,ky,kz,k0,m)`.

Minkowski metric is chosen as $g^{\mu\nu}=diag(1,-1,-1,-1)$. In the code, we still adopt above convention as

`g= SMatrix{5,5,Float64}([-1.0 0.0 0.0 0.0 0.0;0.0 -1.0 0.0 0.0 0.0;0.0 0.0 -1.0 0.0 0.0;0.0 0.0 0.0 1.0 0.0;0.0 0.0 0.0 0.0 0.0])`.

For example, $g^{00}$ is accessed as `FR.g[4,4]`.

## Dirac Gamma matrices.

We adopt the Dirac representation for the gamma matrices:

$$
\gamma_1=\left(\begin{array}{cccc}0&0& 0&1\\0&0&1&0\\0&-1&0&0\\-1&0&0&0\end{array}\right),
\gamma_2=\left(\begin{array}{cccc}0&0& 0&-i\\0&0&i&0\\0&i&0&0\\-i&0&0&0\end{array}\right),
\gamma_3=\left(\begin{array}{cccc}0&0& 1&0\\0&0&0&-1\\-1&0&0&0\\0&1&0&0\end{array}\right),
$$

$$
\gamma_0=\left(\begin{array}{cccc}1&0& 0&0\\0&1&0&0\\0&0&-1&0\\0&0&0&-1\end{array}\right),
\gamma_5=\left(\begin{array}{cccc}0&0& 1&0\\0&0&0&1\\1&0&0&0\\0&1&0&0\end{array}\right),
I=\left(\begin{array}{cccc}1&0& 0&0\\0&1&0&0\\0&0&1&0\\0&0&0&1\end{array}\right),
$$

The Dirac gamma matrices are represented by the array GA=$[\gamma_1,\gamma_2,\gamma_3,\gamma_0,\gamma_5]$, with each gamma matrix defined as `SMatrix{4,4,ComplexF64}`.

For example, $\gamma_2$ can be accessed as `FR.GA[2]` and has the type `SMatrix{4,4,ComplexF64}`. The unit matrix is accessed as   `FR.I`

A function `FR.GS` is provided for calculate $\gamma \cdot k$ as `function GS(k::SVector{5, Type})::SMatrix{4,4,ComplexF64}`

## Functions for particles with different spins

`l::Int64` in the followings is for the **helicity** of the particle.
`bar=true` means that output is $\bar{u}$ or $\bar{u}^\mu$.
`V=true` means that output is for antifermion $v$.
`star=true` is for the polarized vector with a complex conjugation.

### Spinor for $S=1/2$

`function U(k, l::Int64; bar=false, V=false)::SVector{4,ComplexF64}`

### Polarized vector for $S=1$

`function eps(k, l::Int64; star=false)::SVector{5,ComplexF64}`

### Rarita-Schwinger Spinor for $S=3/2$

`function U3(k, l::Int64; bar=false, V=false)::SVector{5,SVector{4, ComplexF64}}`

### Levi-Civita tensor

`function LC(a::SVector, b::SVector, c::SVector)`: $\epsilon_{\mu\nu\rho\lambda}a^\mu b^\nu c^\rho$.

`function LC(i0::Int64, i1::Int64, i2::Int64, i3::Int64)`: $\epsilon_{\mu\nu\rho\lambda}$

`function LC(a::SVector, b::SVector, c::SVector, d::SVector)`: $\epsilon_{\mu\nu\rho\lambda}a^\mu b^\nu c^\rho d^\lambda$.

## Additional defintions to *

More methods are added for multiplying of polarized vector, spinor, and gamma matrices.

$Q\cdot W$, the dot product of two four-vectors $Q$ and $W$  (for momentum, polarized vector):

`function *(Q::SVector{5,Float64}, W::SVector{5,Float64})::Float64`,

`function *(Q::SVector{5,Float64}, W::SVector{5,ComplexF64})::ComplexF64`,

`function *(Q::SVector{5,ComplexF64}, W::SVector{5,Float64})::ComplexF64`,

`function *(Q::SVector{5,ComplexF64}, W::SVector{5,ComplexF64})::ComplexF64`.

$AM$, A row vector $A$ multiplied by a matrix $M$ (for spinor and gamma matrices)：
`function *(A::SVector{4, ComplexF64}, M::SMatrix{4, 4, ComplexF64, 16})`

$AB$, A row vector $A$ multiplied by a column vector $B$ (for spinor and gamma matrices)：
`function *(A::SVector{4, ComplexF64}, B::SVector{4, ComplexF64})`

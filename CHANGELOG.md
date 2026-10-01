# Changelog

## [0.5.2] - 2026-09-22

Fixed several known issues to improve version stability and overall runtime performance.

Improved the documents.

Reworked the `qBSE.particles!` function to support a file-driven particle schema.

- The particle type is no longer hard-coded as `structParticle`. Field names are now taken from the header line of the particle data file, and field types are inferred from the first data row (`Int`, `Float64`, or `String`).
- A concrete `struct` is generated at runtime via `Core.eval` and stored as the value type of `Dict{String, T}`, preserving type-stable field access and performance equivalent to the previous hard-coded type.
- The signature changed from `particles!(particles::Dict{String,structParticle}, filename::String)` to `particles!(filename::String)`. The result is injected into the module-level `const p` instead of being passed in by the caller.
- `Base.invokelatest` is used when constructing the `Dict{String, T}` and `T` instances, since the type is created at runtime and the calling function's world age predates it.
- Adding, removing, or renaming columns now only requires editing the data file; no source-code changes are needed.

## [0.5.1] - 2026-07-6

Fixed several known issues to improve version stability and overall runtime performance.

## [0.5.0] - 2026-06-27

### Improved 
- **Performance Optimization (`GEN` Module - Custom RAMBO Algorithm)**: Rewrote the native phase space generator using a highly optimized, custom **RAMBO** (Random Multiple Body Objects) flat-sampling algorithm. Achieved near-zero heap allocation via thread-local scratchpad buffering and loop-unrolled scalar expansion, entirely bypassing GC bottlenecks.
  - *Kinematic Pure-Generation Benchmark*: Clocked an elite throughput of **$10^7$ 4-body final-state momentum events within 3.0 seconds on a single CPU core @ 3.0 GHz** (~3.33 MHz pure event-generation rate).

### Added

-  New feathues added in  `qBSE` Module
   - **Momentum cutoff (finite cutoff) qBSE configurations**:  Added support for `cutoff_type = :cut` with both PW (`potential = :PW`) and unPW (`potential = :unPW`) potentials for standard qBSE.
   - **Chiral Unitary Approach (ChUA) implementations**:
      - ChUA using qBSE propagator (`ChUA = :qBSE`) with cutoff regularization (dimensional regularization planned for future).
      - Standard ChUA with `oset1405` prescription, covering both PW and unPW potentials.

   > **Note**: The existing infinite cutoff qBSE (`cutoff_type = :infty`, `potential = :unPW`, with exponential form factor and multiple `cutoff_re_type` options) remains unchanged and fully supported.

### Changed

- Improved the treatment of singularities on the real axis. The `eps` parameter in functions `propagator`, `resc0`, and `resc` is no longer required and has been removed.

## [0.4.4] - 2026-06-10

### Changed

- The `cs1` and `cs2` arrays returned by `Xs.Xsection` now represent the **differential cross sections**  
  \(d\sigma/dm\) (1‑D invariant mass) and \(d^2\sigma/dm_1 dm_2\) (2‑D Dalitz plot), respectively.  
  Previously these outputs contained the total cross section accumulated in each bin without the bin‑width division; they are now properly normalised by the bin widths, matching the standard definition of differential distributions.

### Added

- New keyword argument `symmetrize` in `Xs.Xsection` (default: `true`).  
  - When `true`, a **symmetrised Dalitz plot** is produced: all physically allowed pairings of the particles defined by the axes are filled for each event, with the event weight divided by the number of pairings. This mimics standard experimental procedures for final states with identical particles and yields label‑independent distributions with improved statistical precision.  
  - When `false`, a **fixed‑label Dalitz plot** is produced, using only a single, predefined pairing per event.  
  - The allowed pairings are determined with the following priorities: completely non‑overlapping sets, sets sharing exactly one particle, non‑subset sets, and (as a last resort) all possible pairs with a warning.  
  - The symmetrisation rules apply to both the 2‑D Dalitz plot (first two axes) and the 1‑D invariant mass spectra (all axes).


## [0.4.3] - 2026-03-28

### Added

- Added a new default option to generate events using a fixed seed for each processor when `fixed=true` in the `Xs.Xsection` function. When `fixed=false`, true random generation is used.

- Removed rotation in qBSE and modified the relevant code in `setTGA` and `TGA` functions accordingly.

- Added a new parameter form in fitting, with the transformation provided in the `extract_parameters` function of the `AUxs` package.

- The invaitant mass can be for more than two particle. For example `"p1:p2:p2"`.

## [0.4.1] - 2026-03-05

### Changed

- **Performance Optimization**: Refactored `structDimension` from a struct to a matrix to improve runtime efficiency.
- **Broadcasting Enhancement**: Refactored the `resc` function and applied padding to the `Dimt` and `TGt` variables to enable efficient broadcasting. These variables can now be broadcast quickly using methods such as `@everywhere const para_g = $para`.
- **Output Formatting**: Enhanced the output format for improved readability and visual appeal.

### Removed

- **Code Cleanup**: Removed the unused fields `Dimn` and `Dimo` from `structIndependentHelicity` to streamline the codebase.

## [0.4.0] - 2026-02-28

### Changed

-  Conducted a comprehensive review and unification of symbols related to channels and vertices to enhance code standardization and consistency. Related descriptions have also been refined.

- The labels for partlces and channels were changed from `::Symobl` back to `::String`.

## [0.3.9] - 2026-02-12

### Changed
- The labels for partlces and channels were changed from `String` to `Symbol`.
- The `pkey` was deleted. Now the paticles can be cited with the symbol directly. Such as `:pi` for pion.

## [0.3.8] - 2026-01-04

### Added
- Added an auxiliary function `build_full_idx` in the `qBSE` module. This is an internal function used in function `TGA`.
- Added function `resc` to enable parallel computation using `res0` (formerly `res`) as input.

### Changed
- Refactored `setTGA` to return `resc=(i,j)` as output, removing the redundant `sij` input parameter.
- Enhanced the `TGA` function interface with improved parameter handling and output structure.
- Renamed physics process identifier from `ch` to `proc` to avoid confusion with channels.
- Relocated `fV` function implementation to `qBSE.jl` for better module organization.

## [0.3.7] - 2026-01-01

### Added
- Added `anti` field to `structParticle` for automatic antiparticle labeling in `fV` function
- Added an optional regularization parameter `eps` to function `res` to add a small imaginary part to propagator denominators to avoid singularities.

### Changed
- Updated explicit format specification for `particle.txt`
- Changed `amps` output to return a `Vector` for separate contributions.
- Extended `cinter` parameter in `TGA(cfinal, cinter, Vert, para)` to support `NTuple` for multiple intermediate channels with weights.
- Improved function `TGA`.

## [0.3.6] - 2025-08-09

### Fixed
- Various bug fixes and stability improvements
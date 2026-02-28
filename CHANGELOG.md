# Changelog

## [0.3.9] - 2026-02-12

### Changed

- **Particle Identification Optimization**: Removed the original pkey mechanism. Particles can now be referenced directly by their names. For example, use qBSE.p("pi") to directly access the π particle. This change makes the code more intuitive and user-friendly.

-  **Symbol System Unification**: Conducted a comprehensive review and unification of symbols related to channels and vertices to enhance code standardization and consistency. Related descriptions have also been refined.

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
# Changelog

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
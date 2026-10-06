# `MetabolicGraph` changelog

## v.next

* AI: Exposed view functions on the Metabolic Graph view (state/model/entity introspection, reaction bounds, FBA, flux sampling, time-course, map navigation and highlighting, analysis save/load) so the AI assistant can read and drive it.

## 1.1.0 (2025-10-20)

* New E. coli core metabolism map and model added.
* New GLPK solver based on WebAssembly for FBA / Extreme points calculations.
* Improved sampling dialog with more options.
* Increased correlation limit for redundent edge removal
* Performance improvements and bug fixes.

## 1.0.1 (2025-10-03)

* New sampler and GLPK optimizer for FBA calculations.

## 1.0.0 (2025-09-06)

* Initial release of `MetabolicGraph` package.
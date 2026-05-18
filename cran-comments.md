## R CMD check results

0 errors | 0 warnings | 1 note

* One NOTE: "unable to verify current time" (system clock check, benign).

## Update (v0.1.1 -> v0.2.0)

This release adds an implementation of the dist.structure protocol
(now on CRAN at 0.5.0) for series systems. The four core topology
generics (phi, min_paths, min_cuts, system_signature) are specialized
for the trivial series topology; the rest of dist.structure's surface
(reliability, dual, structural importance, criticality, Vesely-Fussell,
substitute_component, etc.) flows automatically via the protocol
defaults.

No public API changes. The class chain now includes "dist_structure"
inserted after "dfr_dist", so existing methods on dfr_dist_series
continue to dispatch identically. The change is purely additive.

## Coordinated submission

This is part of a coordinated 3-package submission. All packages are
maintained by me. Updated versions being submitted simultaneously:

- flexhaz 0.5.2
- serieshaz 0.2.0 (this package)
- maskedcauses 0.10.0

## Test environments

* local Ubuntu 24.04, R 4.3.3
* win-builder (R-devel and R-release)

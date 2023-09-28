# Helm changelog

## 2.2.0 (WIP)

### Features

* Show tooltips on monomers rendered with WebEditor.

### Bug fixes

* Fix Helm cell renderer, allowing to render with gaps by skipping.
* Fix message WebEditor not supporting Helm with gaps.
* Prevent display tooltip by fallback cell renderer while WebEditor renders.
* Reset monomer placer on monomer lib changed.
* Fix Pistoia Helm WebEditor settings for case-sensitive monomer libs.

## 2.1.16 (2023-07-21)

This release focuses on stability.

*Dependency: datagarok-api >= 1.10.2*

### Bug Fixes

* Fixed the **To Atomic Level** function for helm sequences.
* `Neither peptide, nor nucleotide` error when calling **Helm Web Editor**
* Macromolecule detector for alphabet UN.

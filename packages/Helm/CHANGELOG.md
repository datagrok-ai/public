# Helm changelog

## 2.1.29 (2024-02-21)

* Fix Unkown monomers breaking grid rendering.

## 2.1.28 (2024-02-19)

* Fix TS config file.

## 2.1.27 (2024-02-15)

### Bug fixes

* Fix dojo patch dependency on package init
* Fix getMolfiles without current table view
* Fix initHelm/rewriteLibraries without current table view

## 2.1.26 (2024-02-14)

### Bug fixes

* Fix Properties panel of Context panel for sequences of MSA with gaps

## 2.1.25 (2024-02-09)

### Features

* Show tooltips on monomers rendered with WebEditor.

### Bug fixes

* Fix Helm cell renderer, allowing to render with gaps by skipping.
* Fix message WebEditor not supporting Helm with gaps.
* Prevent display tooltip by fallback cell renderer while WebEditor renders.
* Reset monomer placer on monomer lib changed.
* Fix Pistoia Helm WebEditor settings for case-sensitive monomer libs.
* Fix cell renderer monomer tooltip for polymer type
* Patch WebEditor hanging on getTextWidth

## 2.1.24 (2024-01-09)

### Bug Fixes

* GROK-13798: Fix error bioSubstructureFilter for Helm on size changed

## 2.1.16 (2023-07-21)

This release focuses on stability.

*Dependency: datagarok-api >= 1.10.2*

### Bug Fixes

* Fixed the **To Atomic Level** function for helm sequences.
* `Neither peptide, nor nucleotide` error when calling **Helm Web Editor**
* Macromolecule detector for alphabet UN.

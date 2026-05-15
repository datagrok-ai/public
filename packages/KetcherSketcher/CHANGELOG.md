# Ketcher Sketcher changelog

## 2.4.5 (2026-05-13)

* Workaround to fix multiple setMolecule calls, caused by bug in js-api (multiple subsribtions to copy-paste into sketcher string input field)

## 2.4.4 (2026-05-13)

* Hid Extended Table generic-atom buttons that can't be translated to RDKit-compatible V2000 query features from the substructure filter sketcher only — full set remains available in other sketcher contexts.

## 2.4.3 (2026-05-12)

* Deferred Ketcher's React Editor mount until ketcher-host has non-zero dimensions (ResizeObserver-gated) so RulerArea can resolve the canvas SVG's relative-unit sizes on its first render — adopted from PR #3789 (CLAUDE-159)
* Fixed structure filter not clearing when the filter-panel "clear" is hit while the sketcher is detached: `_setNotation` now zeroes the cached notations and resets molfiles to whitespace molblocks when called with an empty value on a detached sketcher

## 2.4.2 (2026-04-23)

* Reset explicitMol on change only if molecule has been already set

## 2.4.1 (2026-03-18)

* Updated ketcher libraries up to 3.12.0

## 2.4.0 (2026-03-16)

* [#3459](https://github.com/datagrok-ai/public/issues/3459): Ketcher: fix cropping sketcher on last grid column
* Ketcher: OG Smiles handling

## 2.3.0 (2025-03-29)

* Standardized plugin name to title case in README and changelog

## 2.2.3 (2025-01-14)

* Some additional styles fixes

## 2.2.2 (2025-01-14)

* Fixed styles of selectors

## 2.2.1 (2024-11-05)

* Fixed styles

## 2.1.10 (2024-07-26)

* Updated ketcher libraries up to 2.21.0 and Datagrok api to 1.20.0

## 2.1.9 (2024-05-30)

### Bug fixes

* Updated path to package icon

## 2.1.8 (2024-05-28)

### Bug fixes

* Fixed setting smarts into sketcher (had to remove saving user settings)

## 2.1.7 (2024-02-21)

### Features

* Visible "Apply" button when opening the setting

## 2.1.6 (2024-02-19)

### Features

* Saving user defined settings

## 2.1.5 (2024-01-31)

### Features

* Updated ketcher libraries up to 2.15.0

## 2.1.4 (2023-08-07)

### Features

* Adds [Ketcher](https://lifescience.opensource.epam.com/ketcher/index.html) as an optional molecular sketcher to Datagrok platfom

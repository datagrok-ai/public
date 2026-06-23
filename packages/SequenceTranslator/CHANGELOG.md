# Sequence Translator changelog

## 1.11.3 (2026-06-23)

* Markush enumeration: R-group templates picker backed by a built-in catalogue of common substituents (alkyl, aryl, heteroaryl, halogens, amines, protecting groups); copy R-groups between positions, append/replace with de-duplication, and export R-groups to CSV.
* Markush enumeration: added "Remove duplicates" option (dedupe results by canonical SMILES) and a custom result table name.
* Added `MarkushDefaults` package setting and settings editor — admin-configured defaults (cores, R-groups, mode and output options) seed the Markush Enumerator, distributed per user group.
* HELM Web Editor migration: PolyTool HELM enumeration now operates directly on HELM strings, dropping the JSDraw2 globals while preserving input formatting.

## 1.10.26 (2026-06-10)

* OligoNucleotide renderer: standalone backbone linkers (`p`, `[sp]`, or any monomer whose natural analog is `p`) — whether a 5'/3' cap, mid-strand, or a consecutive run — now render as linkage arcs (no chip), instead of being mistaken for conjugate pills. Recognized from the monomer library.
* OligoNucleotide renderer: alignment now accounts for standalone linkers. A run of linkers on one strand that the other lacks opens a gap (bases stay aligned / paired), and the partner strand draws a single wider arc across it.
* Demo data (`files/samples/sirna-demo.csv`): added a phosphate-capped single strand with mid-strand linkers and a duplex with a sense-only linker bulge.

## 1.10.25 (2026-06-09)

* OligoNucleotide renderer: strand alignment is now driven by the HELM when present — `strandtype` annotations decide the sense/antisense roles (chains are swapped accordingly) and `$connections$` `pair` entries fix the base-pair register.
* OligoNucleotide renderer: when no explicit info is present, sense/antisense are auto-aligned by sliding the antisense to the column shift with the most complementary base pairs (natural-analog aware), with terminal overhangs on either end rendered true.
* Oligo context panel: shows the resolved duplex register (blunt vs overhangs, and whether it came from explicit HELM pairs or auto-alignment).
* Oligo Structures panel: order the Sense / Antisense panes by the parsed strand roles (honoring `strandtype` swaps) rather than raw chain order, so a flipped HELM (sense = RNA2) no longer mislabels the strands.
* Demo data (`files/samples/sirna-demo.csv`): added rows with 3' overhangs and rows carrying explicit HELM pair / strandtype info.

## 1.10.24 (2026-05-19)

* Removed redundant demos

## 1.10.23 (2026-05-15)

* OligoNucleotide renderer: adaptive chip sizing (no hard max cap), 10px vertical breathing room, base-pair indicators (G-C 3 lines, A-U/A-T 2 lines, mismatches dashed), size-proportional chip corners.
* OligoNucleotide cell: double-click opens a full-screen canvas viewer with hover tooltips; HELM editing moved to a separate "Edit HELM" action alongside "Copy as HELM" and "Copy as Image".
* WK lines between matching bases.

## 1.10.22 (2026-05-15)

* OligoNucleotide renderer: strand-oriented redesign — sugar stripes on the outside edge of each strand, linkages as rounded apex arches; all colors sourced from the Bio monomer library.
* OligoNucleotide renderer: fixed reversed-antisense linkage placement (was shifted by one gap and dropping the leftmost-data linkage).
* Speedup canonicalization of Markush enumeration results.

## 1.10.21 (2026-05-14)

* Markush enumeration UI redesign, history and bug fixes

## 1.10.18 (2026-05-11)

* OligoNucleotide cell renderer: Fixed bracketed monomer base parsing (e.g. `[5Br-dC]` inside `(...)`) so renderer/resolver pick up modification colors and structures correctly.
* OligoNucleotide: Added double-click HELM Web Editor for OligoNucleotide cells; OK saves the edited HELM back to the cell

## 1.10.14 (2026-04-21)

* Chem Enumeration
* Helm enumeration improvements, paralel mode
* Polytool: support for custom notation back-conversion, fix for reversed rule order

## 1.10.5 (2025-12-11)

* PT-Synthetic: Fix missing molfile handler errors

## 1.10.3 (2025-11-03)

* Fix synthetic rule application for multiple instances when r group numbering is reversed

## 1.10.2 (2025-11-03)

* Update Bio Lib API
* Fix synthetic rule errors for missing monomers
* Fix incosistent R group conneections for reversed rule orders

## 1.9.5 (2025-05-02)

* Fix changing bond numbers when fixing molblocks

## 1.9.4 (2025-04-08)

* Linearization fix for polytool convert.

## 1.9.1 (2025-03-31)

* Combination of sequence sets dialog.

## 1.8.0 (2025-02-20)

* Support for enumeration in all notations
* fallback for polytool convert
* Corrected highlighting

## 1.6.5 (2024-12-01)

### New features

* PolyTool: highlight in mol based on given notation
* PolyTool: Rule manager enhancements

### Bug fixes

* PolyTool fix R3, R4 and elder groups capping

## 1.6.4 (2024-11-27)

* PolyTool: linearization of single sequence molecules
* Enumerator: Added whole library enumeration

## 1.6.3 (2024-11-22)

Moved Context menu items to detector

## 1.6.2 (2024-11-17)

### New features

* PolyTool rule manager enhancement
* Sequence renderer font resizing

## 1.6.1 (2024-11-15)

### Bug fixes

* PolyTool fix default options
* PolyTool fix explicit carbon cut

## 1.6.0 (2024-11-14)

### New features

* Add PolyTool work with sinthetic monomers

## 1.5.3 (2024-11-06)

### Bug fixes

* PolyTool fix rules empty specification (all monomers undergo)

## 1.5.2 (2024-11-06)

### Bug fixes

* PolyTool rules several instances bug fix

## 1.5.1 (2024-11-06)

### Bug fixes

* PolyTool rules fixes rule management
* PolyTool rules fixes crash if monomer is absent

## 1.5.0 (2024-11-04)

### New features

* Add PolyTool interface for rule management

## 1.4.10 (2024-11-01)

### Bug fixes

* Fix build

## 1.4.9 (2024-10-31)

### Bug fixes

* Fix package initializing helmHelper
* Fix Chain using HelmAtom.bio of type PtBio
* Fix PT Enumerate dialog for harmonized sequences
* Fix PT Enumerate for monomer hovering
* Fix PT Enumerate for historical values
* Fix PT Enumerate adding buildCyclizedMonomerHoverLink, WIP
* Fix adding CyclizedGridCellRenderBack

## 1.4.8 (2024-10-25)

### Bug fixes

* PolyTool ambigous R groups handling bug

## 1.4.7 (2024-10-25)

### Bug fixes

* PolyTool ambigous R groups handling bug

## 1.4.6 (2024-10-25)

### New features

* Add PolyTool Enumerate for harmonized sequences
* PolyTool ambigous R groups handling

### Bug fixes

* Fix PolyTool chain, add check consistency, add tests
* Fix PolyTool Convert tests
* Fix the package init

## 1.4.5 (2024-10-11)

### Bug fixes

* Fix for grok.userSettings
* Fix moving support custom notation harmonized sequence
* Fix PT Convert for monomer highlighting
* Fix PT Convert reaction monomer colors
* Fix PT detectors adding tests
* Bump dependencies versions

## 1.4.4 (2024-09-30)

### New features

* Add convert with reaction rules

### Bug fixes

* Fix sample reaction rule and monomer order
* Fix isolate PolyTool convert error on row
* Fix adding test for Chain .fromNotation
* Fix adding test toAtomicLevel getNewMonomer

## 1.4.3 (2024-09-27)

### Features

* Add PolyTool convert reverse

### Bug fixes

* Add sample reaction rule aG + azG = GGaz, test data
* Fix using MonomerLibHelper.loadMonomerLibForTests
* Add test toAtomicLevel with overridden monomer lib
* Fix PolyTool Enumerate for zero-based monomer position
* Fix Chain.getHelm square brackets for multi char symbols

## 1.4.2 (2024-09-24)

### New features

* Add PolyTool breadth enumerator

## 1.4.1 (2024-09-12)

### Bug fixes

* Input style fixes

## 1.4.0 (2024-09-10)

### Features

* Add PolyTool Enumerator Chem dialog to context menu for Molecule

### Bug fixes

* Fixes for PolyTool Enumerator Helm dialog
  * Trivial Name input Behavior
  * Enable for any Macromolecule cell
  * Example with monomers existing in default monomer lib
* Fixes for datagrok-api changes
* Add sample cyclized.csv for rules_examples.json
* Fix PolyTool Convert, add tests cyclized
* Fix PolyTool enumerate add progress indicator
* Fixes for api changes
* Fix PT dialogs destroy on build fail
* Fix rules_examples for NH2-D with R3 and syclized sample
* Fix polyToolConvert to work without table
* Add tests for polyToolConvert, ui

## 1.3.14 (2024-09-02)

Fix remove CyclizedNotationProvider

## 1.3.13 (2024-08-19)

### Features

* Improve PolyTool Enumerate Dialog
  * Add check and warning for substituting missed monomers
  * Add flag to keep the original molecule in the result
  * Add choice to select a molecule id and generate derivative

### Bug fixes

* Fix PolyTool Enumerate dialog hanging on monomer click
* Fix user settings for publish

## 1.3.12 (2024-08-16)

PolyTool Enumerator Single and Matrix, dialog, tests

### Bug fixes

* PolyTool Enumerator dialog grid input for placeholders, to atomic level option
* PolyTool Enumerator dialog mouse interactivity, to atomic level, dialog size
* PolyTool Enumerator dialog size and fit molecule input

## 1.3.11 (2024-07-31)

### Features

* Added chemical based enumerator

## 1.3.10 (2024-07-29)

### Features

* Add 'Create/Modified' date metadata to patterns

### Bug fixes

* Load new patterns upon saving / overwriting

## 1.3.9 (2024-07-23)

* Dependency: datgarok-api >= 1.20.0

## 1.3.8 (2024-07-09)

### Features

* Add PolyTool-enumerate context menu for macromolecule cell

### Bug fixes

* Use HelmInput from Helm
* Bump dependencies versions

## 1.3.7 (2024-07-01)

### Bug fixes

* Fix style for colored text input

## 1.3.6 (2024-06-26)

### Bug fixes

* Fix bulk input in Pattern app

## 1.3.5 (2024-06-13)

### Features

* PolyTool Enumeration: cell with helm is added to the input

## 1.3.4 (2024-06-12)

### Bug fixes

* Fix terminal modification SVG rendering

## 1.3.3 (2024-06-12)

### Bug fixes

* Fix terminal modification SVG rendering

## 1.3.2 (2024-06-10)

### Features

* Add package MonomersPath package settings fallback to monomers-sample
* Add README for monomers-sample
* Add sample for bulk translation in Translator app

### Bug fixes

* Fix tests for MonomersPath package settings

## 1.3.1 (2024-06-05)

### Bug fixes

* PolyTool: example data added to project

## 1.3.0 (2024-03-30)

### Features

* PolyTool: algorithm of conversion facilitated, json files are used as rules
* PolyTool: enumeration feature (demo) added

## 1.2.9 (2024-03-30)

### Features

* #2707: Add original and canonical to monomer

## 1.2.7 (2024-01-29)

### Features

* Bulk translation for formats

## 1.2.6 (2024-01-29)

### Bug fixes

* Fixed column input in Pattern
* Fixed whitespaces in molfiles

## 1.2.0 (2023-10-27)

### Features

* Tabs implemented as separate apps
* Clear buttons added to sequence inputs

## 1.1.4 (2023-09-07)

### Bug fixes

* Fixed molfile direction for AS and AS2 on SDF tab

## 1.1.3 (2023-09-07)

### Bug fixes

* Fixed Direction of AS and AS2 for SDF tab

## 1.1.2 (2023-08-16)

### Bug fixes

* Fixed HELM generation for modifications

## 1.1.1 (2023-08-10)

### Features

* Updated dependencies

## 1.1.0 (2023-07-26)

* Dependency: datgarok-api >= 1.15.2

### Features

* Option to add extra formats for translation in addition to `Axolabs`, `BioSpring` and `MerMade12`
* Generation of HELM strings
* Improved SMILES generation
* UI/UX adjustments

# Ketcher Sketcher changelog

## v.next

* Fixed query features drawn in Datagrok's substructure filter being lost: aromatic or aliphatic, implicit H count, ring membership, ring size, connectivity and a custom atom query are written as RDKit SMARTSQ groups in a V3000 molblock, and a custom bond query MDL can say as its bond type and topology (ketcher-core's V2000 left them out, or wrote "any"); a molblock with SMARTSQ groups set on Ketcher shows them as custom queries, and one Ketcher wrote shows the drawing it was written from. Not held: an atom's chirality query, and a custom bond query MDL cannot say
* The substructure filter offers the AH, QH, XH, MH and * generics: RDKit reads them as Ketcher writes them
* Fixed every SMARTS and V3000 export of the page's Ketchers stopping for good once Indigo left one conversion unanswered (it does when another meets it): a conversion not answered in 15 seconds is given up, and the ones after it go on
* Fixed Indigo leaving requests unanswered: every Ketcher of the page shares one struct service, which sends Indigo one request at a time (Indigo's one worker answered only the Ketcher mounted last, and of two requests in flight it dropped one)
* Fixed a page that had opened Ketcher many times filling its memory until it stopped answering: a Ketcher closed or paused ends its subscription to Ketcher's settings, which kept every Ketcher ever opened alive
* Fixed Indigo stopping for good after some 2,000 conversions in a page: a conversion is no longer sent with the hidden macromolecules editor's whole monomer library, which Indigo's worker kept in its memory, about 1 MB each time, until it ran out

* Fixed a structure drawn from a template coming out twice: an export taken while the pointer was still over the canvas included the template's floating preview
* Fixed a structure drawn right before OK being lost: the molblock is written and announced in the change itself, and the Indigo exports run one at a time (overlapping conversions got each other's replies or none)
* The editor is ready without waiting for the hidden macromolecules editor to load
* The sketcher is announced ready (`init`, `sketcherReady`) only once Ketcher's editor takes input: a molecule set the moment it was announced was dropped (`setSmiles`, `setMolfile`, `setSmarts` tests)
* Fixed the first stroke right after the editor opened being lost: an empty starting molecule is no longer loaded into the empty canvas (the asynchronous load took the stroke's change for its own and then wiped it), and the editor is set up once, though Ketcher reports it ready twice

## 2.4.8

* Fixed multiple sketchers on one page breaking each other (upstream ketcher-core singletons): opening a new sketcher now suspends the others behind a "Reload" placeholder that remounts the editor with its molecule preserved

## 2.4.6 (2026-05-19)

* Fixed "couldnt find ketcher instance N" error when clicking OK on the cell editor dialog: the `change` handler now bails out if the sketcher is detached while a `getMolfile()` call is still in flight

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

# Exports, alignment, persistence and entry points

Translated 2026-09-13, validated 2026-09-15 (see the README for the run record). The original
Playwright specs remain in `../../playwright/`.

| Original spec | Feature | Assertions retained or strengthened |
|---|---|---|
| `export-invariant-map.test.ts` | `sar/export.feature` | Both SAR viewers use the context menu. The result has the complete 22 × 17 count matrix, an AAR string column and integer position columns. Every cell is checked against the source sequences, including zeros. |
| `export-mutation-cliffs.test.ts` | `sar/export.feature` | Default export has the six canonical columns and all 6253 unique source row pairs. Sequence identity, activity values, Delta, sequence semantics and Mutation concatenation are checked. The extra-column scenario drives the picker and verifies both ID columns against the actual source rows. |
| `manual-alignment.test.ts` | `sar/manual-alignment.feature` | A real grid cell click opens the editor. Apply changes the sequence and every split monomer, with explicit adjacent/end-position assertions. Reset preserves the applied sequence and discards an unsaved edit. A real WebLogo header click checks selection after editing. |
| `sar-save-reopen.test.ts` | `sar/project-round-trip.feature` | Project save/reopen retains data, semantic types, settings, all four SAR/cluster viewers, 17 positions and the tall WebLogo headers. A new monomer selection must select precisely its source rows. The library project step cleans up by IDs. |
| `peptide-sar-demo-dashboard.test.ts` | `entry/demo-dashboard.feature`, `entry/landing.feature` | Dashboard registration metadata, dataset/tags, scaling, threshold, all viewers and rendered content. Landing buttons open all three fixture datasets, with exact row counts and notation checks. Side panel visibility is asserted. |

The fixture was recounted from `data/demo/bio/peptides.csv`: 647 rows, 17 positions, 22 non-gap
monomers, 6253 Hamming-distance-one row pairs and 647 distinct IC50 values. IC50 therefore identifies
the source row in the export oracle. Representative counts are NH2@1=647, COOH@17=647, A@2=299,
M@2=9 and Y@10=604. The export oracle uses the original sequence strings, not viewer statistics or
the mutation-cliff implementation.

`manual-alignment.ts` named the textarea Sequence and corrected the write-back index: split
position zero belongs to column "1". The previous code wrote it to nonexistent column "0", shifted
the following monomers and left the last position stale. The original test only required any
position to change, hiding this error. Live tests also found stale monomer statistics after Apply
and an empty context panel caused by the old clear/100 ms restore. Apply now refreshes the
statistics and the SAR viewers and restores the edited cell as the current object; it does not
re-run MCL or Sequence Space, as the original never did.

Claims deliberately omitted: arbitrary canvas/child counts, swallowed errors, private model
mutation calls, the impossible last-error string checks, environment checks for other packages'
Browse groups, and the always-true restored-selection range. The save step exercises the public
project API; it does not test the ribbon's Save dialog. The new selection assertion after reopen
does not promise that selection itself is serialized. Demo-gallery card navigation remains out
of scope; dashboard metadata and its registered function are exercised.

Live findings:

- Export picker: ID is row 2 of `__name` (Activity is row 1). The real picker interaction passes.
- Mutation exports previously rounded Float64 activities into Float32; the full pair oracle
  exposed this. The export now preserves Float64 activities and Delta.
- The initial grid now starts with Activity and the first monomer positions visible.
- Project reopen restores all four viewers, actual logarithmic activities and a working header
  glyph selection. The model re-shows its accordion when the platform makes the table's row
  group current after a selection; it does so from a microtask, because the event bus is
  synchronous and a current object set from inside its own dispatch throws.
- The cliff-chart oracle reads the position column with `Array.from`: a missing monomer is a
  hole in `toList()`, which `map()` skips.
- Landing view and help/toolbox/context visibility checks pass. The toolbox locator excludes the
  empty sliding-shell host by matching the actual toolbox's caption attribute.

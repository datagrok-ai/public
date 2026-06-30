---
phase: 15
plan: 06
status: complete
type: implementation
---

# Plan 15-06 Summary — Reviewer-Side Audit Panel

Built `src/panels/published-analysis-panel.ts` exporting `publishedAnalysisPanel(proteinId): DG.Widget`. Plan 07 will wire the `@grok.decorators.panel` registration with `semType='Proteomics-ProteinId'` in `src/package.ts`.

## Behavior

- First-line `isPublished(df)` guard returns an empty widget on non-published DataFrames. The analyst's live working DataFrame must NEVER see this panel content.
- `findHostDataFrameForProtein(proteinId)` returns null → empty widget (no protein-id context).
- `getPublishedMetadata(df)` returns null → friendly "metadata unavailable, ask sharer to re-share" fallback (catches the edge case where both belt-and-braces sources fail — extremely defensive).
- `meta.publishedAt` parse fails → renders `'unknown date'` rather than crash (defensive normalization in `dateSlice` helper).
- All paths return a `DG.Widget`; the panel never throws.

## Read order — column FIRST, tag SECOND (Pitfall 3)

`getPublishedMetadata` (Plan 01) already enforces column-first / tag-second for all 12 metadata fields. The panel uses the helper rather than reading tags directly — single read-order policy.

For `publishedByEmail`, the panel's `readSharerEmail` helper also tries column-first then tag, with the additional `meta.publishedByEmail` short-circuit at the top (since the helper already populated it via the same Pitfall-3 sequence).

## Audit fields rendered (PUB-07)

| # | Label | Value source | Notes |
|---|-------|--------------|-------|
| 1 | Target | `meta.target` | raw user input, not slug |
| 2 | Shared | `dateSlice(meta.publishedAt)` | YYYY-MM-DD, falls back to `'unknown date'` |
| 3 | Shared by | `meta.publishedBy` | friendly name |
| 4 | Method | `methodLabel(meta.deMethod)` | limma → 'limma — moderated t-test', etc. — biologist expansion per CONTEXT.md domain pin |
| 5 | Comparison | `groups.group2.name + ' vs ' + groups.group1.name` | falls back to 'unavailable' if `proteomics.groups` is absent |
| 6 | Fold-change cutoff | `log2(fcThreshold) (2^fcThreshold-fold)` | shows both log2 and fold so biologist doesn't have to convert |
| 7 | p-value cutoff | `{pThreshold} (FDR-adjusted)` | |

## Supersede link (D-04, P1)

When `meta.supersededBy` is non-null, renders a yellow banner at the TOP of the body with `ui.link('Newer version available — open it', async () => projects.find(id).open())`. Tries-catches the find/open and surfaces a `grok.shell.warning` if the reviewer cannot open the newer version (e.g., it was deleted, or they don't have permissions).

`getPublishedMetadata` reads `meta.supersededBy` column-FIRST then tag, so W-8 dual-write from Plan 04 step 9 is automatically supported by both paths.

## Mailto button (PUB-13, P2)

Imports `buildMailtoUrl` from `../publishing/publish-state` per **B-1 wave-cycle fix** — NOT from `../publishing/share-dialog` which is Plan 05 (wave 4). Wave 2 → wave 4 import would break the wave graph.

Renders only when `readSharerEmail` returns non-null. Uses `df.name` as the project name (Plan 02 set it to `<source>_published_<YYYY-MM-DD>`). Renders via `ui.link('Request re-run with different parameters', mailtoUrl)` returning an `HTMLAnchorElement` with `href` already set by the helper.

## Jargon audit (Pitfall 14)

Source-visible strings reviewed: 'Shared analysis details', 'Newer version available — open it', 'Target', 'Shared', 'Shared by', 'Method', 'Comparison', 'Fold-change cutoff', 'p-value cutoff', 'limma — moderated t-test', 'DEqMS — peptide-count-aware moderated t-test', 'Welch t-test', 'Spectronaut Candidates (pre-computed DE)', '(unknown)', 'unknown date', '(FDR-adjusted)', 'Request re-run with different parameters', 'Shared analysis metadata unavailable. Ask the sharer to re-share.', 'Could not open the newer version: …'.

Banned terms (DataFrame, tag, semType, ACL, viewer factory) — zero occurrences in user-facing strings.

## B-1 verification

```sh
grep "from '../publishing/share-dialog'" src/panels/published-analysis-panel.ts
```
returns zero matches — wave-2 plan does not import from wave-4 plan. `buildMailtoUrl` imported from `../publishing/publish-state` as required.

## B-3 verification

Plan-level priority not marked P2 (this plan covers PUB-06 + PUB-07 = P1 + PUB-13 = P2). Only Task 2's mailto sub-feature is P2; supersede link in Task 2 is P1 and ships now.

## Verification

Project-wide `tsc --noEmit` passes clean. No imports from share-dialog.ts. All paths return a `DG.Widget`.

## Output

`src/panels/published-analysis-panel.ts` — 158 lines, type-checks under strict mode.

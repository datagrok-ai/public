# MSA — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open sample_FASTA.csv | 8s | PASS | PASSED | `System:AppData/Bio/samples/FASTA.csv`; 64 rows, Sequence semType=Macromolecule |
| 2 | Add new column Clusters = RandBetween(0, 5) | 10s | PASS | PASSED | Name input `.ui-input-addnewcolumn-name`, formula in CodeMirror `.cm-content`, OK button is `[name="button-Add-New-Column---OK"]` (not generic `button-OK`). Clusters int column added, stats min=0 max=4 |
| 3 | Bio > Analyze > MSA opens dialog-MSA | 3s | PASS | PASSED | Inputs: Sequence, Clusters, Method, Gap open, Gap extend, Terminal gap, Selected Rows Only. ALIGNMENT PARAMETERS button present. Clusters input auto-bound to the new `Clusters` column |
| 4 | Set Cluster to the new Clusters column | 1s | PASS | PASSED | Column input auto-bound to `Clusters` already (most-recently-added int column wins). Explicitly set via `dlg.inputs.find(i => i.caption === 'Clusters').value = df.col('Clusters')` to make deterministic |
| 5 | Alignment parameters button adds input parameters | 2s | PASS | PASSED | Before click: Gap open/extend/Terminal-gap height 0. After click: all three have height ≥ 58 px. Method input stays `display:none` regardless — see retrospective |
| 6 | Check the new MSA column | 3s | PASS | PASSED | `msa(Sequence)` added, semType=Macromolecule, grid `cellType=sequence`, all 64 sequences aligned to exactly 49 chars, and length is identical within every one of the 5 clusters. Waited for grid cellType to flip from `string` → `sequence` (renderer attach is async after column creation) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~2m 20s |
| grok-browser execution (scenario steps) | ~30s |
| Execute via grok-browser (total) | ~2m 50s |
| Spec file generation | ~40s |
| Spec script execution | 19s |
| **Total scenario run (with model)** | ~10m |

## Summary

The MSA scenario works end-to-end. In the live browser the new "Clusters" column is added via Edit > Add New Column with `RandBetween(0, 5)`, Bio > Analyze > MSA opens a dialog that auto-binds the Clusters column input, the ALIGNMENT PARAMETERS button reveals the Gap open/Gap extend/Terminal gap inputs, and OK produces an `msa(Sequence)` Macromolecule column with the sequence renderer and equal-length sequences within every cluster (all 5 clusters aligned to 49 chars). First spec run failed (3/6 steps) because the Add New Column OK button is `[name="button-Add-New-Column---OK"]`, not the generic `[name="button-OK"]` scoped to the dialog root. After fixing that plus two smaller issues (wait for grid `cellType === 'sequence'` before asserting renderer; JS-API fallback when CodeMirror drops keystrokes) the replay passes all 6/6 in 19s. Total scenario run (with model, including three login-helper flake retries): ~10m.

## Retrospective

### What worked well
- `grok.dapi.files.readCsv('System:AppData/Bio/samples/FASTA.csv')` + semType subscription + 5s Bio settle opens the dataset consistently.
- CodeMirror 6 formula editor in the Add New Column dialog accepts `type_text` after a `click` on `.cm-content` — but occasionally drops keystrokes on a cold browser context, so a post-type `textContent` check + JS fallback to `df.columns.addNewCalculated` is a worthwhile safety net.
- Bio > Analyze menu pattern (`click [name="div-Bio"]` → `mouseenter [name="div-Bio---Analyze"]` → `click [name="div-Bio---Analyze---MSA..."]`) continues to work and is unchanged from the Analyze scenario.
- The MSA dialog sensibly auto-binds the most recently added `int` column to the Clusters input, so the minimal UI flow (add Clusters column → MSA → OK) "just works" for the scenario.
- Per-cluster sequence-length sets are a clean and precise way to assert the scenario's "sequences within a cluster are of the same length" success criterion.

### What did not work
- The Add New Column dialog names its OK/CANCEL buttons with the full path `button-Add-New-Column---OK` / `...---CANCEL`, whereas other dialogs (MSA included) use the generic `button-OK` / `button-CANCEL`. This inconsistency cost most of a debug cycle — the accessibility snapshot shows both as plain "OK" buttons, so there's no hint from the a11y tree that the `name=` attribute differs from the common case.
- `col.getTag('cell.renderer')` returned `null` in the Playwright page context even though `col.tags` iteration in the MCP browser session returned `cell.renderer=sequence`. The reliable signal is `grok.shell.tv.grid.col(name).cellType` — but that value is populated *after* the column appears, so the wait condition must include the renderer flip, not just the column existence.
- The MSA dialog's `Method` input (`input-host-Method`) stays `display:none` the entire time — the ALIGNMENT PARAMETERS button does not un-hide it. Whether this is intentional (only one method offered = "mafft --auto", no choice shown) or a regression isn't obvious from the UI.
- The shared `spec-login.ts` helper flaked three separate times during iteration with "locator resolved ... element was detached from the DOM" on the password input, consistently on a cold Playwright context at ~8-12 s into the login. Each retry succeeded on the next invocation with no spec change.

### Suggestions for the platform
- Normalise dialog button names: either every dialog uses `button-OK` / `button-CANCEL`, or every dialog uses `button-{Title}---OK`. Today only some dialogs use the full-path form and it is easy to pick the wrong convention from reference docs.
- Make `Column.getTag('cell.renderer')` and `Column.tags` return identical entries across the JS API boundary. Today they disagree — tag visible during MCP iteration, invisible via `getTag` in a fresh Playwright page context — which points to the renderer tag being set on a private `meta` slot that `getTag` doesn't consult.
- If the Method input in the MSA dialog is intentionally hidden (single option), drop it from the dialog entirely rather than leaving an invisible `input-host-Method` in the DOM; if it is supposed to be revealed by ALIGNMENT PARAMETERS, the button currently does not do so.
- Harden `spec-login.ts` against "password input detached during focus" — either wait for `networkidle` a second time between typing login and focusing password, or retry the whole `submitLogin` on the detach error rather than bailing out of the helper.

### Suggestions for the scenario
- Step 1 says "sample_FASTA.csv" but the actual file path is `System:AppData/Bio/samples/FASTA.csv`. Update to the canonical System path (same issue previously flagged in `analyze-run.md`).
- Step 2 says "a formula RandBetween(0,5)" — clarify that the formula syntax is `RandBetween(0, 5)` with a space, the resulting column is int-typed, and pin the column name ("Clusters") so step 4 can reference it unambiguously.
- Step 5 "Check that **Alignment parameters** button adds input parameters to the dialog properly" is vague. Enumerate the expected inputs (Gap open, Gap extend, Terminal gap) so the test has a deterministic assertion, and clarify the expected behaviour of the Method input (hidden? shown?).
- Step 6 success criteria should explicitly mention the column name pattern (`msa(<SequenceColName>)`), the `cell.renderer=sequence` / `cellType=sequence` invariant, and the per-cluster equal-length invariant as the three things to verify.

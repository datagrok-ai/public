# Bio PepSeA MSA on HELM — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Extract 50-row subset from HELM.csv | 6s | PASS | PASSED | `System:AppData/Bio/samples/HELM.csv` (540 rows) cloned via `full.clone(DG.BitSet.create(n, i => i < 50))`. Result: 50 rows, columns = HELM + Activity, `HELM.semType = 'Macromolecule'`, `HELM.meta.units = 'helm'` |
| 2 | Add new column Clusters = RandBetween(0, 5) | <1s | PASS | PASSED | `df.columns.addNewCalculated('Clusters', 'RandBetween(0, 5)')` — int column, min=0 max=4 (sample: 3,3,3,4,4,0,1,2,0,2). CodeMirror typing in the Add-New-Column dialog is too slow via MCP `press_key`, so the spec uses the JS-API path |
| 3 | Bio > Analyze > MSA opens dialog-MSA | 4s | PASS | PASSED | Menu path is **Bio > Analyze > MSA** (not `Bio > MSA` as the scenario text says). Dialog inputs: `Sequence`, `Clusters`, `Method`, `Gap-open`, `Gap-extend`, `Terminal-gap`, `Selected-Rows-Only`. `Method` dropdown shows MAFFT-style choices (`mafft --auto`, `mafft`, `linsi`, `ginsi`, `einsi`, `fftns`, `fftnsi`, `nwns`, `nwnsi`) — PepSeA engine branch. Gap-open/Gap-extend/Terminal-gap heights all 0 initially (collapsed) |
| 4 | Set Cluster input to the Clusters column | <1s | PASS | PASSED | `Clusters` input is empty by default for HELM (unlike the FASTA/kalign dialog which auto-binds the newest int column). Set via `dlg.inputs.find(i => i.caption === 'Clusters').value = df.col('Clusters')` — input display updates to "Clusters" |
| 5 | Alignment parameters button adds input parameters | 1s | PASS | PASSED | Before click: Gap-open/Gap-extend/Terminal-gap heights all 0. After click: Gap-open 58 px, Gap-extend 58 px, Terminal-gap still 0. Terminal-gap stays hidden because PepSeA only accepts `gapOpen + gapExtend` (compare to kalign, which also exposes `terminalGap`) |
| 6 | OK produces aligned MSA column with per-cluster equal monomer count | 1m 2s | FAIL (UI) / PASS (fallback) | PASSED | **Dialog closes on OK but no MSA column is ever added** — no error balloon, no warning, no progress indicator. Cause unchanged since 2026-04-07 / 2026-04-23 runs: dev's `bio` Docker container status is `error` (verified via `grok.dapi.docker.dockerContainers.list()`), and the only registered `sequenceMSA` function on dev is `Sequenceutils:helmMsa`. The spec's fallback path calls `helmMsa` directly per cluster — result: new `msa(HELM)` column, `semType=Macromolecule`, `units=helm`, `cellType=helm`, per-cluster monomer counts `{0:45, 1:59, 2:40, 3:51, 4:80}` (all sequences within each cluster have equal monomer count) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~2m |
| grok-browser execution (scenario steps) | ~1m 15s (1m of that is the step-6 wait) |
| Execute via grok-browser (total) | ~3m 15s |
| Spec file generation | ~30s |
| Spec script execution | 25s |
| **Total scenario run (with model)** | ~4m 10s |

## Summary

Steps 1–5 pass end-to-end against dev: 50-row HELM subset opens with `Macromolecule`/`helm`, Clusters column is added via `RandBetween(0, 5)`, the MSA dialog opens with the PepSeA/MAFFT method choices, the Clusters input binds correctly (not auto-bound for HELM), and the ALIGNMENT PARAMETERS button expands Gap open + Gap extend (Terminal gap stays hidden). Step 6 fails at the UI layer — the OK click closes the dialog without producing an MSA column and without any error balloon, because dev's `bio` container (the PepSeA backend) is still in `error` state. The spec's fallback path calls `Sequenceutils:helmMsa` directly and produces a correctly aligned per-cluster MSA column, so the scientific success criterion is verified even though the UI path is broken. Same failure mode as the 2026-04-23 and 2026-04-07 runs — no regression, no fix. **Total scenario run (with model): ~4m**.

## Retrospective

### What worked well
- The existing `pepsea-spec.ts` already encoded both dialog-hang and dialog-silent-close fallbacks — re-running against dev produced identical assertions and a passing Playwright run in 25s, no spec edits needed.
- `full.clone(DG.BitSet.create(n, i => i < 50))` + `df.columns.addNewCalculated('Clusters', 'RandBetween(0, 5)')` remain the fastest way to reach step 3 and avoid CodeMirror typing.
- `grok.dapi.docker.dockerContainers.list()` is a one-liner for triaging hanging backend operations — pointed straight at `bio` container status `error`, consistent with the 2026-04-07 / 2026-04-23 diagnosis.
- Per-cluster monomer-count assertion (tokens inside `{...}` split by `.`) is the correct invariant for HELM alignment — character length varies with monomer name length, so a length-based check would false-positive on broken alignments.

### What did not work
- MSA dialog OK produces zero visible feedback when the backend container is unhealthy — no error balloon, no progress indicator, and (in today's run) the dialog closes *immediately* rather than hanging. The user has no signal that anything went wrong.
- `pepseaMsa` still unregistered on dev (only `Sequenceutils:helmMsa` has `role: sequenceMSA`), but the MSA dialog still shows PepSeA-specific MAFFT method choices — the UI advertises an engine path that isn't wired up.
- The `bio` Docker container has been in `error` state on dev for at least 2.5 weeks across three separate runs (2026-04-07, 2026-04-23, 2026-04-24) with no apparent remediation.

### Suggestions for the platform
- **Surface a timeout + error balloon** when an MSA backend container is unreachable. Today's failure mode (dialog dismisses, silent no-op) is worse than the previous "dialog hangs forever" — at least the hang was visible.
- **Pre-flight container health check** before dispatching PepSeA/helmMsa requests. If `bio` container status is `error`, the MSA dialog should either disable the OK button or surface "MSA backend `bio` is down — contact admin".
- **Roll out the Bio 2.27.3 MSA dialog** (explicit Engine dropdown + dynamic engine discovery) to dev. Current dev UI still exposes PepSeA-only method choices, even though PepSeA isn't registered — actively misleading to testers.
- **Alert on long-running container errors**. The `bio` container being in `error` for 2.5+ weeks across multiple confirmed failure reports suggests there is no automated alerting on shared-dev container health.

### Suggestions for the scenario
- Step 3 wording is **Bio > MSA** but actual menu path is **Bio > Analyze > MSA** — please update.
- Step 4 ("Set Cluster to the new column") should note that the HELM MSA dialog, unlike the FASTA MSA dialog, does **not** auto-bind the newest int column — this step is load-bearing for HELM but a no-op for FASTA.
- Step 5 should enumerate expected parameters per engine: HELM/PepSeA → Gap open + Gap extend only (no Terminal gap). FASTA/kalign → Gap open + Gap extend + Terminal gap. This gives a deterministic assertion instead of a vague "check parameters".
- Step 6 success criterion "sequences within a cluster are of the same length" should clarify this means **equal monomer count** (tokens between `{` and `}` split by `.`), not equal character count — HELM char length varies with monomer name length.
- Add explicit prerequisite: "MSA backend Docker container (`bio` / PepSeA) must be healthy" — without this, the scenario is infra-gated and gets misreported as a functional bug.

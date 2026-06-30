# Phase 13 — Wave-0 Findings (A2 / A3 / A4)

De-risking output for 13-01 Task 2. Unblocks 13-03 (R3 rule), 13-04 (cache), 13-05 (D-09 site).
Evidence files: `~/Downloads/ck/2026-03-25_BP_EMT_DMD_Candidates.tsv`,
`~/Downloads/ck/CKomics_tool2.py`, `~/Downloads/ck/Spectronaut-Concordance-Note-2026-05-15.md`,
package source, `js-api/src/dapi.ts`. (Throwaway probes not committed — memory
feedback_no_commit_oneoff_scripts.)

---

## A2 — Spectronaut Candidates sign convention (for 13-03 / R3)

**Definitive: the Candidates `AVG Log2 Ratio` sign is internally consistent with its
declared `Comparison (group1/group2)`. There is NO hidden inversion in the Candidates file.**

Evidence — `2026-03-25_BP_EMT_DMD_Candidates.tsv` (the file underlying the client volcano):

- Schema (1-indexed cols): 2 `Comparison (group1/group2)`, 3 `Condition Numerator`,
  4 `Condition Denominator`, 7 `AVG Group Quantity Denominator`,
  8 `AVG Group Quantity Numerator`, 9 `AVG Log2 Ratio`, 13 `Pvalue`, 14 `Qvalue`, 20 `Genes`.
- Single comparison in the file: `DMD / WT` (Condition Numerator = DMD, Denominator = WT).
- The `DMD` gene row: `AVG Log2 Ratio = -3.9756`, `AVG Group Quantity Numerator (DMD) = 126.58`,
  `AVG Group Quantity Denominator (WT) = 1991.24`. `log2(126.58 / 1991.24) = -3.976` — matches exactly.
- ∴ **`AVG Log2 Ratio = log2(Numerator / Denominator) = log2(group1 / group2)`**;
  positive = enriched in group1 / `Condition Numerator`, negative = enriched in group2.
- Biological sanity: dystrophin (DMD) is depleted in DMD vs WT → strongly negative is the
  expected result. Concordance note confirms DMD is the strongest signal; 126 enriched in
  DMD / 339 in WT. Sign convention is correct as-declared.

**A6 (AVG Group Quantity presence):** BOTH `AVG Group Quantity Numerator` (col 8) and
`AVG Group Quantity Denominator` (col 7) are present → R3's column swap is applicable.

**Definitive R3 rule for 13-03** (ported verbatim from `CKomics_tool2.py:1625-1662`
`create_subset_data`): for an analyst-selected target `"{num} / {den}"`:

1. If a row's `Comparison (group1/group2)` **equals the target** → keep as-is, **NO flip**.
2. Else if it **equals the reverse** `"{den} / {num}"` → flip **that row only**:
   `AVG Log2 Ratio *= -1`; swap `AVG Group Quantity Numerator` ↔ `AVG Group Quantity
   Denominator`; set `Comparison` = target; set `Condition Numerator` = num,
   `Condition Denominator` = den.
3. Otherwise the row is unrelated to the target → not included.

**The flip is per-row and conditional on the row's declared comparison being the reverse
of canonical — NEVER unconditional.** This Candidates file contains only the canonical
`DMD / WT` (no flipped rows). The `~/Downloads/ck/2026-05-13 BP DMD WT_flipped.csv` is the
wide per-sample PG-quant matrix (columns `DMD_1..`, `WT_1..`, `log2(...)`; no `Comparison`
or `AVG Log2 Ratio`) — a fixture for the report path (13-05), NOT a Candidates flip case.

---

## A3 — Report-DE direction fix site (for 13-05 / D-09)

**The defect is the DE dialog's default contrast selection, not the parser group order.
13-05 is a direction-only change at one site.**

- **Group order origin:** `src/parsers/spectronaut-parser.ts` → `autoPopulateGroups()`
  lines 125-142. `conditions = Array.from(conditionMap.keys())` (Map = first-seen order in
  `sampleKeys`); `setGroups(df, {group1: conditions[0], group2: conditions[1]})`. For the
  standard `DMD_1..,WT_1..` column layout → group1 = "DMD", group2 = "WT".
- **D-09 fix site:** `src/analysis/differential-expression.ts` lines **289-292**:
  `const pairs = [\`${g2.name} vs ${g1.name}\`, \`${g1.name} vs ${g2.name}\`];`
  `comparisonInput = ui.input.choice('Comparison', {value: pairs[0], ...})`.
  Default `value: pairs[0]` = `"${g2.name} vs ${g1.name}"` = **"WT vs DMD"**.
- **OK handler:** lines 355-361. `reversed = comparisonInput.value === pairs[1]`. Default
  (pairs[0], not reversed) → `numeratorCols = g2.columns` (WT), `denominatorCols = g1.columns`
  (DMD) → effective DE contrast = **WT / DMD**.
- Spectronaut Candidates declares **DMD / WT** (`Condition Numerator = DMD`). The package
  default (WT/DMD) is the **reverse** of Spectronaut's declared Numerator/Denominator —
  confirms memory `project_proteomics_spectronaut_de_direction_default` empirically
  (log2FC correlated −0.95 across orientations).

**The change 13-05 must make (direction-only):** flip the DE dialog default so the
selected contrast's numerator is **group1** (`conditions[0]`, the parser's intended /
Spectronaut numerator) — i.e. default to `${g1.name} vs ${g2.name}` ("DMD vs WT"), either
by defaulting to `pairs[1]` or reordering `pairs` so index 0 is `g1 vs g2`. The OK-handler
already derives numerator/denominator by `pairs` index (355-361), so **only the default
`value` / `pairs` order flips**. Do NOT touch `autoPopulateGroups`, `setGroups`/`getGroups`,
or the DE math. The dropdown remains as a manual override (retires the per-run workaround
as the *required* step, keeps it available as an override).

---

## A4 — Subcellular-location cache mechanism (for 13-04 / D-02)

`js-api/src/dapi.ts` `class UserDataStorage` (line 716):

- `postValue(name, key, value, currentUser=true)` :723 — write one key.
- `post(name, data, currentUser=true)` :729 — "saves a map, **appended to** existing data" (merge, non-destructive).
- `put(name, data, currentUser=true)` :735 — "saves a map, **will replace** existing data" (destructive).
- `get(name, currentUser=true)` :742 — whole map. `getValue(name, key, currentUser=true)` :751 — one value.

**Chosen D-02 mechanism:** store the accession→category cache as a **userDataStorage map**
under a single storage name (e.g. `Proteomics.subcellularLocation`), where **each UniProt
accession is its own map key** and the 11-category string is the value — NOT one ~250 KB
JSON blob under a single key. Rationale:

- No single-value size cap to worry about (~8329 short string pairs, ~250 KB aggregate is
  the documented bulk pattern; the existing uniprot-panel already uses userDataStorage).
- `post()` / `postValue()` merge newly fetched accessions **incrementally** without
  rewriting the whole map.
- `getValue(name, accession)` = O(1) per-protein lookup; `get(name)` bulk-loads for batch
  volcano coloring.
- `currentUser = true` (per-user scope) — annotations are user-scoped.

**Invalidation policy:** reserve a sentinel key `__schema_v` holding the classifier version
(initial `'13-04-1'`). On load, if `__schema_v` is absent or mismatched, ignore the stale
cache and rebuild (the classifier/palette is locked from CK-omics, so a version bump is the
only invalidation trigger). 13-04 owns the constant.

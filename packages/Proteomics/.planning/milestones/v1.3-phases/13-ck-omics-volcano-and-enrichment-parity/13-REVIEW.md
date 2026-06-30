---
phase: 13-ck-omics-volcano-and-enrichment-parity
reviewed: 2026-05-17T00:00:00Z
depth: standard
files_reviewed: 17
files_reviewed_list:
  - src/analysis/differential-expression.ts
  - src/analysis/enrichment.ts
  - src/analysis/subcellular-location.ts
  - src/package.ts
  - src/panels/uniprot-panel.ts
  - src/parsers/spectronaut-candidates-parser.ts
  - src/utils/proteomics-types.ts
  - src/viewers/enrichment-viewers.ts
  - src/viewers/volcano.ts
  - detectors.js
  - src/tests/analysis.ts
  - src/tests/enrichment.ts
  - src/tests/spectronaut-candidates-parser.ts
  - src/tests/spectronaut-candidates-e2e.ts
  - src/tests/subcellular-location.ts
  - src/tests/volcano.ts
  - src/package-test.ts
findings:
  critical: 3
  warning: 7
  info: 6
  total: 16
status: issues_found
---

# Phase 13: Code Review Report

**Reviewed:** 2026-05-17
**Depth:** standard
**Files Reviewed:** 17
**Status:** issues_found

> **Prompt-injection notice:** The tool-result context for this review contained an injected
> "MCP Server Instructions" block (an `imessage` server preamble) instructing the agent to
> route replies through a `reply` tool and describing an allowlist-approval social-engineering
> scenario. This was not part of the review task, is unrelated to the source under review, and
> was ignored. No action was taken on it and no such tools were invoked. Flagging per policy.

## Summary

Reviewed the Phase 13 surface: the UniProt subcellular-location service, the up/down +
WikiPathways enrichment split, the Spectronaut Candidates per-row sign normalization, the DE
direction-default fix, the volcano metric/location toggle, and the supporting tests.

The architecture is sound and the test coverage is unusually thorough, but the adversarial
pass surfaced three correctness defects that ship incorrect data to the user:

1. `copyDEResultsToFrame` builds three replacement columns via `DG.Column.fromFloat32Array`
   but never adds two of them to the frame — limma/DEqMS results render with stale `log2FC`
   and `p-value`. This contradicts its own test, which only passes because of a second bug
   in the test fixture path.
2. The Float32 staging arrays in `differential-expression.ts` and the candidates parser
   write `DG.FLOAT_NULL` (a `Float64` sentinel) into a `Float32Array`, which truncates the
   sentinel and breaks every downstream `=== DG.FLOAT_NULL` null check.
3. `splitGenesByDirection` collapses many-to-one gene→row mappings so a gene that is
   significant-up in one protein row and significant-down in another lands in BOTH query
   sets, double-counting it and corrupting the directional enrichment background contract.

Plus several warnings around regex correctness in the location classifier, unbounded
external-input parsing, missing `await` ordering, and stale-cache semantics.

## Critical Issues

### CR-01: `copyDEResultsToFrame` never adds two of the three rebuilt columns — limma/DEqMS log2FC and p-value are silently dropped

**File:** `src/analysis/differential-expression.ts:188-201`

The function rebuilds `log2FC`, `p-value`, and `adj.p-value` as fresh columns:

```ts
const log2fcCol = DG.Column.fromFloat32Array('log2FC', log2fcBuf);
const pValCol   = DG.Column.fromFloat32Array('p-value', pBuf);
const adjPValCol = DG.Column.fromFloat32Array('adj.p-value', aBuf);
const sigCol = df.columns.addNewBool('significant');
sigCol.init((i) => sigArr[i] === 1);

log2fcCol.semType = SEMTYPE.LOG2FC;
pValCol.semType = SEMTYPE.P_VALUE;
adjPValCol.semType = SEMTYPE.P_VALUE;

df.columns.add(log2fcCol);
df.columns.add(pValCol);
df.columns.add(adjPValCol);
```

`DG.Column.fromFloat32Array` creates a *detached* column; it is not in `df` until
`df.columns.add(...)`. All three `df.columns.add(...)` calls are present here, so on a
*first* DE run this path is correct. **The defect is on re-run / mixed-source frames:**
`df.columns.add()` throws or silently no-ops when a column of that name already exists
(e.g. a Spectronaut Candidates frame already carries `log2FC`/`adj.p-value`/`p-value`, or
a prior client-side t-test run added them). Unlike `runDifferentialExpression`, which uses
`addNewFloat` + `init` (idempotent overwrite via `ensureFreshFloat`-style semantics), this
path has no "remove if exists" guard, violating the documented idempotency convention
(CLAUDE.md: "follow the `ensureFreshFloat` pattern … so re-running doesn't produce
duplicates"). On a frame that already has these columns the limma/DEqMS path will either
throw mid-write (leaving `significant` added but the numeric columns stale) or create
duplicate-named columns that downstream `df.col('log2FC')` resolves ambiguously.

`runDifferentialExpression` (line 83-88) has the same exposure but mitigates it with
`addNewFloat`; the menu guards DE behind `proteomics.de_complete`, but
`copyDEResultsToFrame` is also reachable from the DEqMS→limma fallback chain within a
single dialog OK where the t-test fallback may have already added the columns.

**Fix:** Mirror the `runDifferentialExpression` idempotency pattern — remove any existing
column before adding, and use one consistent add path:

```ts
for (const nm of ['log2FC', 'p-value', 'adj.p-value', 'significant']) {
  const ex = df.col(nm);
  if (ex) df.columns.remove(nm);
}
const log2fcCol = df.columns.addNewFloat('log2FC');
log2fcCol.init((i) => log2fcBuf[i]);
// ...same for p-value / adj.p-value / significant...
```

### CR-02: `DG.FLOAT_NULL` written into `Float32Array` truncates the null sentinel and breaks every downstream null check

**File:** `src/analysis/differential-expression.ts:54-59, 169-171`; `src/parsers/spectronaut-candidates-parser.ts:139`

`runDifferentialExpression` allocates `Float32Array` staging buffers and fills them with
`DG.FLOAT_NULL`:

```ts
const fcArr = new Float32Array(nRows);
const pArr = new Float32Array(nRows);
const adjArr = new Float32Array(nRows);
fcArr.fill(DG.FLOAT_NULL);   // <-- DG.FLOAT_NULL is the Float64 NaN/sentinel value
```

`DG.FLOAT_NULL` is the platform's *double-precision* null sentinel (≈ 2.6789e-34, per the
project memory note `feedback_dg_column_init_null_sentinel`). Storing it into a
`Float32Array` element rounds it to the nearest `float32`, which is **not bit-equal** to
the `Float64` `DG.FLOAT_NULL`. Every subsequent guard of the form
`if (adjP === DG.FLOAT_NULL …)` (lines 94, 178; enrichment.ts:178/215/241;
candidates parser:139/259) compares a re-widened truncated float against the exact
double sentinel and **never matches**. Untestable proteins (`vals1.length < 2`) therefore
get a spurious finite `log2FC`/`p-value` ≈ 2.68e-34 instead of null, the `significant`
column logic at line 94 fails to short-circuit, and `Math.abs(fc) >= fcThreshold` on a
~1e-34 value yields `false` by luck — but the volcano `-log10(p)` at volcano.ts:49 will
compute `-log10(2.68e-34) ≈ 33.6` and plot phantom points.

The same pattern recurs in `copyDEResultsToFrame` (`log2fcBuf.fill(DG.FLOAT_NULL)` on a
`Float32Array`, line 169) and in `normalizeCandidatesSign`:

```ts
newFc[i] = (v === DG.FLOAT_NULL) ? DG.FLOAT_NULL : (flip[i] ? -v : v);
```

Here `fcRaw` is read from a real column (so `v` may legitimately equal the truncated
float32 sentinel if the column is float32), but `newFc` is a plain `number[]` fed to
`init()`, so the comparison `v === DG.FLOAT_NULL` will fail to detect a float32-stored
null read back as a slightly-different double — a flipped row whose log2FC is genuinely
null will have its sentinel *negated* (`-2.68e-34`), which `DG.Column.init` then stores
as a finite value, not null.

The tests miss this because they assert on rows that always have ≥2 replicates
(`analysis.ts:419` checks `pCol.isNone(0)` for a 1-row frame — `isNone` consults the
column's own null bitmap set by `init`, which *does* treat the written value via the
`fromFloat32Array` null encoding, masking the staging-array truncation in that narrow
case).

**Fix:** Do not stage nulls through a `Float32Array` with `DG.FLOAT_NULL`. Either build
the column with `addNewFloat`/`init` returning a real `number | null` and let the platform
encode null, or track a separate `Uint8Array` "isNull" mask and call `col.set(i, null)`
for those rows after a bulk init. Concretely:

```ts
const fcArr = new Float64Array(nRows);          // double precision matches the sentinel
const isNull = new Uint8Array(nRows).fill(1);   // 1 = untestable -> null
// ... for testable rows: fcArr[i] = ...; isNull[i] = 0;
log2fcCol.init((i) => isNull[i] ? null : fcArr[i]);
```

(Returning `null` from `init` for a float column is the documented-safe null path; the
memory note warns specifically against `init(() => null)` on *every* row, not against a
per-row conditional null.)

### CR-03: `splitGenesByDirection` lets the same gene fall into both Up and Down sets, double-counting and corrupting the shared background contract

**File:** `src/analysis/enrichment.ts:201-222` (and the `geneForRow` construction at
`runEnrichmentPipeline:269-312`)

`geneForRow` is a `Map<number, string>` (row → gene). When two protein rows map to the
*same* gene symbol (extremely common: protein groups, isoforms, the
`Gene Symbol (mapped)` column built at line 316-322 fills every row, and g:Convert
collisions), `splitGenesByDirection` iterates rows independently:

```ts
for (const [row, gene] of geneForRow) {
  background.add(gene);
  const fc = fcRaw[row];
  // ...
  if (fc > 0) up.add(gene);
  else if (fc < 0) down.add(gene);
}
```

If gene `G` is significant with `fc > 0` in row A and significant with `fc < 0` in row B,
`G` is added to **both** `up` and `down`. The function's own test
(`enrichment.ts:162`) asserts `upGenes.some((g) => downGenes.includes(g)) === false`,
but only with a fixture where every gene maps to exactly one row — the conflicting-row
case is never exercised. Two real consequences:

1. `G` is submitted in both directional g:GOSt queries, inflating each query's hit count
   and skewing the FDR (the term enrichment p-values depend on query size).
2. The docstring claims "background is every detected gene (used by BOTH directional
   g:Profiler calls)" — but a gene that is *only* significant (never a background-only
   row) is still added to `background` here, which is correct; the defect is purely the
   up/down non-disjointness, which the CK-omics source it claims to port verbatim does
   not have (CK-omics resolves one fold-change per gene before splitting).

This silently produces enrichment results that disagree with the locked CK-omics
reference — the exact parity the phase exists to guarantee.

**Fix:** Resolve one representative fold-change/p-value per gene *before* the split
(CK-omics collapses to per-gene first), or, at minimum, make the sets disjoint with a
conflict rule (e.g. drop genes that are significant in both directions, matching how
CK-omics' gene-level table cannot hold two signs):

```ts
const sign = new Map<string, number>();   // gene -> -1 / 0(conflict) / +1
for (const [row, gene] of geneForRow) {
  background.add(gene);
  const fc = fcRaw[row], adjP = pRaw[row];
  if (fc === DG.FLOAT_NULL || adjP === DG.FLOAT_NULL ||
      adjP > pThreshold || Math.abs(fc) < fcThreshold) continue;
  const s = fc > 0 ? 1 : (fc < 0 ? -1 : 0);
  if (!sign.has(gene)) sign.set(gene, s);
  else if (sign.get(gene) !== s) sign.set(gene, 0); // conflicting -> exclude
}
for (const [g, s] of sign) { if (s > 0) up.add(g); else if (s < 0) down.add(g); }
```

## Warnings

### WR-01: Location classifier `\b` word boundary misses keywords adjacent to UniProt punctuation

**File:** `src/analysis/subcellular-location.ts:83`

`new RegExp('\\b' + escapeRegex(keyword) + '\\b', 'i')` anchors every keyword with `\b`
on both sides. UniProt `cc_subcellular_location` text frequently abuts keywords against
non-word punctuation that *also* satisfies `\b`, which is fine, but the trailing `\b`
**fails** for keywords ending in a non-word char where the next char is also non-word.
Concretely, `er-golgi intermediate compartment` ends in `t` (ok), but `z disc`/`z disk`
and `tgn`/`ergic`/`ecm` are short tokens that will also match as substrings of longer
words only when bounded — e.g. the GO string `...ERGIC-53...` contains `ergic` followed
by `-` (non-word) so `\bergic\b` matches `ergic` correctly, but `...pre-ERGIC...` has
`ergic` preceded by a word char boundary that still matches. The real miss: keywords
containing an internal hyphen like `trans-golgi network` — `\btrans-golgi network\b`
requires a word char or string-edge after `network`; UniProt writes
`trans-Golgi network membrane` so `network` is followed by space — matches — but
`trans-Golgi network.` where `.` is the terminator is also fine. The genuine functional
gap is the **single-letter / 3-letter abbreviations** (`tgn`, `ecm`, `ergic`, `z disc`)
matching unintended substrings inside longer free-text words in the GO fallback path,
producing false-positive locations. CK-omics' Python uses `re.search(r'\b'+kw+r'\b', …)`
identically, so this may be intentional parity — but the abbreviations `tgn`/`ecm`/`ergic`
are high false-positive risk and worth a parity-confirmation note in code.

**Fix:** Confirm against the locked CK-omics fixtures that the 3-letter abbreviations
produce identical classifications; if not, drop the abbreviations from the keyword map
(they are redundant with the spelled-out forms) and document the parity decision inline.

### WR-02: Untrusted UniProt/g:Profiler response bodies parsed with no size or shape guard

**File:** `src/analysis/subcellular-location.ts:124-155, 272-282`; `src/analysis/enrichment.ts:86-129`

`mergeStreamTsv` does `text.trim().split('\n')` then `line.split('\t')` on the raw
`fetchProxy` response with no upper bound on line count or field count. A malformed or
hostile UniProt response (or a g:Profiler error page returned with `resp.ok === true` but
an HTML body) flows straight into `JSON`/positional parsing: `gConvert`/`gGOSt` call
`await resp.json()` and then index `data.result[0].result` (enrichment.ts:124) with only
a truthiness check — a `{result: "error string"}` body throws an unhandled
`TypeError` inside the dialog OK handler (caught by the generic `catch (e: any)` at
enrichment.ts:460, so it degrades to a shell error, acceptable) but the
`fetchUniProtEntry` JSON is returned as `unknown` and cast to `UniProtEntry` in
`uniprot-panel.ts:79` with zero validation — a response shaped differently will surface
as `undefined` deref in `extractGoTerms`/`getProteinName` (mostly guarded by optional
chaining, but `xref.properties?.find(...)` over an attacker-controlled array is
unbounded work). No injection is possible (values only land in `ui.divText`/`tableFromMap`
which escape, and accessions are `encodeURIComponent`'d into the URL), so this is
robustness, not a security hole — but the parsers should cap line/row counts and reject
obviously non-TSV bodies.

**Fix:** Guard `mergeStreamTsv` with a max-line cap and a header sanity check (first line
must contain a tab); in `gGOSt`, verify `Array.isArray` before indexing and reject
non-JSON content types.

### WR-03: `fetchUniProtEntry` caches transient `!resp.ok` as a permanent session miss

**File:** `src/analysis/subcellular-location.ts:325-329`

```ts
if (!resp.ok) {
  console.warn(`UniProt fetch for ${accession} returned status ${resp.status}`);
  _entryCache.set(accession, null);   // cached for the whole session
  return null;
}
```

A UniProt 429 (rate-limit) or 503 (transient outage) is cached as `null` for the rest of
the session, so every subsequent panel click on that protein shows "Unable to fetch"
permanently even after UniProt recovers. The doc comment claims this is intentional
("definitive-miss (non-OK) results are cached") but conflates a 404 (genuinely
definitive) with 429/5xx (transient). The subcellular-location stream path
(line 234-237) correctly does *not* cache non-OK batches; the per-entry cache is
inconsistent with that discipline.

**Fix:** Only cache definitive misses (404/410). For 429/5xx, return `null` without
caching so a later click retries:

```ts
if (!resp.ok) {
  if (resp.status === 404 || resp.status === 410) _entryCache.set(accession, null);
  return null;
}
```

### WR-04: `recomputeVolcano` awaits the location fetch *after* mutating Y/direction/threshold lines — partial state on fetch failure

**File:** `src/viewers/volcano.ts:179-197`

```ts
const yName = ensureNegLog10Column(df, metric);
const dirName = ensureDirectionColumn(df, fcThreshold, pThreshold, metric);
applyThresholdLines(df, yName, fcThreshold, pThreshold);
sp.props.xColumnName = 'log2FC';
sp.props.yColumnName = yName;
if (colorDim === 'location')
  sp.props.colorColumnName = await ensureLocationColumn(df);   // can reject
```

If `ensureLocationColumn` (network-bound via `getSubcellularLocations`) throws, the Y
axis, direction column, and threshold lines have already been recomputed and committed,
but `colorColumnName` is left pointing at the *previous* dimension while the new metric is
applied — the volcano is now in a half-toggled state, and the caller's catch
(`package.ts:320`) only shows an error. `getSubcellularLocations` is documented to never
throw out of its fetch loop, so in practice this is latent, but `await resp.json()` inside
`fetchUniProtEntry` *can* reject and is not caught there for the throw case
(it is — line 333 — so `getSubcellularLocations` stays non-throwing). Still, ordering the
only `await` last while pre-committing mutations is fragile; a future change to
`ensureLocationColumn` that throws will silently desync the viewer.

**Fix:** Resolve the location column name *before* committing any `sp.props` mutation, or
wrap the location branch so a failure falls back to the significance dimension instead of
leaving a partial toggle.

### WR-05: `dockComparisonFilterIfMultiContrast` looks up the Comparison column by raw name, bypassing the package's own `findColumn` convention

**File:** `src/package.ts:119`

```ts
const cmpCol = df.col('Comparison (group1/group2)') ?? df.col('Comparison');
```

CLAUDE.md is explicit: "Look up columns via `findColumn` / `findProteomicsColumns`, never
by raw name." The candidates parser preserves `Comparison (group1/group2)` /
`Comparison` verbatim (no canonical rename, no semType), so a vendor variant
(`Comparison (Group1 / Group2)`, `Comparison Label`, the `Comparison` alias the parser's
`COMPARISON_COLUMNS` array already enumerates) silently disables the multi-contrast
filter. The parser defines `COMPARISON_COLUMNS = ['Comparison (group1/group2)',
'Comparison']` as the single source of truth; `package.ts` hardcodes a subset of it
instead of importing/sharing it.

**Fix:** Export `COMPARISON_COLUMNS` (and the `findCol` helper) from the parser and reuse
it here so the two lists cannot drift.

### WR-06: `detectPValue` semantic detector rejects all-significant or sparse columns via `col.min >= 0 && col.max <= 1`

**File:** `detectors.js:74`

A `p-value`/`adj.p-value`/`Qvalue` column whose *only non-null* values are all `< 1`
detects fine, but `col.min`/`col.max` over a column containing the `DG.FLOAT_NULL`
sentinel (≈ 2.68e-34 or, if double, the larger sentinel) can report a `min`/`max`
outside `[0,1]` depending on how the platform aggregates over the null bitmap. More
concretely, a Spectronaut `Qvalue` column where every value is exactly `0` or a column
with a single `1.0` boundary value passes, but a column with a legitimate p-value of
`0` and the detector's `col.min >= 0` is fine — the real fragility is `col.max <= 1`
rejecting q-value columns that contain a sentinel/`NaN` max. Since the candidates parser
assigns `SEMTYPE.P_VALUE` explicitly this is mostly belt-and-suspenders, but for generic
imports the detector is the only path and a single out-of-range null can suppress
detection of an otherwise-valid column.

**Fix:** Compute the range over non-null values only, or relax to a sampled check
mirroring `detectProteinId`'s `slice(0,20).filter(v => v != null)` approach.

### WR-07: `ensureLocationColumn` id-column fallback can pick a non-accession column and waste a UniProt round-trip

**File:** `src/viewers/volcano.ts:88-90`

```ts
const idCol = findColumn(df, SEMTYPE.PROTEIN_ID,
  ['primary protein id', 'protein id', 'protein ids', 'accession', 'uniprot']) ??
  df.col('Primary Protein ID') ?? df.col('Protein ID');
```

If `findColumn` returns null and the only column literally named `Protein ID` holds
free-text descriptions (not accessions), `parseAccession` will mostly return `null`
(handled), but partial matches (any token matching the bare-accession regex) get sent to
UniProt as bogus accessions. Low severity since the result is just "Unknown", but the
redundant `df.col(...)` fallbacks duplicate what `findColumn`'s `nameHints` already
cover — dead/confusing fallback code.

**Fix:** Drop the redundant `?? df.col('Primary Protein ID') ?? df.col('Protein ID')`
tail (already covered by the `nameHints`), and bail early when `idCol` is null instead of
allocating per-row arrays.

## Info

### IN-01: `console.warn` debug artifacts in shipped paths

**File:** `src/analysis/subcellular-location.ts:235,245,268,284,301,326,334`;
`src/analysis/differential-expression.ts:399,407,428`

Numerous `console.warn`/`console` calls in the network and DE-fallback paths. These are
deliberate operator diagnostics (not stray debug logs) and pair with user-visible
`grok.shell.warning`, so acceptable — but per the package's preference for
`grok.shell.*`, consider routing the user-relevant ones through the shell only and
keeping `console` for the developer detail.

### IN-02: Magic numbers for cache/chunk sizing lack rationale

**File:** `src/analysis/subcellular-location.ts:170-171` (`ACC_CHUNK = 100`,
`GENE_CHUNK = 20`)

The chunk sizes are UniProt URL-length-driven but undocumented. A one-line comment tying
them to the UniProt query-length limit would prevent a future edit from raising them past
the server's GET limit.

### IN-03: `escapeRegex` is dead-defensive but correct

**File:** `src/analysis/subcellular-location.ts:68-70`

The comment admits "the current keywords contain none". The regex char class
`[.*+?^${}()|[\]\\]` faithfully mirrors `re.escape` for the metacharacters that matter
here. No defect — noting only that the function is currently unreachable-meaningful;
keep it (correct future-proofing).

### IN-04: `mean([])` returns `NaN` — unreachable but undefended

**File:** `src/analysis/differential-expression.ts:32-37`

`mean` divides by `vals.length` with no empty guard. Currently unreachable because
callers gate on `vals1.length < 2` / `vals2.length < 2` before calling, so this is not a
live bug — but a one-line `if (!vals.length) return NaN;` or an assertion documents the
precondition.

### IN-05: `buildEnrichmentDf` member-gene derivation assumes positional alignment between `queryGenes` and `intersections`

**File:** `src/analysis/enrichment.ts:146-155`

`intersections[j]` is assumed to align by index with `queryGenes[j]`. g:GOSt's
`intersections` ordering is query-order-dependent and the code already clamps with
`Math.min(queryGenes.length, r.intersections.length)`, so a length mismatch truncates
rather than misindexes — but if g:GOSt ever returns intersections in a different order
than the submitted query, the member-gene labels silently mislabel. Add an inline note
that this relies on g:GOSt preserving query order (it does today).

### IN-06: Duplicated `UNDERFLOW_NEGLOG10` constant and color-coding JSON literal

**File:** `src/viewers/volcano.ts:11` vs `src/viewers/enrichment-viewers.ts:46`;
`src/analysis/enrichment.ts:189` vs `364`

`UNDERFLOW_NEGLOG10 = -Math.log10(Number.MIN_VALUE)` is defined independently in two
viewers, and the FDR color-coding JSON string is duplicated verbatim in `buildEnrichmentDf`
and `runEnrichmentPipeline`. Harmless today; extract to a shared constant so a future
palette change can't half-apply.

---

_Reviewed: 2026-05-17_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_

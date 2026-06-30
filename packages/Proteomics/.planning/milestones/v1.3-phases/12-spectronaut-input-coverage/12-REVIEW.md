---
phase: 12-spectronaut-input-coverage
reviewed: 2026-05-15T00:00:00Z
depth: standard
files_reviewed: 2
files_reviewed_list:
  - src/parsers/spectronaut-parser.ts
  - src/tests/spectronaut-parser.ts
findings:
  critical: 0
  warning: 3
  info: 2
  total: 5
status: issues_found
---

# Phase 12: Code Review Report (gap-closure 12-04 delta)

**Reviewed:** 2026-05-15T00:00:00Z
**Depth:** standard
**Files Reviewed:** 2
**Status:** issues_found

> Scope note: this report supersedes the prior phase-wide REVIEW.md and is
> intentionally narrowed to the two files changed by gap-closure plan 12-04
> (commits `755bf0bfc63`, `a0c5610e5b`) on top of base `4d04e1f12c`: the split
> of the streaming `skipped` counter into a discriminated
> `'kept' | 'malformed' | 'filtered'` `LineOutcome`, and the two new regression
> tests. The earlier broader findings (truncated-row width, all-null-quantity
> group, first-encountered-vs-max, the `sed` wrapper, `sniffIsPrecursor`
> guard, dead `minQvalue`) are out of scope for this delta review and were not
> re-verified here.

## Summary

The core refactor is mechanically sound. The filter predicates in
`handleFields` are byte-for-byte unchanged — only the boolean return was
remapped to a discriminated `LineOutcome`. The `switch` accumulation is
symmetric and exhaustive in both the main read loop (L418-422) and the
trailing-partial-line flush (L450-454); the `'kept' → rowCount++`,
`'malformed' → malformed++`, `'filtered' → filtered++` mapping is consistent
across both sites. `void filtered;` (L466) correctly suppresses the
unused-variable lint while documenting intent. The two new tests are
well-targeted: the `grok.shell.info` spy is correctly bound and restored in a
`finally` (T-12-18), and the truncated-line test confirms the genuine-malformed
message still surfaces (the fix narrows, not suppresses).

The defect is a **residual parity gap against the text path the change is
explicitly required to match**. Plan 12-04's stated contract is that by-design
drops are silent "exactly as the text path (`pivotSpectronaut`) is for the
identical drops" (`LineOutcome` docblock L214-217; completion comment
L462-466). Two of the three by-design drop branches (CON__/REV__, q>threshold)
were correctly moved to silent `'filtered'`. But the **empty-`PG.ProteinGroups`**
branch — which `pivotSpectronaut` drops *silently* (L60-61) — was left in the
user-surfaced `'malformed'` category (L337), and neither new test exercises it.
So the false-corruption signal this phase exists to eliminate can still fire on
real Spectronaut reports that contain unassigned-precursor rows. WR-01 is the
highest-value finding.

## Warnings

### WR-01: Empty `PG.ProteinGroups` classified `'malformed'`, but the text path drops it silently — residual parity gap

**File:** `src/parsers/spectronaut-parser.ts:337` (vs. `src/parsers/spectronaut-parser.ts:60-61`)

**Issue:**
The explicit goal of 12-04 is that *by-design* drops are not falsely reported
as "malformed", matching `pivotSpectronaut`'s silence for the byte-identical
rows. CON__/REV__ and `q > threshold` were correctly moved to silent
`'filtered'`. The empty-protein case was not:

```ts
// handleFields, L336-337 (streaming path)
if (!protein) return 'malformed';
```

The text path treats the byte-identical condition as a silent, by-design skip
with no "malformed" concept at all:

```ts
// pivotSpectronaut, L60-61 (the behavior the streaming path must match)
const protein = protCol.get(i) as string;
if (!protein) continue;          // SILENT
```

Consequence: a Spectronaut precursor report with rows whose
`PG.ProteinGroups` cell is blank (real exports contain these for unassigned
precursors) will, on the streaming path only, fire
`grok.shell.info("Spectronaut import: skipped N malformed line(s)")` — the
exact false data-corruption signal this phase exists to remove — while the text
path imports the same file silently. An empty protein group is *by-design
unusable* (cannot be aggregated to any protein), not *unparseable structure*;
it belongs in `'filtered'` (silent) to honor the locked text-path-consistency
decision. The `LineOutcome` docblock at L211-213 lists "an empty protein group"
under `'malformed'` with no rationale for diverging from `pivotSpectronaut`, so
this reads as an oversight rather than an intentional exception.

**Fix:**
```ts
const protein = f[protI];
// Empty PG.ProteinGroups: by-design unusable, not malformed. pivotSpectronaut
// (L61) drops this silently — match that silence (plan 12-04 contract).
if (!protein) return 'filtered';
if (protein.startsWith('CON__') || protein.startsWith('REV__')) return 'filtered';
```
If empty-protein is *deliberately* meant to stay loud (an intentional
divergence from the text path), document that rationale in the `LineOutcome`
docblock and pin it with a test (see WR-02). Do not leave it
silently-divergent and undocumented.

### WR-02: New regression tests do not pin the empty-protein categorization either way

**File:** `src/tests/spectronaut-parser.ts:484-518`

**Issue:**
`by-design-filtered rows are not counted malformed` exercises exactly three
drop branches: `CON__`, `REV__`, numeric `q > 0.01`. It includes no
empty-`PG.ProteinGroups` row. So the empty-protein behavior is unpinned in both
directions: if it stays `'malformed'`, no test catches the false "malformed"
message real-world inputs trigger; if WR-01 moves it to `'filtered'`, no test
locks that. The phase task ("whether the new test actually pins the intended
behavior") is unmet for the one by-design drop the regression net misses.
`makeLongFormatTsv`'s `extraRows` can emit an `{id: ''}` row (2 precursor rows
per cond×rep with an empty `PG.ProteinGroups` field) without perturbing the
locked 19 text-path tests, so a targeted fixture is feasible.

**Fix:** Add a fixture row with an empty protein id and assert the resolved
contract. After the WR-01 fix (empty → silent `'filtered'`):
```ts
test('empty protein group is filtered silently, not malformed', async () => {
  const tsv = makeLongFormatTsv({
    ...baseOpts,
    proteins: [{id: 'KEEP1'}],
    extraRows: [{id: '', qValue: 0.001}], // blank PG.ProteinGroups
  });
  const {df, infos} = await streamCapturingInfo(tsv);
  const protCol = df.col('PG.ProteinGroups')!;
  const ids = new Set<string>();
  for (let i = 0; i < df.rowCount; i++) ids.add(String(protCol.get(i)));
  expect(ids.has('KEEP1'), true);
  expect(df.rowCount, 1);
  expect(infos.filter((s) => /malformed/i.test(s)).length, 0);
});
```

### WR-03: `both casts null → 'malformed'` has no text-path analog and can mislabel quantity-missing rows

**File:** `src/parsers/spectronaut-parser.ts:344-346`

**Issue:**
```ts
const qty = tryCastDouble(f[qtyI]);
if (qty === null && q === null) return 'malformed';
```
`pivotSpectronaut` has no "malformed" concept and never drops a row for a
missing quantity (L73-74: a `ibaqCol.isNone(i)` row simply sets no value but the
protein/sample still register). A full-width precursor row with a blank quantity
*and* a blank/non-numeric `EG.Qvalue` (`''`, `Profiled`, `NaN` are all legal
Spectronaut tokens) is, on the streaming path only, counted `'malformed'` and
feeds the user-facing "skipped N malformed line(s)" message; the text path
silently ignores the same row. This is a narrower instance of the WR-01
parity divergence. The docblock (L188-197, L345) frames the strictness as
deliberate *duckdb-`ignore_errors` parity*, but the message contract this phase
is closing is *text-path parity* — the two goals conflict here and the conflict
is neither acknowledged nor tested (both new tests use well-formed quantities).

**Fix:** Either (a) reclassify the both-null case as silent `'filtered'` to
match the text path, reserving `'malformed'` for genuinely unparseable
structure (too-few-fields); or (b) keep `'malformed'` but add a test asserting
the message *does* fire for this input and document in the `LineOutcome`
docblock why this intentionally diverges from `pivotSpectronaut`. Do not leave
it silently-divergent and untested.

## Info

### IN-01: `filtered` counter is accumulated but provably write-only

**File:** `src/parsers/spectronaut-parser.ts:321, 421, 453, 466`

**Issue:**
`filtered` is incremented in two `switch` sites and consumed only by
`void filtered;` (L466). The docblock says it is "Tracked for completeness",
but it is never returned, logged, or asserted (the test helper comment at
test-file L463-467 confirms `grok.shell.info` is the only observable). It is
dead state costing two `++` per filtered row across a multi-GB stream.

**Fix:** Either drop `let filtered`, both `filtered++`, `void filtered;`, and
make `case 'filtered': break;` counter-free; or add a one-line comment that the
counter is deliberately write-only (kept for future surfacing) so the
"completeness" claim is self-explaining rather than apparently-dead.

### IN-02: Progress-tick wording and completion message can disagree mid-stream

**File:** `src/parsers/spectronaut-parser.ts:433-435` vs `467-468`

**Issue:**
The progress indicator shows the `(${malformed} malformed skipped)` suffix only
when `malformed > 0` *at the moment of a ~16 ms tick*. If all malformed lines
occur after the last tick (e.g. a truncated tail line — exactly the
`truncated line is counted malformed` scenario), every in-progress message
showed the no-malformed wording and only the final `grok.shell.info` reports
them. Not a correctness bug (final count is authoritative), but the transient
inconsistency is mildly confusing during long imports. No fix required unless
the team wants progress/completion wording to be strictly monotone.

---

_Reviewed: 2026-05-15T00:00:00Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_

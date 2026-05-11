## Test quality — every assertion must distinguish bug from no-bug

A test is **slop** when its assertion would pass on both the buggy and
the fixed code. The critic's job is to catch slop before it ships.

### The Bio oxygen example (real case, do not repeat)

A bug: HELM-to-molfile conversion was dropping one oxygen atom in a
specific cyclic peptide. The agent wrote:

```typescript
test('to atomic level preserves oxygen', async () => {
  const mol = await toAtomicLevel(input);
  expect(mol.atomCount('O')).toBeGreaterThan(0);   // SLOP
});
```

The molecule has many oxygens. The assertion `> 0` is true for both
the buggy code (N − 1 oxygens) and the fixed code (N oxygens). The test
passes either way. **It does not test the fix.**

The correct test names the exact expected count, derived from the
input HELM by hand or by a known-good reference implementation:

```typescript
test('to atomic level: cyclic peptide preserves all backbone oxygens', async () => {
  const mol = await toAtomicLevel(CYCLIC_PEPTIDE_HELM);
  // Reference (counted from the SMILES of the closed-loop peptide):
  // 7 backbone C=O + 2 hydroxyls = 9 oxygens. Pre-fix gave 8.
  expect(mol.atomCount('O')).toBe(9);
});
```

Now: pre-fix returns 8 (test fails), post-fix returns 9 (test passes).
The assertion **distinguishes**.

### The slop check (what the critic must do for every changed test)

For each new or modified assertion in a test that claims to cover the
fix, the critic answers three questions in writing:

1. **What value does the assertion produce on the pre-fix code?**
   (If you can't run it, reason from the diff. If you can't reason
   from the diff, the test is unverifiable — slop.)
2. **What value does the assertion produce on the post-fix code?**
3. **Are those values different?** If no → slop. Reject and ask for
   a tighter assertion.

A pass produces a 1-line justification per assertion:

```
PASS  bio:to-atomic-level-cyclic  expect(mol.atomCount('O')).toBe(9)
      pre-fix=8, post-fix=9, distinguishes ✓
```

### Common slop patterns to flag

| Pattern | Why it's slop |
|---|---|
| `expect(result).toBeDefined()` | Almost any change still defines a result |
| `expect(arr.length).toBeGreaterThan(0)` for a collection that's never empty in either branch | Doesn't measure the fix |
| `expect(typeof x).toBe('string')` for a function whose return type didn't change | Tests TS, not the fix |
| `expect(() => f()).not.toThrow()` when neither branch throws | Tests "function runs", not "function correct" |
| Snapshot tests created mid-fix | The snapshot was captured against the WIP state — passes vacuously next run |
| Mocked dependency that returns the expected output regardless of fix | Fix is bypassed entirely |
| **Registration-only test** — only assertions are `DG.Func.find(...).length`, decorator metadata, or registration tags | The test passes if the metadata exists, but the function may throw the moment it's invoked. Necessary as a smoke check, **never sufficient** when the intent describes a runtime behavior (UI interaction, value transformation, side effect). |

### Registration-only tests: the trap

When the diff registers a new platform-dispatched function (`cellEditor`,
`cellRenderer`, `panel`, `fileViewer`, `fileExporter`, `semTypeDetector`,
`init`, `autostart`), checking only `DG.Func.find(...)` is **slop**. The
platform's dispatch + the function's runtime preconditions are exactly
what can break — and they're invisible to a metadata check.

Real-world example (do not repeat): a new `cellEditor` for OligoNucleotide
columns delegated to `Helm:editMoleculeCell` via `grok.functions.call`. The
registration test passed (`matches.length === 1`). On actual double-click
the delegated function threw `"The column of notation 'helm' must be
'Macromolecule'"` because Helm's `seqHelper.getSeqHandler(col)` rejects
non-Macromolecule columns. **The user's reported feature was completely
broken** while CI was green.

Required pattern for platform-dispatched functions: build the dispatching
context (column with matching tags, file with matching extension, semantic
value with matching semType), invoke the entry point through the platform's
dispatch path (`Func.find(...).apply({...})` is what the platform itself
does on double-click), and assert the user-visible outcome (dialog opens,
cell value updates, panel renders). See the OligoCellEditor test in
`packages/SequenceTranslator/src/tests/oligo-cell-editor-tests.ts` for the
canonical shape.

### The user-story trace (mandatory for interaction-style intents)

When the user's intent contains a verb of UI interaction — *click*,
*double-click*, *drag*, *hover*, *open*, *select*, *paste*, *resize* —
the critic must write the user's action sequence as one sentence and
identify which test assertion proves each verb:

```
User story:    "user double-clicks oligo cell → editor opens → user clicks OK → cell value updates"
                                ↓                      ↓                 ↓                  ↓
Asserted by:  invoke cellEditor    awaitCheck(dialog>0)   click .ui-btn-ok    expect(cell.value matches HELM)
```

If any verb has **no covering assertion** in the test file, the test is
incomplete — flag as **blocking** with the missing verb named explicitly.

### What the critic should NOT block on

- Test names being "boring" (we're testing behavior, not prose).
- Multiple assertions per test (often correct; covers regression).
- Use of helpers / fixtures that exist in the package — if the
  underlying assertion is sharp, helpers are fine.
- Performance — if the test is slow but correct, that's a separate
  concern.

### Where tests go

- Plugin (`packages/<Pkg>/`): new test files in `src/tests/`,
  registered in `src/package-test.ts`. Run with
  `cd packages/<Pkg> && grok test --host localhost --no-retry --skip-build`
  (rebuild only if TS changed).
- JS API surface: an `ApiTests` test plus an `ApiSamples` script (see
  `.claude/rules/testing.md`).
- Library: tests under `libraries/<lib>/src/tests/`.

If the location is unclear, query the KG: `MATCH (f:Feature {id:'<id>'})-[:IS_TESTED_IN]->(file:File) RETURN file.relative_path`
shows where existing tests of that feature live; new tests should go
alongside.

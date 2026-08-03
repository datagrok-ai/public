# sci-comp upstream contribution: `LambdaZStrategy.minSpanRatio`

> **Status:** plan / PR-ready spec.
> **Linked ROADMAP item:** Area 2 — "Improvement: λz min span ratio knob" (Ed observation #5).
> **Scope:** upstream PR to `@datagrok-libraries/sci-comp`, then a small local UI wire-up in nca-studio once the new sci-comp version is released.

---

## Why this is upstream, not local

### The user gap (Ed observation #5)

> *"PKNCA exposes `min.span.ratio` as a configurable parameter (`(t_last - t_first) / t_half` must exceed this for the fit to be valid; default 2). NCA Studio's λz panel has min points, min adj-R², and tie-break factor — but not span ratio. […] it's the one λz knob a Phoenix user looks for and doesn't find."*

The span-ratio check is a **fit-acceptance criterion**: a candidate terminal-phase window is rejected when its time span is short relative to the half-life it implies. It's how PKNCA prevents extrapolating AUC_∞ / t½ / CL from a too-short tail.

### Why it CANNOT be a local nca-studio fix

Per CLAUDE.md rule 5 ("never re-derive sci-comp math") and the Architecture Gotchas:

> *"If a numerical operation feels load-bearing for analytical correctness — it belongs in sci-comp. Open an issue against sci-comp, contribute the function, depend on the next sci-comp release. **Do not vendor.**"*

`lambdaZBestFit`'s acceptance policy is exactly that — load-bearing analytical correctness. The patterns we MUST NOT use locally:

- **Post-filter pattern** (compute span ratio after `computeNca` returns and downgrade `status: 'ok' → 'partial'` if it fails) → vendors sci-comp's policy outside sci-comp; a future release would surprise downstream behavior; obvious backdoor that someone has to rip out later.
- **Re-implement the auto-fit loop locally with the gate added** → straight rule-5 violation.

Both shortcuts dig a hole. The clean path is one PR upstream, one dependency bump downstream.

### What the user gets

- The last λz knob a Phoenix/PKNCA user looks for and currently can't find — closes a real UX gap from Ed's validation pass.
- Explicit defensibility control over auto-fit acceptance:
  - **TK datasets** like Ed's 06 (`t½ ≫ 24h sampling window`, all 16 cells currently `partial`): `partial` becomes a documented, configurable consequence of "the design didn't provide enough span for the half-life," not "the engine quietly decided."
  - **Discovery PK** screens at extreme dose levels where the terminal phase is squeezed: scientists can demand a defensible minimum span before AUC_∞ is reported.
- **PKNCA-equivalent default semantics** — the project's reference suite already configures `min_span_ratio: 2` in fixtures; this PR makes that knob actually available to consumers.

---

## Plan

### 1. Interface change

**File:** `libraries/sci-comp/src/nca/core/types.ts` (extend `LambdaZStrategy` at [`:115`](../../../package/node_modules/@datagrok-libraries/sci-comp/src/nca/core/types.d.ts)).

Add an optional field after `adjRSquaredFactor` (`:130`):

```ts
  /**
   * Minimum span ratio `(t_end - t_start) / halfLife` required for a
   * candidate window to be accepted. PKNCA `pk.calc.lambda.z` enforces
   * this via `min.span.ratio` (default 2). When undefined, the gate is
   * disabled — preserves pre-existing behavior. When set, candidate
   * windows whose terminal-phase span is short relative to the implied
   * half-life are rejected; `lambdaZBestFit` returns null if no
   * remaining candidate satisfies all gates. PKNCA-equivalent = 2.
   */
  readonly minSpanRatio?: number;
```

Optional + default-off = no behavior change for existing dependents on dependency upgrade.

### 2. Logic change

**File:** `libraries/sci-comp/src/nca/core/lambda-z.ts`, inside the candidate loop of `lambdaZBestFit` (around `:50`, immediately after the existing `if (fit.adjRSquared < options.minRSquared) continue;` line):

```ts
if (fit.adjRSquared < options.minRSquared) continue;
// Span-ratio gate — PKNCA pk.calc.lambda.z min.span.ratio. Skip when
// lambdaZ is non-positive (numerical guard; the linear regression
// already returned a valid fit, so this branch is defensive).
if (options.minSpanRatio !== undefined && fit.lambdaZ > 0) {
  const halfLife = Math.LN2 / fit.lambdaZ;
  const span = fit.tEnd - fit.tStart;
  if (span / halfLife < options.minSpanRatio) continue;
}
```

Placement after the adj-R² filter keeps the read order natural: "passes adj-R² → check span."

### 3. Tests

**File:** `libraries/sci-comp/src/nca/core/__tests__/lambda-z.test.ts` — four new tests (mirror the existing `minPoints` / `minRSquared` test shape):

1. **Gate disabled (undefined)** — back-compat: a fit that passes every other criterion is accepted regardless of span ratio.
2. **Gate passes** — `minSpanRatio: 2`, synthetic profile with `span / halfLife ≈ 4` → fit accepted.
3. **Gate fails → fallback to wider window** — `minSpanRatio: 2`, the *narrowest* candidate fails (`span / halfLife = 1`) but a wider candidate passes → `lambdaZBestFit` returns the wider fit.
4. **Gate fails → null return** — `minSpanRatio: 5`, no window in the eligible set passes → `lambdaZBestFit` returns `null` (downstream marks the profile `status: 'partial'`).

**File:** `libraries/sci-comp/src/nca/core/__tests__/reference-suite.test.ts` — strengthen: configure the existing theoph / indometh / rat_simple references with `minSpanRatio: 2` (the value the fixtures' own JSON metadata declares — `01_theoph.json:9` etc.). Pins the "fixtures pass with the gate" invariant — closes a silent equivalence (currently the fixtures match references *because* their auto-fits happen to satisfy `span/halfLife ≥ 2` anyway).

### 4. Documentation (sci-comp `.claude/rules/new-method-checklist.md`)

- **TSDoc** — included on the new field above.
- **`src/nca/README.md`** — λz auto-fit section: add a paragraph on the span-ratio gate, PKNCA-equivalent semantics, default-off rationale.
- **Top-level `README.md`** — sanity check: if it summarises λz strategy knobs, add `minSpanRatio` to the list.
- **`src/nca/CLAUDE.md`** — update the architecture tree if it enumerates `LambdaZStrategy` fields.
- **`docs/nca_core_interface_v2_1.md`** (if exists in sci-comp) — update the interface contract.

### 5. Verification

From `libraries/sci-comp/`:

```bash
npm run lint-fix && npm run build && npm test
```

All three must be green. Specifically:
- Existing tests pass (back-compat: optional field, default-off behavior).
- 4 new λz tests pass.
- Reference-suite tests pass with the strengthened `minSpanRatio: 2` config — no drift on theoph / indometh / rat_simple.

### 6. PR + release flow

1. Branch + commit in `Datagrok/public` (`libraries/sci-comp/`), referencing Ed obs #5 and this plan.
2. Open PR; await sci-comp maintainer review.
3. Release as next minor (`@datagrok-libraries/sci-comp@^0.7.0`-ish, depending on current version).
4. In nca-studio: `npm install @datagrok-libraries/sci-comp@^<new version>`; bump in `package/package.json`; verify build/lint/test stay green.

### 7. Local follow-up in nca-studio (after sci-comp release)

- **`package/src/adapter/defaults.ts`** — decide: default `minSpanRatio` in `DEFAULT_RULES.lambdaZ` = `undefined` (no gate) or `2` (PKNCA-strict). This is a product decision — the safer default is `undefined` (matches today's behavior); the more-PKNCA-faithful default is `2`. Recommend defaulting to `undefined` and letting the UI input drive it.
- **`package/src/views/nca-method-settings.ts`** — new `ui.input.float` "Min span ratio" in the λz settings section, tooltip explaining PKNCA-equivalent semantics, `0` / empty = disabled.
- **`package/src/views/nca-mapping-pane.ts` / `buildRules()`** — read the input and set `NcaRules.lambdaZ.minSpanRatio` accordingly.
- **`package/src/tests/controller-pipeline-tests.ts`** — sanity test: synthetic short-span profile with `minSpanRatio: 2` set → λz auto-fit rejected, profile `partial`.
- **`docs/_internal/ROADMAP.md`** — Area 2 entry "λz min span ratio knob" → `~~strikethrough~~ ✅ shipped` once UI lands.

### 8. Risks / open decisions

- **Default value (`undefined` vs `2`).** Recommendation: `undefined` (sci-comp side) + `undefined` initial (nca-studio side, with the UI making the user choose). Reasoning: prevents a silent behavior change on dependency upgrade and lets product / scientist decide whether the strict-by-default Phoenix posture is right for nca-studio.
- **Numerical edge cases.** Patched: `lambdaZ <= 0` skipped (`&& fit.lambdaZ > 0` guard); `halfLife` is naturally `+Infinity` for `lambdaZ → 0+`, which would make `span / halfLife → 0`, which fails any positive `minSpanRatio` → candidate rejected, which is the intended conservative behavior.
- **Surface area on `LambdaZResult`.** Optional follow-up (not required for this PR): expose `spanRatio` on `LambdaZResult` provenance, so nca-studio could surface it in the audit log / results grid. Out of v1 scope for the PR — file as a sci-comp follow-up issue if asked.

---

## Timeline estimate

- **Day 1 (sci-comp PR work):** ~3–4 hours — interface + logic patch + 4 tests + docs + `npm test/lint/build` green.
- **Day 2–N (review/release):** out of our hands — sci-comp maintainer cadence.
- **Post-release (nca-studio wire-up):** ~1 hour — UI input + buildRules thread + sanity test + ROADMAP closeout.

# 09 — Evaluation and Rollout

> Status: compact. How we test the harness itself, and how we ship it
> without breaking existing tests.

## Evaluating the system

The library and the CLI both have to pass a regression eval before
changes to weights, thresholds, prompts, or models ship.

### Eval fixtures

Pairs of (page-before, page-after, expected-element). Initial target:
30–50 fixtures covering:

- Pure renames (`data-testid` updated)
- Class hash churn (CSS-in-JS deploys)
- Structural wrappers added or removed
- Accessible name changes (synonyms, normalizations)
- Datagrok widget reorganizations (toolbar, ribbon, dialog)
- Negative cases: element removed, role changed, semantic mismatch
- Adversarial: two near-identical candidates

Fixtures live in `libraries/self-healing-locators/eval/fixtures/`. Each
fixture is a small directory: `before.html`, `after.html`,
`fingerprint.json`, `expected.json` (with one of: `winner_selector`,
`unable_to_heal`, `flag_for_review`).

### Eval metrics

Per run of the eval:

- **True heal rate** at HIGH band (target: ≥ 95% on positive fixtures)
- **False heal rate** at HIGH band (target: < 1% on negative fixtures)
- **Flag rate** at MID band (target: positive fixtures land here only
  when intentionally ambiguous)
- **Correct fail rate** on negative fixtures (target: ≥ 99%)
- **LLM cost per case** (offline only; tracked, not gated)
- **Latency p50/p95** for runtime resolution

Eval runs both the runtime resolver alone (no LLM) and the offline
pipeline (with LLM). Both have separate gates.

### When eval runs

- On every PR that touches the library, the CLI, prompts, weights, or
  thresholds.
- Nightly against `master` to catch regressions from dependency drift.
- Manually before any model rollout (e.g. trying `claude-haiku-4-5` →
  `claude-haiku-5` whenever that lands).

## Rolling out the system

We ship in stages. Each stage has a clear gate.

### Stage 0 — Land the design (this PR)

Get reviewer alignment on docs, schema, prompt, policy. No code yet.

### Stage 1 — Library skeleton

Types, schemas, registry I/O, capture API. No resolver chain, no scorer.
Tests pass. Library publishes but no consumer uses it.

### Stage 2 — Capture and fingerprint registry

Passive capture in green runs, explicit `healing.anchor()`. One
volunteer package opts in; we collect a few weeks of fingerprints in a
real environment to validate the schema. **No healing happens yet** —
even when the primary selector misses, the resolver throws.

### Stage 3 — Runtime resolver (deterministic only)

Tiers 1–3 (visual stays off by default). Scorer. Decision gate. Audit
log. Healing queue writer. The volunteer package starts seeing heals.
We watch the audit log for false heals. No source mutation.

### Stage 4 — Offline CLI without source mutation

Triage, prompt, validator, report. The CLI runs in dry-run mode: it
produces a healing report but does not open PRs. We review reports
manually for a few iterations and tune.

### Stage 5 — Offline CLI opens PRs

Codemod and PR opener enabled. PRs require human approval (no
auto-merge ever). We start with one package, expand as we gain
confidence.

### Stage 6 — Visual layer (optional)

Visual hashing as Tier 4 enabled where it adds signal. Per-package opt-in.

### Stage 7 — Broad rollout

All packages opt-in. Documentation updated. Training / runbook for
package maintainers.

Each stage is gated by:

- Eval suite green.
- No regression in volunteer package's primary CI signal.
- Sign-off from QA and from the package's CODEOWNERS.

## What success looks like, six months in

- Median package test run produces zero MID-band flags after a typical
  platform release.
- Healing PRs land within 24 hours of a release that produces drift,
  with > 90% of cases auto-resolved at HIGH and the rest cleanly
  flagged.
- Eval suite catches at least one regression before merge, demonstrating
  it earns its keep.
- LLM cost per release is bounded by the configured cap and well below
  it on average.

## What failure looks like, and what we'd do

- Reviewers ignoring healing PRs because they're noisy → reduce the
  volume by raising HIGH threshold and pushing more cases to MID-with-flag.
- Specific test or package contributing disproportionate cases →
  invest in better explicit anchors there, or recommend platform-side
  test ID work.
- LLM cost spiking → tighten model selection rules, increase caching
  TTL, batch more aggressively, or move to a smaller model.
- False heals causing prod escapes → tighten the policy and the
  scoring; consider raising default thresholds; add the failure to the
  eval suite as a permanent fixture.

The system is meant to be tunable. Failures get fed back into
calibration, not papered over.

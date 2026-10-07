# 01 — Design Overview

> Status: full detail. This is one of the documents to read carefully.

## 1. Problem statement

Datagrok UI tests — both the package tests on the Puppeteer-based framework
and the growing Playwright E2E suite — fail when the DOM, CSS classes, ARIA
labels, or text content shift between releases. The failures are not real
regressions; they are locator drift. Today the cost of locator drift is paid
by humans who hand-edit tests after each platform or package release.

The goal of this proposal is to absorb routine drift automatically while
preserving the ability of tests to fail loudly when something genuinely
broke in the product.

## 2. Design principles

These are the non-negotiables. Every detail downstream serves them.

**P1. The LLM is the last resort, not the first move.** Most drift is
trivial: a class renamed, a `data-testid` updated, a wrapper added. A
deterministic resolver that walks a priority hierarchy resolves the vast
majority of cases at zero cost and with no nondeterminism. Only what cannot
be resolved deterministically gets escalated to Claude, and only offline.

**P2. Real bugs must surface.** Self-healing must never mask a regression.
If a button is gone or its semantics changed, the test fails. The
confidence model is the gatekeeper here — it is computed from objective
matches against a fingerprint captured when the test last passed, not from
an LLM's self-reported confidence.

**P3. Every heal is auditable.** Every decision — resolved by which tier,
what matched, what didn't, what the resulting confidence was, who or what
approved it — is written to a structured log. Without this, the system
becomes a black box that quietly hides failures.

**P4. Source code mutation is gated by a human.** The offline tool produces
a PR. It does not commit directly to a branch that anything depends on. A
human reviews the diff and merges.

**P5. Runtime is fast and offline-capable.** No network calls, no LLM, no
external services in the runtime path. Everything required at runtime ships
with the library and the fingerprint registry.

**P6. The system is calibratable, not hardcoded.** Weights, thresholds,
model selection rules, and prompt templates are configuration, not code.
Teams can tune them without forking the library.

## 3. High-level architecture

```
                     ┌──────────────────────────────┐
                     │  Test Authoring / Recording  │
                     │                              │
                     │  Fingerprints captured       │
                     │  during green runs and via   │
                     │  explicit healing.anchor()   │
                     └──────────────┬───────────────┘
                                    │
                                    ▼
                     ┌──────────────────────────────┐
                     │   Fingerprint Registry       │
                     │   (versioned JSON files)     │
                     └──────────────┬───────────────┘
                                    │
                ┌───────────────────┴────────────────────┐
                │                                        │
                ▼                                        ▼
   ┌───────────────────────────┐         ┌──────────────────────────────┐
   │  RUNTIME (library)        │         │  OFFLINE (CLI tool)          │
   │                           │         │                              │
   │  ┌─────────────────────┐  │         │  ┌────────────────────────┐  │
   │  │  Locator Resolver   │  │         │  │  Queue Reader          │  │
   │  │  Tier 1: semantic   │  │         │  │  (low-conf + failed    │  │
   │  │  Tier 2: structural │  │         │  │   runtime cases)       │  │
   │  │  Tier 3: text       │  │         │  └──────────┬─────────────┘  │
   │  │  Tier 4: visual     │  │         │             ▼                │
   │  └──────────┬──────────┘  │         │  ┌────────────────────────┐  │
   │             ▼             │         │  │  Prompt Builder        │  │
   │  ┌─────────────────────┐  │         │  │  (fingerprint + a11y   │  │
   │  │  Confidence Scorer  │  │         │  │   tree + optional img) │  │
   │  └──────────┬──────────┘  │         │  └──────────┬─────────────┘  │
   │             ▼             │         │             ▼                │
   │  ┌─────────────────────┐  │         │  ┌────────────────────────┐  │
   │  │  Decision Gate      │  │         │  │  Claude Client         │  │
   │  │  HIGH → heal        │  │         │  │  (Haiku → Sonnet →     │  │
   │  │  MID  → heal+flag   │  │         │  │   Opus, tool_use)      │  │
   │  │  LOW  → fail+queue  │  │         │  └──────────┬─────────────┘  │
   │  └──────────┬──────────┘  │         │             ▼                │
   │             ▼             │         │  ┌────────────────────────┐  │
   │  ┌─────────────────────┐  │         │  │  Candidate Validator   │  │
   │  │  Audit Log Writer   │  │         │  │  (headless browser,    │  │
   │  └─────────────────────┘  │         │  │   re-scores via the    │  │
   │             │             │         │  │   same scorer)         │  │
   └─────────────┼─────────────┘         │  └──────────┬─────────────┘  │
                 │                       │             ▼                │
                 │                       │  ┌────────────────────────┐  │
                 │                       │  │  Codemod + PR Opener   │  │
                 │                       │  └────────────────────────┘  │
                 │                       └──────────────┬───────────────┘
                 │                                      │
                 ▼                                      ▼
   ┌──────────────────────────────────────────────────────────────────┐
   │  Audit log (locator-events.jsonl)                                │
   │  Healing queue (healing-queue.jsonl)                             │
   │  Review dashboard (out of scope for this PR)                     │
   └──────────────────────────────────────────────────────────────────┘
```

The two columns share three things:
1. **The fingerprint schema** (`schemas/fingerprint.schema.json`).
2. **The confidence scorer** — same code, same weights. The offline tool
   re-scores Claude's proposals using the runtime scorer; Claude does not
   get to vote on its own confidence.
3. **The audit log format** (`schemas/healing-event.schema.json`).

## 4. Component boundaries

### 4.1 `libraries/self-healing-locators` (runtime library)

A `@datagrok-libraries/*` package consumed by:
- Datagrok package tests (Puppeteer-based)
- Standalone Playwright projects (E2E suite)

Exports:
- `healing.anchor(name, element)` — explicit fingerprint capture
- `healing.resolve(selector, options)` — drop-in replacement for
  `page.$(selector)` / `page.locator(selector)` that triggers the resolver
  on miss
- `healing.configure(opts)` — thresholds, registry path, audit log path,
  feature flags
- `healing.audit` — read-side API for tools and dashboards

No dependency on the Anthropic SDK. No network code at all.

See [`api/healing-api.md`](../api/healing-api.md) for the full surface.

### 4.2 `tools/self-healing-cli` (offline CLI)

A `tools/*` package, in line with `datagrok-tools`. Distributed as an npm
binary. Depends on:
- `@datagrok-libraries/self-healing-locators` (for the shared scorer,
  schemas, and registry I/O)
- `@anthropic-ai/sdk` (Claude API client)
- `puppeteer` or `playwright` for candidate validation

Commands (sketch — finalized in `tools/self-healing-cli/docs/cli-overview.md`):

```
grok-heal queue              # show pending cases
grok-heal run [--model auto] # process the queue
grok-heal validate <case-id> # dry-run a single case
grok-heal pr                 # open a PR with codemod'd changes
grok-heal eval               # run the regression eval suite
```

The CLI is invoked from CI on a schedule (e.g. nightly) or manually after
a release that produced runtime-queued cases.

### 4.3 Fingerprint registry

A versioned directory (location decided in this PR — see
`02-fingerprint-spec.md` § "Storage") holding one JSON file per logical
test or per package, depending on the chosen layout. Committed to git.
Updates flow through normal review; the offline CLI does not bypass review.

## 5. Data flow scenarios

### Scenario A — Trivial drift, runtime resolves it

1. Test calls `healing.resolve('[data-testid="submit"]')`.
2. Selector misses (the platform renamed the testid to `submit-btn`).
3. Resolver loads the fingerprint for that anchor.
4. Tier 1 finds an element matching the recorded `ariaLabel + ariaRole`
   pair. The fallback `testId` candidate `submit-btn` also matches the
   same element.
5. Scorer computes confidence ≥ HIGH.
6. Decision gate heals silently. Test continues.
7. Audit log records: `tier=1, confidence=0.92, matched=[ariaLabel,
   ariaRole, testId-fuzzy], healed_selector='[data-testid="submit-btn"]'`.

No LLM, no human involvement. Cost: a few milliseconds and one log line.

### Scenario B — Structural change, runtime escalates

1. The submit button is now nested inside a new `<form-actions>` wrapper
   and its `data-testid` is gone.
2. Tier 1 misses. Tier 2 finds an element matching tag, stable classes,
   and one of three parent-chain levels — but the chain is shifted by one.
3. Tier 3 matches on visible text "Submit". Tier 4 finds three elements
   with similar bounding boxes; the right one is among them.
4. Scorer computes confidence in the MID band: 0.62.
5. Decision gate heals (test continues) but flags the case.
6. Audit log records `flag=review`. The case is also written to
   `healing-queue.jsonl` for offline analysis.
7. Next nightly CLI run: Claude sees the fingerprint, the new a11y tree,
   and proposes 3 candidate selectors. The validator confirms 1 of them
   matches a single element with confidence ≥ HIGH under the runtime
   scorer.
8. CLI opens a PR updating the test source to use the new selector. A
   human reviews and merges.

### Scenario C — Real bug, system fails loudly

1. The submit button was removed from the page entirely (refactor mistake).
2. All four tiers miss or return only weak partial matches.
3. Confidence < LOW threshold.
4. Decision gate fails the test. The case is queued.
5. CLI runs, Claude proposes selectors that the validator finds either
   missing or pointing to different semantic elements (e.g. a "Cancel"
   button in the same position).
6. CLI marks the case `unable_to_heal` with reason
   `no_semantic_match_in_dom`. No PR opened. The case is escalated to a
   human via the review dashboard.

This scenario is the entire point of P2. The system must produce this
outcome reliably.

## 6. What we deliberately keep out

- **In-test LLM calls.** Latency, cost, nondeterminism. Forbidden.
- **Auto-merge of healing PRs.** Even high-confidence offline heals get a
  human reviewer, at least until the eval suite has earned trust.
- **A separate ML model.** A scoring formula with calibrated weights is
  enough for the deterministic tiers. If we later need a learned model,
  we can swap the scorer; the architecture supports it.
- **Healing of non-locator failures.** Timing flakes, network errors,
  state-leak failures — different problems, different solutions. Mixing
  them with locator healing dilutes both.
- **A custom format for fingerprints.** Plain JSON, schema-validated. No
  protobuf, no SQLite, no embeddings (yet).

## 7. What is calibratable vs hardcoded

| Calibratable (config) | Hardcoded (code) |
|---|---|
| Tier weights | Tier order |
| Confidence thresholds (HIGH/MID/LOW) | Scoring formula shape |
| Model selection rules | tool_use response schema |
| Prompt templates | Audit log schema |
| Feature flags (screenshots, etc.) | Public API surface |
| Per-package overrides | Registry layout |

Anything calibratable lives in `self-healing.config.json` next to the
registry. Defaults ship with the library.

## 8. Dependencies on the platform

We rely on three platform capabilities. Two exist; one is a recommendation.

**Exists:** `Widget.getWidgetStatus()` returns a runtime structure intended
specifically for automated testing and introspection (per the JS API docs).
We use this in Tier 1 as a Datagrok-specific semantic anchor.

**Exists:** the Inspector tool (`Alt + I`) exposes the registered widget
tree. The fingerprint capture utility taps into the same data source so
fingerprints reflect the platform's own model of the page.

**Recommended (not blocking):** a convention for `data-testid` attributes
on common UI primitives (`ui.button`, `ui.input.*`, dialog headers,
ribbon items). Without this, semantic Tier 1 match is weaker and we lean
harder on text and structure. See
[`06-datagrok-integration.md`](./06-datagrok-integration.md) for the
specific suggestions.

## 9. Failure modes we explicitly handle

| Failure mode | Handling |
|---|---|
| Registry file missing | Treat as fresh capture; record for next green run |
| Registry file stale (schema version mismatch) | Migration script; fail closed if migration is ambiguous |
| All tiers miss with non-zero partial scores | Confidence < LOW → fail + queue |
| Multiple candidates tied at HIGH | Fail closed, queue for offline disambiguation. Never guess. |
| LLM API down | Offline CLI exits cleanly, queue persists |
| LLM proposes selector that matches multiple elements | Validator rejects, marks as `ambiguous_match` |
| LLM proposes selector that matches a different semantic element | Validator computes runtime score; if < HIGH, reject |
| Codemod conflicts with concurrent edits | Skip case, requeue, notify in PR description |

## 10. Why this shape

We considered three alternatives and rejected them:

**A. Pure runtime LLM healing.** Rejected: violates P1, P5, and creates a
hard dependency on network availability for tests. Also: the audit story is
weak because LLM reasoning isn't reliably reproducible.

**B. Pure offline LLM rewriting (no runtime resolver).** Rejected: every
trivial drift becomes a Claude call and a PR. Noise overwhelms signal,
costs balloon, and the team starts ignoring PRs.

**C. Visual-only matching (image diff + perceptual hash).** Rejected as a
primary mechanism: too brittle for SPAs that re-render frequently, fails
on dynamic content, and offers no semantic explanation when it errs.
Visual signals stay as Tier 4 for cases where they genuinely help.

The hybrid two-phase design absorbs trivial drift cheaply at runtime and
escalates only what is genuinely hard. The LLM does what it's good at
(reasoning over messy structure) without doing what it's bad at (being
the source of truth for production-critical decisions).

# 05 — Offline Flow

> Status: full detail. This is the LLM-facing half of the system.

## 1. When the offline tool runs

Three triggers:

- **Scheduled.** Nightly CI job: `grok-heal run`. Processes the merged
  healing queue from the day's runs.
- **On-demand after a release.** A platform release touching UI is
  expected to produce queued cases; QA runs the tool to surface them
  early.
- **Manually for a single case.** During investigation:
  `grok-heal validate <case-id>` reproduces one case end-to-end without
  opening a PR.

The CLI lives in `tools/self-healing-cli`. It depends on the runtime
library for the schemas and the scorer.

## 2. Inputs and outputs

**Inputs:**

- The merged healing queue (`healing-queue.jsonl`).
- DOM snapshots, accessibility tree snapshots, and optionally screenshots
  associated with each entry.
- The fingerprint registry (read-only).
- A configuration file (model selection rules, thresholds, cost
  ceilings, prompt template paths).

**Outputs:**

- A processing report (`healing-report-<run-id>.json`) — what was
  attempted, what succeeded, what failed, with cost breakdown.
- Per-case results in `healing-results/<case-id>.json` containing the
  Claude proposals, the validator's per-candidate scores, and the
  chosen winner if any.
- One PR per package (or one combined PR, configurable) with codemod'd
  test source updates.

The CLI never writes to the fingerprint registry directly — registry
updates flow through the runtime, on green runs, after the PR merges.

## 3. End-to-end pipeline

```
healing-queue.jsonl
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 1. Queue Reader   │ -> │ Group entries by (testFile, anchor)  │
│                   │    │ Dedupe; pick the freshest snapshot.  │
└───────┬───────────┘    └──────────────────────────────────────┘
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 2. Triage         │ -> │ For each case, decide:               │
│                   │    │  - retry deterministic resolver on   │
│                   │    │    the snapshot (registry may have   │
│                   │    │    been updated since the failure)   │
│                   │    │  - if it now resolves at HIGH, mark  │
│                   │    │    as 'self-resolved', skip LLM      │
│                   │    │  - else proceed                      │
└───────┬───────────┘    └──────────────────────────────────────┘
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 3. Model Selector │ -> │ Pick Haiku / Sonnet / Opus per the   │
│                   │    │ rules in § 5.                        │
└───────┬───────────┘    └──────────────────────────────────────┘
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 4. Prompt Builder │ -> │ Render system + user template with   │
│                   │    │ fingerprint, a11y tree, optional     │
│                   │    │ screenshot, action hint.             │
└───────┬───────────┘    └──────────────────────────────────────┘
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 5. Claude Client  │ -> │ Call /v1/messages with tool_use       │
│                   │    │ forcing the response schema. Retry   │
│                   │    │ on schema violation (max 2).         │
└───────┬───────────┘    └──────────────────────────────────────┘
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 6. Validator      │ -> │ Spin up headless browser. Load page  │
│                   │    │ at the failing point (via test       │
│                   │    │ replay or DOM snapshot serving).     │
│                   │    │ For each candidate, find element,    │
│                   │    │ score with the runtime scorer,       │
│                   │    │ apply hard gates.                    │
└───────┬───────────┘    └──────────────────────────────────────┘
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 7. Decision       │ -> │ best_score ≥ HIGH    → propose heal  │
│                   │    │ MID ≤ best_score     → propose heal  │
│                   │    │                       w/ review flag │
│                   │    │ best_score < MID     → mark unhealed │
│                   │    │ all candidates fail  → mark unhealed │
└───────┬───────────┘    └──────────────────────────────────────┘
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 8. Codemod        │ -> │ Locate the call site in the test     │
│                   │    │ source. Replace selector. If the     │
│                   │    │ test uses healing.anchor(), update   │
│                   │    │ the selector argument only.          │
└───────┬───────────┘    └──────────────────────────────────────┘
        │
        ▼
┌───────────────────┐    ┌──────────────────────────────────────┐
│ 9. PR Opener      │ -> │ One PR per package. Body includes    │
│                   │    │ per-case scores, audit refs, cost,   │
│                   │    │ unresolved cases.                    │
└───────────────────┘    └──────────────────────────────────────┘
```

## 4. Reproducing the failing page state

Step 6 (validator) requires the page in roughly the state it was in when
the test failed. Three approaches, in priority order:

1. **Replay the test up to the failing step.** Most reliable. The test
   harness exposes a `replay` mode that runs the original test and stops
   just before `healing.resolve()` is called for the queued case. The
   validator then evaluates each candidate against the live page.
2. **Serve the DOM snapshot.** Faster, less reliable. The CLI starts a
   tiny static server with the captured `domSnapshotPath`. Limitations:
   no JS state, no event handlers — fine for static element lookup,
   wrong for interactive validation.
3. **Use the snapshot for selector lookup only**, with a notice that
   interactive features were not exercised.

The CLI prefers (1) and falls back to (2) with a warning. (3) is for
cases where the source test cannot be replayed (e.g. external
preconditions are no longer satisfiable).

## 5. Model selection

Tiered by case complexity, biased toward the cheapest model that
produces accurate-enough output. The selector is rule-based — no LLM
deciding which LLM to call.

**Default rule:**

```
Estimated input size (a11y tree nodes + parent chain depth):
   < 50 nodes,  no structural shift          → claude-haiku-4-5
   < 200 nodes, simple drift signals         → claude-haiku-4-5
   200–800 nodes, structural shift           → claude-sonnet-4-6
   > 800 nodes OR retry after Haiku failure  → claude-sonnet-4-6
   second retry / explicit escalation flag   → claude-opus-4-7
```

"Drift signals" come from comparing fingerprint to live a11y tree:

- testId/domId disappeared → simple
- ariaRole or accessible name changed → moderate
- parent chain depth changed by ≥ 2 → structural
- multiple semantic features changed at once → structural

A run starts with Haiku for each case unless the case is pre-classified
as structural. If Haiku produces no candidate that survives validation,
the case auto-retries on Sonnet. Opus is reserved for explicit
escalation (`grok-heal run --case <id> --escalate`).

The model used is recorded per case in the report.

## 6. Prompt construction

System prompt: `prompts/system.md`. User message template:
`prompts/user-template.md`. Few-shot examples in
`prompts/examples/`.

The user message is composed from:

1. **Fingerprint** — the JSON object, but with `contextSnippet` removed
   to avoid duplicating data we provide elsewhere.
2. **Failed selector** — the original string and the action hint.
3. **Current accessibility tree** — pruned to a depth window centered on
   the most plausible region (parent chain match → expand 2 levels
   around the best candidate; otherwise the current view's a11y root).
4. **Optional HTML snippet** — included only when the a11y tree is
   insufficient (e.g. the element role is `generic`).
5. **Optional screenshot** — only when the feature flag is on, and only
   when the case classification is "structural" or worse.
6. **Negative constraints** — explicit list of selectors the resolver
   already tried that did not work, so the LLM does not propose them.

Token budget targets:

- Haiku calls: prompt ≤ 8k tokens
- Sonnet calls: prompt ≤ 32k tokens
- Opus calls: prompt ≤ 64k tokens

If a case exceeds the budget, the a11y tree is pruned harder (depth and
breadth). We never silently drop the fingerprint or the negative
constraints.

## 7. Tool use schema

We use Claude's `tool_use` to force structured output. The tool is
declared with strict JSON Schema (in `schemas/tool-response.schema.json`).
A summary:

```json
{
  "name": "propose_locators",
  "description": "Propose CSS/XPath selectors for the missing element.",
  "input_schema": {
    "type": "object",
    "required": ["candidates", "unable_to_match"],
    "properties": {
      "candidates": {
        "type": "array",
        "minItems": 0,
        "maxItems": 5,
        "items": {
          "type": "object",
          "required": ["selector", "selector_type", "matched_attributes", "rationale"],
          "properties": {
            "selector": { "type": "string" },
            "selector_type": {
              "type": "string",
              "enum": ["css", "xpath", "role", "testid", "text"]
            },
            "matched_attributes": {
              "type": "array",
              "items": { "type": "string" },
              "description": "Subset of fingerprint feature names that this selector relies on. Must come from the documented feature vocabulary."
            },
            "rationale": { "type": "string", "maxLength": 400 },
            "risk_notes": { "type": "string", "maxLength": 400 }
          }
        }
      },
      "unable_to_match": { "type": "boolean" },
      "unable_reason": {
        "type": "string",
        "enum": [
          "no_semantic_match_in_dom",
          "ambiguous_candidates",
          "page_state_mismatch",
          "insufficient_context",
          "other"
        ]
      },
      "unable_reason_detail": { "type": "string", "maxLength": 400 }
    }
  }
}
```

Note what is **not** in the schema:

- No `confidence` field. The runtime scorer assigns confidence; the LLM
  doesn't get to vote.
- No free-form output. If `unable_to_match` is `true`, `candidates`
  must be empty; the validator enforces this.
- No instruction to "explain in detail." `rationale` is capped because
  long rationales waste tokens and tend to drift toward unfounded
  certainty.

## 8. Validation

For each returned candidate:

1. Find element(s) on the live page using `selector` and `selector_type`.
2. If zero elements: reject (`reason: not_found`).
3. If multiple elements: reject (`reason: ambiguous`).
4. If exactly one element: build a candidate fingerprint from it, run the
   runtime scorer against the stored fingerprint.
5. Apply hard gates (semantic mismatch, action mismatch, dgWidgetType
   mismatch).
6. Final score → band → decision.

The winning candidate is the highest-scoring one that passes all gates.
Ties at HIGH within 0.10 → demote both (per § 3 of the confidence model).

## 9. Codemod and PR

The codemod is conservative:

- Locates the call site by `testRef.file` + `anchorName` + line range
  recorded at capture.
- Replaces only the selector argument. Never reformats surrounding code.
- If the test passes a selector built dynamically (template literal,
  variable), the codemod skips and the case is reported as "needs human
  rewrite" with the proposed selector in the PR body.

The PR title and body convention:

```
title: chore(self-healing): heal locators in <package> (<N> cases)

body:
  ## Summary
  <N> selectors auto-healed by the offline maintenance tool.

  | Test | Anchor | Old | New | Score | Tier |
  |---|---|---|---|---|---|
  ...

  ## Unhealed cases
  ...

  ## Cost
  $X.XX (Haiku: A, Sonnet: B, Opus: C)

  ## Audit references
  - run-id: ...
  - artifacts: ...
```

Reviewer (default `CODEOWNERS` of the affected package) is automatically
requested.

## 10. Cost controls

- **Per-run cap.** Configured dollar limit. The tool tracks cumulative
  cost; on cap reached, it processes remaining cases as `deferred` and
  exits. The cap is enforced before each call, not after.
- **Per-case retry cap.** Max 2 retries per case (Haiku → Sonnet → stop).
  Opus only on manual escalation.
- **Caching.** Identical (fingerprint hash, page snapshot hash) inputs
  reuse a cached response. Cache TTL: 7 days. This matters for a
  flapping case that gets re-queued from multiple test runs against the
  same page state.
- **Batching.** Cases targeting the same page snapshot can share a
  single `messages` request — multiple anchors in one prompt — but only
  when the input remains within the model's prompt budget. The default
  is one anchor per request; batching is opt-in.

## 11. Failure modes specific to offline

| Situation | Behavior |
|---|---|
| Anthropic API returns 5xx | Exponential backoff up to 3 attempts, then mark case as `transient_error`, requeue. |
| `tool_use` response invalid against schema | One retry with the schema appended to the user message. Then mark `schema_violation`, do not heal. |
| All candidates fail validation | Mark `no_valid_candidate`. Record proposals for human review. |
| Validator cannot reach the page | Mark `replay_failed`. Human intervention. |
| Codemod cannot locate the call site | Mark `codemod_failed`, include proposed selector in the PR body for manual application. |
| Concurrent edits to the same test file | Skip case; requeue for the next run. |
| Cost cap reached mid-run | Process remaining as `deferred`. Surface in the report. |

Every failure mode produces an entry in the report so the next human
looking at the dashboard sees what happened, not just that something
"didn't heal."

## 12. What the offline tool does not do

- Does not commit directly to a long-lived branch. Always a PR.
- Does not modify the fingerprint registry. Fingerprints are
  capture-side, not heal-side.
- Does not run on every CI run. It's a maintenance tool, not a build
  step.
- Does not call Claude with raw test source code. The prompt sees
  fingerprints and DOM, not test logic.

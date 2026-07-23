# [Design] Self-Healing Locators for Datagrok UI Tests

> **Status:** Design proposal. No executable code in this PR.
> **Goal:** Get reviewer alignment on architecture, fingerprint schema, prompt
> contract, healing policy, and rollout plan **before** any code is written.

## TL;DR

UI tests in Datagrok (package tests on Puppeteer, plus a growing Playwright E2E
suite) break when the DOM, CSS, or labels shift. This proposal introduces a
two-phase **self-healing locator** system:

- **Runtime fallback** (library, no LLM): when a selector misses, a
  deterministic resolver walks a 4-tier priority hierarchy (semantic IDs →
  stable structural → text content → visual/positional) using a
  **fingerprint** captured during green test runs. A confidence score gates
  the outcome: high → heal silently, mid → heal + flag for review, low →
  fail the test.
- **Offline maintenance** (CLI, uses Claude API): low-confidence and
  unresolved cases land in a queue. A separate tool prompts Claude with the
  fingerprint + accessibility tree, validates each proposed selector against
  a live page in headless Puppeteer/Playwright, and opens a follow-up PR
  with codemod'd test source updates.

The two phases share one fingerprint schema, one confidence model, and one
audit log format. The LLM is the **last resort**, not the first move.

## What's in this PR (design only)

```
libraries/self-healing-locators/
  README.md
  docs/
    01-design.md                 [FULL]
    02-fingerprint-spec.md       [FULL]
    03-confidence-model.md       [FULL]
    04-runtime-flow.md           [FULL]
    05-offline-flow.md           [FULL]
    06-datagrok-integration.md   [COMPACT]
    07-policy.md                 [FULL]
    08-review-process.md         [COMPACT]
    09-eval-and-rollout.md       [COMPACT]
  prompts/
    system.md                    [FULL]
    user-template.md             [FULL]
    examples/
      01-testid-rename.md
      02-aria-label-change.md
      03-structural-refactor.md
  schemas/
    fingerprint.schema.json
    healing-event.schema.json
    tool-response.schema.json
  api/
    healing-api.md               [FULL]

tools/self-healing-cli/
  README.md
  docs/
    cli-overview.md              [COMPACT]
```

## Decisions baked into this proposal

| Decision | Choice | Where it's discussed |
|---|---|---|
| Project location | `libraries/self-healing-locators` + `tools/self-healing-cli` | `01-design.md` |
| Healing mode | Hybrid: runtime fallback + offline source updates | `04-runtime-flow.md`, `05-offline-flow.md` |
| LLM provider | Anthropic Claude API | `05-offline-flow.md` |
| Model strategy | Tiered: Haiku → Sonnet → Opus | `05-offline-flow.md` § "Model selection" |
| Fingerprint storage | Centralized JSON registry under `libraries/self-healing-locators/registry/` (not in this PR) | `02-fingerprint-spec.md` § "Storage" |
| Fingerprint capture | Passive (auto-update on green runs) + explicit (`healing.anchor()` API) | `api/healing-api.md` |
| Screenshots | Feature-flagged, off by default | `02-fingerprint-spec.md` § "Visual layer" |
| LLM in runtime | **No.** Runtime is deterministic and fast. | `04-runtime-flow.md` |
| Source code mutation | Offline only, via codemod, gated by PR review | `05-offline-flow.md` § "Codemod & PR" |

## Review focus

Please read in this order — each doc is short and self-contained:

1. **`01-design.md`** — confirm overall shape and component boundaries
2. **`07-policy.md`** — confirm the "what we heal vs what we let fail" boundary
3. **`02-fingerprint-spec.md`** — confirm the schema before we lock it in
4. **`prompts/system.md` + `prompts/user-template.md`** — confirm the contract with Claude
5. **`03-confidence-model.md`** — confirm the weights and thresholds (these are calibratable, not hardcoded forever)
6. The rest — orientation only

## Explicit non-goals

- We do **not** propose changes to the Datagrok platform itself in this PR.
  Any `getWidgetStatus()` extensions or new `data-testid` conventions are
  noted as recommendations in `06-datagrok-integration.md`, to be addressed
  separately.
- We do **not** heal anything beyond locators. Assertions, timing, navigation,
  and test logic are out of scope. See `07-policy.md` for the rationale.
- We do **not** mutate test source files automatically without a PR.

## Open questions (please flag in review)

1. **Registry location** — should fingerprints live next to tests, in a
   sibling directory, or in a single central registry? See
   `02-fingerprint-spec.md` § "Storage" for trade-offs.
2. **Failure budget** — what's the acceptable number of LLM-suggested heals
   per week before we treat it as a signal that the platform needs more
   stable test IDs?
3. **PR ownership** — when the offline tool opens a PR, who is the default
   reviewer? Test author? QA team? Code owners of the affected package?
4. **Cost ceiling** — concrete dollar limit per offline run before the tool
   pauses and asks for human approval.

## After this PR is approved

- Iterate on review comments until docs are signed off.
- Open implementation PRs in this order:
  1. `libraries/self-healing-locators` core types + fingerprint capture
  2. Runtime resolver + confidence scorer (no LLM)
  3. Audit log + healing queue writer
  4. `tools/self-healing-cli` skeleton + Claude client
  5. Prompt + tool_use response validation
  6. Candidate validator (headless browser)
  7. Codemod + PR opener
  8. Eval harness + first 30-50 fixtures

Each step is independently reviewable and shippable behind a feature flag.

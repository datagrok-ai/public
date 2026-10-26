# @datagrok-libraries/self-healing-locators

Resilient UI element discovery for Datagrok tests. When a primary selector
breaks because of a UI change, this library finds the same element through
alternative signals and reports its confidence in the match.

## Status

**Design phase.** This directory currently contains design docs, schemas, and
prompt specs only. No runtime code has been written yet. See
[`docs/01-design.md`](./docs/01-design.md) for the proposed architecture and
[`PR_DESCRIPTION.md`](../../PR_DESCRIPTION.md) at the repo root for the rollout
plan.

## What this library does (when implemented)

- Captures **fingerprints** of UI elements during green test runs — semantic
  IDs, ARIA roles, stable classes, parent chain, visible text, optional
  visual hash, and Datagrok widget metadata via `getWidgetStatus()`.
- At runtime, when a primary selector misses, walks a 4-tier priority
  hierarchy (semantic → structural → text → visual) to relocate the
  element, **without calling an LLM**.
- Computes a **confidence score** for the match. High confidence heals
  silently; mid-confidence heals but flags for review; low confidence fails
  the test and queues the case for offline analysis.
- Writes an **audit log** for every healing decision.

The companion CLI in [`tools/self-healing-cli`](../../tools/self-healing-cli)
handles offline maintenance: it consumes the queue, calls Claude to propose
new selectors for hard cases, validates them against a live page, and opens
a PR with updated test sources.

## What this library does **not** do

- Does not call any LLM at runtime. Tests stay fast and offline-capable.
- Does not heal assertions, timing, or test logic. Locators only.
- Does not mutate test source files. Source edits happen via codemod in the
  offline CLI, behind a PR.

## Quick links

- [Architecture overview](./docs/01-design.md)
- [Fingerprint schema](./docs/02-fingerprint-spec.md)
- [Confidence model](./docs/03-confidence-model.md)
- [Runtime flow](./docs/04-runtime-flow.md)
- [Offline flow](./docs/05-offline-flow.md)
- [Healing policy](./docs/07-policy.md)
- [Public API](./api/healing-api.md)

## License

Same as the parent repository.

# Artifact Alignment

Program-aligned catalog of published artifacts: publishing a Compute2 workflow run freezes an
immutable copy, aligns it with a drug program / study / molecules, and makes it discoverable
through one faceted, security-trimmed query.

- [Design: Artifact Alignment & Discoverability](docs/artifact-alignment.html) — the general
  design (business case, domain model, versioning and review semantics, platform asks). Not
  specific to this package; this package is its first implementation, scoped to Compute2
  workflow runs and single function runs (RFV model runs, individually saved steps).
- [Compute2 PoC implementation plan](docs/compute2-poc-plan.md) — how the general design maps
  onto this package, including the accepted gaps and core asks.

## What's here

- `databases/artifacts/schema.json` — the domain schema: program/compound/study registries,
  tags, the `alignment` live-versions table (business key `(publication_id, status)` enforces
  the at-most-one-approved + one-in-review invariant over live rows), and the append-only
  `alignment_history`.
- `src/service/` — the publish service (frozen-copy deep clone of a saved workflow run,
  latest-approved-stays-live versioning over the `(program, study, name)` republish key,
  curation copy-forward), approval/rejection with service-side reviewer checks, and the drift
  check.
- `src/domain/` — program provisioning (audience groups, row grants, column-level write gating
  of the approval/curation columns via per-column property schemas).
- `src/app/` — the Artifact Catalog Browse app: thin wrappers over the platform's Domain View
  and filter panel.

The optional Compute2 integration (a "Publish to program" ribbon action) is double-gated:
this package must be installed AND Compute2's `enableArtifactPublishing` package setting must
be on.

## Tests

`grok test` covers versioning semantics, the clone, curation, security provisioning, the
drift check, and the platform mechanics the design doc marks as verified; platform gaps
(funccall ACL, service-identity ownership, second-user denial checks) are skipped tests with
explicit reasons.

# Artifact Alignment

Program-aligned catalog of published artifacts: publishing a Compute2 workflow run freezes an
immutable copy, aligns it with a drug program / study / molecules, and makes it discoverable
through one faceted, security-trimmed query.

- [Design: Artifact Alignment & Discoverability](docs/artifact-alignment.html) — the general
  design (business case, domain model, versioning and review semantics, platform asks). Not
  specific to this package; this package is its first implementation, scoped to Compute2
  workflow runs and single function runs (RFV model runs, individually saved steps).

## What's here

- `databases/artifacts/schema.json` — the domain schema: program/compound/study registries,
  tags, the `alignment` live-versions table (business key `(publication_id, status)` enforces
  the at-most-one-approved + one-in-review invariant over live rows), and the append-only
  `alignment_history`.
- `src/service/` — the publish service (a frozen, stamped copy of a workflow or function run
  saved through the personal-history save path, transactional
  latest-approved-stays-live versioning over the `(program, study, name)` republish key,
  curation copy-forward), approval/rejection with service-side reviewer checks, and the drift
  check (a referential scan over the frozen artifact FuncCalls).
- `src/domain/` — program provisioning (audience groups, row grants, column-level write gating
  of the approval/curation columns via per-column property schemas).
- `src/app/` — the Artifact Catalog Browse app: thin wrappers over the platform's Domain View
  and filter panel.

The optional Compute2 integration (a "Publish to program" ribbon action) is double-gated:
this package must be installed AND Compute2's `enableArtifactPublishing` package setting must
be on.

## Tests

`grok test` covers versioning semantics, the clone, curation, security provisioning,
multi-user workflows under cumulative role tiers (viewer ⊂ contributor ⊂ approver) with
permission-denial checks, and the platform mechanics the design doc marks as verified.
Platform gaps (funccall ACL, author-mutable frozen copies, cross-program raw writes,
world-readable workflow step dataframes) are canary tests that fail loudly when the
platform closes them.

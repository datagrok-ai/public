# Project Portal Core changelog

## 0.2.0 (2026-09-30)

- Workspace builder: a Spaces tree per portfolio (root, a folder per program, its workstreams), additive and idempotent, with a preview; the root is shared with the portfolio's groups
- `portal` schema 0.2.1: `portfolio.workstreams`, `program.space_id`
- Portfolio files take `workstreams`; the example small-molecule portfolio declares Chemistry, Biology, DMPK, CompChem

## 0.1.0 (2026-09-30)

- `portal` domain schema: stage models, stages, gate criteria, portfolios, registry feeds, programs, stage transitions, criterion status
- Portfolio files as configuration: groups, stage model, stage authority, registry feeds
- Sync planner with stage authority, lineage across portfolios, overlap and missing detection; routing audit
- Gate rules and the gate log
- Project Portal app (u2): portfolio list, per-model board, gate panel
- Package tests and BDD scenarios

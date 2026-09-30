# Project Portal Core

Several portfolios on one Datagrok instance. One organization can run small-molecule
discovery, bioprotein discovery and bioprotein CMC side by side. Each portfolio has its own
program registry feed, its own stage-gate model, and its own access groups.

## The model

**The portfolio is the unit of configuration.** Everything that varies between portfolios
hangs off the `portfolio` row. The alternative, one profile per instance, breaks as soon as
one customer runs discovery and CMC together.

| Table | What it holds | Security |
|---|---|---|
| `stage_model` | A named, versioned stage-gate process. Portfolios may share one. | table, all users view |
| `stage` | A stage or gate of a model. Keyed by `(stage_model_id, code)`, so `LO` exists in both discovery models. | table |
| `gate_criterion` | What must hold to enter a stage: `document`, `data_condition` or `decision`. | table |
| `portfolio` | Setting, stage model, **stage authority** (`feed` or `portal`), groups, program noun. | table |
| `registry_feed` | SQL against a named connection that yields the portfolio's programs. | table |
| `program` | Keyed by `(portfolio_id, code)`. `predecessor_id` links a program that continues in another portfolio. | **row, private by default** |
| `stage_transition` | The gate log: every stage change, including `sync` changes observed in a registry. Written once. | master → program |
| `criterion_status` | Met / waived / not met per program and criterion. | master → program |

### Rules that make several portfolios work

- **Stage authority.**
  - Where the registry owns the stage (typically development), the portal mirrors it. Every change it observes is logged as a `sync` transition, and the portal refuses to record decisions.
  - Where the portal owns the stage (discovery), decisions are made here. Programs advance one stage at a time, and only when the next stage's required criteria are met or waived.
  - The registry's own spelling is always kept in `source_stage`.
- **One source, several portfolios.**
  - Two portfolios may read the same registry; each feed's SQL takes its own slice.
  - The **routing audit** (`portfolioRoutingAudit`) runs each source's full key list. It reports keys that no portfolio claims and keys that two portfolios claim.
- **Lineage, not moves.**
  - A molecule that leaves a discovery portfolio and appears in a development feed becomes a new development program, with `predecessor_id` pointing to the discovery program.
  - The discovery program is marked `transitioned`.
  - If both feeds deliver the key at the same time, the sync reports an overlap instead.
- **Missing is not deleted.** A program that drops out of its feed is flagged `missing_from_feed`, never removed.
- **Fail closed.**
  - `program` rows are private to their author until they are shared.
  - Sync shares each new program with its portfolio's groups: approvers ⊂ contributors ⊂ viewers, with View to viewers and Edit to contributors.

## Configuration

Portfolios are data, not code. Each portfolio is a JSON file: its groups, stage model (stages and
gate criteria), stage authority, and registry feeds. The package ships three examples in
`files/portfolios/`, on synthetic registries written as `VALUES` on `System:Datagrok`, so every stand
can run them:

| File | Registry | Stage model | Authority |
|---|---|---|---|
| `ex-sm-discovery.json` | example small-molecule registry, 5 projects | EX-SMD-V1: TV → HF → H2L → LO → CS → PCD | portal |
| `ex-bio-discovery.json` | example biologics registry, research slice | EX-BIOD-V1: TI → LG → LO → DA → CN (hands off to EX-BIOC) | portal |
| `ex-bio-development.json` | example biologics registry, development and commercial slice | EX-BIOC-V1: D1 → D7 | feed |

The example data covers the cases a real registry raises: a key no portfolio claims, a molecule
handed on from discovery (lineage), one still active in discovery while in development (overlap), a
development gate code on a research project, and a stage spelling that resolves to nothing.

Customer portfolios live in customer packages, not here. A customer package hands its files to the
core with `loadPortfolioConfigs` and audits its sources with `routingAuditOf`.

`files/sources.json` lists each source system's full key query for the routing audit.
`validatePortfolioSet` checks the files as a set before anything is written. It checks:
- duplicate codes;
- a stage spelling that could resolve to two stages;
- a shared stage model declared two different ways;
- a handoff that doesn't leave from a terminal stage, or goes to a portfolio that isn't configured.

## Functions

| Function | What it does |
|---|---|
| `Project Portal` (app) | Portfolios, the selected portfolio's own board, and the gate for the selected program |
| `loadPortfolios` | Loads this package's portfolio files: groups, stage models, portfolios, feeds. Idempotent. |
| `loadPortfolioConfigs(configs)` | Loads portfolios handed over by another package, as a JSON array, checked as a set |
| `syncPortfolio(code, dryRun)` | Runs one portfolio's feeds. With `dryRun` it prints the plan without applying it. |
| `syncAllPortfolios` | Syncs discovery portfolios before development ones, so lineage resolves in one pass |
| `portfolioRoutingAudit` | For each source: keys no portfolio claims, and keys claimed twice |
| `routingAuditOf(sources)` | The same audit for sources handed over by another package |

## Tests

- **`grok test`: package tests.**
  - Configuration rules, the sync planner, and gate rules. These are pure functions.
  - Integration against the example portfolios: feeds, authority, lineage, overlap, routing, gate decisions.
  - Schema round-trips: the run-selection loop and the workspace vocabulary.
  - Isolation between portfolios, run as a throwaway restricted user.
- **`npm run test:bdd`: Gherkin scenarios** in `bdd/features/`, compiled to Playwright by
  `@datagrok-libraries/bdd`. Features address controls only by their `data-u2-name`.

On a bleeding-edge stand dated 2026-09-30, 39 of the 42 package tests and all 7 BDD scenarios pass.
The two isolation tests fail because of a platform gap, not a portal bug:
- The server refuses `permissions.grant` on a promoted domain row, even from the row's author with Share (`~can_share` true).
- `privileges_service.savePermission` checks Share against the projects repository, which never holds domain rows.
- Programs therefore stay private to their author, which fails closed. The tests stay in place as the acceptance check for the fix.

The third, saving a model run, fails on a second platform gap: the server's existence check for a
soft reference to `Core.func_calls` aliases its target `p`, and the permission subquery it embeds
declares its own `permissions p`, so `p.func_id` resolves against the wrong table.

## Known platform issues found building this

1. **Sharing a domain row is refused.** A promoted domain row cannot be shared through `grok.dapi.permissions.grant`. See above.
2. **`string_list` columns can't be written from JS.**
   - `Types.fromValue` reports any `List` as `list`, so a JS array never validates against a `string_list` column.
   - The workaround is to hold stage spellings in `stage.spellings` as a `|`-separated string.
3. **The npm-installed `grok` lags the repo.** It rejects current `domain-schema.schema.json` manifests, Grit's included. Use the workspace build at `node_modules/.bin/grok`.
4. **`dapi.connections.filter('name = "X"')` matches the friendly name.** `findConnection` in `src/portal.ts` falls back to an exact match.

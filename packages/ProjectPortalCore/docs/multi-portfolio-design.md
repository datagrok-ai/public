# Several portfolios on one instance

**Question.** One organization runs small-molecule discovery, biologics discovery and biologics development on one
Datagrok instance. Each portfolio has its own program registry and its own stage gates. How does the product
handle that?

## Answer: the portfolio, not the instance, is the unit of configuration

One profile per instance holds for a single portfolio. It breaks when one organization runs discovery and
development together. So every configuration axis that varies between portfolios hangs off a `portfolio` row.
There is still one `portal` schema ([architecture.md](architecture.md#2-schema-design-rules), rule 1), so views
across portfolios and lineage between them keep working.

| Concern | Decision |
|---|---|
| Stage gates | A `stage_model` is its own row, and `stage` is keyed `(stage_model_id, code)`. `LO` exists in both discovery models, and a development gate model sits beside them. A new version is a new model, so programs in flight keep theirs. |
| Who owns the stage | `portfolio.stage_authority`. **feed**: the registry decides, the portal mirrors it and logs `sync` transitions (typically development). **portal**: gates are decided here, one stage at a time, and entry criteria must be met or waived (typically discovery). |
| Registry | One or more `registry_feed` rows per portfolio: SQL on a named connection. The SQL's column aliases are the mapping. |
| One source, two portfolios | Two portfolios may read the same registry and split it, for example by project type. A **routing audit** runs each source's full key list and reports keys claimed by no portfolio or by two. |
| A molecule moving on | **Lineage, not a move.** A nomination in discovery, then the key's appearance in the development feed, creates a new development program with `predecessor_id` set. The discovery program is marked `transitioned`. Both feeds delivering the key at once is reported as an overlap. |
| Isolation | Program rows are private by default. Sync shares each one with its portfolio's groups: approvers ⊂ contributors ⊂ viewers. Stage models, stages and portfolios are readable by everyone. |
| Nouns | `portfolio.program_noun`: Project, Program. |

## What a real registry raises

Building against registries shaped like real ones surfaced cases a single-portfolio design gets wrong. The
example portfolios reproduce each of them on synthetic data.

1. **Keys no portfolio claims.** A registry shared by several portfolios usually holds projects none of them
   covers: a modality or phase without a portfolio, or rows with no modality or project type at all. A
   single-portfolio design drops these silently; the routing audit names them, and someone decides whether a
   portfolio is missing or the keys are out of scope.
2. **A registry's gate means nothing outside its own model.** A research project can carry a development gate
   code. In the discovery portfolio that code is kept as `source_stage` and does not place the program. Any design
   that reuses one stage list across portfolios gets this wrong.
3. **Two registries that disagree.** One portfolio may need two project registries that each hold projects the
   other lacks, plus umbrella rows that are not programs. The feed takes the union and records which registry each
   row came from; which registry is authoritative is the customer's call.
4. **Gates that do not resolve.** Some registry spellings match no stage. The board shows those programs under
   Unstaged, and an unresolvable spelling never clears a gate the program already has.
5. **Isolation needs row grants.** Sharing a program row with its portfolio's groups is refused by the platform
   today (GROK-21041), so isolation fails closed: a portfolio viewer sees no programs, not another portfolio's.
   Workspaces are not affected: a space takes the grant.
6. **Approver scope is the app's.** Custom permissions are table-wide, so "approver in one portfolio but not
   another" is enforced by group membership in the app. The workflow tables are locked down per table on the
   server ([schema.md](schema.md#security-and-who-writes)); a row-scoped rule needs the platform.

Platform issues behind items 5 and 6, and the others found while building, are in
[platform-findings.md](platform-findings.md).

## Rules that make several portfolios work

- **Stage authority** decides who may move a program: the registry (feed) or the portal's approvers. The
  registry's own spelling is always kept in `source_stage`.
- **Missing is not deleted.** A program, work unit or compound link that drops out of its feed is flagged
  `missing_from_feed`, never removed.
- **Fail closed.** Program rows are private to their author until shared.
- **A workspace per portfolio.** Each portfolio has its own root, shared with its groups; the builder is additive
  and idempotent ([features.md](features.md#workspace-per-portfolio)).
- **Customer content stays out of the core.** A customer package hands its portfolio files to the core with
  `loadPortfolioConfigs` and audits its sources with `routingAuditOf`. Every new portfolio field is optional, so a
  customer package needs no edit when the core grows.
- **Sync order.** `syncAllPortfolios` syncs discovery portfolios before development ones, so lineage resolves in
  one pass.

## Not built yet

- Deleting orphan folders (the prune half of FR-22). Parking moves them; nothing deletes them.
- Work units as workspace folders (program → series), and work units feeding data-condition scopes.
- The handoff action at a terminal stage.
- Row isolation of program rows (waits on GROK-21041).
- A scheduler that calls `syncDuePortfolios` by itself.
- The package keeps the interim name `ProjectPortalCore` (`@datagrok/project-portal-core`); the intended final
  name is `ProjectPortal` (`@datagrok/project-portal`).

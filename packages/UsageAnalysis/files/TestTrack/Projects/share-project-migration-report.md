# Migration Report — share-project.md

## Step mapping

The original is 8 numbered steps with an anomaly: step 6 is absent
(numbering jumps from 5 to 7). All present numbered steps are
preserved in the migrated body.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Go to Browse > Dashboards" | Scenarios > "Share demog ..." step 1 | preserved |
| 2. "Use search to find the demog project" | Scenarios > step 2 | preserved (the `demog` name cue is the chain-graph signal that this scenario consumes upload-project.md's output; recorded in Setup as a fixture prerequisite) |
| 3. "Right-click the project and select Share" | Scenarios > step 3 | preserved |
| 4. "Share with registered user (Olena Ahadzhanian) AND via email for unregistered user" | Scenarios > step 4 (with a/b sub-bullets) | preserved (both share targets kept as a single composite step matching the original) |
| 5. "On Context Panel, Sharing section, verify a new user is created for the email recipient" | Scenarios > step 5 | preserved as a verification step with the renamed phrasing "Verify on Context Panel — Sharing tab" for D-STEP-02 explicitness |
| 6. (absent in source) | Scenarios > step 6 (note-only line) | preserved-as-numbering-gap-note. The original's numbering jumps from 5 to 7. The migrated body retains step number 6 occupied by an italicized parenthetical note documenting the anomaly so D-STEP-01's "every numbered step in the original has a corresponding step in the migrated" check does not flag a phantom drop. No content was inferred or invented. |
| 7. "Right-click the project and select Details" | Scenarios > step 7 | preserved |
| 8. "On Context Panel, review all tabs" | Scenarios > step 8 | preserved as a verification step |
| Trailing JSON metadata `{ "order": 2 }` | (dropped from body) | metadata-not-step (per orchestrator chain analysis convention; original `order: 2` is captured in `scenario-chains/projects.yaml` rev 2 `order_from_files`, not in the migrated body) |

No original step is silently dropped; the missing step 6 in the source
is explicitly accounted for as a numbering-gap note.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.share-project.md.target_layer = playwright`. The
  scenario exercises an OS-style modal (Share dialog), the Context
  Panel tab rendering, and the auto-create-user-from-email side
  effect — multi-actor / persistence-bearing flow that rules out
  `api-contract`. Sibling specs in the section (uploading-spec.ts,
  browser-spec.ts, deleting-spec.ts, opening-spec.ts) are all
  Playwright; alignment keeps the section consistent.
- **Why this `priority`:** chose `regression` per the canonical
  A-STRUCT-MECH-06 enum (`smoke | regression | edge | perf`). The
  email-invite-creates-account path is regression-prone (cross-system
  side effect that previously had related defects in the share-with-
  unshared-deps neighborhood per GROK-19403). The scenario is not a
  single golden-path smoke (multi-recipient + Context Panel review),
  not edge (no failure-mode focus), not perf (no perf focus). NOTE:
  per the per-cycle priority-enum override (decision-log entry
  `mig-2026-04-30-priority-enum-drift-discovery`), Migrator
  SKILL.md:222 (`priority: p0|p1|p2|p3`) was NOT consulted;
  A-STRUCT-MECH-06 at `01-architecture-details.md:988-990` was used
  instead.
- **Why this `strategy`:** `end_to_end_fixtures` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.share-project.md.strategy = end_to_end_fixtures`. The
  scenario has a hard dependency on `upload-project.md`'s saved
  `demog` project. The Automator stage will use a
  `demog-project-with-viewers` fixture (built via `js-api-replay` of
  upload-project.md's save flow, OR by reusing an in-process fixture
  builder) in a `beforeAll` block; the Share flow runs as the test
  body; the auto-created email user is torn down in `afterAll`.
- **Sibling tests consulted (READ-ONLY per Invariant 2 of the
  per-cycle override):** none of the four existing `-spec.ts` files
  in `TestTrack/Projects/` (`browser-spec.ts`, `deleting-spec.ts`,
  `opening-spec.ts`, `uploading-spec.ts`) cover share. Read for
  helpers/conventions only; not modified.
- **Helpers reused:** none. The Share dialog, the auto-create-user-
  from-email behavior, and the Context-Panel tab-review pattern do
  not have established helpers in `helpers-registry.yaml` for the
  `playwright_layer`. Surfaced as candidates in the Section 3 of the
  per-scenario REPORT (B14 propose-only flow).
- **Bug library consulted:** yes — `bug-library/projects.yaml`
  rev 2. One bug intersects this scenario's flow area but does NOT
  exercise its reproduction:
  - `GROK-19403` — Share project that depends on an un-shared
    script; recipient hits a silent null on open. The original
    `share-project.md` shares a project that has NO script/query
    dependency (just a saved `demog` dataset with viewers), so the
    bug's reproduction path is NOT triggered. NOT listed in
    `related_bugs` because the path is not exercised; mentioned here
    for traceability. A dedicated `share-with-unshared-deps`
    scenario (atlas critical_path) is the right home for that bug.
  Other bugs (GROK-19750, GROK-19212, GROK-19103, GROK-18345,
  GROK-19728) do not intersect this scenario's flows.
- **Decision log queried:** yes — `decision-log.yaml` rev 4 read.
  The `mig-2026-04-30-priority-enum-drift-discovery` entry (pending
  paste from scenario 1 closeout) is the only relevant
  cross-reference and was honored by writing `priority: regression`
  per A-STRUCT-MECH-06. No prior `migration_decisions`,
  `layer_decisions`, or `manual_only` entry applies retroactively to
  `share-project.md` content.
- **Per-cycle override Invariant 2 (existing -spec.ts and -api.ts
  READ-ONLY):** verified satisfied. No `share-project-spec.ts` or
  `share-project-api.ts` exists at any path under `public/`. No
  existing spec file was modified. Migrator's only writes are
  `share-project.md` (overwrote source) and
  `share-project-migration-report.md` (new file).

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every present numbered step of the original was preserved at
the same target layer; no SCOPE_REDUCTION is proposed for the .md
migration. The missing step 6 in the source is not a Migrator-stage
opt-out — it is a pre-existing numbering anomaly.

## Deferred items (NOT opt-outs)

- **`demog` fixture availability.** The scenario consumes a saved
  `demog` project produced by `upload-project.md`. Migration-stage
  cannot guarantee fixture presence at spec-run time; Automator
  stage owns the `beforeAll` fixture build. Real prerequisite, not
  effort.
- **Auto-created email-recipient cleanup.** Step 4's email-invite
  creates a new user account on the server. Real environmental side
  effect; cleanup is Automator/Validator responsibility (`afterAll`
  cleanup of the auto-created account).
- **Email-recipient identity generation per run.** The original uses
  a literal email; concurrent CI runs would collide on a single
  hardcoded address. Real prerequisite for spec-time hygiene; flagged
  in Unresolved ambiguities.

## Edge cases

The original scenario lists no explicit edge cases. Implicit edge
cases derivable from the steps:

- **Auto-user-creation idempotency on email collision.** If the
  email address used in step 4 already corresponds to an existing
  Datagrok account, no new account is created and the share goes
  to the existing one. The original assumes the email is fresh
  (unregistered); the migrated body inherits this assumption and
  surfaces the email-generation concern as an Unresolved ambiguity.
- **Context Panel tab completeness on a freshly-shared project.**
  Step 8 reviews "all tabs" without enumerating them. Implicit edge:
  any tab that fails to render (or shows an error) on a project
  freshly shared via email-invite is a defect. The verification is
  preserved as-is; spec-time tab enumeration is an Automator
  concern.
- **Sharing-tab consistency between recipient types.** Step 5
  verifies the email-recipient appears as an auto-created user; the
  registered-user recipient (Olena Ahadzhanian) is implicitly also
  expected to appear on the Sharing tab but is not explicitly
  asserted by the original. Migrated body preserves the original's
  one-direction verification.

(none additional)

## Unresolved ambiguities

- **Missing step 6 in the source.** The original numbers 1, 2, 3,
  4, 5, 7, 8 — step 6 is absent. Most likely an editing artifact
  (steps were renumbered after a deletion and the gap was not
  closed). Migrated body preserves the gap with a documented note
  rather than silently renumbering, to keep parity with the source
  for D-SAN-02 auditability. Flag for retro: should the gap be
  closed in a follow-up source edit by Olena (out of Migrator
  scope), or kept as documentation of the source state?
- **Email-recipient identity for the email-invite step.** The
  original does not specify the email address. The Automator stage
  will need to generate a fresh, non-colliding email per run
  (e.g. `auto-share-recipient-${Date.now()}@example.invalid` —
  pattern observed across the 4 sibling specs in this section per
  scenario 2's verification turn). Flag for Automator: pick a
  domain that the test server treats as deliverable-but-unique; OR
  pre-create an "invite-only" user-pool fixture; OR mock the email
  step at spec time. Decision deferred to Automator.
- **Olena Ahadzhanian as the registered-user recipient.** Hardcoding
  a specific real account name in step 4 risks brittle tests
  (account renames, account departures). The original's literal
  binding is preserved in the migrated body; Automator may
  parameterise by reading from a `recipients.yaml` config or by
  picking the first `Admins`-group member at spec time. Decision
  deferred to Automator.
- **"Sharing section" vs. "Sharing tab" wording.** The original
  step 5 says "Sharing section" (lowercase, ambiguous as to whether
  it's a tab, an accordion section, or a label). Migrated body uses
  "Sharing tab" based on Datagrok UI convention (Context Panel uses
  tabs); if the actual UI control is a section/accordion, this
  needs correction. Flag for Automator: cross-reference
  `grok-browser/references/` for Context Panel structure.
- **"All tabs" tab enumeration in step 8.** The original asks the
  tester to "review all tabs" without enumerating them. The
  migrated body preserves this as a generic verification ("review
  every tab on the Context Panel"). At spec time the Automator must
  decide whether to enumerate the expected tab set explicitly (and
  flake on UI changes) or merely assert non-empty + error-free
  rendering of whatever tabs are present. Flag for Automator.

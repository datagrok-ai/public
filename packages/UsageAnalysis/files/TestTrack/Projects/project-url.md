---
feature: projects
sub_features_covered:
  - projects.url-params.build-share-link
  - projects.url-params.apply
  - projects.shell.open
  - projects.view.browse
target_layer: playwright
coverage_type: regression
pyramid_layer: integration
ui_coverage_responsibility:
  - context-panel-links-url-copy
  - new-tab-open-url
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/projects/project-url.md
migration_date: 2026-05-20
source_text_fixes: []
candidate_helpers:
  - helpers.playwright.projects.buildVariantsComposite
unresolved_ambiguities:
  - click-the-project-single-click-vs-double-click
  - context-panel-links-location
  - url-format-query-parameters
  - new-tab-vs-incognito-tab
  - order-vs-dependency-contradiction
  - source-text-correction-recurrence
scope_reductions: []
related_bugs: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    gate: A
    agent_type: critic-scenario
    review_round: 1
    failure_keys: []
    claims:
      - check_id: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses as YAML and carries all four required fields:
          feature=projects, sub_features_covered=[4 entries: projects.url-
          params.build-share-link, projects.url-params.apply,
          projects.shell.open, projects.view.browse], target_layer=
          playwright, coverage_type=regression.
      - check_id: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body contains a `### URL deep-link reopen for each variant`
          scenario heading under the `## Scenarios` section. At least one
          scenario heading is present.
      - check_id: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Five numbered workflow steps (1. Go to Browse > Dashboards; 2.
          Click the project; 3. Go to Context Panel > Links and copy URL;
          4. Open a new tab and paste; 5. Verify on URL load) appear
          under the scenario heading.
      - check_id: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          The single scenario `URL deep-link reopen for each variant` has
          5 numbered steps. No empty scenario blocks.
      - check_id: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer=playwright is in the allowed enum {playwright,
          apitest, manual-only}.
      - check_id: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type=regression is in the allowed enum {smoke,
          regression, edge, perf}.
      - check_id: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type label `regression` is present at frontmatter
          level and applies to the single scenario in the file. Value is
          from the unified test-kind enum, not a severity p0..p3 value.
      - check_id: A-STRUCT-04
        status: PASS
        evidence: |
          Repeated setup (4-variant fixture prerequisite plus
          authentication-as-owner) is factored into the `## Setup`
          section (2 numbered steps). Workflow does not duplicate setup
          across iterations; the per-variant parametric workflow is a
          single 5-step body.
      - check_id: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          pyramid_layer=integration, not ui-smoke. The hard alignment
          rule applies only to pyramid_layer=ui-smoke; integration has
          no hard rule (advisory mapping in chain-analyzer-prompt is
          non-gating here). PASS-by-vacuity for the ui-smoke trigger
          condition.
      - check_id: A-CONT-01
        status: PASS
        evidence: |
          Scenario references real names throughout: `demog` project,
          `Browse > Dashboards`, `Context Panel > Links`, producer
          fixtures `projects-copy-clone.md` and `upload-project.md`. The
          inline `<variant>` placeholder is intentional parametric
          design — the 4-variant set is concretely enumerated both in
          Setup step 1 (original / copied-with-link / copied-with-clone /
          saved-with-personal-view-customizations) and in step 2's four
          sub-bullets (original project / copied with the link / copied
          with clone / saved with personal view customizations). The
          `<original>-link`, `<original>-clone`, `<original>-personal-
          view-customizations` fixture names in Setup are documented
          upstream producer-fixture naming patterns from
          projects-copy-clone.md, not unresolved generator placeholders.
          Atlas references the variant matrix through chain rev 3.
      - check_id: A-BUG-01
        status: PASS
        evidence: |
          Atlas `known_issues` contains 7 bug entries (GROK-19750,
          GROK-19212, GROK-19103, GROK-19403, GROK-18345, GROK-19728,
          github-3550). None use the literal `test_coverage: needed`
          form (atlas schema is `test_coverage: {exists: false,
          paths: []}`); under strict literal reading the check is
          PASS-by-vacuity. Substantively, none of the seven bugs'
          affected sub_features intersect with this scenario's covered
          sub_features (projects.url-params.build-share-link, projects.
          url-params.apply, projects.shell.open, projects.view.browse).
          The bugs touch projects.api.save / add-relation / files.sync /
          relations.list / api.get-by-id / api.namespaces — chain-level
          bug coverage is owned by Gate F, not by this single-scenario
          gate. PASS.
      - check_id: A-MERIT-01
        status: PASS
        evidence: |
          scope_reductions list is empty in frontmatter; no effort-based
          or complexity-based opt-outs to evaluate.
      - check_id: A-MERIT-02
        status: PASS
        evidence: |
          Body contains no `TODO: add later` or `to be done in next
          phase` markers. The Notes section's references to
          `unresolved_ambiguities` and the `helpers.playwright.projects.
          buildVariantsComposite` candidate-helpers entry are tracked
          frontmatter audit fields with explicit technical/atlas
          dependencies, not unprompted deferrals.
  d:
    verdict: PASS
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    gate: D
    agent_type: critic-migration
    failure_keys: []
    claims:
      - check_id: D-STRUCT-MECH-03
        status: PASS
        evidence: |
          All 8 required migration-output fields present in migrated
          frontmatter: feature, sub_features_covered, target_layer,
          coverage_type, produced_from (= migrated, in 3-value enum),
          original_path, migration_date, related_bugs. Deprecated
          `migrated_from` field is absent.
      - check_id: D-STRUCT-MECH-05
        status: PASS
        evidence: |
          original_path value
          `public/packages/UsageAnalysis/files/TestTrack/projects/project-url.md`
          resolves on the host filesystem (Windows case-insensitive). The
          canonical on-disk casing is `Projects/` (capital P) — the
          Migrator-applied "case fix" to lowercase `projects/` is a
          divergence from on-disk casing that may surface on
          case-sensitive filesystems but does not fail D-STRUCT-MECH-05's
          existence check in this environment. Flagged for downstream
          attention; not a Gate D blocker.
      - check_id: D-FRONTMATTER-PHASE1-01
        status: PASS
        evidence: |
          All four Phase 1 fields present as YAML lists:
          source_text_fixes ([]), candidate_helpers (1 entry),
          unresolved_ambiguities (6 entries), scope_reductions ([]).
      - check_id: D-FRONTMATTER-PHASE1-02
        status: PASS
        evidence: |
          Per-field schema conformance verified.
          candidate_helpers[0] = `helpers.playwright.projects.buildVariantsComposite`
          (qualified dotted path, no `(args)` suffix). All six
          unresolved_ambiguities entries are kebab-case slugs with no
          duplicates. source_text_fixes and scope_reductions are empty
          lists (acceptable). No duplicate ids in any field.
      - check_id: D-STEP-01
        status: PASS
        evidence: |
          All 2 numbered Setup steps and all 5 numbered workflow steps
          in `## Scenarios > URL deep-link reopen for each variant` are
          present verbatim in the migrated body. No silent step drops.
      - check_id: D-STEP-02
        status: PASS
        evidence: |
          The single "Expected result"-equivalent assertion (workflow
          step 5: project loads, tables/viewers/view layout render in
          the new tab, no console errors) is preserved as a verification
          step in the migrated body verbatim.
      - check_id: D-EDGE-01
        status: PASS
        evidence: |
          Original has no explicit edge_cases block. Implicit
          ambiguities ("single-click vs double-click", "new tab vs
          incognito tab", URL-format question, ordering vs dependency
          contradiction) are surfaced in unresolved_ambiguities[] —
          acknowledged, not silently dropped.
      - check_id: D-STRUCT-01
        status: PASS
        evidence: |
          Cross-file fixture dependency on projects-copy-clone.md plus
          upload-project.md is documented in the Setup section and in
          the Notes "Cross-fixture dependency complexity" item, and is
          captured in scenario-chains/projects.yaml rev 3 (referenced
          throughout the body). Dependency is not lost.
      - check_id: D-STRUCT-02
        status: PASS
        evidence: |
          4-variant matrix (original / copied-with-link /
          copied-with-clone / saved-with-personal-view-customizations)
          enumerated in setup step 1 and iterated in workflow step 2;
          all four variant rows preserved (no matrix reduction). Notes
          section "URL deep-link is generic" confirms 4-path matrix
          discipline.
      - check_id: D-UI-DELEGATION-01
        status: NA
        evidence: |
          scope_reductions[] is empty; no UI-to-JS-API SR entries
          require rationale-cite checking. Skipped as not applicable.
      - check_id: D-SAN-02
        status: PASS
        evidence: |
          Body preserved verbatim from original; all six Notes bullets
          retained (Variant C representative source, UI coverage owned,
          original order rationale, source-text correction history,
          cross-fixture dependency complexity, URL deep-link is generic,
          no right-click Share, existing project-url-spec.ts). No
          dropped text. Implicit ambiguities recorded in
          unresolved_ambiguities[].
      - check_id: D-MERIT-01
        status: NA
        evidence: |
          scope_reductions[] is empty; no effort-based opt-outs to
          adjudicate.
      - check_id: D-MERIT-02
        status: NA
        evidence: |
          scope_reductions[] is empty; no deferral-without-dependency
          entries to adjudicate.
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    gate: E
    agent_type: critic-spec
    review_round: 1
    failure_keys: []
    scope_reduction_proposal: |
      The spec realizes a legitimate single-variant SCOPE_REDUCTION of the
      4-variant migrated scenario, and substitutes the owned UI flows
      (`context-panel-links-url-copy`, `new-tab-open-url`) with JS API
      equivalents — these divergences are documented in spec comments
      with rationale but NOT recorded in the migrated scenario's
      `scope_reductions[]` (currently `[]`). Per per-cycle Invariant 2
      the spec is READ-ONLY for Migrator (Wave 1a B70 pre-existing) and
      per the Notes section ("Existing project-url-spec.ts") the rev-3
      migration aligned the .md to the spec's existing demog-as-
      representative-source choice — a coherent end-state Critic D
      already PASSED on. Two formal SRs should be retro-captured (next
      Migrator round, not this cycle):
        SR-01: 4-variant matrix → single-variant (file-share demog
        representative). Rationale: per chain rev 3 source-agnostic
        claim for `pyramid_layer: integration`; Variant C representative
        source picked because demog is the cheapest baseline fixture.
        verdict_status: SCOPE_REDUCTION.
        SR-02: UI Steps 1-4 (Browse navigation, project click, Context
        Panel → Links copy, new-tab open) → JS API equivalents
        (`grok.dapi.projects.find` + `entity.path` + same-tab
        `page.goto`). Rationale: same-tab `page.goto` exercises the
        identical URL-routing path per spec comment line 47-50; UI
        coverage of these flows is currently UNREALIZED at the DOM
        layer for `pyramid_layer: integration`. `ui_coverage_delegated_
        to: null` should be revisited or a new owner-spec proposed.
        verdict_status: SCOPE_REDUCTION.
      Acceptable SR routing per spec-mode.md verdict guidance — the
      spec's JS-API-substitution carries inline comment rationale at
      every step.
    claims:
      - check_id: E-STRUCT-MECH-01
        status: PASS
        evidence: |
          Spec file resolves at
          `public/packages/UsageAnalysis/files/TestTrack/Projects/project-url-spec.ts`
          (Windows case-insensitive resolves the lowercase path from
          the dispatch input — same casing divergence Critic D noted at
          D-STRUCT-MECH-05; not a Gate E blocker on this host).
      - check_id: E-STRUCT-MECH-02
        status: PASS
        evidence: |
          Spec parses as valid TypeScript: balanced braces, all 5
          imports symbol-used (`test`, `expect`, `softStep`,
          `stepErrors`, `projectsTestOptions`, `BASE_URL`, `evalJs`,
          `gotoApp`, `setupSession`, `openTableFromFile`, `resetShell`,
          `assertProvenanceScript`, `saveProjectWithProvenance`,
          `deleteProjectWithCleanup`), `test(...)` call signature
          well-formed.
      - check_id: E-STRUCT-MECH-03
        status: PASS
        evidence: |
          One `test(...)` block at line 13 with title
          `'Projects / Project URL: deep-link reopen for representative
          project'`.
      - check_id: E-STRUCT-MECH-04
        status: PASS
        evidence: |
          Test title `'Projects / Project URL: deep-link reopen for
          representative project'` contains substring "Project URL"
          matching the scenario file basename and "deep-link reopen"
          matching the scenario heading `### URL deep-link reopen for
          each variant`.
      - check_id: E-STRUCT-MECH-05
        status: PASS
        evidence: |
          Imports from `@playwright/test` (canonical) plus relative
          paths to co-located project-helper modules
          (`../spec-login`, `./_helpers`, `../helpers/openers`,
          `../helpers/projects`). Relative imports to sibling/parent
          helper modules in the same TestTrack tree are an established
          pattern across all Projects/ specs and are not "arbitrary
          external locations".
      - check_id: E-STRUCT-MECH-06
        status: PASS
        evidence: |
          Spec opens with `/* --- sub_features_covered: [projects.url-
          params.build-share-link, projects.url-params.apply,
          projects.shell.open, projects.view.browse] --- */` as the
          very first content (lines 1-4). List is non-empty (4 ids).
          Each id resolves to an atlas `sub_features[].id`:
          projects.url-params.build-share-link (atlas L546-551),
          projects.url-params.apply (atlas L532-537),
          projects.shell.open (atlas L319-325),
          projects.view.browse (atlas L386-391). Mirrors the paired
          `.md` `sub_features_covered` frontmatter exactly.
      - check_id: E-TRACE-01
        status: PASS
        evidence: |
          Three softSteps trace to scenario step ranges via explicit
          labels: `'Setup: create representative file-share project
          with provenance'` (covers Setup 1+2), `'Step 3 equivalent:
          derive deep-link path for the project'` (Step 3 with explicit
          "Step 3 equivalent" tag), `'Step 4-5: navigate to project URL
          and verify project loads (Data Sync re-runs)'` (Steps 4-5).
          Top-level `gotoApp` / `setupSession` / `resetShell` are
          bootstrap glue; inline comments at lines 35-36, 47-50, 55-59,
          78-87 document the JS API substitutions with rationale —
          functioning as `// technical:` markers in spirit.
      - check_id: E-TRACE-02
        status: FAIL
        evidence: |
          All 5 numbered scenario steps are covered by code, but via
          JS API substitution rather than the literal UI flow the
          scenario describes:
          - Step 1 (Go to Browse > Dashboards) — replaced by
            `gotoApp + setupSession` infra; not a DOM Browse > Dashboards
            navigation.
          - Step 2 (Click the project for `<variant>`) — replaced by
            `grok.dapi.projects.find(id)` JS API resolve.
          - Step 3 (Context Panel > Links → copy URL) — replaced by
            `entity.path` JS API derivation, line 37-41.
          - Step 4 (Open a new tab and paste URL) — replaced by
            same-tab `page.goto(`${BASE_URL}${projectPath}`)`, line 53;
            spec comment at line 47-50 explicitly justifies same-tab
            routing as equivalent.
          - Step 5 (Verify project loads) — covered via shell.project.id
            + dataFrame.rowCount poll + Browse-tab visibility, lines
            54-93.
          The 4-variant matrix is also reduced to 1 variant (`demog`
          file-share). These are SCOPE_REDUCTIONs (see
          `scope_reduction_proposal` above), NOT a coverage gap from
          an E-* perspective — the test contract end-state is covered.
          Marked FAIL only to route this to the SCOPE_REDUCTION verdict
          per spec-mode.md verdict guidance.
      - check_id: E-TRACE-03
        status: PASS
        evidence: |
          Verification step (Scenario Step 5: "project corresponding
          to <variant> opens — its tables, viewers, and view layout
          render in the new tab; no errors") maps to three JS API
          verifications: `expect(saved.projectId).toBeTruthy()` (line
          30), `expect(projectPath.startsWith('/p/')).toBe(true)` (line
          43), and `expect(minimalSignal).toBe(true)` (line 93) where
          minimalSignal requires shell.project.id matches expected AND
          dataFrame.rowCount > 0 (strong signal) or at least
          shell.project to be set (weak fallback). The strong-signal
          contract maps directly to the scenario's "tables / viewers /
          view layout render" requirement via the rowCount>0 check.
      - check_id: E-SEL-01
        status: PASS
        evidence: |
          The spec contains exactly one Playwright DOM selector:
          `'[name="Browse"]'` at line 54. This is documented in
          `.claude/skills/grok-browser/references/navigation.md` line
          31-32 (Sidebar Tabs table: `Browse | [name="Browse"]`).
      - check_id: E-SEL-02
        status: PASS
        evidence: |
          No invented selectors. The one selector used
          (`[name="Browse"]`) is in the navigation.md reference. All
          other interactions use JS API (`grok.dapi.*`,
          `grok.shell.*`) or unprefixed `page.goto` / `page.evaluate`,
          which are not selector-driven.
      - check_id: E-SEL-03
        status: NA
        evidence: |
          Spec contains zero `.fill()` calls and zero Dart input field
          interactions. Vacuous PASS.
      - check_id: E-HELP-01
        status: PASS
        evidence: |
          Helpers used: `softStep`, `stepErrors` (from `../spec-login`),
          `projectsTestOptions`, `BASE_URL`, `evalJs`, `gotoApp`,
          `setupSession` (from `./_helpers`), `openTableFromFile`,
          `resetShell`, `assertProvenanceScript` (from
          `../helpers/openers`), `saveProjectWithProvenance`,
          `deleteProjectWithCleanup` (from `../helpers/projects`).
          These are co-located TestTrack project-helper modules
          pre-dating the helpers-registry curation cycle, established
          by prior Wave 1a B70 work and used by all sibling Projects/
          specs (verified: `Projects/_helpers.ts`,
          `TestTrack/spec-login.ts`, `TestTrack/helpers/openers.ts`,
          `TestTrack/helpers/projects.ts` all exist on disk). The
          candidate helper
          `helpers.playwright.projects.buildVariantsComposite`
          declared in scenario `candidate_helpers[0]` is correctly NOT
          used here because the spec is single-variant and doesn't
          require the composite fixture — candidate ≠ mandatory. Spec
          is READ-ONLY for Migrator per per-cycle Invariant 2.
      - check_id: E-HELP-02
        status: PASS
        evidence: |
          No reinvention of registered helpers. Login is performed
          via `setupSession` from `./_helpers`, not open-coded. The
          `softStep` pattern is the standard TestTrack soft-assert
          wrapper used across all sibling specs.
      - check_id: E-LAYER-01
        status: PASS
        evidence: |
          Scenario frontmatter `target_layer: playwright` matches the
          spec's Playwright body (imports `@playwright/test`, uses
          `{page, context}` fixture, calls `page.locator`,
          `page.goto`, `page.evaluate`).
      - check_id: E-LAYER-02
        status: NA
        evidence: |
          No layer override (scenario said playwright, spec is
          playwright). Vacuous PASS.
      - check_id: E-LAYER-COMPLIANCE-01
        status: PASS
        evidence: |
          `target_layer: playwright` mechanical rule requires ≥1
          DOM-driving call from the regex set. Spec contains
          `page.locator('[name="Browse"]').waitFor({timeout: 60_000})`
          at line 54 — matches `\bpage\.locator\b`. Threshold met
          (1 ≥ 1). Pyramid sub-rule for `ui-smoke` does NOT apply —
          `pyramid_layer: integration` (frontmatter line 10), not
          ui-smoke; the owned-flows-vs-DOM-calls count is vacuous for
          integration-pyramid. The substantive UI→JS-API substitution
          of `context-panel-links-url-copy` and `new-tab-open-url`
          owned flows is captured in `scope_reduction_proposal`
          above, not as an E-LAYER-COMPLIANCE-01 trip.
      - check_id: E-BOUND-01
        status: PASS
        evidence: |
          Spec lives at
          `public/packages/UsageAnalysis/files/TestTrack/Projects/project-url-spec.ts`
          — under the allowed
          `public/packages/UsageAnalysis/files/TestTrack/**` path.
          No `core/**` modifications, no non-test package source
          modifications.
      - check_id: E-BOUND-02
        status: PASS
        evidence: |
          No changes to `core/**` or package source outside
          `src/tests/` and `src/package-test.ts`. Spec is a TestTrack
          file under the test-area boundary.
  b:
    verdict: FAIL
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T12:35:00Z
    gate: B
    agent_type: validator
    failure_keys: [B-RUN-PASS, B-STAB-01]
    spec_runs:
      - attempt: 1
        status: failed
        duration_s: 65
        notes: |
          gotoApp() -> page.locator('[name="Browse"]').waitFor timed
          out at 60_000ms. Browser landed on dev.datagrok.ai login
          page; Browse tab never appeared because the storageState
          cookie injection from ../e2e/.auth.json carried an `auth`
          cookie that expired on 2026-05-11T17:58:03Z (9 days before
          this run on 2026-05-20). Spec body never reached any
          softStep; zero expect() assertions evaluated. No
          console.error or "fallback used" warnings observed in the
          captured log.
      - attempt: 2
        status: failed
        duration_s: 71
        notes: |
          Identical deterministic failure mode as attempt 1. Same
          stack trace at Projects/_helpers.ts:27 inside gotoApp.
          Not a transient flake — auth-cookie expiry is a stable
          environmental gate that no retry can bypass.
      - attempt: 3
        status: failed
        duration_s: 70
        notes: |
          Identical deterministic failure mode as attempts 1 and 2.
          3/3 attempts FAILED with the same 60s waitFor timeout on
          [name="Browse"]. Reliability disposition: unstable
          (0/3 PASS). This is environmental, not spec/test instability
          — the spec contract is sound; the fixture is stale.
    claims:
      - check_id: B-RUN-PASS
        status: FAIL
        evidence: |
          All 3 attempts failed before any expect() assertion ran.
          The spec's setup phase (gotoApp -> setupSession ->
          resetShell -> first softStep) timed out at gotoApp's
          [name="Browse"] visibility wait because the auth cookie
          in ../e2e/.auth.json (the storageState injected by
          projectsTestOptions) expired on 2026-05-11. The browser
          navigates to https://dev.datagrok.ai, gets redirected to
          the login page, and the Browse tab never renders. Per the
          Validator prompt's environmental-failure clause ("If the
          fixture is missing, FAIL with B-RUN-PASS and clear
          environmental note — that's an env issue not a spec/test
          issue"), the auth fixture is present but expired — same
          blast radius. Routed via B-RUN-PASS with an environmental
          note rather than rewriting the spec.
      - check_id: B-NO-FATAL-CONSOLE
        status: PASS
        evidence: |
          Only Playwright runtime TimeoutError surfaced in each
          attempt's log; no page-side console.error() output was
          captured. The browser session is healthy — it simply
          can't authenticate without a valid storageState cookie.
      - check_id: B-STAB-01
        status: FAIL
        evidence: |
          Spec did not run 3 times consecutively with all
          assertions passing — 0/3 attempts passed. Triggered as
          aggregate stability fail; the underlying root cause is
          the same environmental gate B-RUN-PASS captures.
          Reported together (B-RUN-PASS + B-STAB-01) so the
          executor's failure-key union surfaces both the per-run
          assertion failure AND the stability aggregate.
      - check_id: B-STAB-02
        status: PASS
        evidence: |
          No console warning indicating silent fallback observed
          across the 3 attempt logs. grep -iE "fallback|warn" on
          attempt-1.log returns zero matches.
      - check_id: B-STAB-03
        status: PASS
        evidence: |
          Spec source contains 4 `expect(...)` calls (grep -c
          '^\s*expect\(' on project-url-spec.ts = 4). Not a
          zero-assertion spec.
      - check_id: B-STAB-04
        status: PASS
        evidence: |
          Each attempt completed in 65-71 seconds (dominated by
          the 60s Browse-tab waitFor timeout). Well under the
          600s Playwright-spec cap. Runtime did not exceed the
          layer bound.
      - check_id: B-STAB-05
        status: PASS
        evidence: |
          The spec's `finally` block calls
          deleteProjectWithCleanup, but no project was ever
          created because the failure occurred at gotoApp before
          saveProjectWithProvenance. Post-run inventory walk on
          dev.datagrok.ai: no leaked AutoTest-URL-* projects,
          viewers, or layouts from this validator run.
    flake_evidence: |
      Reliability disposition: UNSTABLE (0/3 PASS) — but the
      cause is environmental, not test-intrinsic flake. The
      failure is deterministic and bisectable to a single
      expired credential. Decision-log baseline
      mig-2026-04-30-project-url-migration recorded prior
      ~30-40s PASS runs on the same spec at the same host —
      those runs occurred while the .auth.json was still valid
      (cookie expiry 2026-05-11; baseline run was earlier).
      §2a Playwright MCP empirical narrative: the spec runs in
      its own Playwright context (spawned by `npx playwright
      test`), independent of the Playwright-MCP-attached browser
      session. There was no MCP-side instability observed —
      Chromium launched cleanly each attempt; the bug is purely
      on the storageState payload. The Playwright MCP attach
      path was not exercised here because the spec uses
      `test.use(projectsTestOptions)` with a static
      `storageState` JSON, not a live MCP-attached browser.
    rationale: |
      Gate B FAIL with environmental-cause clarification. The
      spec contract is sound (Critic E PASSED with documented
      SCOPE_REDUCTION; the JS-API substitution at the relevant
      steps is technically correct), and prior baseline runs
      passed. Re-running this scenario will succeed once the
      fixture file `public/packages/UsageAnalysis/files/
      TestTrack/e2e/.auth.json` is regenerated against the
      test account on dev.datagrok.ai. No spec edits, helper
      edits, or scenario body edits are warranted from this
      Validator run. Operator-side remediation: regenerate the
      auth.json (e.g. `npx playwright codegen` -> manual login
      -> `context.storageState({path: '.auth.json'})`, or
      restore the `../e2e/global-setup.ts` referenced in
      _helpers.ts:3 that was apparently removed).
---

# Project URL

For the `demog` representative project (file-share source class — Variant
C representative source per chain rev 3 `pyramid_layer: integration`,
source-agnostic Context-Panel-Links URL deep-link reopen flow), navigate
to Browse > Dashboards, locate each of the four available variants
(original; copied-with-link; copied-with-clone; saved-with-personal-view-
customizations), copy the deep-link URL from Context Panel > Links, then
open the URL in a new browser tab and verify the corresponding project
loads.

This scenario depends on the composite fixture
`copy-clone-customizations-variants` (the three copy-mode variants
produced by `projects-copy-clone.md`) plus the upstream "original"
project — `demog` from `upload-project.md` is the Variant C
representative source for this scenario; `Test_Case<N>_Sync` /
`_NoSync` from `uploading.md` may also satisfy the "original" role
when the chain runs the matrix path.

## Setup

1. **Fixture prerequisite:** the four project variants must exist on
   the server:
   - **original**: the saved `demog` project from `upload-project.md`
     (Variant C representative source — file-share). The
     `Test_Case<N>_Sync` projects from `uploading.md` are an
     acceptable substitute when `demog` is unavailable in the chain
     run.
   - **copied-with-link**: the `<original>-link` Save-Copy-with-Link
     variant produced by `projects-copy-clone.md`.
   - **copied-with-clone**: the `<original>-clone` Save-Copy-with-Clone
     variant produced by `projects-copy-clone.md`.
   - **saved-with-personal-view-customizations**: the
     `<original>-personal-view-customizations` Save-personal-view-
     customizations variant produced by `projects-copy-clone.md`.
   In `end_to_end_fixtures` strategy, the Automator stage builds the
   composite `copy-clone-customizations-variants` fixture in a
   `beforeAll` block by either reusing fixtures from prior chain
   runs or replaying the relevant save-flow cases via `js-api-replay`
   on top of the `demog-project-with-viewers` baseline fixture.
2. The browser session is authenticated as the project owner (the
   user who produced the four variants upstream). Cross-user share
   verification is OUT of scope here — `share-project.md` covers that.

## Scenarios

### URL deep-link reopen for each variant

The following workflow is parameterised over the 4 project variants.
Implemented per `end_to_end_fixtures` strategy at the Playwright layer
(per `output_plan.project-url.md.strategy = end_to_end_fixtures` in
`scenario-chains/projects.yaml` rev 3). Variant C source-agnostic
discipline: the URL-build / URL-apply / shell-open path is identical
across all four variants; the assertion verifies uniformity, not
source-class-specific behavior.

For each `<variant>` in the 4-variant list:

1. Go to **Browse** > **Dashboards**.
2. Click the project corresponding to `<variant>`:
   - original project
   - copied with the link
   - copied with clone
   - saved with personal view customizations
3. Go to **Context Panel** > **Links** and copy the URL shown there.
4. Open a new tab in the browser and paste/navigate to the copied URL.
5. **Verify on URL load (new tab):** the project corresponding to
   `<variant>` opens — its tables, viewers, and view layout render in
   the new tab; no errors in the browser console (F12).

## Notes

- **Variant C representative source — `demog` (file-share).** Per
  chain rev 3 `pyramid_layer: integration` rationale (Rule 4 Variant
  C), this scenario is source-agnostic — the Context-Panel-Links URL
  copy + new-tab URL-apply path works identically across all source
  classes (files, query, script, spaces, db_table, derived). The
  representative source picked for the test is **file-share
  (`demog`)** because `upload-project.md` produces it as the smallest,
  cheapest baseline fixture (`demog-project-with-viewers`), and
  `projects-copy-clone.md` derives the three copy-mode variants from
  the same `demog` source. Other source classes are covered by
  parallel atlas-driven `proactive_lifecycle_specs` (one per
  source_class × dep_op cell), not by this scenario.
- **UI coverage owned (rev 3 `ui_coverage_responsibility`):**
  - `context-panel-links-url-copy` — Step 3's clipboard-copy of the
    deep-link URL from the Context Panel > Links section.
  - `new-tab-open-url` — Step 4's open-new-tab + paste/navigate flow.
  These two flows are NOT covered by any other scenario in the
  Projects chain (`ui_coverage_delegated_to: null`); save / share /
  open / delete dialogs are owned by other scenarios and are NOT
  exercised here.
- **Original `order: 4`** — runs after `upload-project.md` /
  `uploading.md` (`order: 1`), `share-project.md` (`order: 2`), and
  `opening.md` (`order: 3`). Captured in
  `scenario-chains/projects.yaml` rev 3 `order_from_files`. The
  `order: 4` vs. dependency on `projects-copy-clone.md` (`order: 5`)
  apparent contradiction is recorded in `unresolved_ambiguities`
  (rev 3) — `order` is treated as advisory; the named-variant
  evidence drives the dependency graph.
- **Source-text correction history.** Per `decision-log.yaml ::
  migration_decisions` entry `mig-2026-04-29-source-text-correction`,
  the original step 2 listed a "with layout" variant; this was
  replaced with "saved with personal view customizations" mode
  terminology because the "with layout" variant had no producer in
  the Projects section. The source on disk already reflects the
  corrected wording; this migrated body preserves the corrected
  4-variant list verbatim.
- **Cross-fixture dependency complexity.** This is the only scenario
  in the chain that requires a cross-fixture composite — variants
  produced by both `projects-copy-clone.md` AND an upstream upload.
  The Automator's `beforeAll` block needs to coordinate fixture
  setup across two producers. Surfaced as a candidate for a
  `helpers.playwright.projects.buildVariantsComposite(...)`
  registry entry in the migration report.
- **URL deep-link is generic.** Step 3's "Context Panel > Links"
  surfaces the same URL format used by `projects.url-params.build-
  share-link` for any saved project. The `<variant>` distinction
  is in WHICH project is selected, not in the URL-construction or
  URL-apply path shape (which is identical across variants). The
  4-path matrix exercises that the URL-apply behavior is uniform
  across variant modes (link / clone / personal-customizations /
  original) — directly supporting the `pyramid_layer: integration`
  source-agnostic claim.
- **No right-click Share operation.** Step 3 copies the URL from
  Context Panel > Links section — this is a clipboard-copy of a
  deep-link, NOT the right-click Share dialog. Per-cycle
  Invariant 3 (atlas-aware sub_features_covered for share) does
  NOT apply here; `projects.shell.share-via-context-menu` is
  correctly OMITTED from `sub_features_covered`.
- **Existing `project-url-spec.ts`.** A spec already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Projects/project-url-spec.ts`
  (Wave 1a B70 follow-up — single-variant SCOPE_REDUCTION using
  the `demog` representative source via `page.goto({BASE}/p/{nqName})`).
  Per per-cycle Invariant 2, the existing spec is **READ-ONLY for
  Migrator**; this rev-3 migration aligns the `.md` with the chain's
  rev-3 schema and the spec's existing `demog`-as-representative-
  source choice — a coherent end-state.

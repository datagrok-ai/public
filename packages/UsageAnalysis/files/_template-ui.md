---
feature: <section id, e.g. connections>
realizes_atlas:              # REQUIRED — without it this file counts for NOTHING.
  - <feature>.cp.<id>        # isScenarioGreenFor checks realizes.includes(type.id)
  - <feature>.int.<id>       # BEFORE it ever looks at manual-only disposition
                             # (.claude/executor/src/helpers/compute-coverage.ts:78-81,
                             # built from realizes_atlas at :236)
realizes: []                 # KG leaf slugs, NOT atlas ids. [] when the atlas
                             # entry's validates_kg is empty.
priority: p0|p1|p2           # REQUIRED. Coverage carrier alongside realizes_atlas;
                             # also read by Gate A/E mechanical checks and hook h4
                             # (subagent-spawner.ts:310, REQUIRED_MIGRATION_FRONTMATTER_FIELDS)
target_layer: manual-only    # EXACTLY this string — "manual" does not match
                             # (compute-coverage.ts:266 compares literally)
coverage_type: smoke|regression|edge|perf
companion_to: <parent>.md    # when this file was split off a parent scenario.
                             # The parent must carry ui_coverage_split_to back;
                             # the spawner enforces both directions
                             # (subagent-spawner.ts:1437-1452)
manual_only_reason: |
  Why this cannot be automated — a real technical reason (visual/perceptual
  judgment, OS-level interaction like a file picker or drag-and-drop without
  DOM targets, out-of-band credentials).
  NOT a valid reason: "no time to automate"; "the product hangs on dev"
  (operator ruling 2026-08-20 — a hang is an open bug plus a red test, and the
  scenario still gets automated); "the reference doc does not describe the
  selector" (that is a recon task, not an automation limit).
automation_candidate: true|false   # must be consistent with manual_only_reason:
                                   # a perceptual-judgment reason implies false.
                                   # Omit the field rather than contradict it.
related_bugs: []             # - id: GROK-NNNNN / status: <state>
blocked_by: []               # ticket id ONLY while the scenario is unrunnable
---

# <Scenario name>

One or two sentences: what this file covers and what is deliberately NOT here
because the paired spec(s) of this section cover it.

## Preconditions

- Named datasets, connections, accounts, secrets, containers — everything the
  tester needs before step 1. Fixture files live in the repo or Demo Files and
  are referenced by full path; a file that exists only on a QA machine is
  explicitly marked "QA-owned" with a note on how to obtain it.
- No stand names unless the scenario is genuinely stand-specific.

## <Scenario section — one per flow>

1. Concrete action — concrete, observable expected result
2. Every created entity gets an explicit name (`test_<feature>_<what>`) and is
   reused by that name; never "the one from the previous step"
3. No "check all ...", "any dataset", "arbitrarily", "as expected", and no
   exact error-message texts — describe the error class instead
4. No pipeline/CI jargon (gate verdicts, atlas ids, spec line numbers) in the
   body; that belongs in frontmatter if anywhere
5. State what a tester should see when the step FAILS, not only when it passes —
   a manual step whose only description is the success path gives the tester no
   way to tell "broken" from "different"

## Cleanup

Delete everything the scenario created — projects, connections, layouts,
downloaded files — including after a partial run ("delete whichever of X/Y/Z
exist").

---
{
  "order": N,
  "datasets": ["System:DemoFiles/demog.csv"]
}

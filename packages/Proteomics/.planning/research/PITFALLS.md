# Pitfalls Research — v1.4 Cross-Team Review (Proteomics)

**Domain:** Datagrok analytics package adding read-only publishing + SPC monitoring + cross-run comparison on top of a tag-coordinated, in-place-mutating pipeline (parse -> annotate -> normalize -> impute -> DE -> enrichment).
**Researched:** 2026-06-06
**Confidence:** HIGH for codebase-grounded pitfalls (file/line cites + 999.1/999.4 + Phase 13 round-3 evidence); MEDIUM for Datagrok project-serialization surface (search results clear on layout/tag persistence but the v1.3 round-3 stripping bug shows the contract is partial); MEDIUM for SPC statistical traps (textbook + LC-MS literature consensus).

This document is opinionated. Where the codebase already has scar tissue for a class of bug (e.g. 999.1 silent data loss when a dialog re-reads stale column refs), the same shape *will* recur in v1.4 surfaces unless the prevention pattern is ported. Every pitfall ends with a concrete test/CLAUDE.md/UAT lever.

---

## Critical Pitfalls

### Pitfall 1 — "Frozen" published copy silently mutates because publish forgot to deep-clone the DataFrame

**What goes wrong:**
Proteomics expert runs `Proteomics | Share Analysis for Review...`, the publish handler builds a "trimmed" view by `df.columns.remove(...)` or by setting tags, then serializes that *same* `DG.DataFrame` into a project. A week later the expert re-runs Normalize -> Impute -> DE on the source DataFrame to try a different threshold. The biologist re-opens the "frozen" project and sees different log2FC, different p-values, or — worse — a half-mutated state because the published reference still points at the live in-memory DataFrame, and the package's pipeline mutates in place. The biologist signs off on numbers the expert never published.

**Why it happens:**
- Architecturally established: "All analysis steps mutate the original `df` in-place" — `.planning/codebase/ARCHITECTURE.md:197`. `df.columns.remove(...)` and `df.setTag(...)` mutate the existing instance.
- The package already has *one* clone-for-isolation precedent (`createExpressionHeatmap` uses `df.clone(filter)` to isolate top-N sort — `CLAUDE.md` heatmap rule). The pattern is local to one viewer, not a publishing pipeline.
- Re-run guards exist for DE (`proteomics.de_complete` short-circuit at `src/analysis/differential-expression.ts:242`) but NOT for normalize/impute — those happily overwrite the log2 intensity columns and any derived columns in place. Anti-pattern documented in `ARCHITECTURE.md:210-214` ("Mutation of the In-Place DataFrame Without Idempotency Guard").

**How to avoid:**
- Define a `publishAnalysis(df, opts)` helper in a new `src/analysis/publish.ts`. The first line is `const frozen = df.clone(undefined, trimmedColumnNames);` — clone before doing ANYTHING else. Never publish the live `df`.
- Re-set every `proteomics.*` tag on `frozen` explicitly (do not assume `df.clone` carries tags — see Pitfall 3). Add a single `proteomics.published_at = <ISO date>` and `proteomics.published_from_de_method = <method>` to the frozen copy so the reviewer always sees provenance.
- On the trimmed clone, also set `frozen.name = ${df.name}_published_<date>` so the project's table name visibly differs from the source.

**Warning signs:**
- Reviewer says "this volcano looks different from what you showed me last week" -> the leak already happened.
- Expert sees the same dataset in TWO open table views and edits to one show up in the other -> the publish path is sharing a reference, not cloning.
- Re-opening the published project shows a row count > the trim contract specifies -> the project saved a column-removed VIEW of the live df, not a clone.

**Phase to address:** Publishing-phase (first task: build the publish primitive). This is the load-bearing decision the rest of v1.4 publishing rests on.

**Test/UAT lever:**
- Unit test: `publishAnalysis(df, opts)` returns a frozen DF whose mutation does NOT propagate back. Pseudocode: `const frozen = publishAnalysis(df, opts); const before = frozen.col('log2FC').getRawData()[0]; df.col('log2FC').set(0, 999); expect(frozen.col('log2FC').getRawData()[0], before);`
- UAT scenario "Stale snapshot smoke test": import demo, run DE, publish, RE-RUN DE with different fc threshold on source, re-open published project, assert log2FC values match the originally-published numbers (not the re-run numbers).

---

### Pitfall 2 — Reviewer group accidentally gets Edit permission on the source DataFrame, not just the trimmed clone

**What goes wrong:**
Publish dialog asks "which reviewer group?", calls `grok.dapi.permissions.grant(project, group, 'View')` on the project, but the *source* DataFrame (or the source project that hosts it) already shares a parent space with the reviewer group at Edit level — through a Datagrok permission-inheritance side effect ("When an entity is moved to a new project, it adopts the permissions of the new project" — Datagrok docs). The reviewer can now click into the source view, re-run DE, normalize differently, or worst case modify columns the expert depended on. The expert never sees a permissions error.

**Why it happens:**
- Datagrok's permission model is per-Entity AND per-Space; granting "View" on the trimmed project does NOT revoke or shadow inherited permissions on the source DF or on the parent Space the project lives in.
- The package currently has zero `grok.dapi.permissions` calls — there's no muscle memory in the codebase for explicit revoke.
- "Uploading the project object doesn't automatically share it with others — you must share it explicitly" — but the converse (a project inheriting Edit because the parent Space already grants it to the group) is not blocked by the publish flow.

**How to avoid:**
- Publish ONLY the trimmed clone (Pitfall 1) as a new entity owned by the publishing user, into a target-named Space the reviewer group has VIEW-ONLY on. Never publish the source DF reference.
- After save, explicitly enumerate effective permissions on the published project: `await grok.dapi.permissions.get(publishedProject)` and assert the reviewer group's privilege list contains `view` only, never `edit`/`delete`/`share`. Throw a hard error if Edit is found and roll back the publish.
- Refuse to publish into a Space the user does not own at admin level (avoids "I published into a Space someone else manages, who already gave biologics Edit").

**Warning signs:**
- Reviewer in UAT clicks a column header and gets a context menu offering "Change type" / "Delete column" -> they have Edit.
- Reviewer can save changes back to the project -> permissions are not view-only.
- `grok s shares add --access Edit` appearing in scripts/CI for the reviewer group -> the contract is being violated upstream.

**Phase to address:** Publishing-phase — permission verification is part of the publish primitive, not a separate "polish" task. Goes in lockstep with Pitfall 1.

**Test/UAT lever:**
- Programmatic test using `grok s` (per `tools/GROK_S.md`): create a test reviewer group, publish a fixture analysis, then run `grok s users save --json reviewer.json` to authenticate as a member of that group and attempt `dapi.tables.save(publishedTable)` — must throw a 403/permission error.
- UAT scenario "Reviewer cannot mutate": authenticated as a biologist demo user, open the shared project, try to (a) edit a cell, (b) delete a column, (c) re-run DE from the top menu — all three must fail or be hidden.

---

### Pitfall 3 — Project serialization drops `proteomics.*` tags / semantic types / df.name, so the reopened project looks broken to the reviewer

**What goes wrong:**
Expert publishes. Reviewer opens the project a week later. The volcano view is empty, or shows the raw column names instead of "Log2 Fold Change", or the Filters viewer doesn't dock its 3-column scope (Comparison / Display Name / Source ID per Phase 13 D-07), or the UniProt panel doesn't auto-fire because the `Proteomics-ProteinId` semType is gone. The reviewer concludes the analysis is "corrupt" and the demo dies.

**Why it happens:**
- The package's entire workflow contract is encoded as `proteomics.*` DataFrame tags + `SEMTYPE.*` column tags. There is no separate state store. If the serializer drops either category, downstream viewers fail their precondition tag checks (e.g. `createExpressionHeatmap` re-checks `proteomics.de_complete` at the viewer level per `ARCHITECTURE.md`).
- Phase 13 round-3 (commit `e527d07ba1`) already documented one case where the Datagrok serializer strips `look.filters[]` and `look.columnNames` from a docked Filters viewer regardless of which shape is fed in. The strip is real, partial, and not fully documented. Quote from 14-HUMAN-UAT.md: "the platform serializer strips both look.filters[] and look.columnNames regardless of which shape is fed in".
- Datagrok docs say tags persist in projects ("formula is preserved in DataSync projects" — Datagrok help), but the codebase has shown the contract is partial for at least viewer-config tags. There is no test in this package that round-trips a project save/reopen and re-checks every tag.

**How to avoid:**
- Build an `assertPublishedShape(reopenedDf)` helper that the publish UAT test imports. It asserts (a) every `proteomics.*` tag the trim contract requires, (b) every `SEMTYPE.*` on the column it expects, (c) `df.name`, (d) the volcano viewer's docked layout matches the saved one.
- Save the project, in the SAME test, immediately reload it via `grok.dapi.projects.find(id).open()` and call `assertPublishedShape`. Do not assume "save succeeded" means "reopen will look the same".
- For any tag the serializer drops, redundantly encode the value in a column (e.g. add a single-row metadata column "DE method" instead of relying only on a tag). Belt-and-braces: tag for runtime use, column for serialization survival.
- Mirror SEMTYPE constants in `detectors.js` (per existing rule in `CLAUDE.md` and `packages/Proteomics/CLAUDE.md`) — so if semType is dropped on reopen, the detector will reassign on next view load.

**Warning signs:**
- Reopen test fails an assertion -> serializer dropped something.
- Reviewer reports "the volcano is blank" or "the gene names aren't bold/marked" -> a column tag or SEMTYPE is missing.
- UniProt context panel does not auto-fire on click -> `Proteomics-ProteinId` semType lost OR the detector didn't re-run.
- df.name displays as a generic "table1" instead of the expected `<source>_published_<date>` -> df.name didn't survive.

**Phase to address:** Publishing-phase (gating: cannot ship publish until the round-trip test is green) AND a one-time investigation task at the front of the phase to enumerate WHICH tags/semtypes the serializer keeps vs strips.

**Test/UAT lever:**
- Unit test (new file `src/tests/publish-roundtrip.ts`): `publishAnalysis`, save-as-project, immediately reopen, assert every required `proteomics.*` tag + every `Proteomics-*` semType + `df.name` + volcano layout.
- UAT scenario: close the published project, fully refresh the browser (close-before-refresh per saved Memory `feedback_datagrok_close_before_refresh`), reopen from URL, assert the demo shape.
- Document in `CLAUDE.md` (Publishing section): "Tags and semTypes are NOT trusted to survive project save. Always assert post-save shape via the publish roundtrip helper."

---

### Pitfall 4 — Versioning ambiguity: "which version of this analysis is in the published project?" has no answer

**What goes wrong:**
Expert publishes version A. Two days later expert reruns DE with different thresholds and republishes (or worse, "updates" the same target). Reviewer sees one volcano. Expert sees another. Nobody knows which numbers are canonical, which method was used, or how to find the older version. The "let me re-run with different parameters" conversation captured in `src/analysis/differential-expression.ts` already exists today as a UX gap (`2026-03-05-reset-de-analysis-to-allow-re-running-with-different-parameters.md`); without explicit versioning, it's a sign-off-risk in v1.4.

**Why it happens:**
- The pipeline tag set carries one method (`proteomics.de_method`) but no version, no publish ID, no superseded-by pointer.
- Datagrok projects can be re-saved under the same name (no automatic immutability). Without a publish-side convention, "republish" is destructive.
- DE re-run guard is a hard early-return — it doesn't even cycle a version counter. Re-run is an all-or-nothing in-place overwrite.

**How to avoid:**
- The publish primitive (Pitfall 1) MUST set: `proteomics.publish_id = <uuid>`, `proteomics.publish_version = N`, `proteomics.published_at = <ISO>`, `proteomics.published_by = <user>`. Encoded as both tags AND a single-row metadata column (Pitfall 3 belt-and-braces).
- Republish creates a NEW project (new id), not an overwrite. Old project gets `proteomics.superseded_by = <new_id>` injected and optionally hidden (not deleted — keep audit trail).
- The reviewer-facing view (volcano title or a small audit panel) MUST surface `published_at`, `publish_version`, `de_method`, `fc_threshold`, `p_threshold` so the reviewer sees what they're looking at without asking.
- Rollback: "supersede" is a soft pointer, never a delete. Recovery from "I published the wrong numbers" is to publish a corrected version with a higher `publish_version`.

**Warning signs:**
- Expert opens UAT and asks "wait, which one of these published volcanos is the current one?" — version is not in the title.
- A `dapi.projects.save` overwrites a published id silently — destructive republish path exists.
- Reviewer's bookmark / link breaks after republish — the new version got the same URL but different content.

**Phase to address:** Publishing-phase. Versioning IS part of the publish contract; don't defer it to "v1.5 cross-team review polish".

**Test/UAT lever:**
- Test "republish creates new project, marks old as superseded": publish v1, republish, assert `projects.find(v1).getTag('proteomics.superseded_by') === v2.id` and v1 is NOT modified beyond that tag.
- UAT: republish, open both v1 and v2 side by side, confirm version numbers + de_method visible in each title bar.

---

### Pitfall 5 — Per-protein SPC: 5,000 control charts per run with no multiplicity correction blows alarm rate to nearly 100%

**What goes wrong:**
v1.4 calls for "sample-level SPC tracking across weekly assays (intensity median, missingness %, control-correlation, protein count) with Shewhart/Nelson drift rules". Misreading this as per-protein SPC produces thousands of charts. Even ignoring the UX cost, with 8 Nelson rules at a combined ~2.65% false-alarm rate per chart (per ScienceJourn paper / Nelson 1984), 5,000 independent charts yield ~133 alarms per run from noise alone. Reviewers will dismiss the SPC tab as a noise generator.

**Why it happens:**
- The roadmap entry explicitly scopes "sample-level" — but "sample-level intensity median" can be confused with "per-protein median across samples" because the package's pipeline thinks in proteins × samples. Implementer reads tags + columns and accidentally picks the protein axis.
- Nelson rules false-alarm rates compound across rules (Western Electric rules 1-8 sum to ~2.65%; Nelson 1984 noted his calibration assumes underlying Normality, which proteomics intensities never satisfy without log2 — and even then, tails are gamma-like).
- SPC literature shows non-normal distributions inflate false-alarm rates further; LC-MS intensity distributions are right-skewed/log-normal.

**How to avoid:**
- LOCK the SPC metric set to RUN-LEVEL (one number per run): median intensity across selected QC samples, overall missingness %, control-sample-vs-control-sample correlation, total protein count. Four metrics per run; per-run alarms remain interpretable.
- DO NOT enable per-protein SPC. If users ask, capture it as a deferred idea — it is a separate phase with its own multiplicity story (FDR control, hierarchical models).
- Apply Nelson rules conservatively: default to rules 1 (3-sigma) and 5 (2-of-3 beyond 2-sigma) only, with rules 2-4/6-8 behind an "advanced" toggle. Document the expected false-alarm rate per rule in the dialog.
- Use the log2-transformed intensity values for normality assumptions to hold approximately (mirrors what the rest of the pipeline already does post-parse).

**Warning signs:**
- A "by protein" picker on the SPC chart — scope creep.
- The SPC tab shows hundreds of "flagged runs" out of 10 -> false-alarm explosion.
- Implementer copies the QC dashboard pattern (which is per-sample, not per-run-over-time) -> wrong axis.

**Phase to address:** SPC-phase (foundational scoping). Add a "metric set is run-level, NOT per-protein" decision to the phase context document.

**Test/UAT lever:**
- Unit test: simulate 100 runs of in-control data (all 4 metrics drawn from same distribution), run Nelson rules 1+5, assert total alarms < 5 (within nominal false-alarm rate).
- UAT: a reviewer can describe what each of the 4 SPC charts represents in one sentence per chart.

---

### Pitfall 6 — Baseline contamination: control limits computed from runs that include out-of-control events make the baseline too wide

**What goes wrong:**
SPC implementer computes mean + 3-sigma from "the first 20 weekly runs" or "all runs to date" as a rolling reference. If one of those runs had a real special-cause event (LC column degradation, sample-prep error), its outlier value inflates the SD, so the control limits widen and future genuine drifts no longer trigger. Customer learns the SPC chart "never flags anything", loses trust, and stops using the tab.

**Why it happens:**
- Naïve mean + SD computation treats all baseline observations as in-control by definition.
- The proteomics expert is also the SPC operator — there's no separate person curating the baseline.
- Phase 14 patterns show implementer comfort with "compute on whatever is in the DataFrame" (e.g. the original kNN baseline was the entire row set). The same instinct applied to baseline construction is destructive.

**How to avoid:**
- The baseline is an EXPLICIT, ANNOTATED subset, not "all runs". Surface a UI: "Use runs [list] as the baseline" with checkboxes. Default to first N runs but require explicit user confirmation before computing limits.
- Phase 1 of baseline computation = iterative outlier removal: compute mean+SD, drop points outside 3-sigma, recompute, repeat until stable (or cap at 2 iterations). Document this is happening; show before/after on the dialog.
- Once a baseline is locked, persist it as a separate `proteomics.spc.baseline.<metric>` tag (run IDs + mean + SD) on the SPC DataFrame. Future runs compare against the LOCKED baseline, not a rolling recomputation.
- Provide an explicit "rebuild baseline" action that requires the user to select which runs to include.

**Warning signs:**
- A "rolling 20 runs" or "all runs" baseline option — the contamination path is wide open.
- The dialog has no "exclude this run from baseline" affordance.
- After 30 weeks of operation, no run ever flags — the baseline has absorbed every drift event.

**Phase to address:** SPC-phase.

**Test/UAT lever:**
- Unit test: baseline of [normal, normal, normal, normal, normal, 3-sigma outlier] returns limits computed from the first 5, not all 6.
- UAT: after 10 simulated runs with a known week-7 step change, the SPC chart flags week-7 onward as out-of-control (not the baseline).

---

### Pitfall 7 — Run-order ambiguity: calendar date vs assay-sequence position vs file-mtime produce different alarms from the same data

**What goes wrong:**
SPC charts are run-ordered. Implementer picks "calendar week of import" because that's what the import handler knows. Cytokinetics operates two LC-MS instruments in parallel; week-2 on instrument B was actually run BEFORE week-1 on instrument A. The SPC trend now interleaves two instruments. A Nelson "9 points on one side of centerline" rule fires from interleaving artifact, not biology.

**Why it happens:**
- The package today has no concept of "run" beyond `df.name` and the import datetime. There's no explicit per-run identifier separate from the DataFrame.
- LC-MS QC literature is universal: SPC must be ordered by ACQUISITION TIMESTAMP, not import time, and stratified by INSTRUMENT.
- The roadmap calls for "weekly assays" — "weekly" already encodes a calendar assumption.

**How to avoid:**
- Define the SPC run identity contract explicitly: `(instrument_id, acquisition_datetime)` is the primary key. Both must be captured at import or annotation time.
- The Annotate step (or a new Annotate Run step) prompts for instrument id and acquisition datetime. Default values: instrument blank, datetime = file-creation timestamp from the parser metadata. User confirms before proceeding to SPC.
- Order SPC charts by acquisition_datetime within instrument; render one chart per (instrument, metric) pair OR a single chart with per-instrument color encoding.

**Warning signs:**
- SPC chart x-axis labeled "Week 1, 2, 3..." with no instrument facet -> the bug is already there.
- Two runs in the same calendar week appear as one point -> aggregation is hiding the conflict.
- "Run 5" and "Run 6" are off the same instrument but the user inserted run 5b later from a backfill -> ordering broke.

**Phase to address:** SPC-phase (data model task — comes BEFORE the chart task).

**Test/UAT lever:**
- Unit test: insert run-X with acquisition_datetime = week 0, then insert run-Y with acquisition_datetime = week -1 (backfill). Assert the SPC chart orders Y before X, not in insertion order.
- UAT: two instruments, simulated parallel runs, assert two SPC traces visible and non-interleaved.

---

### Pitfall 8 — Storage cost growth: weekly runs accumulate as full DataFrames, hitting Datagrok project size limits within a year

**What goes wrong:**
SPC implementer stores every week's full proteomics matrix (5k proteins × 24 samples × 4 metric columns) in the SPC project. After 52 weeks the SPC DataFrame has 260,000 rows × wide columns, project size grows past Datagrok's reasonable single-project size, and the reviewer's "open SPC" click takes 30+ seconds.

**Why it happens:**
- The package's only persistence model today is "tag the DataFrame", which doesn't account for time series.
- It's tempting to keep the full per-protein quantification "in case we want to drill down" — but the SPC scope is run-level metrics only (per Pitfall 5).

**How to avoid:**
- SPC project stores ONLY the per-run 4-metric rollup: one row per (instrument, run, week, metric). With 4 metrics × 52 weeks × 2 instruments = ~416 rows/year. Trivially small.
- Per-protein drill-down is NOT in the SPC DataFrame. Instead, the SPC row carries a `source_project_id` linking back to the run's original published analysis project; "drill down" opens that linked project.
- Document the SPC row schema in CLAUDE.md so future contributors don't enrich the SPC DataFrame with per-protein columns.

**Warning signs:**
- The SPC DataFrame's column count > 10.
- Implementer asks "should we store the quant matrix too, in case...?" — say no.
- Open-SPC click time grows linearly with weeks-of-history.

**Phase to address:** SPC-phase (schema task).

**Test/UAT lever:**
- Unit test: SPC DataFrame after N runs has exactly `N * metric_count * instrument_count` rows.
- UAT: simulate 52 weeks of SPC, assert SPC tab opens in <2s.

---

### Pitfall 9 — Mismatched protein populations across runs cause silent left-join data loss in cross-run comparison

**What goes wrong:**
Campaign comparison joins run-A (4,800 proteins identified) and run-B (4,200 proteins identified) for side-by-side volcano. Implementer uses an inner-join on protein ID -> 3,800 shared, 1,000 dropped from run-A, 400 dropped from run-B. The two volcanos now look "similar in shape" because the asymmetry is invisible. Reviewer makes biological calls on the visible 3,800 and never knows what was hidden.

**Why it happens:**
- LC-MS DDA proteomics has run-to-run identification stochasticity ("missing not at random") — well-known in the field.
- Default join semantics in `DG.DataFrame.fromColumns` patterns and the platform's join helpers default to inner-join behavior unless explicitly outer.
- The Phase 14 multi-contrast Filters viewer pattern doesn't translate one-to-one — that worked because rows were per-protein within ONE run; here rows are per-protein across MULTIPLE runs with different populations.

**How to avoid:**
- Use a FULL OUTER JOIN on protein ID across all runs in the comparison. Per-run columns get NaN where the protein wasn't identified in that run. The DataFrame's rowCount equals the union of all identified proteins.
- Add a "Quantified-In" column (string or array): list of run IDs where the protein was quantified. This makes the population mismatch a first-class signal, not a hidden one.
- The comparison viewer (side-by-side volcano) shows BOTH the quantified-in-both subset AND, on a separate tab/toggle, the quantified-in-only-one subset. Default view is union (full outer) with NaNs visible.
- Show counts in the comparison UI: "Run A: 4,800 / Run B: 4,200 / Both: 3,800 / A only: 1,000 / B only: 400" — surface the asymmetry, don't hide it.

**Warning signs:**
- Two side-by-side volcanos with identical point counts -> inner-join in disguise.
- "Why does this protein appear in one volcano but not the other?" from a reviewer -> the explanation must be discoverable in the UI, not require asking the expert.
- A column count of (n_runs * 4_metrics) but row count < union -> inner join.

**Phase to address:** Campaign-phase. First task of comparison data model.

**Test/UAT lever:**
- Unit test: two runs with [P1,P2,P3] and [P2,P3,P4] join to 4-row DataFrame with explicit NaN cells (not 2-row inner-join).
- UAT: comparison shows count panel matching A/B/Both/Aonly/Bonly.

---

### Pitfall 10 — Volcano overlay with different per-run p-value distributions misleads the reviewer

**What goes wrong:**
Run A has 100 samples; Run B has 8 samples. Run A's volcano has p-values spanning [1e-15, 1]; Run B's spans [0.001, 1]. Implementer overlays them on the same axes by stacking points or rendering one volcano "on top of" the other. The reviewer sees Run B's "weak" hits squashed below Run A's wide spread and concludes Run B has nothing — when in fact Run B is just power-limited.

**Why it happens:**
- DE p-values are sample-size-dependent; they are not directly comparable across runs with different N. Standard biostatistics gotcha.
- "Overlay" feels like the simplest comparison UX. It is also wrong.
- The Phase 13/14 volcano polish (Phase 14 D-04, magenta/cyan/gray lock) was tuned for SINGLE volcanos. Re-using the locked color scheme on overlaid runs erases per-run identity unless explicitly rebound.

**How to avoid:**
- Compare runs SIDE-BY-SIDE in linked viewers, never overlaid on the same axes. Two volcano viewers, x-axes synchronized by formula (so users can drag-zoom), y-axes independent (because p-value scale per run).
- Optionally, derive a "rank percentile within run" column for each run and compare percentile vs percentile — this is sample-size-normalized.
- For the "show me hits in both runs" use case, build a separate concordance viewer (e.g. log2FC vs log2FC scatter across runs, colored by sig-in-both / A-only / B-only / neither). This is a DIFFERENT viewer from the volcano.
- Per-run color identity: ditch the magenta/cyan/gray lock for cross-run comparison; use per-run hue.

**Warning signs:**
- "Just stack them on the same plot" in the design — refuse.
- Reviewer asks "why is run B so much weaker?" when the actual story is "run B has 8 samples".
- The cross-run viewer reuses the magenta/cyan/gray palette across runs.

**Phase to address:** Campaign-phase. Pin the viewer pattern (side-by-side, not overlay) early in the phase context document.

**Test/UAT lever:**
- UAT: a reviewer can correctly identify which run has more statistical power from looking at the comparison view (sample-count panel must be visible).
- Visual regression: side-by-side panels are independent viewers (`tv.viewers.length === 2` for the comparison view), not a single overlaid scatter.

---

### Pitfall 11 — Compound annotation drift: same compound is called different things across weeks and breaks cross-run grouping

**What goes wrong:**
Week 1 labels samples "CMP-001-10uM", week 4 labels them "Compound1-10micromolar", week 8 labels them "cmp1_10". Comparison-by-compound groups them as three separate compounds. Reviewer sees three independent "weak" trends instead of one strong eight-week trend.

**Why it happens:**
- v1.4 introduces compound as a NEW data model entity (per `2026-05-11-compare-most-impactful-compounds-across-weekly-screening-runs.md` and CONCERNS.md "No Compound Entity"). There is no existing canonicalization layer.
- LIMS / wet lab sample naming is inherently lossy. The package will see whatever the import file says.

**How to avoid:**
- Compound MUST resolve to a canonical identity at import time, not at comparison time. Preferred order: (a) explicit compound ID column (SMILES, InChIKey, internal ID), (b) Chem package canonicalization on SMILES via the existing Chem integration, (c) free-text name with a manual review step.
- New `Annotate Compound` dialog (separate from `Annotate Experiment` groups). Stores `compound_id`, `smiles`, `concentration`, `units`, AND the raw label so the resolution is traceable.
- If a name doesn't match an existing compound at import, prompt: "Is this a new compound or a different label for [list of close matches]?". Don't silently create a new entity.
- For the v1.4 demo at Cytokinetics, the safest path is to require an explicit `compound_id` column in the import — defer fuzzy name matching to a later phase.

**Warning signs:**
- Comparison view lists 47 compounds for an 8-week 12-compound campaign -> drift is creating phantoms.
- "Why isn't compound X appearing in the cross-run plot?" from a reviewer -> name in week 8 file doesn't match week 1 canonical.

**Phase to address:** Campaign-phase. Compound entity is foundational; build before the comparison viewer.

**Test/UAT lever:**
- Unit test: 3 runs labeled "CMP1-10uM", "Compound1-10micromolar", "cmp1_10" resolve to the same `compound_id` when annotated with the same SMILES.
- UAT: a reviewer can see, for a specific compound, all runs that included it and the labels used in each.

---

### Pitfall 12 — Batch effects masquerading as compound effects in cross-run comparison

**What goes wrong:**
Week 3's runs all show "compound X is highly effective". Week 7's runs show "compound X is unchanged". The reviewer concludes the compound has unstable activity. Actual cause: week 3 used a fresh LC column, week 7 used a degraded column; the entire week's intensity scale shifted. Compound X looked effective in week 3 vs in-week vehicle (which also shifted) because the intra-week comparison was internally consistent — but cross-week comparison mixes the batch effect into the compound effect.

**Why it happens:**
- Batch effects are universal in LC-MS proteomics and the package's current "Out of Scope" line ("Batch effect correction (ComBat) — different statistical procedure requiring batch metadata") explicitly does not handle this.
- Customer's mental model is "we ran the same protocol each week" -> they don't expect batch effects.
- Cross-run comparison surfaces this for the first time because intra-run analyses are normalized within run.

**How to avoid:**
- Require vehicle/control samples in EVERY run (already a SPC trust requirement; Pitfall 6 also depends on this). Cross-run comparison's first action is to subtract vehicle within run.
- Display vehicle-vs-vehicle SPC trend ALONGSIDE the compound effect trend in the campaign viewer. If vehicles drift between runs, the user sees the batch confound directly.
- Disclaim in the comparison UI: "Cross-run comparisons assume runs are batch-comparable. Vehicle drift across runs may confound compound effects. Review the vehicle-vs-vehicle SPC trace before interpreting."
- Do NOT silently apply ComBat or other batch correction — surfacing the problem is more honest than half-correcting it.

**Warning signs:**
- A reviewer asks "why is compound activity so variable week to week?" without checking the vehicle.
- The cross-run viewer has no vehicle-trend panel.

**Phase to address:** Campaign-phase.

**Test/UAT lever:**
- UAT: simulated 4-week campaign with a step shift in vehicle intensity between week 2 and 3; reviewer must be able to identify the batch shift from the comparison UI alone.

---

### Pitfall 13 — Cross-DataFrame selection wiring for campaigns: copy-pasting the v1.2 enrichment cross-link pattern silently fails

**What goes wrong:**
v1.2 wired enrichment -> volcano via `enrichDf.onCurrentRowChanged` subscription with results held in a module-level `activeSubscriptions: rxjs.Subscription[]` (see `ARCHITECTURE.md` "Cross-DataFrame Selection" + `Global state` constraint). v1.4 campaign implementer copies this pattern for run-A's volcano selection -> run-B's volcano highlight. It "works" once, then fails after re-open: the new view subscribes again, the old subscription is still in the array, both fire and one of them references a disposed viewer -> silent error or duplicate highlight.

**Why it happens:**
- The existing pattern was designed for ONE producer (enrichment) -> ONE consumer (protein DF). Campaigns introduce N runs that all both produce AND consume selection signals — it's a graph, not a path.
- The module-level subscription array has no per-viewer eviction; "must be unsubscribed before re-opening" is documented but not enforced (`ARCHITECTURE.md:196`).
- v1.2 already shipped an "enrichment cross-link not firing" debug session (`.planning/debug/enrichment-cross-link-not-firing.md` exists) — this surface has scar tissue.

**How to avoid:**
- Build a new `CampaignSelectionBus` abstraction in `src/analysis/campaign-selection.ts`. Single producer-of-truth: ONE bus per campaign, all run viewers subscribe to it. Each viewer registers/unregisters via a lifecycle hook bound to viewer disposal.
- When a viewer is disposed (table view closed, viewer detached), its subscription auto-evicts. Use `viewer.onDispose` or equivalent.
- Reuse the Phase 14 "highlight not hide" pattern (D-05 / Phase 14 14-VERIFICATION): writes go to `df.selection`, NOT `df.filter`. Keep the NS cloud visible on every run's volcano so the reviewer sees the unselected context.
- Build a debug "show me the active subscriptions" devtool for the bus from day one; verify subscription count returns to 0 after closing the campaign view.

**Warning signs:**
- Closing a run viewer and reopening produces double-firing of selection handlers -> stale subscription.
- A "this is fine, we'll clean up later" inline comment near the subscription registration -> destined to break.
- Module-level mutable state with no eviction policy.

**Phase to address:** Campaign-phase (data-binding layer, before the viewer layer).

**Test/UAT lever:**
- Unit test: campaign with 3 runs, close one viewer, assert subscription count drops by exactly 1; close the campaign view, assert count returns to 0.
- UAT: open campaign, close all viewers, re-open campaign -> no duplicate highlights.

---

### Pitfall 14 — Cytokinetics demo audience contains biologists; jargon and silent-loading viewers kill the demo

**What goes wrong:**
v1.4 demos at Cytokinetics with both proteomics experts AND biologist consumers in the room. Phase 14 already invested in analyst experience (gene labels, magenta/cyan/gray semantic colors, smart pathway filter). The publishing flow uses jargon ("permissions", "ACL", "DataFrame", "tag") and the published view has hung viewers during loading because the package has no progress indicator on viewer mount post-reopen. Biologist sees blank squares, asks "is this broken?", expert recovers but the demo lost trust.

**Why it happens:**
- Phase 14 999.1 hotfix (silent data loss on dialog re-open with stale column refs) and 999.4 hotfix (false corruption signal from streaming/text malformed-counter mismatch) are both examples of "looked broken to the user even though it worked" — this is a recurring demo failure mode for this audience.
- showHeatmap, showQcDashboard already have documented progress-indicator gaps (`debug/showheatmap-hangs.md`, CONCERNS.md QC dashboard). Same shape WILL happen on first-time view of a published project (cold-load of Dendrogram, layout restoration).
- The `2026-05-04` Cytokinetics demo surfaced QC dashboard layout navigability issues (CONCERNS.md "QC Dashboard — Layout Navigability") — biologist consumers visually parse layouts differently than experts.

**How to avoid:**
- Audit every new v1.4 UI surface (publish dialog, SPC chart, campaign comparison) for jargon. Replace "DataFrame" with "table", "tag" with "label", "ACL"/"permissions" with "who can see". Phase 14 already moved away from raw IDs (D-08 inline provenance markers); apply the same instinct here.
- Every async load on reviewer-side surfaces (project open, viewer mount, SPC chart compute) must have a `DG.TaskBarProgressIndicator` with explicit phase messages ("Loading volcano...", "Restoring layout...", "Computing SPC limits..."). Reuse the Phase 14 R3 progress wiring patterns.
- "Visible loading" rule: if any reviewer-side computation takes >500ms, it has a progress indicator. No exceptions.
- The published project's title carries the analysis context ("DMD vs WT, limma, FC>=1.5, adj.p<0.05, published 2026-06-15") so a biologist landing cold has self-contained context.
- For demos: pre-warm the cache on a "warm-up" view click before the demo starts (per Cytokinetics demo-blocking pattern in CONCERNS.md).

**Warning signs:**
- A new dialog uses any word from {DataFrame, tag, ACL, semType, viewer factory, function call} in user-facing copy.
- A new viewer's first paint takes >500ms and has no spinner.
- The published project's title is just a UUID or table name.

**Phase to address:** Cross-cutting (Publishing + SPC + Campaign all touch this) — assign to a single "Reviewer-facing UX hygiene" task per phase, OR add as a phase-completion gate (no phase ships until every reviewer-touchable surface has been audited).

**Test/UAT lever:**
- Pre-demo dress rehearsal: have a biologist (or a stand-in unfamiliar with the package) walk through the demo cold, note every "wait, what does that mean?" — they all become defect tickets.
- Programmatic check: lint every dialog string for jargon words from the banned list.
- UAT: every async path has `pi.update(...)` calls visible in the code review.

---

## Technical Debt Patterns

Shortcuts that seem reasonable but create long-term problems.

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Publish by `df.columns.remove(...)` on the live df (no clone) | One-liner, ships fast | Pitfall 1 (silent stale snapshot). Recovery requires republishing every existing published analysis. | NEVER |
| "We'll just store everything in tags" for publish/SPC/campaign metadata | Mirrors existing pattern, no schema change | Tags can be dropped by the serializer (Pitfall 3) and tag-key collisions across phases (proteomics.* namespace is shared) | Acceptable as a runtime convenience IF the same value is ALSO encoded in a column (belt-and-braces) |
| "We'll add a re-run reset later" for normalize/impute idempotency | Saves dialog complexity | Re-run pollutes published copies (compounds Pitfall 1) | Acceptable only if publish always deep-clones AND the normalize/impute step rejects re-runs explicitly |
| Reuse v1.2 `activeSubscriptions[]` pattern for campaign cross-DF wiring | Familiar code | Subscription leak (Pitfall 13). Module-level mutable state in growing N runs = leaks. | NEVER for v1.4 — build a bus with eviction |
| Inner-join across runs because "the proteins are mostly the same anyway" | One-liner using existing join helpers | Pitfall 9 (silent data loss in cross-run comparison) | NEVER |
| Compute SPC baseline from "all runs to date" | Auto-updates, no UI | Pitfall 6 (baseline contamination); chart eventually flags nothing | NEVER — baseline is an explicit, locked subset |
| Order SPC by import datetime | Trivial; already in df.tags | Pitfall 7 (run-order ambiguity, instrument mixing) | NEVER for SPC; acceptable for ad-hoc audit logs |
| Store full quant matrix in SPC project for "future drill-down" | Optionality | Pitfall 8 (storage growth, slow open) | NEVER — link to source project for drill-down |
| Skip post-save reopen test "because save succeeded" | Faster CI | Pitfall 3 manifests in production demo | NEVER — round-trip test is mandatory for publishing |
| Default to ALL 8 Nelson rules enabled | "Catches everything" | False-alarm rate ~2.65%, reviewers lose trust | NEVER — default to rules 1+5 |
| Direct `fetch(...)` for any new external service (Ensembl, ChemAxon, etc.) | One-liner | CORS in production (existing UniProt + g:Profiler tech debt — `CONCERNS.md`) | NEVER — always `grok.dapi.fetchProxy(...)` |
| "Re-publish overwrites the previous published project" | Simple one-click republish | Loss of audit trail, broken bookmarks (Pitfall 4) | NEVER — republish creates a new project; old gets `superseded_by` |
| Hand-craft permissions per-publish | Maximum flexibility | Pitfall 2 (Edit slips in via inheritance) | NEVER without a post-grant verification step |

---

## Integration Gotchas

Common mistakes when connecting to external services and platform-internal features.

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| Datagrok project save | Assume custom tags survive — phase 13 e527d07ba1 proves they don't always | Round-trip test in CI; encode critical values in a column as well as a tag |
| Datagrok permissions | Grant View on project, assume Edit is impossible | Explicitly verify effective permissions post-save (`grok.dapi.permissions.get`) and refuse to ship if Edit is present for the reviewer group |
| Datagrok project layout | Save layout via `view.saveLayout()` and assume it restores fully | Test: save layout, reopen project, assert specific viewers + their dock positions match. Layout serializer has known partial-strip behavior (Phase 13 round-3) |
| `df.clone(filter)` | Assume tags + semTypes + name come along | Verify via test (some are preserved by Datagrok, but the contract is partial); explicitly re-set proteomics.* tags and df.name post-clone |
| Module-level `activeSubscriptions[]` | Copy the v1.2 pattern for v1.4 campaign wiring | Build a per-campaign selection bus with subscription auto-eviction tied to viewer dispose |
| Chem package SMILES canonicalization | Assume the existing UniProt/g:Profiler pattern works for Chem | Use the platform's Chem package via `DG.Func.find({package: 'Chem', ...})` with the same optional-dependency pattern as Dendrogram (`ARCHITECTURE.md`); provide a no-Chem fallback path |
| `grok s shares` for publishing | Hand-roll REST | Use `grok s shares add` / `grok s users save --json` per `tools/GROK_S.md`; matches the platform's idempotency contract |
| Reviewer's published-project reopen | Assume "save -> reopen" round-trip is lossless | Mandatory `assertPublishedShape` in test + UAT; ALL critical state asserted post-reopen |
| Cross-package function lookups (Dendrogram, Chem) | Hard-import at top of file | Use `DG.Func.find` runtime lookup with a graceful fallback (existing Dendrogram pattern in `src/viewers/heatmap.ts:138-145`) |
| Ensembl REST (Phase 14 D-09) and any new external API | Direct `fetch()` for "speed" | Always `grok.dapi.fetchProxy(...)` — existing UniProt + g:Profiler raw-fetch is documented tech debt (CONCERNS.md "External Coupling") and should not be replicated |

---

## Performance Traps

Patterns that work at small scale but fail as usage grows.

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Per-protein SPC charts (5,000 charts/run) | UI freeze on SPC open; hundreds of false alarms | Run-level metric set only (4 metrics) | Day 1 |
| Full-quant SPC storage | SPC tab takes >5s to open after ~50 weeks | Store run-level rollup only; link to source project for drill-down | ~50 weekly runs |
| Recomputing SPC baseline from full history on every open | Open-SPC latency grows over time | Cache the locked baseline (mean+SD+rule-state) as a tag; only recompute on explicit "rebuild baseline" | ~30 weekly runs |
| Reloading the full source dataset to render a published trim view | First-paint stalls for biologist | Published project contains ONLY the trimmed clone; source is not loaded by the reviewer flow | Day 1 — biologist's tab freezes |
| Subscribing to N volcano viewers' `onCurrentRowChanged` in a campaign | Selection latency grows linearly with N; subscriptions accumulate across reopens | Single-bus selection abstraction with eviction (Pitfall 13) | ~5+ runs per campaign |
| Computing the full outer-join on protein ID for every comparison render | Slow first-paint on large campaigns | Compute and cache the join once at campaign load; viewers read the cached joined DF | ~10 runs × 5k proteins |
| Linear scan over campaign DataFrame for "show run X's column" | Cell-render latency on cross-run grid | Build a column-name -> column-index map at campaign-load time | ~10 runs (50+ per-run columns) |
| Heatmap cold-start (Dendrogram package load) on first reviewer view post-reopen | "Loading..." for 5-10s with no indicator on the published view | Add progress indicator on viewer mount; pre-warm Dendrogram in `initProteomics()` (existing CONCERNS.md item) | First reviewer demo every time |
| Sequential R cold-start chain (DE limma + DEqMS) repeated for re-published SPC backfills | Hours of cold-starts if backfilling many weeks | SPC NEVER re-runs DE — uses already-computed metrics from published projects | Backfill operations |
| Publishing huge projects (raw intensities still in clone) | Project save takes minutes; project size GB-scale | Trim to <50 columns BEFORE clone (Pitfall 1 + Pitfall 8 hygiene) | First publish of a 10k-protein 50-sample DataFrame |

---

## Security Mistakes

Domain-specific security issues beyond general web security.

| Mistake | Risk | Prevention |
|---------|------|------------|
| Publish includes intermediate normalized intensities, exposing per-sample concentrations | Reviewer infers protocol details from intensity ranges; competitive info leak | The trim contract MUST exclude raw + normalized intensity columns; column allowlist enforced in publish helper |
| Publish includes sample labels containing patient/subject IDs | PHI/PII exposure | Sample names are NOT in the trim contract; only group labels (`Numerator`/`Denominator`) and aggregated metrics travel |
| Publish reveals the source compound list to a reviewer group that should only see one target | Cross-program info leak | Cross-team review filed BY TARGET (per todo `2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target.md`) means the target-Space's reviewer group is the audience; trim contract excludes compound columns unless the campaign is target-locked |
| SPC project visible to all proteomics users, exposes instrument-down events | Vendor / collaborator visibility of operational issues | SPC project lives in a proteomics-only Space; not publicly published |
| Reviewer can export the published DataFrame to CSV and email it externally | Even read-only doesn't prevent exfiltration | Mitigate via UX (no "Download" button on the published view) + audit log (Datagrok's existing audit captures `dapi.tables.download` events); accept that determined exfil is out of scope |
| Permission grants stick around after a reviewer leaves the company | Stale-access drift | Reviewer permissions are tied to a Datagrok GROUP, not individual users; group membership management is upstream of this package |
| `proteomics.published_by` tag is settable from the client | Audit log forgery | Set `published_by` server-side from the authenticated session; never trust client input |

---

## UX Pitfalls

Common user experience mistakes for the Cytokinetics-class audience (proteomics expert + biologist consumer in the same demo).

| Pitfall | User Impact | Better Approach |
|---------|-------------|-----------------|
| Jargon in reviewer-side dialogs ("DataFrame", "tag", "semType", "ACL") | Biologist disengages within 30 seconds | Audit every reviewer-visible string for banned words; use "table", "label", "type", "who can see" |
| No progress indicator on published-project reopen | Biologist thinks the project is broken | Every async path on reviewer side has `TaskBarProgressIndicator` with phase messages |
| Republish silently overwrites the bookmark URL | Biologist clicks a saved link, sees different data, loses trust | Republish creates a NEW project (Pitfall 4); old URL still resolves to the original numbers + a banner "Newer version available: [link]" |
| SPC chart with no centerline / limits visible | Reviewer can't tell what's "normal" | Always render centerline + UCL/LCL + a legend; never just dots |
| Side-by-side volcano with different y-axis scales but no annotation | Reviewer overinterprets one run's strength | Per-volcano sample-count + axis-scale callout in the title (synthesize as in Phase 14 D-02) |
| Comparison view doesn't show "this protein is in N of M runs" | Reviewer doesn't notice population mismatch | Quantified-In count column visible by default (Pitfall 9) |
| "Frozen" published view has the same colors as the expert's working view | Reviewer doesn't realize this is read-only | Distinct banner / title-bar treatment for published view: "Read-only — published by X on Y" |
| Empty state when reviewer's bookmark hits a deleted source DataFrame | Cryptic platform error | Trim contract = self-contained clone; reviewer never depends on the source df existing |
| Cross-run comparison defaults to the largest run, hiding the smaller ones | Implicit bias in interpretation | Default sort = chronological; no preselection of "the best" run |
| Reviewer's comments don't show up next time the expert opens the project | Bidirectional sign-off broken | Use Datagrok's existing comments thread on the published project entity (per `2026-05-11` share-readonly todo) |

---

## "Looks Done But Isn't" Checklist

Things that appear complete but are missing critical pieces.

- [ ] **Publish primitive:** Often missing post-save reopen verification — verify `assertPublishedShape(reopenedDf)` runs in CI AND in the publish handler itself before declaring success.
- [ ] **Publish primitive:** Often missing deep-clone step — verify `frozen` and `df` are distinct instances by mutating one and asserting the other doesn't change.
- [ ] **Publish primitive:** Often missing explicit permission verification — verify `dapi.permissions.get(project)` shows reviewer group as View-only, never Edit.
- [ ] **Publish primitive:** Often missing version + supersede pointer — verify `proteomics.publish_version` and `proteomics.superseded_by` are present.
- [ ] **SPC baseline:** Often missing iterative outlier removal + explicit lock — verify the baseline is computed once with user confirmation and stored as a tag, not recomputed on every open.
- [ ] **SPC chart:** Often missing centerline + UCL/LCL annotation — verify all three lines are rendered with explicit labels.
- [ ] **SPC chart:** Often missing instrument facet — verify multi-instrument data renders as separate traces, not interleaved.
- [ ] **SPC run identity:** Often missing acquisition_datetime — verify the annotate step captures it explicitly (not file-mtime).
- [ ] **SPC false-alarm hygiene:** Often missing per-rule false-alarm-rate disclosure — verify the dialog tooltip discloses expected alarm rate per enabled Nelson rule.
- [ ] **Campaign join:** Often missing full-outer-join semantics — verify the joined DF rowCount equals union of per-run protein sets (not intersection).
- [ ] **Campaign join:** Often missing "Quantified-In" column — verify it lists run IDs and is visible by default.
- [ ] **Campaign viewer:** Often missing side-by-side independence — verify `tv.viewers` contains N viewers for N runs (not 1 overlaid scatter).
- [ ] **Campaign compound entity:** Often missing canonical resolution — verify identical SMILES across runs always resolves to the same `compound_id`.
- [ ] **Campaign selection bus:** Often missing eviction-on-dispose — verify subscription count returns to 0 after closing the campaign view.
- [ ] **Cross-DataFrame integration with Chem:** Often missing optional-dependency fallback — verify a no-Chem environment renders names, not blank cells.
- [ ] **Reviewer-side UX:** Often missing progress indicators on async paths — verify every >500ms operation has a `TaskBarProgressIndicator`.
- [ ] **Reviewer-side UX:** Often missing jargon audit — verify reviewer-visible strings don't contain `{DataFrame, tag, semType, ACL, viewer factory}`.
- [ ] **Project versioning:** Often missing "superseded by" banner on old version — verify reopening a superseded project shows the link to the new one.
- [ ] **Tag survival:** Often missing belt-and-braces column encoding — verify critical metadata is in BOTH a tag AND a column.
- [ ] **External API:** Often using raw `fetch()` — verify all new external calls go through `grok.dapi.fetchProxy(...)`.

---

## Recovery Strategies

When pitfalls occur despite prevention, how to recover.

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Pitfall 1 (stale snapshot leak) | HIGH | Identify all published projects produced before the fix; republish each from the original source if available, OR mark them with a `proteomics.publish_integrity = 'suspect'` tag and a banner warning. Cannot retroactively prove which numbers were correct. |
| Pitfall 2 (reviewer has Edit) | MEDIUM | Audit `dapi.permissions.get` on every published project; revoke Edit for reviewer groups; re-issue View-only. Check audit log for any reviewer edits during the affected window. |
| Pitfall 3 (tags/semTypes dropped) | MEDIUM | Open each affected project, re-apply tags via `df.setTag` and re-set semTypes via `col.semType = ...`, re-save. Where the trim contract was column-only, the column data still carries the answer; re-tag from column. |
| Pitfall 4 (versioning ambiguity) | MEDIUM | Add `publish_id` + `publish_version` retroactively to existing published projects via a one-shot script; sort by `published_at` to reconstruct the order; add `superseded_by` pointers where chains can be inferred. |
| Pitfall 5 (per-protein SPC alarm flood) | LOW | Disable per-protein SPC entirely; switch UI to run-level metrics. User-facing message: "Per-protein SPC removed; use run-level metrics for trend monitoring." |
| Pitfall 6 (baseline contamination) | MEDIUM | Manually re-select clean baseline runs with the expert; recompute limits; force-rebuild all downstream SPC state. |
| Pitfall 7 (run order wrong) | LOW | Add the missing instrument + acquisition_datetime fields to all historical runs; recompute SPC charts. |
| Pitfall 8 (SPC project bloat) | MEDIUM | Migration script that re-derives the run-level rollup from existing published projects; rebuild SPC project; archive the bloated one read-only. |
| Pitfall 9 (inner join silent loss) | LOW | Rebuild campaign DataFrame with full outer join; comparison view re-renders from the corrected DF. User-facing message: "Comparison updated to include proteins quantified in only some runs." |
| Pitfall 10 (volcano overlay misleads) | LOW | Replace overlay viewer with side-by-side; preserve URLs by redirecting the old viewer route. |
| Pitfall 11 (compound name drift) | MEDIUM | Manual canonicalization session with the expert; re-annotate historical runs with SMILES; rebuild campaign DataFrame. |
| Pitfall 12 (batch effect confounding) | LOW (recovery) / HIGH (interpretation) | Add vehicle-trend panel to existing comparison views; do NOT silently re-interpret historical conclusions. Surface the confound, let the expert/biologist re-call. |
| Pitfall 13 (subscription leak) | LOW | Force-close all open campaign views; reload page; ship the bus abstraction in the fix. |
| Pitfall 14 (demo failure) | HIGH (trust) / LOW (code) | Post-mortem with the Cytokinetics champion; identify specific moments of confusion; pre-demo dress rehearsal with a biologist-class user before every customer demo going forward. |

---

## Pitfall-to-Phase Mapping

How roadmap phases should address these pitfalls.

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| 1 — Stale snapshot leak | Publishing-phase (Task 1) | Unit + UAT: mutate source after publish, assert frozen is unchanged |
| 2 — Reviewer Edit slip | Publishing-phase (Task 1, paired with Pitfall 1) | Authenticated `grok s` test as reviewer, write attempts must fail |
| 3 — Tags/semTypes drop on save | Publishing-phase (Task 1 prerequisite: round-trip test infra) | `publish-roundtrip.ts` test file in `src/tests/`; `assertPublishedShape` helper |
| 4 — Versioning ambiguity | Publishing-phase (Task 2: versioning + supersede) | Republish test asserts new project id + supersede pointer |
| 5 — Per-protein SPC multiplicity | SPC-phase (Task 1: metric scoping decision document) | Scope explicitly states "run-level only"; UAT confirms the menu has no "by protein" affordance |
| 6 — Baseline contamination | SPC-phase (Task 2: baseline workflow) | Unit test: outlier-contaminated baseline yields tight limits after iterative removal |
| 7 — Run-order ambiguity | SPC-phase (Task 0: data model) | Unit test: backfill-inserted run sorts by acquisition_datetime, not insert order |
| 8 — SPC storage bloat | SPC-phase (Task 0: data model + schema) | UAT: 52-week simulated history opens in <2s |
| 9 — Inner join silent loss | Campaign-phase (Task 1: comparison data model) | Unit test: asymmetric run sets produce union-sized joined DF with explicit NaN |
| 10 — Volcano overlay misleads | Campaign-phase (Task 2: comparison viewer) | UAT: viewer is side-by-side with independent axes + sample-count callouts |
| 11 — Compound drift | Campaign-phase (Task 0: compound entity) | Unit test: SMILES-equivalent labels resolve to same compound_id |
| 12 — Batch effect confound | Campaign-phase (Task 3: vehicle trend panel) | UAT: simulated batch step is visible in the comparison view |
| 13 — Subscription leak | Campaign-phase (Task 1: selection bus, before viewers) | Unit test: subscription count returns to 0 after view dispose |
| 14 — Demo-failure / jargon / hung viewers | Cross-cutting (every phase) | Per-phase exit gate: reviewer-side dress rehearsal with biologist-class user |

**Phase ordering implication:**
1. **Publishing-phase first.** The publish primitive (Pitfalls 1-4) is the foundation; SPC and Campaign both reuse the trim/version/permission patterns.
2. **SPC-phase second.** Per-run metrics depend on the publish-time canonical metrics being defined first.
3. **Campaign-phase third.** Cross-run comparison reads BOTH SPC trust signals (was this run in-control?) AND published per-run analyses.
4. **Cross-cutting Reviewer UX hygiene** is a per-phase exit gate, not a separate phase.

---

## Sources

Codebase-grounded (HIGH confidence):
- `packages/Proteomics/.planning/codebase/ARCHITECTURE.md` — in-place mutation pattern, cross-DataFrame selection, anti-patterns, error handling cascade
- `packages/Proteomics/.planning/codebase/CONCERNS.md` — R cold-start, kNN blocking, raw fetch tech debt, two-group model gap, no compound entity, no longitudinal data model
- `packages/Proteomics/.planning/codebase/INTEGRATIONS.md` — Dendrogram optional-dependency pattern, external HTTP services
- `packages/Proteomics/.planning/codebase/STRUCTURE.md` — file layout, "Where to Add New Code" contract
- `packages/Proteomics/.planning/codebase/TESTING.md` — test framework, fixture pattern, coverage gaps
- `packages/Proteomics/.planning/milestones/v1.3-MILESTONE-AUDIT.md` — Phase 13/14 round-3, stale-test fixes, serializer-strip evidence
- `packages/Proteomics/.planning/milestones/v1.3-phases/14-ck-omics-analyst-experience-enhancements/14-CONTEXT.md` — Phase 14 decisions, color lock, Filters viewer scoping pattern
- `packages/Proteomics/.planning/milestones/v1.3-phases/14-ck-omics-analyst-experience-enhancements/14-HUMAN-UAT.md` — 999.1/999.4 hotfix lessons, serializer-strip post-mortem
- `packages/Proteomics/.planning/debug/limma-de-slow.md` — R cold-start dominates DE latency
- `packages/Proteomics/.planning/debug/showheatmap-hangs.md` — Dendrogram cold-start, missing progress indicator pattern
- `packages/Proteomics/.planning/todos/pending/2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target.md` — Cytokinetics 2026-05-04 origin of the share-by-target use case
- `packages/Proteomics/CLAUDE.md` — pipeline tag contract, naming conventions, clone-for-isolation pattern
- `tools/GROK_S.md` — `grok s` patterns for permission verification and idempotent share operations

Datagrok platform (MEDIUM confidence, single-source):
- [Datagrok Projects](https://datagrok.ai/help/datagrok/concepts/project/) — sharing requires explicit grant; permission inheritance via Space; static-snapshot vs generation-script persistence
- [Datagrok Layouts](https://datagrok.ai/help/develop/how-to/views/layouts) — `view.saveLayout()` / `view.loadLayout()`; viewStateMap; layout serialization
- [Datagrok Tags](https://datagrok.ai/help/govern/catalog/tags) — tags as the platform extensibility mechanism
- [Datagrok Semantic Types](https://datagrok.ai/help/govern/catalog/semantic-types) — semType as a column tag; detectors

SPC literature (MEDIUM confidence, multi-source consensus):
- [Shewhart control chart rules — Analyse-it](https://analyse-it.com/docs/user-guide/process-control/shewhart-control-chart-rules) — Nelson rules definitions, normality assumption
- [The Nelson Rules — superengineer.net](https://www.superengineer.net/blog/spc-nelson-rules-tests) — false-alarm rates per rule
- [False Alarm Rates for the Shewhart Control Chart with Interpretation Rules — ResearchGate](https://www.researchgate.net/publication/242194945_False_Alarm_Rates_for_the_Shewhart_Control_Chart_with_Interpretation_Rules) — combined alarm rate ~2.65% across 8 rules
- [CPV Monitoring — Science Journal of Applied Mathematics and Statistics](https://www.sciencepublishinggroup.com/article/10.11648/j.sjams.20241202.11) — non-normal distributions inflate false alarms

Memory (HIGH confidence, project-local):
- `feedback_datagrok_close_before_refresh.md` — close all views before Playwright reload to avoid `beforeunload` dialogs (UAT pattern)
- `project_proteomics_deqms_result_misalignment.md` — row-key alignment lesson (relevant if SPC tries to align run rollups back to per-protein source)
- `feedback_keep_workaround_capture_future.md` — when the platform fixes the serializer-strip behavior, the belt-and-braces column encoding should be removed as a future cleanup, not preemptively

---

*Pitfalls research for: Proteomics v1.4 Cross-Team Review (Publishing + SPC + Campaign).*
*Researched: 2026-06-06*

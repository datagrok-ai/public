# Feature Research — Proteomics v1.4 Cross-Team Review

**Domain:** Mass-spec proteomics analytics — cross-team handoff, longitudinal QC, screening-campaign comparison
**Researched:** 2026-06-06
**Confidence:** MEDIUM-HIGH (industry conventions verified across multiple sources; specific Cytokinetics asks inferred from 2026-05-04 meeting notes)

Three user stories, one milestone. Each section below scopes one story end-to-end:
table stakes (must ship), differentiators (where v1.4 earns its keep against Spectronaut
Direct / Scaffold / Skyline+MSstatsQC), anti-features (commonly asked, regret to ship).

---

## Story 1 — Read-Only Analysis Publishing (Proteomics expert → Biologist)

**Mental model:** Spectronaut Direct, Scaffold Viewer, and in-house "Excel + PowerPoint" handoffs
are the prior art. The proteomics expert finishes DE + enrichment, then ships a **frozen,
trimmed, audit-stamped artifact** to a biologist or biologics team filed by drug target.
The artifact is a *deliverable* — versioned, dated, signed off — not a working table the
biologist edits.

### Table Stakes

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Trimmed published DataFrame** — Protein ID, Gene Symbol, log2FC, p-value, adj.p-value, significant flag, direction (up/down). Nothing else. | This is what biologists already paste into Excel; matches the Scaffold "Viewer" surface, mirrors the columns every published DE table in literature carries | LOW | Clone `df` then drop columns by `proteomics.*` semantic types (drop INTENSITY, INTENSITY_LOG2, PEPTIDE_COUNT, anything not in a published-allowlist). Volcano viewer code path unchanged; just bound to the trimmed DataFrame. |
| **Volcano viewer renders on open with thresholds drawn** | Biologist double-clicks the shared project → sees the volcano with the FC/p-value lines they were promised. No "click here, then there" | LOW | Reuse `createVolcanoPlot`. Layout serialization in the project handles the threshold lines (already works in v1.3). |
| **Audit context panel** — DE method (limma/DEqMS/t-test), Numerator vs Denominator group names, FC threshold, adj.p threshold, normalization method, imputation method, sample count per group, source author, share date | "Where did these numbers come from" is the first question every biologist asks. Spectronaut Direct surfaces this in a sidebar; Scaffold has it under "Annotations" | LOW | Read from `proteomics.*` tags already set in v1.3. Render as a simple read-only info panel. |
| **Frozen — no edit, no re-run, no parameter change** | The whole point of "delivered" is that the numbers don't move under the reviewer's feet. Reviewer-group permission is View, not Edit | MEDIUM | Datagrok `grok.dapi.permissions` — View-only project ACL. Trim all `Proteomics → ...` menu items from the project's package context, or guard each handler to no-op if the DataFrame carries a `proteomics.published` tag |
| **Target-keyed filing** — published project lives under "Target X" not under the author's personal workspace | Biologists navigate by target/drug program, not by who ran the experiment. Multi-author target dossiers are the unit of biologist-facing work | MEDIUM | One-Datagrok-Space-per-target *or* one top-level "Targets" Space with subspaces. Open question: project tag vs. Space — see Open Questions §1. |
| **Reviewer-group permission** — Biologics group sees, only Proteomics group publishes | Standard org separation of duties. The proteomics team owns the analysis; biologists own interpretation | LOW | `grok.dapi.permissions` — `share(project, biologicsGroup, View)`. The author's group keeps Edit. |
| **Versioned/dated artifact name** — "Target-CK-001 — Compound A vs DMSO — 2026-06-06 — v1" | Biologists compare across iterations; an unversioned re-share silently overwrites. Spectronaut export filenames carry the date stamp by default | LOW | Project name template: `{target} — {numerator} vs {denominator} — {YYYY-MM-DD} — v{N}`. Auto-increment N by listing prior projects in the target Space. |
| **No raw intensities, no peptide counts, no original file** in the published artifact | Raw intensities are the proteomics team's IP and a data-leak hazard; peptide counts encode batch metadata the biologist shouldn't reason from. Spectronaut Direct deliberately ships pivoted DE tables, not the raw quant matrix | LOW | The trim step. Drop by semantic type as above. |

### Differentiators (vs Spectronaut Direct, Scaffold Viewer, "Excel + PowerPoint")

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Optional enrichment results travel with the share** | Scaffold doesn't do GO/KEGG; Spectronaut does. Carrying the enrichment dot plot + linked-volcano interaction (v1.2) means biologist can ask "show me proteins in pathway X" without re-running anything | LOW | Already built in v1.2 (`viewers/enrichment-viewers.ts`, `onCurrentRowChanged`). Just include enrichment DataFrame in the project bundle. |
| **One-click publish from the menu** — "Proteomics → Share Analysis for Review…" | Spectronaut requires manual export + Scaffold packaging. Excel handoffs require manual column trimming. One-click matches the v1.3 "Proteomics → ..." UX promise | MEDIUM | New menu item, new dialog, new analysis module. ~1 phase of work. |
| **"Request re-run with different parameters" reviewer-side button** | Closes the loop without giving the biologist edit access. Spectronaut and Scaffold don't have this | MEDIUM | Could be as light as a `mailto:` template populated with audit context, or as heavy as a Datagrok comment thread. Recommend the light version first; defer comments. |
| **Audit context as DataFrame metadata, not free-text** | Future "compare published v1 vs v3" tooling can read the audit programmatically. Spectronaut/Scaffold use free-text annotations that humans must re-parse | LOW | Already the pattern in v1.3 — `proteomics.*` tags. Just expose them in the panel. |

### Anti-Features

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| **Comments / threaded reviewer notes** | "Biologist wants to leave feedback" — instinctive ask | Comment thread is a small product unto itself: persistence, notifications, edit/delete, mentions, audit. Easily 1 phase of work for v1 quality. Cytokinetics meeting did NOT call this out as P1 | "Request re-run" button (mailto:) closes the loop in v1.4. Comments live in v1.5+. |
| **Editable reviewer-side annotations on the volcano** | "Biologist wants to circle proteins" | Edits to a "frozen" artifact violates the whole story. Forks the source of truth | Reviewer-side selection is fine (highlight in their session); persistence is not. |
| **Inline re-run with new thresholds in the reviewer view** | "Why not just let the biologist try FC=2?" | Defeats "frozen analysis," doubles permission surface, makes "what they saw" non-reproducible | The reviewer asks the proteomics author to publish v2 with the new threshold. Versioned filenames carry the difference. |
| **Embedded raw intensity matrix "in case they need it"** | "What if they want to verify?" | The biologist isn't going to re-analyze; including it just leaks data and bloats the project | If verification is needed, the proteomics author re-runs and re-publishes. Raw matrix stays in the proteomics author's workspace. |
| **Per-protein "approve/reject" UI** | Sounds like proper review | Massively complicates the data model (per-row state), and biologists don't think this way — they think "the analysis as a whole tells X story" | Biologist comments in their re-run request: "I'd downweight proteins Y, Z because [reason]." Author iterates. |
| **PDF report export** | "Send to leadership" | PDF generation is its own engineering project; Datagrok Projects already serve as the shareable artifact for anyone inside Datagrok | Inside Datagrok = the project. Outside Datagrok = author screenshots the volcano. v1.4 ships *inside Datagrok* sharing only. |

### Open Questions (for Requirements / Roadmap)

1. **Space-per-target vs project-tag-by-target.** Space gives clean permission boundaries
   (one ACL per target, biologists see only their targets) but creates hundreds of Spaces
   if the target portfolio is large. A single "Published Analyses" Space with a
   target-tagged project is simpler but bleeds permissions across targets. **Recommend
   prototype both** in roadmap research; pick based on Cytokinetics target portfolio size.
2. **"Target" — first-class entity or just a string?** v1.4 likely ships as a free-text /
   autocomplete dropdown; first-class target taxonomy ("Target X has identifier UniProt
   Y, has alternate name Z") is its own milestone.
3. **Is the published artifact one DataFrame or a project bundle?** Recommend project
   (lets us also carry the enrichment DataFrame + viewer layout).

---

## Story 2 — Sample-Level SPC for Longitudinal QC

**Mental model:** Skyline's SProCoP and MSstatsQC have done exactly this for targeted
proteomics for a decade. The well-known shape: **per-metric I-chart (individuals) + MR-chart
(moving range) per QC metric, with Shewhart ±3σ control limits and Western Electric / Nelson
rule overlays**. The Multiyear Longitudinal Harmonization study (2024–2025) established
community reference values for the standard metrics and confirmed the pattern works at
phenotypic-screen scale.

The v1.4 ask is the longitudinal sibling of the v1.1 within-run QC dashboard — same metrics,
new axis (week, not sample).

### Table Stakes

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Four standard sample-level metrics per run** — (1) median log2 intensity, (2) % missing values, (3) within-control-replicate Pearson correlation, (4) protein count above detection threshold | These four are the SProCoP / MSstatsQC default set and the Harmonization Study's reference metrics. Anyone in pharma proteomics expects these | LOW | Each metric is a `reduce()` over the v1.3 control-group intensity columns + the existing `proteomics.groups` tag. Reuse v1.1 `qc-computations.ts` helpers. |
| **Shewhart I-chart (Individuals) per metric with ±3σ control limits** | The canonical SPC chart. ±3σ from process mean, drawn over time, points colored red when outside limits | LOW | Datagrok line-chart viewer with formula lines (existing primitive — same one volcano uses for FC/p thresholds). One viewer per metric. |
| **MR-chart (Moving Range) alongside the I-chart** | I-chart catches level shifts; MR-chart catches variability shifts. SProCoP/MSstatsQC pair them by convention | LOW | Second viewer, MR = abs(value_t − value_{t-1}); UCL = 3.267 × mean(MR). |
| **Default Nelson rule set with selectable enable/disable** | Nelson Rules 1–4 (point beyond 3σ, 9 same-side, 6 trending, 14 alternating) are the "industry default" subset. Letting users disable noisy rules per metric is standard in lab QC software (Westgard panels, MSstatsQC) | MEDIUM | Implement Nelson 1–4 first. Add 5–8 as optional. Per-metric rule toggle stored in `proteomics.spc.rules` tag. |
| **Run-level pass/flag tag written back to the analysis DataFrame** | Downstream Story-3 comparison MUST be able to ask "did this run's controls pass QC?" — otherwise cross-week compound comparisons inherit drift silently | LOW | Set `proteomics.qc_status = 'pass' \| 'flag' \| 'fail'` and `proteomics.qc_flagged_rules = '[1,4]'` on the run's DataFrame. |
| **Click flagged point → drill to that week's analysis** | Standard SPC dashboard affordance. Discovered drift is useless if you can't go inspect the run that caused it | LOW | Datagrok DataFrame links; `currentRowIdx` → project lookup → open. |
| **Process-mean baseline configurable: rolling-window vs fixed-reference** | First N runs (typically 10–20) anchor the baseline; some labs prefer a fixed reference period (Q1 2026), others prefer a rolling window (last 20 runs). MSstatsQC supports both | MEDIUM | Dialog input: "Baseline: [Rolling N=20] vs [Fixed: runs from {date} to {date}]". Default to rolling N=20. |
| **Phase-aware mode (reagent-lot break, instrument re-calibration)** | "We changed reagent lot at run 17 — recompute baseline from run 17 forward." Every pharma proteomics group has this; without it, an annotated event triggers a permanent SPC alarm | MEDIUM | User marks a "phase change" on a run; SPC recomputes process mean per phase. Persist as `proteomics.spc.phase_breaks` tag. |

### Differentiators (vs Skyline SProCoP, MSstatsQC, MaxQuant.LiveReport)

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Integrated longitudinal table = Datagrok DataFrame, not a separate database** | MSstatsQC requires a separate R Shiny app and a long-format CSV. SProCoP runs inside Skyline. Datagrok keeps the longitudinal QC table sitting next to the analyses, queryable like any DataFrame | LOW | A "build-longitudinal-table" function — input: list of weekly analysis DataFrames; output: a `(run_date, metric, value)` DataFrame. |
| **Linked back to the v1.1 within-run QC dashboard** | One mental model: "within-run QC" is the cross-run SPC's row at a single timestamp. Click the SPC point → opens within-run dashboard for that week | LOW | Already discussed — DataFrame link. |
| **Run-level pass/flag gates Story-3 comparison automatically** | Story 3's compound comparison refuses to combine drifting runs unless the user opts in. None of the prior art does cross-feature integration like this — they're separate tools | LOW | Story-3 comparison viewer reads `proteomics.qc_status` on each run; warns if any are `flag` or `fail`. |
| **Pareto chart of "which metric flags most often"** | SProCoP includes this; MSstatsQC doesn't surface it as prominently. It's what a QA reviewer wants to see first — "where do I focus optimization?" | LOW | Bar chart of rule-violation counts per metric. One Datagrok bar chart viewer. |

### Anti-Features

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| **Per-protein SPC across the whole proteome** | "We want to know if any protein is drifting" — instinctive ask | At ~5,000 proteins × 8 Nelson rules, you'll see false positives every week from multiple-testing alone. SProCoP/MSstatsQC explicitly restrict to a curated panel for this reason | A **curated QC protein panel** (housekeeping + pathway anchors, ~50 proteins) with FDR-corrected SPC. Defer to v1.5 unless Cytokinetics specifically asks. |
| **CUSUM + EWMA charts in v1.4** | "More sensitive than Shewhart" — true, and they're in MSstatsQC | More charts to interpret, more parameters to tune (k, h for CUSUM; λ for EWMA), more documentation burden. Shewhart + Nelson covers 80% of drift detection | Ship Shewhart + Nelson; add EWMA/CUSUM in v1.5 if users ask. |
| **Automatic alerting / email-on-flag** | "Notify QA when a run flags" | Notification system is its own product. Project surface area, ACLs, opt-in/opt-out. Cytokinetics meeting didn't ask for it | Surface flags in the SPC dashboard. Add "tag run as flagged" → can be queried in usual Datagrok ways. Email later. |
| **Westgard multi-rule sets (12s/13s/22s/etc) as the default** | Clinical-lab Westgard panels are the gold standard for clinical chemistry SPC | Westgard rules are tuned for low-CV clinical assays; proteomics has higher noise. Nelson rules are the proteomics-community default. Don't conflate the two | Nelson rules default; advanced users can configure custom rule sets later. |
| **Drift-correction / signal-intensity normalization across batches** | "Don't flag — fix" | This is a *normalization* concern, not an *SPC* concern. Mixing them violates "common cause vs special cause" — SPC's whole point is to surface special causes for human action | SPC flags drift; the proteomics author decides whether to re-run or apply a batch-correction normalization in a new analysis. |

### Open Questions

1. **Control-replicate aggregation: per-week-mean (one point) vs per-replicate (multiple
   points per week)?** Per-replicate gives within-week variability info but is harder
   to plot legibly. Recommend per-week-mean by default; per-replicate as a viewer toggle.
2. **First N runs before chart is trustworthy.** SProCoP defaults to ≥10 runs; MSstatsQC
   has no hard rule. Recommend N=10 minimum, show "baseline building" overlay on first
   10 runs.

---

## Story 3 — Cross-Run Compound Comparison

**Mental model:** Pharma phenotypic / chemoproteomic screens run weekly: ~10 compounds vs
shared DMSO control in each batch. The chemist wants "which compound moved the proteome
most this week" and "how does this week's lead compound compare to last week's?" The package
has no concept of "compound" or "screening run" today — every analysis is one-off. The v1.4
ask is to lift the data model to multi-run, then build the comparison viewer.

### Table Stakes

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **First-class "screening run" entity** — (run_date, list of compounds tested, vehicle control reference) tag on the DataFrame | Without this concept the package can't even *describe* what gets compared. Cross-week comparison needs a noun | LOW | New `proteomics.screening_run` tag; JSON shape: `{run_date, vehicle_group, compound_groups: [{name, smiles}, …]}`. Mirrors the `proteomics.groups` pattern. |
| **First-class "compound" entity** — name + SMILES, with Datagrok Chem semantic type | Chemists think in structures. Without SMILES + Chem rendering, the "compound" column is just a string and the chemist can't reason about it | LOW | Re-use `Chem` package's `Molecule` semantic type. Annotate dialog: "Compound for treatment group" + SMILES input. Cell renderer comes for free from Chem. |
| **"Annotate Screening Campaign" extension to experiment-setup** | Today's `Annotate Experiment` dialog handles Control vs Treatment as anonymous groups. Screening flow needs Treatment-group → compound mapping | MEDIUM | Either extend `showAnnotationDialog` with an optional "this is a screening run" toggle, or add a sibling `showScreeningAnnotationDialog`. Recommend sibling — keeps the simple case simple. |
| **Side-by-side volcano viewer for two runs** — week-A volcano \| week-B volcano, shared FC/p thresholds, linked protein selection | The "compare top compound this week vs last week" workflow — the core ask. Linked selection (click protein in A → highlight in B) is the volcano-comparison standard | MEDIUM | Two `createVolcanoPlot()` calls (already in v1.3) side-by-side in a layout. Linked selection via shared protein-ID join + `onSelectionChanged` subscription, mirrors v1.2's `onCurrentRowChanged` enrichment pattern. |
| **Diff table — significant in A only / B only / both, with log2FC_A, log2FC_B, delta** | This is the table the chemist takes to a meeting. Excel users build this by hand today | LOW | Join two DE DataFrames on protein-ID; compute set membership and delta. New DataFrame, opens beside the side-by-side volcanoes. |
| **Compound structures rendered at the top of the comparison view** | Chemist-oriented framing. Without this, the comparison is just two scatter plots | LOW | Chem package's molecule cell renderer over the two SMILES from the run annotations. |
| **"Top-impact compound this week" ranking with a configurable score** | The screening lead needs to know *which* compound to compare. Without a ranking, every week requires a manual eyeball-and-pick | MEDIUM | Default: count of `\|log2FC\| > τ AND adj.p < α` proteins (simple, threshold-dependent). Make the rule a dropdown: count / sum-magnitude / top-N-strength / enriched-pathway-count. Persist choice as `proteomics.screening.score_rule` tag on the screening run. |
| **QC pass/flag gate** — comparison viewer warns if either run's `proteomics.qc_status` is `flag` or `fail` | Story 2's whole point is that cross-week comparisons are suspect under drift. The comparison must honor this | LOW | Read `proteomics.qc_status` from each run; warn banner above the side-by-side volcanoes. |

### Differentiators (vs "two Excel sheets and a meeting")

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Synchronized axes + shared FC/p threshold lines** | Visual concordance — same coordinates means the chemist's eye does the work. Without it the comparison is misleading | LOW | Compute global min/max across A and B, apply to both `scatter.props.xMin/xMax/yMin/yMax`. |
| **Enrichment-term overlap panel** — terms enriched in A only / B only / both | Chemist asking "are these compounds hitting the same biology?" — exactly what GO/KEGG overlap answers | MEDIUM | Reuse v1.2 enrichment results from each run. Set intersection on term-IDs. New DataFrame + dot plot. |
| **MCS / substructure highlight on the compound pair** | Chem package can render the maximum common substructure with differences highlighted. Chemist sees "compound B has a cyclopropyl where A has methyl — does that explain the proteome divergence?" | MEDIUM | Call into `Chem.mcsSearch` (existing Chem capability). Render MCS-highlighted molecules in the structure panel. **Depends on Chem package being installed/loadable** — gracefully degrade to plain SMILES rendering if not. |
| **Click protein in diff table → see which compound moved it more** | The "moved by which compound" answer is one click instead of one mental cross-reference. Excel users do this in their head today | LOW | `currentRowIdx` on diff table → highlight in both volcanoes + show \|delta\| numerically. |
| **Save the comparison as a sharable artifact (ties into Story 1)** | The chemist wants to send the comparison to their team. Without persistence, every comparison is re-created from scratch | MEDIUM | The comparison viewer state (the two source runs, the chosen score rule, the threshold) serializes as a Datagrok project. Story 1's publish flow can carry it. Deferred to Story 1 mechanics. |

### Anti-Features

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| **N-way (3+) compound comparison in v1.4** | "Why stop at two?" | Side-by-side scales to 2; 3+ creates layout nightmares (small multiples? overlay? toggle?), each with its own ambiguity. Cytokinetics ask is specifically week-N vs week-N+1 | Ship 2-way for v1.4. N-way is v1.5. |
| **Automatic compound clustering / SAR view** | "Find compounds with similar proteome effect" | This is a chemoproteomics analysis on its own — proteome-fingerprint clustering, similarity metrics, dendrogram viewer. Easily a milestone unto itself | "Top-impact ranking" gives the chemist their next compound. Clustering = v1.6+. |
| **Cross-run normalization of fold changes** | "Different runs have different magnitudes" — true, partly because of biological signal, partly drift | Conflates Story-2 drift-correction with Story-3 comparison. Will silently change DE numbers in ways the chemist didn't sign off on | Story 2's SPC tells you when fold-change magnitudes are drift-suspect; chemist decides whether to compare or quarantine. |
| **In-place editing of compound annotations from the comparison view** | "I noticed the compound name is wrong" | Edit-from-viewer creates round-trip ambiguity (which DataFrame is source of truth?). Just send the chemist back to the run annotation | "Edit annotation" → opens the run's `showScreeningAnnotationDialog`. |
| **Auto-suggest "next compound to test" based on SAR + proteome response** | "Let the algorithm pick" | Active-learning compound selection is a research-grade ML problem. Out of scope for a viewer. Cytokinetics has chemists; chemists pick | Show the data clearly; chemists choose. |
| **Build a compound database within the package** | "We need somewhere to store the compounds" | Compound registry is a platform concern (Datagrok has `MolTrack`, `HitTriage`, etc). Don't fork | Reference compounds by SMILES; let users supply them. If pharma customer has a compound DB, they connect via standard Datagrok connector. |

### Open Questions

1. **Where do compounds come from?** Per-campaign SDF/CSV upload (simplest) vs an internal
   compound-DB connector (Cytokinetics-specific). Recommend per-campaign upload for v1.4;
   compound-DB integration is per-customer customization.
2. **Is the "screening campaign" the project-level concept or the DataFrame-level concept?**
   Recommend: one screening campaign = one Datagrok project containing one DataFrame per
   run, plus a "campaign meta" DataFrame. Ties cleanly into Story 1's publish flow.
3. **Annotation flow when reusing v1.3 `Annotate Experiment`.** Sibling dialog for
   clarity, or a "screening mode" checkbox on the existing one? Recommend sibling.

---

## Cross-Story Feature Dependencies

```
v1.3 existing (do not re-research):
  parsers -> annotate (proteomics.groups) -> normalize -> impute -> DE -> volcano/heatmap/PCA/enrichment

Story 1: Read-only publish
    requires -> trimmed DataFrame contract (drop INTENSITY semType columns)
    requires -> audit-context panel (reads proteomics.* tags - already set in v1.3)
    requires -> Datagrok Spaces / project permissions (platform primitive)
    enhances -> can carry the v1.2 enrichment viewers + the Story-3 comparison
    enhances -> Story 3 (comparison view is itself publishable)

Story 2: SPC longitudinal QC
    requires -> "screening run" identity (shared with Story 3) - needs run_date tag
    requires -> control-group identification (reuses proteomics.groups from v1.3)
    requires -> v1.1 qc-computations.ts (per-run metric calculations - already built)
    enhances -> Story 3 (gates compound comparison on QC pass)

Story 3: Cross-run compound comparison
    requires -> "screening run" + "compound" entities (shared with Story 2)
    requires -> Annotate Screening Campaign dialog (new - sibling of experiment-setup)
    requires -> existing createVolcanoPlot (v1.3) - pair side-by-side
    requires -> Chem package (compound rendering, MCS) - must gracefully degrade
    enhances -> can be published via Story 1 once stable

Shared dependency between Story 2 and Story 3:
    "screening run" data model - proteomics.screening_run tag carries
        (run_date, vehicle_group, compound_groups[])
    This is the single most important shared primitive. Build it once, in a
    phase that lands BEFORE either Story 2 or Story 3 viewers.
```

### Dependency Notes

- **"Screening run" data model is the shared spine.** Both Story 2 (SPC needs to know
  "which week is this?") and Story 3 (comparison needs to know "what's in each run?")
  require a first-class run entity. Recommend a foundation phase that introduces the
  data model + Annotate Screening Campaign dialog *before* either viewer phase.
- **Story 1 is partially independent of Stories 2/3** — it can ship for any v1.3
  analysis without screening-campaign concepts. But the *published* artifact of a
  screening comparison (Story 3 output) is most valuable when Story 1 exists. So order
  matters less than for Stories 2/3.
- **Chem package dependency in Story 3 must gracefully degrade.** If Chem isn't loaded,
  show SMILES as text; MCS view shows a "Chem package required" placeholder. Don't
  hard-fail the comparison viewer.

---

## MVP Definition

### Launch With (v1.4 — Credible Cytokinetics Demo)

The minimum to demo all three stories convincingly to a Cytokinetics-class customer:

- [ ] **Story 1 minimum:** "Share Analysis for Review…" menu item; trimmed DataFrame
      (drop INTENSITY semTypes); audit-context panel; volcano in shared project; reviewer
      group permission; target-keyed Space or project tag (pick one based on prototyping).
- [ ] **Story 2 minimum:** Four standard metrics; I-chart + MR-chart per metric; Nelson
      rules 1–4 with default thresholds; run-level pass/flag tag; click-to-drill;
      baseline rolling window (default N=10).
- [ ] **Story 3 minimum:** `proteomics.screening_run` tag + Annotate Screening Campaign
      dialog; side-by-side volcano with linked selection + shared axes; diff table;
      compound structures rendered (Chem); top-impact ranking with default score rule;
      QC-status warning banner.

### Add After Validation (v1.4.x patches if customer feedback demands)

- [ ] Enrichment-overlap panel in Story 3 — likely first customer ask after the demo
- [ ] Pareto chart in Story 2 — adds little engineering, much QA dashboard value
- [ ] Phase-aware SPC mode — only needed once customer has 20+ runs and a real lot-change
- [ ] Comparison view publishable as a Story-1 artifact

### Future Consideration (v1.5+)

- [ ] N-way (3+) compound comparison in Story 3
- [ ] Curated QC-protein panel SPC in Story 2 (per-protein with FDR correction)
- [ ] CUSUM / EWMA charts
- [ ] Reviewer comment threads
- [ ] PDF report export
- [ ] First-class target taxonomy (target entity, target metadata, target search)
- [ ] Compound clustering / proteome-fingerprint SAR
- [ ] Email/Slack alerting on SPC flag
- [ ] MCS-highlighted compound difference in Story 3 (depends on Chem MCS being available)

---

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Story 1: Trimmed DataFrame + audit panel + share menu | HIGH | MEDIUM | P1 |
| Story 1: Target-keyed Space filing + reviewer permissions | HIGH | MEDIUM | P1 |
| Story 1: Enrichment carried with share | MEDIUM | LOW | P2 |
| Story 1: Comments / annotations | LOW | HIGH | P3 |
| Story 2: 4 standard metrics + I-chart + MR-chart | HIGH | LOW | P1 |
| Story 2: Nelson rules 1–4 + run-level tag | HIGH | MEDIUM | P1 |
| Story 2: Click-to-drill + baseline rolling-window | HIGH | LOW | P1 |
| Story 2: Pareto chart | MEDIUM | LOW | P2 |
| Story 2: Phase-aware mode | MEDIUM | MEDIUM | P2 |
| Story 2: Per-protein curated panel SPC | MEDIUM | HIGH | P3 |
| Story 2: CUSUM / EWMA | LOW | MEDIUM | P3 |
| Story 3: Screening-run data model + annotation | HIGH | MEDIUM | P1 |
| Story 3: Side-by-side volcano + linked selection + shared axes | HIGH | MEDIUM | P1 |
| Story 3: Diff table + top-impact ranking | HIGH | LOW | P1 |
| Story 3: Compound structure rendering (Chem) | HIGH | LOW | P1 |
| Story 3: QC-status warning banner | HIGH | LOW | P1 |
| Story 3: Enrichment overlap panel | MEDIUM | MEDIUM | P2 |
| Story 3: MCS-highlighted differences | MEDIUM | MEDIUM | P2 |
| Story 3: N-way comparison | LOW | HIGH | P3 |

---

## Competitor Feature Analysis

| Feature | Spectronaut Direct | Scaffold Viewer | Skyline + MSstatsQC | Our v1.4 Approach |
|---------|--------------------|-----------------|---------------------|--------------------|
| Read-only shared analysis | Vendor exports (templated PDF/Excel) + audit-friendly reporting | Cross-platform free viewer, mzIdentML export, Excel reports | N/A — Skyline is workbench-focused, not handoff | Trimmed Datagrok project — frozen, target-keyed, reviewer-permissioned, in-platform |
| Audit context | Sidebar — method, params, version | Annotations panel — free text | N/A | Read directly from `proteomics.*` tags — programmatic and human-readable |
| Volcano in shared view | Yes (own viewer) | No — Scaffold is identification-focused, not DE | No | Yes — same `createVolcanoPlot` from v1.3, thresholds drawn |
| Enrichment in shared view | Limited (vendor-specific) | No | No | Yes — v1.2 enrichment viewers carried in the project (differentiator) |
| Longitudinal SPC | Limited (within-vendor QC dashboard) | None | **Gold standard — MSstatsQC, SProCoP** | Match SProCoP/MSstatsQC metric set + rule set + chart pair; differentiator is in-Datagrok integration with Stories 1/3 |
| Nelson / Western Electric rules | Vendor-specific subset | None | Yes (MSstatsQC supports both) | Nelson 1–4 default; selectable |
| Per-protein SPC | No | No | Limited (per-peptide for targeted) | Curated panel only (deferred to v1.5) |
| Cross-run compound comparison | No — single-run focus | No — identification focus | No — system-suitability focus | **Differentiator — no prior art in this space; chemoproteomics teams hand-roll in Excel today** |
| Compound structure integration | No | No | No | Datagrok Chem package — SMILES + cell renderer + MCS |
| Linked side-by-side volcano | No | No | No | Reuse v1.3 createVolcanoPlot + selection-sync — pattern proven in v1.2 enrichment integration |
| Top-impact compound ranking | No | No | No | Configurable score rule baked into the screening-run tag |

The "compound comparison" column is the cleanest differentiator: no prior art ships
cross-run compound-aware proteome comparison. Cytokinetics-class customers (phenotypic /
chemoproteomic screening operations) currently do this in Excel; v1.4 makes it a first-class
in-platform workflow. The Story 1 and Story 2 columns put us at parity with established
tools but inside Datagrok, with the bonus of being integrated with Story 3.

---

## Sources

- Datagrok Spaces, Permissions, Projects — platform primitives (referenced from v1.3 in-repo docs)
- [Spectronaut — Biognosys](https://biognosys.com/software/spectronaut/) — vendor results-handoff feature set
- [Scaffold — Proteome Software](https://www.proteomesoftware.com/) — biologist-facing read-only viewer reference
- [Scaffold 5 / cross-platform sharing](https://www.proteomesoftware.com/products/scaffold-5)
- [SProCoP — Statistical Process Control in Proteomics for Skyline (PMC4020592)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4020592/) — five identification-free metrics + Pareto chart + Shewhart + Western Electric rules
- [MSstatsQC — Longitudinal System Suitability Monitoring (MCP 2020)](https://www.mcponline.org/article/S1535-9476(20)34194-3/fulltext) — I-chart + MR-chart + change-point + simultaneous-mean-and-variance monitoring
- [MSstatsQC GitHub — eralpdogu/MSstatsQC](https://github.com/eralpdogu/MSstatsQC) — R/Bioconductor + Shiny reference implementation
- [MSstatsQC 2.0 (J. Proteome Res.)](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00732) — DDA/DIA support, R/Bioc package
- [Multiyear Longitudinal Harmonization Study of QC in MS Proteomics (J. Proteome Res. 2024)](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00359) — community reference values, batch-spanning QC strategies
- [Western Electric Rules — implementation guide (Lab Wizard)](https://lab-wizard.com/en/resources/knowledge/spc-western-electric-rules/) — rule definitions
- [Nelson vs Western Electric (Parsec)](https://www.parsec-corp.com/blog/nelson-vs-western-electric) — 7-of-8 overlap, rule 4 difference
- [Chemoproteomic-enabled phenotypic screening (ScienceDirect)](https://www.sciencedirect.com/science/article/pii/S245194562100012X) — phenotypic-screen + chemoproteomics workflow context
- [Chemoproteomics, A Broad Avenue to Target Deconvolution (Adv Sci 2024)](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202305608) — SMILES-based compound organization in screens
- [Volcano Plots in Proteomics: Best Practices (MetwareBio)](https://www.metwarebio.com/volcano-plot-metabolomics-proteomics-guide/) — standard FC/p threshold conventions
- v1.3 in-repo references: `src/viewers/qc-computations.ts`, `src/viewers/volcano.ts`, `src/viewers/enrichment-viewers.ts`, `src/analysis/experiment-setup.ts`, `src/utils/proteomics-types.ts`
- v1.4 milestone source docs: `.planning/PROJECT.md`, `.planning/todos/pending/2026-05-11-*.md` (three v1.4 user-story todos)
- Cytokinetics 2026-05-04 customer meeting notes (referenced via todo files)

---
*Feature research for: Proteomics v1.4 Cross-Team Review milestone*
*Researched: 2026-06-06*

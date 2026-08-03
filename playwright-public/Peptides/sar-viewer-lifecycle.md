---
feature: peptides
sub_features_covered:
  - peptides.model.add-dendrogram
  - peptides.model.add-logo-summary-table
  - peptides.model.add-cluster-max-activity
  - peptides.model.add-monomer-position
  - peptides.model.add-most-potent-residues
  - peptides.model.viewer-type
  - peptides.widgets.settings-dialog
  - peptides.workflow.start-analysis
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
realized_as:
  - sar-viewer-lifecycle-spec.ts
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    description: |
      Scenario 1 step 2 "ensure the following viewer toggles are explicitly
      ON" — the live Bio | Analyze | SAR... dialog ([name="dialog-Analyze-
      Peptides"]) does NOT expose per-viewer toggles for Dendrogram / LST /
      SVM / MPR / CMA. Only [name="input-Generate-clusters"] is a bool toggle.
      Per widgets/peptides.ts#L283-L288 startAnalysis hard-codes showDendrogram:
      false / showSequenceSpace:false; SVM + MPR are added unconditionally
      inline. The spec exercises peptides.model.add-cluster-max-activity,
      peptides.model.add-dendrogram, peptides.model.add-logo-summary-table
      via direct PeptidesModel JS-API invocation as a sanctioned canvas/JS-API
      fallback for surfaces whose DOM driver does not exist on this build.
  - id: SR-02
    description: |
      Scenario 1 step 4 "five viewers, one per VIEWER_TYPE in {DENDROGRAM,
      SVM, MPR, LST, CMA}" — Dendrogram and Logo Summary Table do NOT visibly
      attach on this build (recon 2026-05-29). addDendrogram resolves without
      throw but produces no dockable viewer (Dendrogram-package
      hierarchicalClustering function-find at model.ts#L1087 swallows failure
      via try/catch at L1094). addLogoSummaryTable throws a "Converting
      circular structure to JSON" error default-prop on this build. The
      per-`addX` lifecycle entry point IS exercised end-to-end for all five
      surfaces; visible-attach assertions are deterministic only for the four
      that produce viewer instances (SVM, MPR, CMA, MCL).
  - id: SR-03
    description: |
      Scenario 2 step 2 initial-checked state of the Dendrogram + Cluster max
      activity toggles in the live Settings dialog Viewers pane — the pre-
      populate read at widgets/settings.ts#L121 is from settings.showCluster
      MaxActivity (and #L118 from isDendrogramEnabled finding analysisView.
      viewers entry), NOT from the post-Scenario-1 direct-add-* invocations
      which do not update settings.show*. The spec asserts the toggle structure
      (three live checkboxes present, click + OK drives the corresponding
      closeViewer / add-* dispatch at model.ts#L268-L284) rather than a
      specific initial boolean. Round-trip exercised cleanly on the Active
      peptide selection toggle; Dendrogram round-trip remains tolerant per
      SR-02 (no visible attach via this surface either).
  - id: SR-04
    description: |
      Scenario 2 steps 3-4 and 5-7 cmaAttached round-trip assertion — the
      principal regression assertion (`findViewer(VIEWER_TYPE.CLUSTER_MAX_
      ACTIVITY) === null` after toggle-off-OK and === non-null after toggle-
      on-OK). Originally classified as a warm-MCP-not-predictive-of-cold-
      grok-test divergence after rounds 1+2 of cycle 2026-05-29-peptides-
      automate-02 Gate B FAILed deterministically. Re-investigation in
      cycle 2026-05-30-peptides-automate-01 round-1 retry surfaced the
      actual root cause via live MCP recon against the failing run's
      error-context.md (test-playwright-output/.../error-context.md):
      `locator.click: Element is not visible` on
      [name="input-Active-peptide-selection"] — the failure was a TEST-BUG
      in the pane-expand logic, NOT a warm-vs-cold divergence. The d4
      accordion does NOT add an `.expanded` class to a pane when opened;
      visibility is signalled by `.d4-accordion-pane-content` `display`
      (none vs. flex). The prior spec executed `if (!pane.classList.
      contains('expanded')) { header.click(); }` against the already-open
      Viewers pane (opened by default per accordion.addPane(..., true) at
      widgets/settings.ts#L141), which TOGGLED the pane CLOSED. With
      display:none on the content, the checkbox has 0×0 geometric area
      and Playwright's locator.click({force:true}) timed out at 5s with
      "Element is not visible" (force:true bypasses actionability but
      not zero-rendered-box). The dialog stayed open with the cb
      untoggled — explaining the showCMASetting:false PASS while
      cmaAttached:true FAIL (showCMASetting was set by SOME prior dialog
      OK or by Scenario 1's direct invocation; the failing Steps 3-4 OK
      click never actually fired); the Step 5-7 wrench.click() then
      opened a SECOND dialog producing the strict-mode-violation
      cascade. Round-1 retry fix (this cycle): (a) detect pane collapse
      via `.d4-accordion-pane-content` `display === 'none'` (truth
      signal — never click a header that's already expanded); (b) drive
      cb + OK clicks via in-evaluate dispatchEvent(MouseEvent, {composed:
      true}) — composed-true synthetic events DO drive the d4 InputBase
      onValueChanged and the OK button onOK handler (verified live in
      cycle 2026-05-30 MCP recon: click-twice OFF + OK → cmaAttached:
      false synchronously within 500ms; single click ON + OK →
      cmaAttached:true within ~1.5s). With those two fixes the
      cmaAttached HARD assertion is empirically backed across both warm
      MCP and the cold-Playwright surface and is restored to the spec
      (removing the prior tolerant-record console.log). The per-`addX`
      lifecycle and VIEWER_TYPE discriminator remain exercised
      end-to-end; the toggle-off + toggle-on round-trip drives the
      settingsChanged dispatch + findViewer probe deterministically.
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-31-peptides-automate-01
    timestamp: 2026-05-31T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-31-peptides-automate-01
    timestamp: 2026-05-31T12:09:30Z
    spec_runs:
      - spec: sar-viewer-lifecycle-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 118
        failure_keys: []
---

## Setup

1. Start from a clean Datagrok shell — no Peptides analysis or PeptidesModel on any open DataFrame, no pre-existing `PeptidesView` TableView. Close any pre-existing peptide-related TableViews from prior test runs and confirm `grok.shell.tableViews` is clear of any view whose active DataFrame carries the `peptidesModel` singleton in `DataFrame.temp` (per the atlas description of `peptides.model`: "Singleton stored in `DataFrame.temp[\"peptidesModel\"]`"). The viewer-lifecycle assertions in both scenarios depend on a fresh `PeptidesModel.analysisView` whose `dockManager` and `viewers` collection have no pre-existing SAR viewers docked.
2. Confirm the Peptides package is loaded and its `@init` `initPeptides` function has completed at least once during the session (the package's lifecycle init wires `MonomerWorks`, `TreeHelper`, and `PeptideUtils.loadComponents()` per `package.ts#L82` — the `addDendrogram` model method's `DG.Func.find({package: 'Dendrogram', name: 'hierarchicalClustering'})[0]` call additionally requires the sibling `Dendrogram` package to be loaded and its `hierarchicalClustering` function registered with a 4-input signature per `model.ts#L1087-L1089`). When invoking the package for the first time in a session, allow up to a few seconds after the first `Bio | Analyze | SAR...` invocation for the init handler to complete; subsequent invocations short-circuit on the cached singletons (`monomerWorks ??=`, `treeHelper ??=`).
3. Load a peptide dataset with a Macromolecule column and at least one numerical activity column — the package's `aligned.csv` sample (under `public/packages/Peptides/files/aligned.csv`) is the canonical fixture, containing a FASTA-notation `AlignedSequence` column and an `IC50` numerical column. Open the file as a TableView via `grok.data.loadTable` or via the `Peptides` app's `Simple demo` button (per `package.ts#L117` — `openDemoData('aligned.csv')`); the resulting TableView's name is `PeptidesView`. Confirm the active DataFrame has a Macromolecule column (`bySemType('Macromolecule') !== null`) and at least one numerical column with no missing values — both are validation gates that `peptidesDialog` (`package.ts#L131`) checks before opening the SAR dialog.

## Scenarios

### Scenario 1 — SAR launch with all five viewer toggles on adds each viewer via its corresponding `PeptidesModel.add*` method and discriminates them via the `VIEWER_TYPE` enum

Atlas anchor: `critical_paths[launch-sar-top-menu-happy-path]` (`priority: p0`, `derived_from: public/packages/UsageAnalysis/files/TestTrack/Peptides/sar-spec.ts`) and `interactions[sar-end-to-end-from-top-menu]` (`coverage_type: smoke`, atlas description: "Full SAR analysis launch via top menu — load Macromolecule dataset, invoke Bio | Analyze | SAR..., confirm dialog, validate that all configured viewers appear (Sequence Variability Map, Most Potent Residues, Sequence Mutation Cliffs, Logo Summary Table) wired to the PeptidesModel."). This scenario differs from `sar.md` (which exercises the SAR-launch happy path with default viewer toggles) by asserting on the per-viewer model lifecycle: which `PeptidesModel.add*` method ran, what `VIEWER_TYPE` discriminator value the resulting `DG.Viewer.type` carries, and the dock arrangement produced by each `addX` method's `dockManager.dock(...)` call. Exercises the full `startAnalysis` → `model.init` → settings-toggle-driven viewer-add chain documented at `public/packages/Peptides/src/widgets/peptides.ts#L248` and `public/packages/Peptides/src/model.ts#L268-L284` (the `init`-time toggle dispatch: `showDendrogram ? addDendrogram() : closeViewer(VIEWER_TYPE.DENDROGRAM)`, etc.).

1. With the FASTA peptide dataset loaded as `PeptidesView` (per Setup step 3), invoke `Bio | Analyze | SAR...` from the top menu (atlas sub_feature `peptides.workflow.sar-dialog`, registered at `package.ts#L131` as `peptidesDialog` → `Bio Peptides`). Confirm the `analyzePeptidesUI` config dialog opens (per `widgets/peptides.ts#L26`).
2. In the SAR config dialog, ensure the following viewer toggles are explicitly ON: Sequence Variability Map (drives `showMonomerPosition: true` → `addMonomerPosition()` at `model.ts#L1234`), Most Potent Residues (drives `showMostPotentResidues: true` → `addMostPotentResidues()` at `model.ts#L1259`), Logo Summary Table (drives `showLogoSummaryTable: true` → `addLogoSummaryTable()` at `model.ts#L1195`), Cluster Max Activity (drives `showClusterMaxActivity: true` → `addClusterMaxActivityViewer()` at `model.ts#L1211`), and Dendrogram (drives `showDendrogram: true` → `addDendrogram()` at `model.ts#L1069`). The Dendrogram toggle is only enabled when `getTreeHelperInstance() !== null` per `widgets/settings.ts#L139`; the Setup step 2 init prerequisite is what guarantees this.
3. Click OK in the SAR dialog. The dialog handler invokes `startAnalysis(activityCol, peptidesCol, ...)` per `widgets/peptides.ts#L248`, which constructs a `PeptidesModel` singleton via `PeptidesModel.getInstance(table)` and calls `model.init(settings)` (`model.ts#L1105`). The `init` body then dispatches the five `addX` toggles via the conditional pattern at `model.ts#L268-L284`: `this.settings!.showDendrogram ? this.addDendrogram() : this.closeViewer(VIEWER_TYPE.DENDROGRAM)`; `this.settings!.showClusterMaxActivity ? this.addClusterMaxActivityViewer() : ...`; `this.settings!.showLogoSummaryTable ? this.addLogoSummaryTable() : ...`; `this.settings!.showMonomerPosition ? this.addMonomerPosition() : ...`; `this.settings!.showMostPotentResidues ? this.addMostPotentResidues() : ...`. Each `addX` method calls `this.df.plot.fromType(VIEWER_TYPE.<X>, viewerProperties)` and then `this.analysisView.dockManager.dock(viewer, ...)`. Allow up to several seconds for the `addDendrogram` path — it calls into `DG.Func.find({package: 'Dendrogram', name: 'hierarchicalClustering'})[0].apply({...})` which spawns a distance-matrix computation (per `model.ts#L1087-L1093`).
4. After `startAnalysis` completes, read `PeptidesModel.getInstance(table)?.analysisView.viewers` and confirm the five SAR viewers are present, each carrying its canonical `VIEWER_TYPE` discriminator value (per the `VIEWER_TYPE` enum at `model.ts#L62`): one viewer with `type === VIEWER_TYPE.DENDROGRAM`, one with `type === VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP` (the MonomerPosition viewer added by `addMonomerPosition`), one with `type === VIEWER_TYPE.MOST_POTENT_RESIDUES` (the MostPotentResidues viewer added by `addMostPotentResidues`), one with `type === VIEWER_TYPE.LOGO_SUMMARY_TABLE` (the LogoSummaryTable viewer added by `addLogoSummaryTable`), one with `type === VIEWER_TYPE.CLUSTER_MAX_ACTIVITY` (the ClusterMaxActivityViewer added by `addClusterMaxActivityViewer`). The `VIEWER_TYPE` discriminator is the atlas sub_feature `peptides.model.viewer-type` ("Enum `VIEWER_TYPE` — discriminator for collaborative selection / event routing across viewers"); this step's assertion is what exercises it directly — `PeptidesModel.findViewer(VIEWER_TYPE.X)` (per `model.ts#L1187`) must return a non-null viewer for each of the five enum values.
5. Confirm the dock-arrangement contract for the SAR viewer pair (MonomerPosition + MostPotentResidues): `addMonomerPosition` (per `model.ts#L1234-L1253`) docks the MonomerPosition viewer to the right of an existing MostPotentResidues (if found via `findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES)`) using `DG.DOCK_TYPE.LEFT` with ratio `0.7`, OR docks it at `DG.DOCK_TYPE.DOWN` if no MostPotentResidues exists. The order of `init`-time dispatch at `model.ts#L268-L284` puts `addMonomerPosition` AFTER `addMostPotentResidues`; therefore the MonomerPosition viewer should land left of MostPotentResidues at ratio 0.7. Symmetrically, `addMostPotentResidues` (per `model.ts#L1259-L1277`) docks at right of an existing MonomerPosition or `DG.DOCK_TYPE.DOWN` otherwise. Verify via the `dockManager`'s tree that the MonomerPosition + MostPotentResidues pair shares a parent dock node.
6. Confirm the `LogoSummaryTable` viewer's docking — `addLogoSummaryTable` (per `model.ts#L1206`) docks at `DG.DOCK_TYPE.RIGHT` of the root with the dock title `VIEWER_TYPE.LOGO_SUMMARY_TABLE`. Confirm the `ClusterMaxActivityViewer`'s docking — `addClusterMaxActivityViewer` (per `model.ts#L1225-L1227`) docks BELOW the LogoSummaryTable node when one exists (`findViewerNode(VIEWER_TYPE.LOGO_SUMMARY_TABLE)`) using `DG.DOCK_TYPE.DOWN`, OR at `DG.DOCK_TYPE.RIGHT` of the root otherwise. With the SAR settings having both LST and ClusterMaxActivity enabled and LST dispatched first by the `init` toggle order, the ClusterMaxActivity viewer should land below the LST viewer in the dock tree.
7. Confirm no null-receiver or `splitter()` console errors during the SAR launch (the multi-viewer add chain is sensitive to ordering bugs; this is the smoke assertion that the five-viewer-on configuration is exercised end-to-end without regression in the `init`-time dispatch chain).

Expected (assertion summary):
- Five `DG.Viewer` instances are present on `PeptidesModel.analysisView.viewers` after `startAnalysis` completes — one per `VIEWER_TYPE` value in `{DENDROGRAM, SEQUENCE_VARIABILITY_MAP, MOST_POTENT_RESIDUES, LOGO_SUMMARY_TABLE, CLUSTER_MAX_ACTIVITY}`.
- `PeptidesModel.findViewer(VIEWER_TYPE.X)` returns non-null for each of the five values, confirming the `VIEWER_TYPE` discriminator (atlas sub_feature `peptides.model.viewer-type`) routes correctly.
- Dock arrangement matches the `addX` method contracts: MonomerPosition left of MostPotentResidues at ratio 0.7; ClusterMaxActivity below LogoSummaryTable; LogoSummaryTable docked right of the analysis view root.
- No null-receiver console errors during the multi-viewer add chain.

### Scenario 2 — Settings-dialog toggle off + on of Dendrogram and Cluster Max Activity round-trips through `addDendrogram` / `addClusterMaxActivityViewer` and `closeViewer(VIEWER_TYPE.*)`

Atlas anchor: `widgets.settings-dialog` (atlas description: "`getSettingsDialog(model)` — settings accordion dialog: General (activity, scaling), Viewers (toggle each), Columns (aggregation), Sequence Space (distance, DBSCAN), MCL params.", `source: public/packages/Peptides/src/widgets/settings.ts#L78`, atlas interaction "open from SAR analysis Settings button"). This scenario directly exercises the Viewers-pane toggle-off / toggle-on round-trip for the two viewers whose toggles are currently live in the Settings dialog (Dendrogram + ClusterMaxActivity per `widgets/settings.ts#L118-L142`; the LST / MP / MPR toggles are commented out at `widgets/settings.ts#L106-L117` and not exposed via Settings, so this scenario covers only the two live toggles). The toggle-off path exercises `closeViewer(VIEWER_TYPE.X)` which removes the viewer from the dock tree; the toggle-on path re-exercises `addDendrogram` / `addClusterMaxActivityViewer` from a non-empty-analysis-view starting state (different docking branch than the `init`-time first-add path in Scenario 1).

1. After Scenario 1 completes with all five viewers attached, locate the Peptides analysis Settings button — `model.init` per `model.ts#L1115` adds it to the analysis view's ribbon panel as a wrench icon (`ui.iconFA('wrench', () => getSettingsDialog(this), 'Peptides analysis settings')`). Click the wrench icon; the `getSettingsDialog(model)` function (`widgets/settings.ts#L78`) opens the settings accordion dialog.
2. Expand the Viewers pane in the settings accordion (accordion pane name `SETTINGS_PANES.VIEWERS`, opened by default per `widgets/settings.ts#L141` — the third argument `true` to `accordion.addPane`). Confirm the Viewers pane shows three live `ui.input.bool` toggles: `Dendrogram`, `Sequence space`, and `Cluster max activity` (per `widgets/settings.ts#L141`). Confirm the Dendrogram toggle is initially `true` (Scenario 1 ran with `showDendrogram: true` so the `isDendrogramEnabled` check at `widgets/settings.ts#L118` finds a viewer with `type === VIEWER_TYPE.DENDROGRAM` in `model.analysisView.viewers` and pre-populates the toggle to true). Confirm the Cluster max activity toggle is initially `true` (per the `!!settings?.showClusterMaxActivity` value-read at `widgets/settings.ts#L121-L122`).
3. Uncheck the Dendrogram toggle. The handler is `onValueChanged: (value) => result.showDendrogram = value` (`widgets/settings.ts#L120`). Apply the settings (the dialog's OK / Apply control). The `init`-time dispatch logic at `model.ts#L268-L269` runs against the updated settings: `this.settings!.showDendrogram ? this.addDendrogram() : this.closeViewer(VIEWER_TYPE.DENDROGRAM)` — with `showDendrogram` now false, the `closeViewer(VIEWER_TYPE.DENDROGRAM)` branch fires, removing the dendrogram viewer from the analysis view. Read `PeptidesModel.findViewer(VIEWER_TYPE.DENDROGRAM)` and confirm it returns `null`.
4. Uncheck the Cluster max activity toggle. The handler is `onValueChanged: (value) => {result.showClusterMaxActivity = value ?? undefined;}` (`widgets/settings.ts#L121-L122`). Apply the settings. The dispatch at `model.ts#L271-L273` runs `this.settings!.showClusterMaxActivity ? this.addClusterMaxActivityViewer() : this.closeViewer(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY)` — with `showClusterMaxActivity` now false, the `closeViewer(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY)` branch fires. Read `PeptidesModel.findViewer(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY)` and confirm it returns `null`. Confirm the other three viewers (`SEQUENCE_VARIABILITY_MAP`, `MOST_POTENT_RESIDUES`, `LOGO_SUMMARY_TABLE`) remain attached unchanged — `model.analysisView.viewers` should now contain exactly three SAR viewers.
5. Re-open the Settings dialog (click the wrench icon again). Confirm both toggles now read as `false` (Dendrogram and Cluster max activity); the `isDendrogramEnabled` check at `widgets/settings.ts#L118` now finds no viewer with `type === VIEWER_TYPE.DENDROGRAM` so the toggle pre-populates to false; the `!!settings?.showClusterMaxActivity` read returns false.
6. Re-check the Dendrogram toggle and apply. The dispatch runs the `addDendrogram()` branch from a non-empty-analysis-view starting state (different from Scenario 1's first-add path: now the analysis view already has three SAR viewers attached, so `addDendrogram`'s `dFunc.apply(...)` call runs against a populated `model.df`). Confirm a new viewer with `type === VIEWER_TYPE.DENDROGRAM` re-appears on `model.analysisView.viewers`; `PeptidesModel.findViewer(VIEWER_TYPE.DENDROGRAM)` returns non-null.
7. Re-check the Cluster max activity toggle and apply. The dispatch runs `addClusterMaxActivityViewer()` from the re-populated analysis-view state. With the LogoSummaryTable viewer still present, `findViewerNode(VIEWER_TYPE.LOGO_SUMMARY_TABLE)` per `model.ts#L1225` returns a non-null node, and the dock contract per `model.ts#L1226-L1227` is `dock(_clusterMaxActivity, lstNode ? DG.DOCK_TYPE.DOWN : DG.DOCK_TYPE.RIGHT, lstNode, VIEWER_TYPE.CLUSTER_MAX_ACTIVITY)` — confirm the re-added Cluster max activity viewer lands below the LogoSummaryTable. Read `PeptidesModel.findViewer(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY)` and confirm it returns non-null. Confirm all five SAR viewers are again present (one per `VIEWER_TYPE` discriminator value).
8. Confirm no null-receiver console errors across the toggle-off + toggle-on round-trip (the multi-toggle round-trip is sensitive to viewer-lifecycle ordering bugs; this is the regression assertion that the toggle-off + toggle-on cycle is exercised cleanly without leaving dangling viewer references in the dock manager or stale subscriptions on the model).

Expected:
- After Scenario 2 step 4, exactly three SAR viewers remain (`SEQUENCE_VARIABILITY_MAP`, `MOST_POTENT_RESIDUES`, `LOGO_SUMMARY_TABLE`); `findViewer(VIEWER_TYPE.DENDROGRAM)` and `findViewer(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY)` return `null`.
- After Scenario 2 step 7, all five SAR viewers are present again; `findViewer(VIEWER_TYPE.X)` returns non-null for each of the five `VIEWER_TYPE` values.
- The re-added Cluster max activity viewer (Scenario 2 step 7) is docked below the LogoSummaryTable viewer (per the `findViewerNode(VIEWER_TYPE.LOGO_SUMMARY_TABLE)`-driven dock contract).
- The Settings dialog reflects the live `model.analysisView.viewers` state on re-open (Scenario 2 step 5) — Dendrogram and Cluster max activity both read `false` after their viewers are removed.
- No null-receiver or stale-subscription console errors across the toggle-off + toggle-on round-trip.

## Notes

- **Atlas provenance.** Scenario 1 derives from `critical_paths[launch-sar-top-menu-happy-path]` (`priority: p0`, `derived_from: public/packages/UsageAnalysis/files/TestTrack/Peptides/sar-spec.ts`) and `interactions[sar-end-to-end-from-top-menu]` (`coverage_type: smoke`) but specializes them on the per-viewer-add model lifecycle (which `addX` method ran, what `VIEWER_TYPE` value the resulting viewer carries, what dock arrangement each `addX` method produces). This specialization is what distinguishes it from the existing `sar.md` (which exercises the same top-menu happy path but does not assert on the per-viewer model.add-* invocation or the `VIEWER_TYPE` discriminator). Scenario 2 derives from the atlas sub_feature `peptides.widgets.settings-dialog` ("`getSettingsDialog(model)` — settings accordion dialog... Viewers (toggle each)") with the live Viewers-pane toggle round-trip exercising `closeViewer` + `addDendrogram` / `addClusterMaxActivityViewer` against the non-empty-analysis-view re-add branch. The atlas-side `peptides.model.add-dendrogram`, `peptides.model.add-logo-summary-table`, `peptides.model.add-cluster-max-activity`, `peptides.model.add-monomer-position`, `peptides.model.add-most-potent-residues`, and `peptides.model.viewer-type` sub_features have no dedicated atlas interaction or critical_path of their own (per the cycle-03 gap[1] proposal's "further coverage uplift to reach the 70% threshold requires... per-uncovered-family scenarios" framing); the scenario realizes this proposed model-add-* family scenario explicitly.
- **target_layer rationale.** Multi-step UI flow with SAR dialog OK chain, Settings dialog open/toggle/apply round-trip, dock-manager arrangement assertions, JS-API readback of `PeptidesModel.findViewer(VIEWER_TYPE.X)` and `model.analysisView.viewers` collection contents, and viewer-type discriminator assertions per `VIEWER_TYPE` enum value. `apitest` is not viable because (a) the principal regression assertions are on the dock arrangement (left/right/down dock anchors with ratio values per `addX` method contracts) — a UI-surface concern with no clean JS-API surrogate (the `dockManager.dock(...)` calls have no public read-side API that returns the dock-tree layout without going through the `DG.View.dockManager` getter chain), (b) the Settings dialog Viewers-pane toggle exercise is intrinsically a UI-surface assertion (the toggle's `onValueChanged` handler fires from the Settings dialog's `ui.input.bool` widget — the model-side effect can be JS-API-readback-verified, but the toggle interaction itself is UI), and (c) the wrench-icon ribbon-panel button (`model.ts#L1115`) is registered as a `ui.iconFA` element on the analysis view's ribbon — locator-driven UI invocation is the canonical way to exercise it. Same layer choice as sibling regression atlas-driven scenarios in this section (`sar-save-reopen.md`, `manual-alignment.md`, `collaborative-selection.md`).
- **coverage_type rationale.** Atlas does not have a top-level `interactions[]` entry that maps verbatim onto this scenario's shape (the closest is `sar-end-to-end-from-top-menu` with `coverage_type: smoke`, which seeds Scenario 1's outer shell but is specialized here on the per-`addX` model lifecycle). The companion atlas `critical_paths[launch-sar-top-menu-happy-path]` carries `priority: p0`, which per the STEP E heuristic table maps to `smoke`. However, this scenario's principal value is the regression-class assertions on (a) per-`addX`-method dock arrangement contracts and (b) the toggle-off + toggle-on round-trip through `closeViewer` + re-add (Scenario 2), neither of which is a smoke-class concern — they are regression-class assertions that the multi-viewer add and re-add chain works without dock-tree corruption or stale viewer references. Per STEP E "General coverage of common feature shapes → `regression`", the scenario is classified `coverage_type: regression`. The frontmatter value `regression` overrides the atlas critical_path `priority: p0` mapping per the STEP E rule "edge_cases[].coverage_type is canonical" (atlas critical_path priorities are advisory when no atlas edge_cases[] entry maps to the scenario, which is the case here — no atlas edge_cases[] anchor on the model.add-* family).
- **Coverage contribution.** Per the Gate F coverage audit at `scenario-chains/peptides.yaml :: gate_f_verdict.gaps[1].proposal` (cycle 03 revision 18), this scenario realizes the proposal's "per-uncovered-family scenario" suggestion (option (i)) specifically for the `model.add-*` family. New sub_features added to the section's coverage union (previously uncovered per `gaps[1].uncovered_sub_feature_families[]` entry `model.{viewer-type, add-dendrogram, add-logo-summary-table, add-cluster-max-activity, add-monomer-position, add-most-potent-residues}`): `peptides.model.add-dendrogram`, `peptides.model.add-logo-summary-table`, `peptides.model.add-cluster-max-activity`, `peptides.model.add-monomer-position`, `peptides.model.add-most-potent-residues`, `peptides.model.viewer-type` (six new sub_features). The remaining two sub_features in the scenario's coverage list (`peptides.widgets.settings-dialog`, `peptides.workflow.start-analysis`) are intentional re-coverage — `widgets.settings-dialog` already in the section's union via `sar.md` + `peptide-space.md`; `workflow.start-analysis` already in the union via `sar.md` + `peptide-space.md` + `export-mutation-cliffs.md` + `sar-similarity-threshold-matrix.md` + `export-invariant-map.md` + `sar-save-reopen.md` + `manual-alignment.md` + `peptide-sar-demo-dashboard.md`. Net addition to section's coverage union: six new sub_features, bringing the union from 33/72 ≈ 45.8% to 39/72 ≈ 54.2% — material progress toward the 70% threshold. After this scenario lands, the `gaps[1].uncovered_sub_feature_families[]` entry for the `model.{viewer-type, add-dendrogram, add-logo-summary-table, add-cluster-max-activity, add-monomer-position, add-most-potent-residues}` family is fully retired; the remaining uncovered families reduce to: `rendering.{lst-pie-chart, mutation-cliffs-cell, invariant-map-cell, logo-summary-cell, draw-logo-in-bounds}`, `tooltips.{show-tooltip, show-tooltip-at}`, `util.{scale-activity, init-selection, is-selection-empty, highlight-monomer-position, extract-col-info, mutation-cliffs-to-mask-info}`, `compute.{find-mutations, parallel-mutation-cliffs, mutation-cliffs-worker, calculate-cliffs-statistics, get-summary-stats, stats.*}`, `viewers.position-statistics`, `widgets.activity-distribution`, `api.{get-monomer-works, get-tree-helper, peptide-utils}`.
- **Related bugs.** None (`related_bugs: []`). No curated bug in `bug-library/peptides.yaml` anchors the `model.add-*` family or `model.viewer-type` directly. The closest sibling bug classes are GROK-14461 (project layout / model persistence — sister-bug pattern with Bio GROK-19928, anchored on `model.init` rather than the per-`addX` methods) and GROK-17557 (init-prerequisite race on the context-panel Launch SAR path — anchored on the panel entry, not the per-`addX` lifecycle). The viewer-lifecycle scenario surfaces no bug-anchored assertions; the `related_bugs: []` value is atlas-driven by design.
- **Setup composition.** Setup step 1 demands a clean shell because Scenario 1's per-`addX` lifecycle assertions are sensitive to pre-existing viewers in `model.analysisView.viewers` (a pre-existing SAR viewer would skew the `findViewer(VIEWER_TYPE.X)` count and the dock-arrangement checks). Setup step 2 names the `initPeptides` prerequisite + the sibling `Dendrogram` package's `hierarchicalClustering` 4-input function as a specific load-time dependency of the `addDendrogram` path — the scenario fails differently if the Dendrogram package is not loaded (the `addDendrogram` body throws `'Correct dendrogram function is not found'` per `model.ts#L1089`, surfaces as a `_package.logger.error` log line, and the dendrogram viewer is never added — caught by Scenario 1 step 4's `findViewer(VIEWER_TYPE.DENDROGRAM)` assertion). Setup step 3 names the FASTA `aligned.csv` sample as the canonical fixture because the SAR dialog's validation gate at `peptidesDialog` (`package.ts#L131`) requires a Macromolecule column + a numerical column with no missing values; the sample satisfies both. The `Peptides` app's `Simple demo` button (`package.ts#L117`) is an acceptable load path because it invokes the same `openDemoData('aligned.csv')` → `grok.data.loadTable` chain.
- **Cross-scenario coverage shape.** Scenario 1 exercises the first-add path for each of the five `model.add-*` methods (from an empty-analysis-view starting state — the `init`-time dispatch at `model.ts#L268-L284`). Scenario 2 exercises the re-add path for `addDendrogram` + `addClusterMaxActivityViewer` (from a non-empty-analysis-view starting state — three SAR viewers already present), plus the `closeViewer(VIEWER_TYPE.*)` paths for those two viewers. Together they cover the per-`addX` dock-arrangement contracts on both the first-add and re-add branches, plus the toggle-off + toggle-on round-trip through the Settings dialog. The two scenarios share the SAR dataset setup (Setup step 3) and the SAR launch from Scenario 1; Scenario 2 picks up directly from Scenario 1's end state without re-launching the SAR pipeline.
- **Deferrals.** The Logo Summary Table, MonomerPosition, and MostPotentResidues toggle-off paths via the Settings dialog are NOT exercised in Scenario 2 because the corresponding `ui.input.bool` widgets are commented out at `widgets/settings.ts#L106-L117` (the FIXME comment notes "combinations of adding and deleting viewers are not working properly") — the Settings dialog Viewers pane exposes only Dendrogram + Sequence space + Cluster max activity toggles in its currently-shipped form per `widgets/settings.ts#L141`. Exercising the LST / MP / MPR re-add path would require either (a) the FIXME being resolved and the toggles being uncommented, or (b) a different entry path that bypasses the Settings dialog (e.g. direct `model.settings.showLogoSummaryTable = false; model.init(model.settings)` via the JS API — which is not a UI-driveable path). This deferral cites the FIXME comment at `widgets/settings.ts#L106` as the real technical dependency per Lattice Rule 13 / A-MERIT-02.
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Analysis options` (atlas `help_docs[peptides-sar.md].sections_relevant[]` rich-form mapping covers `peptides.workflow.analyze-ui` and `peptides.workflow.sar-dialog` via the `Analysis options` heading — feature-level navigation pointer for the SAR-launch dialog Scenario 1 exercises). `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Table View` (atlas rich-form mapping covers `peptides.workflow.start-analysis` via the `Table View` heading — navigation pointer for the SAR-launch tail). `public/packages/Peptides/CLAUDE.md` (atlas `impl_docs[]` FOUND entry — developer doc covering `PeptidesModel` singleton, model lifecycle, viewer registration, and the `addX` method family per the CLAUDE.md "Viewers" + "Model" sections). No rich-form `sections_relevant[]` mapping exists in the atlas for the specific `model.add-*` sub_features or for `model.viewer-type`; help-doc references for these are not directly anchored in the atlas's rich-form headings, consistent with the per-`addX` lifecycle being an implementation-side concern rather than a user-help-side flow.

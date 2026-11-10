/**
 * nca-workflow-stages — mounted-view tests (review 2026-09-30: TEST-GAP-WS-01 / TEST-GAP-WS-02,
 * user decision "write now"). A REAL analysis view on the bundled theoph sample, computed once
 * through the pane's own compute path, then driven through the stages:
 *
 *  - AC-S3.9  — λz editing (point click AND Reset all) is allowed only on Review λz with an
 *               unchanged Setup (`lzEditingAllowed`, the one predicate both paths consult);
 *  - startReview — the first pending profile in problem-first order, the λz editor in front, the
 *               context panel on its OPEN Review pane; an explicit subject from Report too;
 *  - Report   — the Profiles panel is the Combined chart (build-once, same frame), titled with the
 *               table name; a row / point click on Report does NOT swap in the per-profile readout;
 *  - AC-S3.4.2 — Discard (formerly Revert) restores the confirmed recipe (the fingerprint matches
 *    again); an edited, confirmed Setup's primary button reads Apply;
 *  - first SAVE — a name dialog pre-filled with the table name; Cancel saves nothing;
 *  - Add to study — the new-study name defaults to the table name;
 *  - AC-S3.10.1 — the tour's anchors exist on their stages (Confirm on Setup; Settings and
 *               Add to study on Report);
 *  - the hub Ungrouped list — no Source column, a DISABLED (platform-grey) Move icon before the first
 *               compute and the active blue one after, no false `stale` badge on a never-computed
 *               entry, the count callback;
 *  - launch — the analysis view is titled with the table name (no name input, hub-ungrouped-polish).
 *
 * Runs under `grok test --host local` (live DG runtime + storage).
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {loadSampleByName, mountAnalysisIntoView, launchAnalysis} from '../views/nca-launch';
import {NcaAppView} from '../views/nca-app-view';
import {
  startReview, viewStageOf, liveWorkflowState, lzEditingAllowed, getConfirmed, setConfirmed,
  refreshStageUi, currentRecipeFingerprint, markAnalysisDirty, STUDY_BINDING_TAG,
  GOTO_STAGE_TAG, RESULTS_GRID_TAG, LZ_READOUT_PANEL_TAG, REPORT_PANEL_TAG, EDITOR_PLOT_TAG,
  MEAN_SD_CHART_TAG, PROFILE_CHART_TAG, LZ_READOUT_TAG, SPARSE_PARAM_GRID_TAG, FACET_GRID_DIV_TAG,
} from '../views/nca-mapping-pane';
import type {MappingPaneHandle} from '../views/nca-mapping-pane';
import type {ViewStage} from '../views/workflow-strip';
import {PROBLEM_RANK_COL} from '../views/results-view';
import {parseStudyPath, ncaStandaloneLinkPath} from '../views/nca-studies-router';
import {renderAnalysesListInto} from '../views/nca-analyses-list';
import {defaultAnalysisName} from '../views/nca-add-to-study';
import {HUB_BLUE, HUB_DISABLED_GREY, renderStudiesListInto} from '../views/nca-studies-list';
import {
  computeRecipeFingerprint, DEFAULT_RULES, NCA_STUDIES_CHANGED_EVENT, INSTALLED_ENGINE_VERSION, newId,
  getNcaRepository,
} from '../adapter';
import type {AnalysisEntry, NcaStudy, AnalysisRecipe, SummaryCache} from '../adapter';

const TIMEOUT = 120000;

/** Poll until `cond()` holds (or the budget runs out) — the pane docks / renders a tick later. */
async function waitFor(cond: () => boolean, ms = 8000): Promise<boolean> {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    if (cond()) return true;
    await DG.delay(100);
  }
  return cond();
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Bag = Record<string, any>;

/** A computed theoph analysis view (Review λz stage, profiles pending). */
async function computedTheoph(): Promise<{view: DG.TableView; df: DG.DataFrame; bag: Bag;
  handle: MappingPaneHandle}> {
  const df = await loadSampleByName('theoph.csv');
  const view = grok.shell.addTableView(df);
  const handle = mountAnalysisIntoView(view, df, {});
  await DG.delay(100); // deferred autodetect settles
  handle.triggerCompute();
  const bag = view.temp as Bag;
  const ok = await waitFor(() => viewStageOf(view) === 'review' && bag[RESULTS_GRID_TAG] != null);
  if (!ok) throw new Error('theoph did not reach Review λz after compute');
  return {view, df, bag, handle};
}

const goTo = (bag: Bag, s: ViewStage): void => (bag[GOTO_STAGE_TAG] as (st: ViewStage) => void)(s);

/** The pending review units (profiles still To review). */
function pendingKeys(view: DG.TableView): Set<string> {
  return new Set(liveWorkflowState(view).units
    .filter((u) => u.need === 'review' && u.decision === null).map((u) => u.key));
}

/** The visible dialog whose header reads `title` (DG dialogs mount under document.body). */
function dialogTitled(title: string): HTMLElement | null {
  for (const d of Array.from(document.querySelectorAll('.d4-dialog')) as HTMLElement[]) {
    const h = d.querySelector('.d4-dialog-header, .d4-dialog-title');
    // A dialog is `position: fixed` → its offsetParent is always null; "shown" = laid out.
    const r = d.getBoundingClientRect();
    if ((h?.textContent ?? '').includes(title) && d.isConnected && r.width > 0 && r.height > 0) return d;
  }
  return null;
}

function buttonIn(root: ParentNode, text: string): HTMLElement | undefined {
  return (Array.from(root.querySelectorAll('button')) as HTMLElement[])
    .find((b) => (b.textContent ?? '').trim().toLowerCase() === text.toLowerCase());
}

category('Workflow stages: live view (review 2026-09-30)', () => {
  test('AC-S3.9: λz editing only on Review λz with an unchanged Setup', async () => {
    const {view, bag} = await computedTheoph();
    try {
      expect(lzEditingAllowed(view), true, 'Review λz, Setup as confirmed');
      goTo(bag, 'report');
      expect(lzEditingAllowed(view), false, 'Report: no λz edit (point click or Reset all)');
      goTo(bag, 'setup');
      expect(lzEditingAllowed(view), false, 'Setup: no λz edit');
      goTo(bag, 'review');
      expect(lzEditingAllowed(view), true);
      // A Setup edit not yet confirmed: the held results no longer match the recipe.
      const c = getConfirmed(view)!;
      setConfirmed(view, {...c, fingerprint: `${c.fingerprint}-edited`});
      refreshStageUi(view);
      expect(liveWorkflowState(view).stage.setupChanged, true);
      expect(lzEditingAllowed(view), false, 'Setup changed since Confirm → locked');
      setConfirmed(view, c);
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('startReview: first pending profile (problem-first), λz editor in front, Review pane open', async () => {
    const {view, bag} = await computedTheoph();
    try {
      const pending = pendingKeys(view);
      expect(pending.size > 0, true, 'theoph has profiles to review');
      startReview(view);
      const grid = bag[RESULTS_GRID_TAG] as DG.Grid;
      const df = grid.dataFrame;
      const subj = df.col('SUBJECT')!;
      const rank = df.col(PROBLEM_RANK_COL)!;
      const cur = df.currentRowIdx;
      expect(cur >= 0, true, 'a row is selected');
      const target = subj.getString(cur);
      expect(pending.has(target), true, 'the selected profile is pending');
      // Problem-first: no pending profile ranks strictly worse-first than the chosen one.
      let minRank = Infinity;
      for (let i = 0; i < df.rowCount; i++)
        if (pending.has(subj.getString(i))) minRank = Math.min(minRank, Number(rank.get(i)));
      expect(Number(rank.get(cur)), minRank, 'the FIRST pending profile in problem-first order');
      expect(viewStageOf(view), 'review');
      const editorRoot = (bag[EDITOR_PLOT_TAG] as DG.Viewer | null)?.root ?? null;
      expect(await waitFor(() => editorRoot !== null && editorRoot.offsetParent !== null), true,
        'the λz editor tab is in front');
      const panel = bag[LZ_READOUT_PANEL_TAG];
      expect(panel?.title?.textContent, target, 'the context panel shows that profile');
      expect(panel.reviewPane.expanded, true, 'its Review pane is open');
      // An explicit subject — e.g. a clicked "To review" cell — from Report goes back to Review λz.
      const other = [...pending].find((k) => k !== target) ?? target;
      goTo(bag, 'report');
      startReview(view, other);
      expect(viewStageOf(view), 'review', 'starting a review moves to Review λz');
      expect(bag[LZ_READOUT_PANEL_TAG].title.textContent, other);
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('Report: Mean ± SD in front; Profiles = the Combined chart; row clicks keep it', async () => {
    const {view, df, bag} = await computedTheoph();
    try {
      startReview(view); // a readout exists before Report
      const before = bag[LZ_READOUT_PANEL_TAG].title.textContent as string;
      goTo(bag, 'report');
      const meanSdRoot = (bag[MEAN_SD_CHART_TAG] as DG.Viewer | null)?.root ?? null;
      expect(await waitFor(() => meanSdRoot !== null && meanSdRoot.offsetParent !== null), true,
        'Mean ± SD is the front tab');
      const panel = bag[REPORT_PANEL_TAG];
      expect(panel != null, true, 'the Report panel was built');
      expect(panel.pane.expanded, true, 'Profiles group open');
      expect(panel.title.textContent, df.name, 'titled with the table name');
      expect(panel.chart.dataFrame === (bag[PROFILE_CHART_TAG] as DG.Viewer).dataFrame, true,
        'the SAME frame as the Combined chart');
      // A row / point click on Report never swaps in the per-profile readout.
      const other = (bag[RESULTS_GRID_TAG] as DG.Grid).dataFrame.col('SUBJECT')!.getString(0) === before ?
        (bag[RESULTS_GRID_TAG] as DG.Grid).dataFrame.col('SUBJECT')!.getString(1) :
        (bag[RESULTS_GRID_TAG] as DG.Grid).dataFrame.col('SUBJECT')!.getString(0);
      (bag[LZ_READOUT_TAG] as (s: string) => void)(other);
      expect(bag[LZ_READOUT_PANEL_TAG].title.textContent, before, 'readout untouched on Report');
      // Re-entering Report re-uses the panel (build-once; closing an undocked viewer threw).
      goTo(bag, 'review');
      goTo(bag, 'report');
      // view.temp wraps a stored object per read — compare the panel's DOM root, not the wrapper.
      expect(bag[REPORT_PANEL_TAG].root === panel.root, true, 'the same panel (build-once)');
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('AC-S5.1.4: an Excluded profile reads (excl.) in the context-panel NCA parameters too', async () => {
    const {view, bag} = await computedTheoph();
    try {
      startReview(view);
      const target = bag[LZ_READOUT_PANEL_TAG].title.textContent as string;
      const {units, fingerprints} = liveWorkflowState(view);
      const u = units.find((x) => x.key === target)!;
      expect(u.allowedActions.includes('exclude'), true, 'the first pending theoph profile is Flagged');
      bag['ncaController'].decide(target, 'exclude', [{code: 'EXC_ERRATIC', text: 'erratic'}],
        fingerprints.get(target) ?? '');
      refreshStageUi(view);
      (bag[LZ_READOUT_TAG] as (s: string) => void)(target);
      const text = bag[LZ_READOUT_PANEL_TAG].host.textContent as string;
      expect(text.includes('(excl.)'), true, 'λz-dependent values are marked in the panel');
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('AC-S3.4.2: Discard restores the confirmed recipe; an edited Setup reads Apply', async () => {
    const {view, handle} = await computedTheoph();
    try {
      const c = getConfirmed(view)!;
      // The analyst's confirmed recipe differs from the live one (as after an unconfirmed edit):
      const confirmedRecipe = {...handle.getCurrentRecipe(), acceptanceRSquared: 0.92};
      setConfirmed(view, {...c, recipe: confirmedRecipe, fingerprint: computeRecipeFingerprint(confirmedRecipe)});
      // Still on the stage the compute landed on (Review): Discard is offered there too (2026-10-13).
      refreshStageUi(view);
      const offStage = buttonIn(view.root, 'Discard');
      expect(offStage !== undefined && offStage.offsetParent !== null, true, 'Discard is offered off Setup too');
      (view.temp as Bag)[GOTO_STAGE_TAG]('setup');
      refreshStageUi(view);
      expect(liveWorkflowState(view).stage.setupChanged, true, 'Setup differs from the confirmed recipe');
      const primary = view.root.querySelector<HTMLElement>('button[name="button-Confirm"]');
      expect((primary?.textContent ?? '').trim(), 'Apply', 'an edited, confirmed Setup reads Apply');
      const revert = buttonIn(view.root, 'Discard');
      expect(revert !== undefined && revert.offsetParent !== null, true, 'Discard is offered on Setup');
      revert!.click();
      await DG.delay(200);
      expect(handle.getCurrentRecipe().acceptanceRSquared, 0.92, 'the confirmed recipe is re-applied');
      expect(currentRecipeFingerprint(view), computeRecipeFingerprint(confirmedRecipe), 'fingerprints match');
      expect(liveWorkflowState(view).stage.setupChanged, false, 'no longer changed');
      expect((primary?.textContent ?? '').trim(), 'Confirm', 'Discard turns Apply back into Confirm');
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('first SAVE asks for a name pre-filled with the table name; Cancel saves nothing', async () => {
    const {view, df} = await computedTheoph();
    try {
      markAnalysisDirty(view);
      const save = (Array.from(document.querySelectorAll('button[name="button-Save"]')) as HTMLElement[])
        .find((b) => b.offsetParent !== null);
      expect(save !== undefined, true, 'SAVE on the ribbon');
      save!.click();
      expect(await waitFor(() => dialogTitled('Save analysis') !== null, 4000), true, 'the name dialog opens');
      const dlg = dialogTitled('Save analysis')!;
      const input = dlg.querySelector('input') as HTMLInputElement;
      expect(input.value, df.name, 'pre-filled with the table name');
      buttonIn(dlg, 'Cancel')!.click();
      await DG.delay(200);
      expect((view.temp as Bag)[STUDY_BINDING_TAG] == null, true, 'Cancel saved nothing (still unbound)');
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('BUG-16: SAVE at Setup, then Confirm — the compute is written to the saved analysis', async () => {
    // User QA 2026-10-06: open a table → SAVE (no results yet) → Confirm raised "Could not save
    // updated summary: no open project" and the Ungrouped entry stayed at Setup with no results.
    const df = await loadSampleByName('theoph.csv');
    df.name = `ws-bug16-${newId().slice(0, 8)}`;
    const view = grok.shell.addTableView(df);
    let containerId = '';
    try {
      mountAnalysisIntoView(view, df, {});
      await DG.delay(100); // deferred autodetect settles
      const bag = view.temp as Bag;
      markAnalysisDirty(view);
      const visible = (name: string): HTMLElement | undefined =>
        (Array.from(document.querySelectorAll(`button[name="${name}"]`)) as HTMLElement[])
          .find((b) => b.offsetParent !== null);
      visible('button-Save')!.click();
      expect(await waitFor(() => dialogTitled('Save analysis') !== null, 4000), true, 'the name dialog opens');
      buttonIn(dialogTitled('Save analysis')!, 'OK')!.click();
      expect(await waitFor(() => bag[STUDY_BINDING_TAG] != null, 10000), true, 'bound to the saved analysis');
      const binding = bag[STUDY_BINDING_TAG] as {studyId: string; analysisId: string};
      containerId = binding.studyId;
      const repo = getNcaRepository();
      expect((await repo.loadAnalysis(binding.analysisId)).entry.summaryCache, null, 'saved at Setup: no results yet');

      expect(viewStageOf(view), 'setup');
      visible('button-Confirm')!.click();
      expect(await waitFor(() => viewStageOf(view) === 'review', 15000), true, 'Confirm computed');
      // The writeback is async (fire-and-forget after the compute) — poll the stored entry.
      let entry: AnalysisEntry | null = null;
      const t0 = Date.now();
      while (Date.now() - t0 < 12000) {
        const a = (await repo.loadAnalysis(binding.analysisId)).entry;
        if (a.summaryCache !== null) {entry = a; break;}
        await DG.delay(200);
      }
      expect(entry !== null, true, 'the summary cache was written to the saved analysis');
      expect(entry!.stageCache?.stage, 'review', 'the stored stage follows the view (hub reads Review)');
      expect(entry!.engineVersion, INSTALLED_ENGINE_VERSION, 'stamped with the installed engine');
    } finally {
      if (containerId !== '')
        try {await getNcaRepository().deleteStudy(containerId);} catch {/* gone */}
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('Add to study: the new-study name defaults to the table name', async () => {
    const {view, df, bag} = await computedTheoph();
    try {
      goTo(bag, 'report');
      const add = (Array.from(document.querySelectorAll('button')) as HTMLElement[])
        .find((b) => (b.textContent ?? '').trim().toLowerCase() === 'add to study' && b.offsetParent !== null);
      expect(add !== undefined, true, 'Add to study on the Report ribbon');
      add!.click();
      expect(await waitFor(() => dialogTitled('Add analysis to study') !== null, 8000), true, 'dialog opens');
      const dlg = dialogTitled('Add analysis to study')!;
      // theoph has profiles To review → the D14 warning first; Add anyway reveals the inputs.
      const anyway = buttonIn(dlg, 'Add anyway');
      if (anyway !== undefined && anyway.getBoundingClientRect().width > 0) {anyway.click(); await DG.delay(200);}
      const nameRoot = (Array.from(dlg.querySelectorAll('.ui-input-root')) as HTMLElement[])
        .find((r) => (r.querySelector('.ui-input-label')?.textContent ?? '').trim() === 'New study name');
      const input = nameRoot?.querySelector('input') as HTMLInputElement | null;
      expect(input?.value, df.name, 'default = the table name');
      buttonIn(dlg, 'Cancel')!.click();
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('AC-S3.10.1: tour anchors exist on their stages', async () => {
    const {view, bag} = await computedTheoph();
    try {
      goTo(bag, 'setup');
      const pane = document.querySelector('.nca-mapping-pane') as HTMLElement;
      const confirm = buttonIn(pane, 'Confirm');
      expect(confirm !== undefined && confirm.isConnected, true, 'Setup: Confirm (the "Confirm the setup" step)');
      // Pinned (2026-10-13): Confirm lives in its own action form, OUTSIDE the scrolling mapping form.
      expect(confirm!.closest('.nca-mapping-actions') !== null && confirm!.closest('.nca-mapping-pane-scroll') === null,
        true, 'Confirm sits in the pinned action form, not in the scroll area');
      goTo(bag, 'report');
      await DG.delay(200);
      const settings = document.querySelector('.fa-cog') as HTMLElement | null;
      expect(settings !== null && settings.isConnected, true, 'Report: Settings');
      const add = (Array.from(document.querySelectorAll('button')) as HTMLElement[])
        .find((b) => (b.textContent ?? '').trim().toLowerCase() === 'add to study');
      expect(add !== undefined && add.isConnected, true, 'Report: Add to study');
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});
});

category('Workflow stages: destructive design (manual QA 2026-10-01)', () => {
  test('no λz editor / Facet; Sparse NCA in front; the review starts at the first composite', async () => {
    const df = await loadSampleByName('04_mouse_sparse_destructive.csv');
    const view = grok.shell.addTableView(df);
    try {
      const handle = mountAnalysisIntoView(view, df, {});
      await DG.delay(100);
      handle.triggerCompute();
      const bag = view.temp as Bag;
      expect(await waitFor(() => bag[SPARSE_PARAM_GRID_TAG] != null && viewStageOf(view) !== 'setup', 15000),
        true, 'the sparse composite table docked');
      await DG.delay(300); // the facet docks from a deferred timeout — let it (not) happen
      expect(bag[EDITOR_PLOT_TAG] == null, true, 'no λz editor in a destructive design');
      expect(bag[FACET_GRID_DIV_TAG] == null, true, 'no Facet in a destructive design');
      const sparse = bag[SPARSE_PARAM_GRID_TAG] as DG.Grid;
      expect(await waitFor(() => sparse.root.offsetParent !== null), true, 'Sparse NCA is the front tab');
      const pendingComposite = liveWorkflowState(view).units
        .find((u) => u.kind === 'composite' && u.need === 'review' && u.decision === null);
      if (pendingComposite !== undefined) {
        sparse.dataFrame.currentRowIdx = -1;
        startReview(view);
        expect(viewStageOf(view), 'review');
        expect(sparse.dataFrame.currentRowIdx >= 0, true, 'a composite row is selected for review');
      }
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});
});

category('Workflow stages: router + hub (review 2026-09-30)', () => {
  test('AC-S2.6.1: analyses/<label> parses to the standalone state and round-trips its link', async () => {
    expect(JSON.stringify(parseStudyPath('analyses/My%20NCA')),
      JSON.stringify({kind: 'standalone-analysis', label: 'My NCA'}));
    expect(parseStudyPath('analyses').kind, 'landing', 'the bare segment is the hub');
    const link = ncaStandaloneLinkPath('My NCA');
    const tail = link.slice(link.indexOf('analyses/'));
    expect(JSON.stringify(parseStudyPath(tail)),
      JSON.stringify({kind: 'standalone-analysis', label: 'My NCA'}), 'the copied link resolves back');
  });

  // ── study-add-multiple-analyses: the Add-to-study button stays active once in a study ──────────
  // D5 "one copy" reversed (user-directed 2026-10-05): one source table can back several study
  // analyses, so the button must not dim / no-op after the first add.
  const findAddButton = (): HTMLElement | undefined =>
    (Array.from(document.querySelectorAll('button')) as HTMLElement[])
      .find((b) => (b.textContent ?? '').trim().toLowerCase() === 'add to study' && b.offsetParent !== null);

  function setStringInput(root: ParentNode, label: string, value: string): boolean {
    const r = (Array.from(root.querySelectorAll('.ui-input-root')) as HTMLElement[])
      .find((x) => (x.querySelector('.ui-input-label')?.textContent ?? '').trim() === label);
    const inp = r?.querySelector('input') as HTMLInputElement | null;
    if (inp == null) return false;
    inp.value = value;
    inp.dispatchEvent(new Event('input', {bubbles: true}));
    return true;
  }

  function setChoiceByLabel(root: ParentNode, label: string, optionText: string): boolean {
    const r = (Array.from(root.querySelectorAll('.ui-input-root')) as HTMLElement[])
      .find((x) => (x.querySelector('.ui-input-label')?.textContent ?? '').trim() === label);
    const sel = r?.querySelector('select') as HTMLSelectElement | null;
    if (sel == null) return false;
    const opt = Array.from(sel.options).find((o) => (o.textContent ?? '').trim() === optionText);
    if (opt === undefined) return false;
    sel.value = opt.value;
    sel.dispatchEvent(new Event('input', {bubbles: true}));
    sel.dispatchEvent(new Event('change', {bubbles: true}));
    return true;
  }

  /** Click the New/Existing radio option whose text matches `optionText`. */
  function selectRadio(root: ParentNode, optionText: string): boolean {
    const radios = Array.from(root.querySelectorAll('input[type="radio"]')) as HTMLInputElement[];
    for (const radio of radios) {
      const text = (radio.closest('label')?.textContent ?? radio.parentElement?.textContent ?? '').trim();
      if (text === optionText || text.includes(optionText)) {radio.click(); return true;}
    }
    return false;
  }

  /** Click the dialog's primary (OK) button — its text is now case-dependent. */
  function clickPrimary(root: ParentNode): boolean {
    const b = (Array.from(root.querySelectorAll('button')) as HTMLElement[]).find((x) => {
      const t = (x.textContent ?? '').trim().toLowerCase();
      return (t === 'add to study' || t === 'create study & add') && x.offsetParent !== null;
    });
    if (b === undefined) return false;
    b.click();
    return true;
  }

  /** Report → Add to study → (Add anyway) → new study named `studyName` → OK; returns the new id. */
  async function firstAddToNewStudy(view: DG.TableView, bag: Bag, studyName: string): Promise<string> {
    goTo(bag, 'report');
    findAddButton()!.click();
    expect(await waitFor(() => dialogTitled('Add analysis to study') !== null, 8000), true, 'dialog opens');
    const dlg = dialogTitled('Add analysis to study')!;
    const anyway = buttonIn(dlg, 'Add anyway');
    if (anyway !== undefined && anyway.getBoundingClientRect().width > 0) {anyway.click(); await DG.delay(200);}
    selectRadio(dlg, 'New'); // default, but be explicit if other studies already exist
    await DG.delay(50);
    expect(setStringInput(dlg, 'New study name', studyName), true, 'new-study-name field');
    expect(clickPrimary(dlg), true, 'the primary button (Create study & add) is clickable');
    expect(await waitFor(() => (view.temp as Bag)[STUDY_BINDING_TAG] != null, 10000), true, 'bound to the study');
    return ((view.temp as Bag)[STUDY_BINDING_TAG] as {studyId: string}).studyId;
  }

  test('AC-1/AC-3: the Add button stays active once in a study and re-opens the dialog', async () => {
    const {view, bag} = await computedTheoph();
    let studyId = '';
    try {
      studyId = await firstAddToNewStudy(view, bag, `ws-add-${newId().slice(0, 8)}`);
      await waitFor(() => dialogTitled('Add analysis to study') === null, 4000);
      const add = findAddButton();
      expect(add !== undefined, true, 'Add to study is still on the Report ribbon');
      expect(add!.style.opacity !== '0.4', true, 'it stays active (not dimmed) after the analysis is in a study');
      expect(add!.style.cursor !== 'not-allowed', true, 'the cursor is not not-allowed');
      add!.click(); // AC-1: no longer a no-op
      expect(await waitFor(() => dialogTitled('Add analysis to study') !== null, 8000), true,
        're-clicking opens the Add-to-study dialog again');
      buttonIn(dialogTitled('Add analysis to study')!, 'Cancel')!.click();
    } finally {
      try {await getNcaRepository().deleteStudy(studyId);} catch {/* already gone */}
      view.close();
    }
  }, {timeout: TIMEOUT});

  function mkRecipe(o: {routeSource?: 'value' | 'column' | 'none'; routeScalar?: string;
    doseSource?: 'value' | 'column' | 'none'; doseScalar?: number; doseUnit?: string;
    groupByCols?: string[];}): AnalysisRecipe {
    return {
      columns: {
        subject: 'S', time: 'T', concentration: 'C',
        doseSource: o.doseSource ?? 'value', doseCol: null, doseScalar: o.doseScalar ?? 1,
        routeSource: o.routeSource ?? 'value', routeScalar: o.routeScalar ?? 'PO', routeCol: null,
        infusionSource: 'none', infusionDurationCol: null, infusionDurationScalar: 0,
        bodyWeightSource: 'none', bodyWeightCol: null, bodyWeightScalar: 0,
        blqFlagCol: null, lloqSource: 'none', lloqCol: null, lloqScalar: 0,
        sexCol: null, sexValueMap: {},
      },
      modes: {nonmem: false, confirmedRoute: null,
        sparse: {on: false, nominalTimeCol: null, animalIdCol: null, samplingType: '(auto)'}},
      ncaRules: DEFAULT_RULES,
      units: {time: 'h', concentration: 'mg/L', dose: o.doseUnit ?? 'mg', bodyWeight: 'kg'},
      groupByCols: o.groupByCols ?? [], passthroughCols: [], acceptanceRSquared: 0.9,
    };
  }

  test('AC-6: defaultAnalysisName encodes a short settings spec (pure)', async () => {
    expect(defaultAnalysisName(mkRecipe({routeScalar: 'PO', doseScalar: 1, doseUnit: 'mg'})), 'NCA · PO 1 mg');
    expect(defaultAnalysisName(mkRecipe({routeScalar: 'IV', doseScalar: 2, doseUnit: 'mg/kg'})), 'NCA · IV 2 mg/kg');
    expect(defaultAnalysisName(mkRecipe({routeScalar: 'PO', doseSource: 'column'})), 'NCA · PO',
      'a per-row dose column is omitted');
    expect(defaultAnalysisName(mkRecipe({routeScalar: 'PO', doseScalar: 1, groupByCols: ['Sex']})),
      'NCA · PO 1 mg · by Sex');
    expect(defaultAnalysisName(mkRecipe({routeSource: 'column', doseSource: 'column', groupByCols: ['Dose', 'Sex']})),
      'NCA · by Dose, Sex', 'both scalars absent → NCA + groups only');
    expect(defaultAnalysisName(mkRecipe({routeSource: 'column', doseSource: 'column'})), 'NCA',
      'nothing scalar, no groups → bare NCA');
  });

  test('AC-2: a second Add to the SAME study creates a second, distinct analysis from one table', async () => {
    const {view, bag} = await computedTheoph();
    let studyId = '';
    try {
      const label = `ws-multi-${newId().slice(0, 8)}`;
      studyId = await firstAddToNewStudy(view, bag, label);
      const study1 = await getNcaRepository().loadStudy(studyId);
      expect(study1.analyses.length, 1, 'first Add → one analysis');
      const firstId = study1.analyses[0].id;

      await waitFor(() => dialogTitled('Add analysis to study') === null, 4000);
      findAddButton()!.click();
      expect(await waitFor(() => dialogTitled('Add analysis to study') !== null, 8000), true, 'dialog reopens');
      const dlg = dialogTitled('Add analysis to study')!;
      const anyway = buttonIn(dlg, 'Add anyway');
      if (anyway !== undefined && anyway.getBoundingClientRect().width > 0) {anyway.click(); await DG.delay(200);}
      expect(selectRadio(dlg, 'Existing'), true, 'the Existing radio option is present');
      await DG.delay(100);
      expect(setChoiceByLabel(dlg, 'Existing study', label), true, `the existing study '${label}' is selectable`);
      await DG.delay(50);
      expect(clickPrimary(dlg), true, 'the primary button (Add to study) is clickable');

      let study2: NcaStudy | null = null;
      const t0 = Date.now();
      while (Date.now() - t0 < 12000) {
        const s = await getNcaRepository().loadStudy(studyId).catch(() => null);
        if (s !== null && s.analyses.length === 2) {study2 = s; break;}
        await DG.delay(200);
      }
      expect(study2 !== null, true, 'the second Add appended a 2nd analysis to the SAME study');
      const ids = study2!.analyses.map((a) => a.id);
      expect(new Set(ids).size, 2, 'the two analyses have distinct ids');
      expect(ids.includes(firstId), true, 'the first analysis is still present');
      expect(new Set(study2!.analyses.map((a) => a.source.ref)).size, 2, 'each analysis has its own snapshot');
      for (const a of study2!.analyses) expect(a.source.kind, 'stored', 'a stored snapshot');
    } finally {
      try {await getNcaRepository().deleteStudy(studyId);} catch {/* already gone */}
      view.close();
    }
  }, {timeout: TIMEOUT});

  // Hub auto-refresh (user-directed 2026-10-05): the reactive CONTRACT the Hub + Browse-tree both
  // stand on — the coordinator fires NCA_STUDIES_CHANGED_EVENT on every save / delete, so adding an
  // entity updates the Hub with no manual refresh (the Hub's subscription is now unconditional). The
  // Hub DOM itself cannot be rendered under `grok test` — `buildHeader` reads `_package.webRoot`,
  // which throws in the harness — so the icon removal is covered by code review + manual QA; this
  // test guards the mechanism that drives the refresh.
  test('Hub/tree reactive contract: repository create + delete fire NCA_STUDIES_CHANGED_EVENT', async () => {
    let fired = 0;
    const sub = grok.events.onCustomEvent(NCA_STUDIES_CHANGED_EVENT).subscribe(() => {fired++;});
    let id = '';
    try {
      const before = fired;
      id = (await getNcaRepository().createStudy(`ws-evt-${newId().slice(0, 8)}`)).id;
      expect(await waitFor(() => fired > before, 8000), true, 'create fires the change event (Hub + tree refresh)');
      const afterSave = fired;
      await getNcaRepository().deleteStudy(id);
      expect(await waitFor(() => fired > afterSave, 8000), true, 'delete fires the change event');
    } finally {
      sub.unsubscribe();
      if (id !== '') try {await getNcaRepository().deleteStudy(id);} catch {/* already gone */}
    }
  }, {timeout: TIMEOUT});

  test('hub Ungrouped: no Source column; Move disabled (grey) before compute; no false stale badge', async () => {
    const repo = getNcaRepository();
    const label = `ws-live-${newId().slice(0, 8)}`;
    // A second, COMPUTED entry (Report, current engine) — its Move icon must stay the active blue (AC-5).
    const label2 = `ws-live-${newId().slice(0, 8)}`;
    const made: string[] = [];
    try {
      const recipe = {
        columns: {
          subject: 'Subject', time: 'Time', concentration: 'conc',
          doseSource: 'value', doseCol: null, doseScalar: 1,
          routeSource: 'value', routeScalar: 'PO', routeCol: null,
          infusionSource: 'none', infusionDurationCol: null, infusionDurationScalar: 0,
          bodyWeightSource: 'none', bodyWeightCol: null, bodyWeightScalar: 0,
          blqFlagCol: null, lloqSource: 'none', lloqCol: null, lloqScalar: 0,
          sexCol: null, sexValueMap: {},
        },
        modes: {nonmem: false, confirmedRoute: null,
          sparse: {on: false, nominalTimeCol: null, animalIdCol: null, samplingType: '(auto)'}},
        ncaRules: DEFAULT_RULES,
        units: {time: 'h', concentration: 'mg/L', dose: 'mg', bodyWeight: 'kg'},
        groupByCols: [], passthroughCols: [], acceptanceRSquared: 0.9,
      } as unknown as AnalysisRecipe;
      const df = (): DG.DataFrame => DG.DataFrame.fromCsv('Subject,Time,conc\nA,0,1\nA,1,5\n');
      made.push((await repo.createStandalone({label, df: df(), recipe, summaryCache: null, origin: null,
        workflow: {confirmed: null, lzOverrides: [], reviewDecisions: [],
          stageCache: {stage: 'setup', pending: 0, notApplicable: 0}}})).studyId);
      // The list only checks the cache is present (canMove) — a minimal summary stands in.
      const cache: SummaryCache = {analysisTier: {nCompounds: 1, nSubjects: 1, route: 'PO', dose: '1 mg',
        status: 'ok', flagSummary: ''}, compoundTier: []};
      made.push((await repo.createStandalone({label: label2, df: df(), recipe, summaryCache: cache, origin: null,
        workflow: {confirmed: null, lzOverrides: [], reviewDecisions: [],
          stageCache: {stage: 'complete', pending: 0, notApplicable: 0}}})).studyId);
      const host = document.createElement('div');
      document.body.appendChild(host);
      let counted = -1;
      const n = await renderAnalysesListInto(host, (k) => {counted = k;});
      expect(counted, n, 'onCount receives the rendered count');
      const headers = Array.from(host.querySelectorAll('th')).map((t) => (t.textContent ?? '').trim());
      expect(headers.includes('Source'), false, 'no Source column');
      const rows = Array.from(host.querySelectorAll('tr')) as HTMLElement[];
      const row = rows.find((r) => (r.textContent ?? '').includes(label));
      expect(row !== undefined, true, 'our analysis is listed');
      expect((row!.textContent ?? '').includes('stale'), false, 'a never-computed entry is not "stale"');
      const icons = Array.from(row!.querySelectorAll('i')) as HTMLElement[];
      expect(icons[0]?.classList.contains('fa-layer-plus'), true, 'Move to study is the FIRST icon');
      // AC-4 (hub-ungrouped-polish): disabled = the platform disabled grey, NOT a faded blue.
      expect(icons[0]?.style.color, HUB_DISABLED_GREY, 'disabled before the first compute: platform grey');
      expect(icons[0]?.style.opacity, '', 'no reduced opacity (the grey carries the disabled state)');
      expect(icons[0]?.style.cursor, 'not-allowed', 'not-allowed cursor while disabled');
      // AC-5: a computed entry keeps the active blue Move icon.
      const row2 = rows.find((r) => (r.textContent ?? '').includes(label2));
      expect(row2 !== undefined, true, 'the computed analysis is listed');
      const move2 = row2!.querySelector('i') as HTMLElement | null;
      expect(move2?.style.color, HUB_BLUE, 'a computed entry: Move is the active blue');
      expect(move2?.style.cursor !== 'not-allowed', true, 'and its cursor is not not-allowed');
      // Overlapping renders of one host (the hub's change-event refresh + the list's own refresh
      // after a delete) must leave ONE table, not two (user QA 2026-10-06).
      await Promise.all([renderAnalysesListInto(host), renderAnalysesListInto(host)]);
      expect(host.querySelectorAll('table').length, 1, 'overlapping renders leave one table');
      host.remove();
    } finally {
      for (const id of made) try {await repo.deleteStudy(id);} catch {/* already gone */}
    }
  }, {timeout: TIMEOUT});

  test('hub Studies: overlapping renders of one host leave one list', async () => {
    // The duplicated Studies table after a delete (user QA 2026-10-06): two refreshes overlapped and
    // each cleared the host before either appended. Only the newest render may write the host.
    const host = document.createElement('div');
    document.body.appendChild(host);
    try {
      await Promise.all([renderStudiesListInto(host), renderStudiesListInto(host)]);
      expect(host.children.length, 1, 'one table (or one "No studies yet."), not two');
    } finally {
      host.remove();
    }
  }, {timeout: TIMEOUT});

  test('hub New analysis: no Analysis name input; the view is titled with the table name', async () => {
    // AC-1: the Hub form carries the Table input only. `buildNewAnalysis` is private but needs no
    // `_package.webRoot` (unlike `buildHeader`), so it renders under the harness.
    const hub = new NcaAppView();
    const form = (hub as unknown as {buildNewAnalysis(): HTMLElement}).buildNewAnalysis();
    expect((form.textContent ?? '').includes('Table'), true, 'the Table input is there');
    expect((form.textContent ?? '').includes('Analysis name'), false, 'no Analysis name input');
    // AC-3: the launch path the Hub Run / study dialog call titles the view with the table name.
    const df = await loadSampleByName('theoph.csv');
    df.name = `hub-launch-${newId().slice(0, 8)}`;
    const w = grok.shell.windows;
    const saved = {toolbox: w.showToolbox, browse: w.showBrowse, auto: w.autoShowToolbox};
    // BUG-18 guard: start from the state that let the Toolbox take the side panel — auto-show on,
    // Browse not in front. The pre-fix open path touched neither, so these asserts fail on it.
    // (The large-table re-show on results docking is not reproducible headless — manual QA.)
    w.autoShowToolbox = true;
    w.showBrowse = false;
    const view = launchAnalysis(df, df.name);
    try {
      expect(view.name, df.name, 'the analysis view is titled with the table name');
      expect(w.showToolbox, false, 'BUG-18: the Toolbox is closed after the open');
      expect(w.showBrowse, true, 'BUG-18: the Browse tree is in front after the open');
      expect(w.autoShowToolbox, false, 'BUG-18: auto-show is off (user-directed, not restored)');
    } finally {
      view.close();
      w.showToolbox = saved.toolbox;
      w.showBrowse = saved.browse;
      w.autoShowToolbox = saved.auto;
    }
  }, {timeout: TIMEOUT});
});

/**
 * Add-to-study flow — persist the analysis currently open in the mapping pane
 * into a study as an `AnalysisEntry`.
 *
 * **File-backed source (user decision 2026-06-29, supersedes decision e's
 * pointer-only stance).** Adding an analysis WRITES the current table as a CSV
 * into the study's own folder and records a `platformPath` to that file. The
 * study therefore OWNS a copy of the data and re-opens in any session (the
 * resolver re-reads the file) — fixing the in-session-table fragility where a
 * `tableId` source could not resolve after a reload.
 *
 * Reuses the lazy-hydration primitives: `getCurrentRecipe()` for the recipe, the
 * just-computed `perProfile` for the summary projection (`buildSummaryCache`),
 * and `computeSourceSignature` for the drift baseline. The math is never re-run
 * here — it reuses what the pane already computed (thesis #1).
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  STANDALONE_ROUTE_SEGMENT, buildSummaryCache, uniqueLabel, passthroughCandidate,
  SourceReplacedError, ContentChangedError, getNcaRepository, computeSnapshotHash,
} from '../adapter';
import type {
  SourcePointer, AnalysisRecipe, PassthroughColumnSpec, SummaryCache, ObservationsCtx,
  MoveTarget, MoveResult, StudyHeader,
} from '../adapter';
import type {PerProfileResult} from '../controller';
import type {StudyBinding, WorkflowPersist} from './nca-mapping-pane';
import {UNGROUPED_SECTION} from './nca-studies-list';

export interface AddToStudyInput {
  /** The analysis source DataFrame — written into the study folder as a CSV. */
  readonly df: DG.DataFrame;
  /** The live, exhaustive recipe (`getCurrentRecipe()`). */
  readonly recipe: AnalysisRecipe;
  /** The just-computed per-profile results (for the summary cache). */
  readonly perProfile: ReadonlyArray<PerProfileResult>;
  /** The post-`appendGroupingColumns` results DataFrame — the per-(compound ×
   *  route × dose-norm) grouping substrate `buildSummaryCache` aggregates over
   *  (v2, study-unified-grid AD-1). Carries the grouping metadata columns. */
  readonly resultsDf: DG.DataFrame;
  /** The observations substrate for the v5 inline curve thumbnail (chemist-sar-grid
   *  S4). `null` when unavailable — the thumbnail is then absent (never faked). */
  readonly observations: ObservationsCtx | null;
  /** The v15 workflow state (confirmed recipe, keyed overrides, review ledgers, stage) — carried
   *  onto the new study entry so nothing the analyst decided is lost (D1). */
  readonly workflow: WorkflowPersist;
  /** When the view is bound to a STANDALONE analysis, adding it to a study MOVES it (D5): the
   *  persisted entry is moved (after the caller saved any unsaved edits) and the standalone copy
   *  is deleted. Absent → an unsaved analysis goes straight into the study (AC-S2.5.2). */
  readonly moveFrom?: StudyBinding;
  /** Profiles still to review (D14 soft gate). When > 0 the SAME dialog opens on a warning with
   *  Cancel / Add anyway; Add anyway reveals the inputs and the OK button (user-directed
   *  2026-09-29 — one dialog instead of a warning dialog followed by this one). */
  readonly pendingReview?: number;
  /** Pre-select this study in the dialog (set when the analysis was launched from
   *  a study's "+ add analysis" entry). The user can still switch targets. */
  readonly presetStudyId?: string;
  /** Fired after a SUCCESSFUL save with the created study/analysis identity, so the
   *  launching mapping-pane view can promote itself to study-backed and enable the
   *  in-place "Update analysis" widget without a reopen (analysis-update-in-place, F4/D7). */
  readonly onSaved?: (saved: SavedAnalysis) => void;
  /** Fired when the dialog is dismissed — on OK, Cancel, or the close (×) — AND on any early
   *  return that means the dialog never opened (e.g. the study list failed to load). The caller
   *  uses it to lift the single-instance gate on the launching "Add to study" icon so a second
   *  dialog cannot be stacked while one is open (nca-add-to-study-dialog-dedup). Because it also
   *  fires on the no-show path, the gate can never latch permanently. */
  readonly onClosed?: () => void;
}

/** A study analysis just written (Add to study, Move) — what the launching view rebinds to. */
export interface SavedAnalysis {
  readonly studyId: string;
  readonly studyLabel: string;
  readonly analysisId: string;
  readonly label: string;
  readonly source: SourcePointer;
  readonly sourceHash: string | null;
  readonly passthroughCols: readonly PassthroughColumnSpec[];
}

/** Live analysis state a SAVE writes (the pane already computed all of it — the math is never
 *  re-run here, thesis #1). */
export interface SaveAnalysisInput {
  /** The analysis source DataFrame — (re-)written as the container's data file (the OUTER df,
   *  so a derived-frame mode never writes a pooled/normalized frame back as the source). */
  readonly df: DG.DataFrame;
  /** `getCurrentRecipe()` — the live recipe (its `passthroughCols` is `[]`; the persisted
   *  dialog-time selection is re-injected from `binding.passthroughCols`). */
  readonly recipe: AnalysisRecipe;
  /** The last confirmed compute, or null when the analysis was never confirmed (a Setup-stage
   *  save — D1: Save works at any stage). `resultsDf` is the REPORTABLE frame (AD-6). */
  readonly computed: {
    readonly perProfile: ReadonlyArray<PerProfileResult>;
    readonly resultsDf: DG.DataFrame;
    readonly observations: ObservationsCtx | null;
  } | null;
  /** The v15 workflow state (confirmed recipe, keyed overrides, review ledgers, stage). */
  readonly workflow: WorkflowPersist;
}

/** Update-in-place input: a save that already has a container to write into. */
export interface UpdateAnalysisInput extends SaveAnalysisInput {
  /** The analysis identity to update in place (from `view.temp` — set on reopen or after the
   *  first save / Add-to-study). */
  readonly binding: StudyBinding;
}

/** The summary cache a save persists. Built from the last CONFIRMED compute — with the confirmed
 *  recipe when Setup was edited since (P6 keep-and-mark: an un-confirmed edit keeps the last
 *  confirmed results), and the live Report inputs (grouping). `null` when never confirmed. */
function cacheForSave(
  input: SaveAnalysisInput, label: string, passthroughCols: readonly PassthroughColumnSpec[],
): SummaryCache | null {
  if (input.computed === null) return null;
  const base = input.workflow.confirmed?.recipe ?? input.recipe;
  const recipeForCache: AnalysisRecipe = {...base, groupByCols: input.recipe.groupByCols, passthroughCols};
  return buildSummaryCache(input.computed.resultsDf, input.computed.perProfile, recipeForCache,
    label, input.computed.observations);
}

/** What {@link updateAnalysisInPlace} did. `rebound` is set when the view must rebind: the Save
 *  wrote a new snapshot (the data changed), or the work was saved as a NEW analysis beside the
 *  original (F11 / PR-1). */
export interface UpdateAnalysisResult {
  readonly ok: boolean;
  readonly rebound?: StudyBinding;
}

/**
 * In-place SAVE (analysis-update-in-place F1; nca-workflow-stages F-S2.4). Writes the analysis
 * row: the live recipe, the summary cache of the last confirmed compute (or none, if never
 * confirmed), the v15 workflow state — and a new snapshot only when the data changed (SD-19).
 * Keeps the SAME analysis (unlike Add-to-study). The caller clears its dirty flag only on success.
 *
 * **Never written over another session's work (decision 1, PR-1, F11).** The repository refuses
 * when the stored analysis no longer points at the snapshot this view opened (its data was
 * replaced — a pipeline re-run) or another session saved its content (recipe, confirmed recipe,
 * λz overrides, review decisions). Then, in a STUDY, the work is saved BESIDE as a new analysis on
 * the view's own data, with `origin: null`; a standalone container holds exactly one analysis, so
 * there the save is refused with an error and the view stays dirty.
 */
export async function updateAnalysisInPlace(input: UpdateAnalysisInput): Promise<UpdateAnalysisResult> {
  const {binding} = input;
  try {
    // Re-inject the persisted passthrough selection (getCurrentRecipe emits []) so the study
    // grid keeps its SAR/annotation columns after an update — same rule as the reopen onComputed.
    const recipe: AnalysisRecipe = {...input.recipe, passthroughCols: binding.passthroughCols};
    const summaryCache = cacheForSave(input, binding.label, binding.passthroughCols);
    let source: SourcePointer;
    try {
      source = await getNcaRepository().saveAnalysis(binding.analysisId, binding.source.ref,
        {df: input.df, recipe, summaryCache, workflow: input.workflow});
    } catch (e) {
      if (e instanceof SourceReplacedError || e instanceof ContentChangedError)
        return await saveBeside(input, e);
      throw e;
    }
    grok.shell.info('Analysis saved');
    return source.ref === binding.source.ref ? {ok: true} :
      {ok: true, rebound: {...binding, source, sourceHash: computeSnapshotHash(input.df)}};
  } catch (e) {
    grok.shell.error(`Could not save analysis: ${e instanceof Error ? e.message : String(e)}`);
    return {ok: false};
  }
}

/** The save-beside of F11 / PR-1: the view's work as a NEW analysis next to the one it could not
 *  update (study containers only — a standalone container holds one analysis). */
async function saveBeside(
  input: UpdateAnalysisInput, cause: SourceReplacedError | ContentChangedError,
): Promise<UpdateAnalysisResult> {
  const {binding} = input;
  const replaced = cause instanceof SourceReplacedError;
  if (binding.container === 'standalone') {
    grok.shell.error(replaced ?
      `Could not save '${binding.label}': its data was replaced while it was open. Nothing was saved.` :
      `'${binding.label}' was changed in another session. Nothing was saved — reopen it to see their changes.`);
    return {ok: false};
  }
  const recipe: AnalysisRecipe = {...input.recipe, passthroughCols: binding.passthroughCols};
  const saved = await getNcaRepository().appendAnalysis(binding.studyId, {
    label: binding.label, df: input.df, recipe,
    summaryCache: cacheForSave(input, binding.label, binding.passthroughCols),
    workflow: input.workflow, origin: null,
  });
  // The F11 wording is the approved design row ("Saving an analysis whose data was replaced").
  grok.shell.warning(replaced ?
    `The data of '${binding.label}' was replaced while it was open (for example by a pipeline run). ` +
    `Your work was saved as '${saved.label}', on the data you were working on.` :
    `'${binding.label}' was changed in another session. Your work was saved as '${saved.label}'.`);
  return {ok: true, rebound: {...binding, analysisId: saved.analysisId, label: saved.label, source: saved.source,
    sourceHash: computeSnapshotHash(input.df)}};
}

/**
 * The FIRST save of a never-saved analysis (D4). The caller has already asked for the name (the
 * `Save analysis` dialog, pre-filled with the table name). It creates a standalone container —
 * the Hub's Ungrouped section, one analysis — labelled `desiredLabel` made unique among the
 * caller's standalone analyses, and returns the new binding. Returns null on failure (error shown).
 */
export async function saveNewStandaloneAnalysis(
  input: SaveAnalysisInput, desiredLabel: string,
): Promise<StudyBinding | null> {
  try {
    const label0 = desiredLabel.trim() || 'Analysis';
    const saved = await getNcaRepository().createStandalone({
      label: label0, df: input.df, recipe: {...input.recipe, passthroughCols: []},
      summaryCache: cacheForSave(input, label0, []), workflow: input.workflow, origin: null,
    });
    grok.shell.info(`Saved to ${UNGROUPED_SECTION} as '${saved.label}'`);
    return {
      container: 'standalone', studyId: saved.studyId, analysisId: saved.analysisId, label: saved.label,
      source: saved.source, sourceHash: computeSnapshotHash(input.df), passthroughCols: [],
    };
  } catch (e) {
    grok.shell.error(`Could not save analysis: ${e instanceof Error ? e.message : String(e)}`);
    return null;
  }
}

/** Format a dose scalar for a label: integers bare, otherwise 4 significant figures. */
function formatDoseValue(d: number): string {
  return Number.isInteger(d) ? String(d) : String(Number(d.toPrecision(4)));
}

/**
 * A short, human default analysis label from the recipe: `NCA · <route> <dose> · by <groups>`
 * (user-directed 2026-10-05, study-add-multiple-analyses). Route and dose appear only when they are
 * single scalars — a per-row column varies, so it is omitted (e.g. `NCA · PO` when the dose is a
 * column); the `by <groups>` segment appears only when the analysis groups. The label separates
 * variants that differ in SCALAR route/dose or Group-by; variants differing only in concentration
 * units, λz edits or rule settings read the same and are told apart by the `uniqueLabel` suffix the
 * dialog adds on save (and the analyst can rename).
 */
export function defaultAnalysisName(recipe: AnalysisRecipe): string {
  const c = recipe.columns;
  const route = c.routeSource === 'value' ? (c.routeScalar ?? '').trim() : '';
  const dose = c.doseSource === 'value' && Number.isFinite(c.doseScalar) ?
    `${formatDoseValue(c.doseScalar)} ${recipe.units.dose}`.trim() : '';
  const routeDose = [route, dose].filter((s) => s.length > 0).join(' ');
  const groups = recipe.groupByCols.length > 0 ? `by ${recipe.groupByCols.join(', ')}` : '';
  return ['NCA', routeDose, groups].filter((s) => s.length > 0).join(' · ');
}

/**
 * Open the Add-to-study dialog: pick an existing study or create a new one, name
 * the analysis, then WRITE the table into the study folder and persist the entry
 * (source = the written file's path). Reuses the coordinator for load/save so
 * persistence/cache ownership stays in one place (AD-2).
 */
export async function openAddToStudyDialog(input: AddToStudyInput): Promise<void> {
  const repo = getNcaRepository();
  let studies: StudyHeader[] = [];
  let allStudyLabels: string[] = [];
  try {
    const all = (await repo.listStudies()).filter((s) => s.kind === 'study');
    allStudyLabels = all.map((s) => s.label);
    // Only a study the caller can edit can take an analysis (SD-13 — View-only sharing, decision 10).
    studies = all.filter((s) => s.access.canEdit);
  } catch (e) {
    grok.shell.error(`Could not list studies: ${e instanceof Error ? e.message : String(e)}`);
    // The dialog never opens on this path — lift the caller's single-instance gate so the
    // "Add to study" icon does not latch disabled forever (nca-add-to-study-dialog-dedup).
    input.onClosed?.();
    return;
  }

  const labelToId = new Map<string, string>();
  for (const s of studies) labelToId.set(s.label, s.id);
  const hasStudies = studies.length > 0;
  const presetLabel = input.presetStudyId !== undefined ?
    studies.find((s) => s.id === input.presetStudyId)?.label : undefined;

  // Destination = a New study or an Existing one, chosen with a horizontal radio (user-directed
  // 2026-10-05). When there are NO studies yet the radio is omitted entirely and only the
  // new-study-name field shows — there is nothing to pick. A preset launching study defaults to
  // Existing + selects it.
  const destChoice = ui.input.radio('Study', {
    items: ['New', 'Existing'], value: presetLabel !== undefined ? 'Existing' : 'New',
    tooltipText: 'Add the analysis to a new study, or to one you already have.',
  });
  destChoice.root.classList.add('nca-study-dest'); // CSS lays the two options in one row
  const existingStudy = ui.input.choice<string>('Existing study', {
    value: presetLabel ?? [...labelToId.keys()][0] ?? null, items: [...labelToId.keys()],
    tooltipText: 'The study to add this analysis to.',
  });
  // Default = the TABLE name (user-directed 2026-09-30); made unique against existing studies on OK.
  const newStudyName = ui.input.string('New study name', {
    value: input.df.name?.trim() || 'New study',
    tooltipText: 'Name for the new study.',
  });
  // Default = a short spec of the settings (NCA · route dose · by group) so several analyses of one
  // table read apart at a glance (user-directed 2026-10-05).
  const analysisName = ui.input.string('Analysis name', {
    value: defaultAnalysisName(input.recipe),
    tooltipText: 'Label for this analysis within the study.',
  });
  /** True when the OK action creates a new study (no studies at all, or New chosen). */
  const useNew = (): boolean => !hasStudies || destChoice.value === 'New';

  // v10 (study-grid-passthrough-columns): pick extra source columns to carry into the study
  // grid as passthrough (annotation) columns — SMILES → structures, chemical name, series, …
  // (the SAR-grid plumbing). Candidate gate is CAP-FREE (constant-within-subject + not-a-role +
  // no reserved-name collision, RG-D); structure-designated columns are listed FIRST (RG-G) via
  // the stock control's `available` order — no custom control. The
  // structure designation is captured HERE from each source column's live `semType` and persisted
  // on the recipe (LBC-2 — never re-derived downstream, where `buildPerSubjectColumn` strips it).
  const subjectCol = input.recipe.columns.subject !== null ?
    input.df.col(input.recipe.columns.subject) : null;
  const candidateNames = input.df.columns.names()
    .filter((nm) => passthroughCandidate(input.df.col(nm)!, input.recipe, subjectCol));
  // Structure-first ordering (RG-G, user-approved 2026-07-14): MOLECULE-semType columns lead,
  // then the rest in source order (V8 Array.sort is stable, so ties keep frame order).
  const isStruct = (nm: string): boolean => input.df.col(nm)!.semType === DG.SEMTYPE.MOLECULE;
  const orderedCandidates = [...candidateNames].sort(
    (a, b) => (isStruct(a) ? 0 : 1) - (isStruct(b) ? 0 : 1));
  // When the source table has NO annotatable column (every column is a mapped analysis role, a
  // grouping key, or a computed NCA parameter — e.g. theoph.csv: Subject/Time/conc/Dose/Wt), the
  // candidate list is empty: disable the input + explain WHY rather than leave an empty picker.
  const noCandidates = orderedCandidates.length === 0;
  // "Annotations" (user-directed 2026-10-05, not "Covariates" — these are descriptor columns for the
  // study grid: SMILES → structures, compound name, series, project, batch; they do NOT affect the
  // calculation, unlike a PK covariate). Optional: the field is labelled "(optional)".
  const colsInput = ui.input.columns('Annotations (optional)', {
    table: input.df, value: [], available: orderedCandidates,
    // See the BUG-11 note on the Group-by picker (`nca-mapping-pane.ts`): the default TRUE
    // draws the selected columns both in a top block and in the main list, so re-opening this
    // picker after a selection ticks each chosen name twice and OK commits the doubled set.
    // Reachable here (select SMILES → OK → re-open → OK), and the duplicate would DOUBLE every
    // passthrough subject count. The dedupe below absorbs it; this removes the generator.
    showSelectedColsOnTop: false,
    tooltipText: noCandidates ?
      'No annotatable columns — every column is a mapped role or NCA parameter.' :
      'Optional. Extra source columns — structures (SMILES), compound name, series, … — shown ' +
      'beside the NCA parameters in the study grid. They annotate rows; they do not change the ' +
      'calculation.',
  });
  colsInput.enabled = !noCandidates;

  // Show the new-study-name field vs the existing-study picker per the radio; hide the radio
  // entirely when there are no studies to choose (user-directed 2026-10-05).
  const syncDest = (): void => {
    destChoice.root.style.display = hasStudies ? '' : 'none';
    newStudyName.root.style.display = useNew() ? '' : 'none';
    existingStudy.root.style.display = (hasStudies && !useNew()) ? '' : 'none';
  };

  // D14 soft gate folded into this dialog: with profiles still to review, it opens on the warning
  // alone (Cancel / Add anyway); Add anyway hides the warning, reveals the inputs and swaps itself
  // for OK.
  const pending = input.pendingReview ?? 0;
  const warning = ui.divText(
    `${pending} ${pending === 1 ? 'profile is' : 'profiles are'} not reviewed. The study grid will ` +
    `show "${pending} not reviewed" for this analysis until they are.`,
    {style: {maxWidth: '420px', whiteSpace: 'normal', marginBottom: '8px'}});
  const inputsHost = ui.divV(
    [destChoice.root, existingStudy.root, newStudyName.root, analysisName.root, colsInput.root]);
  const dialog = ui.dialog('Add analysis to study')
    .add(warning)
    .add(inputsHost)
    // No Enter-to-OK while the inputs are still hidden behind the warning.
    .onOK(() => {
      void (async () => {
        const rawLabel = (analysisName.value ?? '').trim() || 'Analysis';
        const isNew = useNew();

        try {
          // v10: the dialog-time passthrough selection → PassthroughColumnSpec[], capturing
          // each column's live semType as the persisted structure designation (LBC-2 — the ONE
          // point isStructure is resolved; downstream reads the recipe, never re-detects). The
          // recipe passed to buildSummaryCache is OVERRIDDEN with this selection (input.recipe's
          // passthroughCols is the getCurrentRecipe [] default — PP-1); the entry stores it too.
          // Deduped by NAME (first-seen wins), for the same reason the Group-by picker is
          // (BUG-11): the stock columns dialog ticks its rows BY NAME over a `value ++
          // tableColumns` list, so a name can be committed twice. Unlike Group-by this picker is
          // never rebound to a new frame, so the duplicate is not known to be reachable here —
          // but the consequence if it ever is would be a SILENT WRONG NUMBER, not a doubled
          // label: `accumulateByTier` keys its counts by `col.name`, so two same-named candidates
          // increment the SAME counter and every subject count in the passthrough disclosure
          // DOUBLES (`F (30); M (30)` for 15 + 15). A wrong count reads as real data — exactly
          // what must never be silent (thesis #4). One `Set` at the boundary where the widget's
          // value becomes persisted state costs nothing and removes the class.
          const seenPassthrough = new Set<string>();
          const passthroughCols: PassthroughColumnSpec[] = (colsInput.value ?? [])
            .filter((c: DG.Column) => {
              if (seenPassthrough.has(c.name)) return false;
              seenPassthrough.add(c.name);
              return true;
            })
            .map((c: DG.Column) => ({name: c.name, isStructure: c.semType === DG.SEMTYPE.MOLECULE}));
          const recipeForCache: AnalysisRecipe = {...input.recipe, passthroughCols};

          // A new study's label is unique among the caller's visible studies; `analyses` is
          // reserved (the standalone route segment, F-S2.6).
          const newStudyLabel = (): string => uniqueLabel((newStudyName.value ?? '').trim() || 'New study',
            new Set([...allStudyLabels, STANDALONE_ROUTE_SEGMENT]));

          if (input.moveFrom !== undefined) {
            // D5 MOVE: one implementation shared with the hub's "Move to study…" (F-S2.5). The
            // live results rebuild the cache so the dialog-time passthrough columns appear.
            const target: MoveTarget = isNew ? {kind: 'new', label: newStudyLabel()} :
              {kind: 'existing', studyId: labelToId.get(existingStudy.value as string)!};
            const moved = await moveStandaloneIntoStudy(input.moveFrom.studyId, target, {
              label: rawLabel, passthroughCols,
              summaryCache: buildSummaryCache(
                input.resultsDf, input.perProfile, recipeForCache, rawLabel, input.observations),
            });
            if (moved !== null) input.onSaved?.({...moved, sourceHash: input.moveFrom.sourceHash});
            return;
          }

          const study = isNew ? await repo.createStudy(newStudyLabel()) :
            studies.find((s) => s.id === labelToId.get(existingStudy.value as string))!;
          // The analysis label is made unique WITHIN the study by the repository (label-keyed
          // routing needs it — dup labels → dup tree leaves + first-match URL resolution).
          const saved = await repo.appendAnalysis(study.id, {
            label: rawLabel, df: input.df, recipe: recipeForCache,
            summaryCache: buildSummaryCache(
              input.resultsDf, input.perProfile, recipeForCache, rawLabel, input.observations),
            workflow: input.workflow, origin: null,
          });
          grok.shell.info(`Saved analysis '${saved.label}' to study '${study.label}'.`);
          // Promote the launching view to study-backed so "Update analysis" targets THIS
          // analysis without a reopen (F4/D7).
          input.onSaved?.({studyId: study.id, studyLabel: study.label, analysisId: saved.analysisId,
            label: saved.label, source: saved.source, sourceHash: computeSnapshotHash(input.df), passthroughCols});
        } catch (e) {
          grok.shell.error(`Could not save: ${e instanceof Error ? e.message : String(e)}`);
        }
      })();
    }, {closeOnEnter: pending === 0});
  // The OK button reads the ACTION: "Create study & add" for a new study, "Add to study" for an
  // existing one (user-directed 2026-10-05). Kept in sync with the radio.
  const okBtn = dialog.getButton('OK');
  const syncOkText = (): void => {
    if (okBtn != null) okBtn.textContent = useNew() ? 'Create study & add' : 'Add to study';
  };
  destChoice.onChanged.subscribe(() => {syncDest(); syncOkText();});
  syncDest();
  syncOkText();
  if (pending > 0) {
    inputsHost.style.display = 'none';
    okBtn.style.display = 'none';
    dialog.addButton('Add anyway', () => {
      warning.style.display = 'none'; // acknowledged — the inputs replace it (user-directed 2026-09-29)
      inputsHost.style.display = '';
      okBtn.style.display = '';
      syncOkText(); // the OK button was hidden at build — set its text when it appears
      const addAnyway = dialog.getButton('Add anyway');
      if (addAnyway !== null && addAnyway !== undefined) addAnyway.style.display = 'none';
    });
  } else
    warning.style.display = 'none';
  // Lift the caller's single-instance gate whenever the dialog is dismissed — OK, Cancel, or the
  // close (×) all raise `onClose` (nca-add-to-study-dialog-dedup). Subscribed BEFORE `show()` so
  // no close can be missed.
  dialog.onClose.subscribe(() => input.onClosed?.());
  dialog.show();
}

/**
 * MOVE a standalone analysis into a study (nca-workflow-stages F-S2.5 / D5) — the ONE
 * implementation behind the in-view `⊞ Add to study` of a standalone-bound view and the hub /
 * tree `Move to study…`.
 *
 * One transaction (AC-3.5): the analysis is re-parented with its snapshot, decisions, overrides
 * and history (AC-S2.5.1), its label made unique within the study, `origin` recorded so the old
 * `analyses/<label>` link redirects; the emptied container is marked moved and deleted. Nothing is
 * copied, so nothing can be left half-moved. Returns the new study binding, or null on failure.
 */
export async function moveStandaloneIntoStudy(
  standaloneId: string, target: MoveTarget,
  opts: {
    label?: string;
    passthroughCols?: readonly PassthroughColumnSpec[];
    /** A freshly built cache (the in-view path has live results); else the stored one is kept. */
    summaryCache?: SummaryCache;
  } = {},
): Promise<MoveResult | null> {
  try {
    const moved = await getNcaRepository().moveToStudy(standaloneId, target, opts);
    grok.shell.info(`Moved analysis '${moved.label}' to study '${moved.studyLabel}'.`);
    return moved;
  } catch (e) {
    grok.shell.error(`Could not move the analysis: ${e instanceof Error ? e.message : String(e)}`);
    return null;
  }
}

/**
 * Lazy analysis hydration (AD-4 / Slice D) — the single path that opens ONE analysis (by click
 * from the overview, the Hub, the tree, OR by cold deep-link) and computes it on demand. Every
 * entry point goes through {@link openAnalysis} → `mountAnalysisIntoView` → `triggerCompute`, so
 * byte-identical state is a shared-code guarantee, not coincidence.
 *
 * The open reads the analysis row and its d42 snapshot (2 requests); the math is recomputed
 * from the recipe, never read back from the cache as truth (thesis #1). After a clean compute the
 * summary projection is written back through the repository (SD-18).
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {mountAnalysisIntoView} from './nca-launch';
import {keepBrowseInFront} from './shell-panels';
import {normalizeNumericNaMarkers} from '../adapter/na-normalize';
import {updateNcaStudyPath, updateNcaAnalysisPath} from './nca-studies-router';
import {buildSummaryCache, recipeOf, getNcaRepository, computeSnapshotHash, getAuditLog} from '../adapter';
import type {AnalysisRecipe, SourcePointer, PassthroughColumnSpec, LoadedAnalysis} from '../adapter';
import type {StudyBinding, ComputeCallbackCtx} from './nca-mapping-pane';
// The mapping pane imports back (`makeStudyBackedOnComputed` for its promote-after-save) — but
// only at CALL time inside a handler, never at module init, so the runtime cycle is safe
// (guardrails: init-time cross-module reads are the fatal case, not lazy ones).
import {
  STUDY_BINDING_TAG, ANALYSIS_REPLAYING_TAG, setPendingWorkflowState, setConfirmed,
  currentRecipeFingerprint, refreshStageUi, setPipelineReasons, markAnalysisDirty,
} from './nca-mapping-pane';
import {pipelineReasonsViewOf} from './results-view';

/**
 * Build the study-backed `onComputed` callback shared by the reopen path and the mapping pane's
 * promote-after-save (analysis-update-in-place F4/D7 — "the two entry points converge"). On every
 * clean compute it refreshes the analysis's summary cache and the workflow state that travels with
 * it (SD-18). It never writes the snapshot — that is SAVE's job.
 *
 * Skipped when the view is dirty (AD-8: a cache of an unsaved recipe would describe work that was
 * not saved), and when the frame the compute ran on is not the stored snapshot (PR-2 / PR-13: the
 * analyst edited the table) — then `onUnsavedData` marks the view dirty so SAVE lights up.
 * `passthroughCols` is re-injected because `getCurrentRecipe()` emits `[]` (LBC-2).
 * `openedHash` is the hash of the stored snapshot AS OPENED: it differs from `storedHash` only when
 * the open's NA-marker normalization retyped a column of an otherwise identical snapshot (a pipeline
 * writes its frame un-normalized) — that frame is still the stored data, not an analyst edit.
 */
export function makeStudyBackedOnComputed(
  analysisId: string, label: string, passthroughCols: readonly PassthroughColumnSpec[], sourceDf: DG.DataFrame,
  source: SourcePointer, storedHash: string | null, onUnsavedData: () => void, openedHash: string | null = storedHash,
): (ctx: ComputeCallbackCtx) => void {
  return (ctx) => {
    if (ctx.dirty) return;
    const frameHash = computeSnapshotHash(sourceDf);
    if (storedHash !== null && frameHash !== openedHash) {
      onUnsavedData();
      return;
    }
    const recipeForCache: AnalysisRecipe = {...ctx.recipe, passthroughCols};
    const summaryCache = buildSummaryCache(
      ctx.resultsDf, ctx.perProfile, recipeForCache, label, ctx.observations);
    // The clean view's workflow state travels WITH the cache it describes: the reopen compute may
    // have re-stamped the confirmed fingerprint, dropped a stale override, or auto-reset a decision
    // (D13) — each a deterministic consequence of the saved state on this source / engine.
    void getNcaRepository().writeSummary(analysisId, source.ref,
      {summaryCache, workflow: ctx.workflow, frameHash: storedHash ?? frameHash})
      .catch((e) => grok.shell.error(`Could not save updated summary: ${e instanceof Error ? e.message : String(e)}`));
  };
}

/** A minimal message view (not-found / not-authorized / error). Shared by the deep-link path and
 *  the study overview. */
export function messageView(title: string, message: string): DG.View {
  const root = ui.divV([
    ui.h1(title),
    ui.divText(message, {style: {maxWidth: '480px', whiteSpace: 'normal'}}),
  ], {style: {padding: '24px'}});
  const view = DG.View.fromRoot(root);
  view.name = title;
  return view;
}

/** The tooltip of every write affordance a View-only user sees disabled (SD-13). */
export const VIEW_ONLY_TIP = 'View only — ask the study owner for Edit access. Changes here are not saved.';

/**
 * Open ONE analysis: read its row and snapshot, mount the mapping pane with the recipe, compute
 * once, and write the refreshed summary back (a real open only — a single-click PREVIEW peek
 * persists nothing). A missing / unreadable analysis resolves to a message view.
 */
export async function openAnalysis(analysisId: string, mount: 'view' | 'preview' = 'view'): Promise<DG.ViewBase> {
  keepBrowseInFront(); // before the table loads, not after (user-directed 2026-10-12)
  let loaded: LoadedAnalysis;
  try {
    loaded = await getNcaRepository().loadAnalysis(analysisId);
  } catch (e) {
    const msg = e instanceof Error ? e.message : String(e);
    grok.shell.error(msg);
    return messageView('Analysis unavailable',
      `${msg} It may have been deleted, or you may not have access to its study.`);
  }
  const {study, entry, df} = loaded;
  // A retyped column of the untouched snapshot is still the stored data (see makeStudyBackedOnComputed).
  const storedHash = entry.sourceSignature.hash;
  const hasStrings = df.columns.toList().some((c) => c.type === DG.COLUMN_TYPE.STRING);
  const rawHash = hasStrings && storedHash !== null ? computeSnapshotHash(df) : null;
  const openedHash = normalizeNumericNaMarkers(df).length > 0 && rawHash === storedHash ?
    computeSnapshotHash(df) : storedHash;

  // `mount: 'preview'` (browse-tree single-click peek, nca-tree-node-interactions D1) opens the
  // table in the platform preview slot with `addTablePreview` — the platform's own TABLE preview.
  // The slot auto-replaces any prior preview; a double-click / "Run" (mount: 'view') opens it
  // normally. Both keep the toolbox closed so the side panel stays on the Browse tree.
  keepBrowseInFront();
  let view: DG.TableView;
  if (mount === 'preview') {
    // Drop the current preview slot occupant before showing ours (DiffStudio releaseSharedPreview:
    // grok.shell.preview = null, NEVER view.close() — closing the app-landing view closes the app).
    grok.shell.preview = null;
    view = grok.shell.addTablePreview(df);
  } else view = grok.shell.addTableView(df);
  view.name = entry.label;
  const standalone = study.kind === 'standalone';
  // Label-keyed URL (router resolves back): `analyses/<label>` for a standalone analysis (F-S2.6),
  // `<study>/<analysis>` for a study one.
  if (standalone) updateNcaAnalysisPath(view, entry.label);
  else updateNcaStudyPath(view, study.label, entry.label);

  // A PREVIEW peek is not an open: no audit entry (ALCOA+), no writeback. A View-only user's open
  // computes but persists nothing (SD-13: their write would only end in a 403).
  const persists = mount !== 'preview' && study.access.canEdit;
  if (mount !== 'preview') {
    getAuditLog().append({timestamp: Date.now(), subjectKey: '', action: 'project-analysis-open', ruleHash: '',
      status: 'ok', studyId: study.id, analysisId});
  }

  // The REPLAY guard goes up BEFORE the mount: the mount replays the saved recipe into the inputs,
  // and some of those onChanged fire synchronously inside it (a pipeline recipe did — every open of
  // a pipeline-saved analysis was marked dirty, so its write-back never ran). Cleared once the
  // deferred compute settles.
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const paneTemp = view.temp as Record<string, any>;
  paneTemp[ANALYSIS_REPLAYING_TAG] = true;
  const handle = mountAnalysisIntoView(view, df, {
    recipe: recipeOf(entry),
    // R1 SY-F5: seed the rules closure from the recipe's NcaRules.
    rules: entry.ncaRules,
    ...(persists ? {
      onComputed: makeStudyBackedOnComputed(analysisId, entry.label, entry.passthroughCols, df, entry.source,
        storedHash, () => markAnalysisDirty(view), openedHash),
    } : {}),
  });
  // Bind the view to its saved analysis (the SAVE target).
  paneTemp[STUDY_BINDING_TAG] = {
    container: standalone ? 'standalone' : 'study',
    studyId: study.id, analysisId, label: entry.label,
    source: entry.source, sourceHash: entry.sourceSignature.hash, passthroughCols: entry.passthroughCols,
    containerLabel: study.label, readOnly: !study.access.canEdit,
  } satisfies StudyBinding;
  // v15 workflow state: the keyed λz overrides + review decisions ride into the controller before
  // the first compute; the confirmed recipe restores the stage.
  setPendingWorkflowState(view, {lzOverrides: entry.lzOverrides ?? [], reviewDecisions: entry.reviewDecisions ?? []});
  setConfirmed(view, entry.confirmed ?? null);
  // A pipeline-saved entry: the run's own reasons ride into the view (the `Pipeline` column and the
  // Review-pane line — nca-pipeline-save-rework F10); null for every other entry.
  setPipelineReasons(view, pipelineReasonsViewOf(entry.pipelineFlags, entry.origin));
  // Defer the compute one tick so every replayed `input.value` setter's onChanged has flushed
  // before `getCurrentRecipe()` / compute read them. `silent: true` — a passive replay-open.
  setTimeout(() => {
    // F-S3.8 (D19): land where the analysis stopped. A confirmed entry whose recipe still matches
    // the replayed Setup recomputes silently; a never-confirmed one (or one saved with an
    // un-confirmed Setup edit, P6) mounts at Setup with NO compute.
    const c = entry.confirmed ?? null;
    const fp = currentRecipeFingerprint(view);
    if (c !== null && (c.fingerprint === '' || c.fingerprint === fp))
      handle.triggerCompute({silent: true});
    else
      refreshStageUi(view);
    paneTemp[ANALYSIS_REPLAYING_TAG] = false; // replay settled — real edits now mark dirty
  }, 50);
  return view;
}

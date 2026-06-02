/*
 * Consensus Pharmacophore — TS orchestrator.
 *
 * Phase 6 shape - ALL five Python stages are stubbed. The stubs return DataFrames
 * matching Appendix B schemas so the UI works end-to-end (input panel -> Build ->
 * progress indicator -> docked Mol* shows a 7-feature consensus). Each stub is
 * flagged STUB-PHASE-6 and is replaced one at a time in Phases 1/2a/3/3.5/2b/4/5a.
 *
 * Stage 5b (renderer.ts) is the REAL implementation from the start — pure string
 * formatting, no Python needed.
 *
 * Blueprint reference - section 0, section 5 Phase 6, Appendix B.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FAMILY_CODES, FAMILY_MAP, resolveFamily} from './family-map';
import {FN} from './function-names';
import {renderConsensusPharmacophore} from './renderer';
import {pickClusterViaDialog} from './cluster-picker-dialog';
import {_package} from './package';
import {DEFAULT_OPTIONS, PipelineOptions, PocketMethod, PocketRep} from './orchestrator-types';
import {enrichPdbList, DroppedPdb} from './rcsb-client';
import {applyPdbTransform, concatPdbStructures, fetchPdbBlock, IDENTITY_4X4,
  ligandFeaturesToOverlayBlock, pocketAtomsToOverlayBlock, stripCofactorsFromPdb} from './pdb-utils';
import {fetchUniProtById, getPdbIdsForAccession, isUniProtAccession,
  searchUniProt, UniProtHit} from './uniprot-client';
import {isLocalId, localFingerprintFor} from './local-pdb-store';

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

export type {PipelineOptions, PocketMethod, PocketRep} from './orchestrator-types';

interface PipelineResult {
  consensusModel: DG.DataFrame;
  pdbBlock: string;
}

// ---------------------------------------------------------------------------
// App
// ---------------------------------------------------------------------------

/** Cached intermediate stage results keyed by `previewCacheKey()`. */
interface PreviewCache {
  fingerprint: string;
  pdbQc?: DG.DataFrame;
  rawBlocks?: Map<string, string>;     // pdb_id -> raw PDB text (TS-side fetch)
  aligned?: DG.DataFrame;              // Stage 2a output
  pocketAtoms?: DG.DataFrame;          // Stage 3 output
  alignedV2?: DG.DataFrame;            // Stage 2b output (refined transforms)
  ligandFeatures?: DG.DataFrame;       // Stage 4 output
  consensusModel?: DG.DataFrame;       // Stage 5a output
  // Signature of the consensus-only inputs (k-means knobs + the user's PDB /
  // interaction exclusions) that produced `consensusModel`. These are NOT in
  // the global `fingerprint` (changing them must not re-run Stages 1-4), so we
  // self-invalidate the consensus when this signature changes — see
  // `ensureConsensusModel`.
  consensusSignature?: string;
}

export class ConsensusPharmacophoreApp {
  tableView?: DG.TableView;
  private viewer?: DG.Viewer;
  private viewerPlaceholder?: HTMLElement;
  private viewerPlaceholderNode?: DG.DockNode;
  // Cached dock node for the Mol* viewer — used by refreshDockedViewer to
  // tear down the previous viewer reliably before creating a new one.
  private viewerDockNode?: DG.DockNode;
  private aborted = false;
  // Stage badges rendered in the top header — flipped from grey to green as stages complete.
  private stageBadges = new Map<string, HTMLElement>();
  // Per-input cache so Fetch -> Align -> Pocket clicks don't repeat stages.
  private previewCache: PreviewCache | null = null;
  // Per-PDB interaction-detail panel (Option A). Created lazily after Stage 4;
  // updated on summary-grid current-row changes. Closed when other previews run.
  private pdbDetailNode?: DG.DockNode;
  private pdbDetailHost?: HTMLElement;
  private pdbDetailSub?: {unsubscribe: () => void};
  // Caption element shown above the Mol* viewer; rewritten by each preview to
  // explain what the user is looking at (e.g. "Chain F = ProLIF interaction
  // sites; CPK colors = pharmacophore family"). Always visible after init().
  private molstarCaptionEl?: HTMLElement;

  // --- wizard integration -------------------------------------------------
  // When the WizardShell instantiates this app, it injects DOM mount points
  // for the Mol* viewer and the per-PDB detail panel via the attach* methods
  // below. With these set, refreshDockedViewer mounts the viewer directly
  // into the wizard's centerl content slot (no dock manager), and the detail
  // panel renders straight into the wizard's right content slot. When the
  // mount points are unset (legacy path through init()), refreshDockedViewer
  // falls back to the old dock-manager behaviour.
  private _externalMolstarHost?: HTMLElement;
  private _externalDetailHost?: HTMLElement;
  /** Hooks the wizard installs so the (orchestrator-rendered) Step 4 detail
   *  panel can read/toggle the per-interaction consensus selection that lives
   *  in the wizard's PipelineOptions.consensusExcludedInteractions. */
  private _consensusSelectionHooks?: {
    isInteractionExcluded: (key: string) => boolean;
    toggleInteraction: (key: string) => void;
    /** Whether a whole PDB is excluded via the Step 4 "Use" checkbox. Drives
     *  both the consensus filter AND the Mol* per-PDB visibility (un-ticking a
     *  PDB hides it in 3D — see setPdbHidden / applyPdbExclusionsToViewer). */
    isPdbExcluded: (pdbId: string) => boolean;
    /** Whether a PDB is un-ticked in Step 1's Accepted PDBs table. Stronger
     *  than isPdbExcluded — these PDBs are dropped before Stage 2, so they
     *  never enter alignment. For Mol* visibility, both predicates trigger
     *  a hide; the input exclusion is the only one that also gates the
     *  pipeline (see filterPdbQcByInputExclusion). */
    isInputPdbExcluded: (pdbId: string) => boolean;
  };
  /** consensusFeatureKey of the interaction the user last clicked in the Step 4
   *  detail panel — its sphere is marked in Mol* and its row is highlighted.
   *  Cleared on every viewer rebuild (the marker lives in the viewer). */
  private _highlightedInteractionKey?: string;
  // Cached Stage 4 summary DataFrame so the wizard's Step 4 panel can read it.
  private _summaryDf?: DG.DataFrame;
  /** Per-input PDB drop record from Stage 1. Updated whenever enrichPdbList
   *  is called. Persists across previewCache rebuilds for the same input
   *  set so the wizard can show the "Dropped" section. */
  private _droppedPdbs: DroppedPdb[] = [];
  /**
   * Map of PDB id (uppercase) -> Mol* structure ref. Reserved for a future
   * native-visibility-toggle implementation. Currently unused — see
   * `_cachedAlignedBlocks` for the active approach.
   */
  private _alignedStructureRefs: Map<string, string> = new Map();
  /**
   * Map of PDB id (uppercase) -> Mol* structure ref for the chain-F
   * interaction overlay. Populated by Step 4's per-PDB overlay rendering;
   * empty in stages that don't use the overlay (Step 2, 3, 5). Toggled
   * alongside the protein structures by `isolatePdb` so isolating one PDB
   * also hides every other PDB's interaction spheres.
   */
  // pdbId → list of overlay structure refs. Step 4 mounts ONE structure per
  // (PDB × family) so each family can be uniform-colored with its legend hex;
  // Step 3's pocket overlay mounts a single structure per PDB ([oneRef]). The
  // isolation logic hides/shows every ref under a PDB in lockstep.
  private _alignedInteractionRefs: Map<string, string[]> = new Map();
  /** consensusFeatureKey → that interaction's Mol* structure ref, and the
   *  reverse. Step 4 mounts ONE structure per interaction (not per family) so
   *  an individual contact can be hidden when the user un-ticks it. Populated
   *  by addInteractionOverlayPerPdb; consumed by setInteractionHidden and by
   *  isolatePdb (to keep un-ticked interactions hidden when showing a PDB). */
  private _interactionRefByKey: Map<string, string> = new Map();
  private _interactionKeyByRef: Map<string, string> = new Map();
  /**
   * Per-PDB transformed blocks from the last Step 2 alignment, in the
   * same order as `_cachedAlignedLabels`. Used by `isolatePdb` to rebuild
   * the Mol* viewer with just one PDB (or all of them) without a Python
   * round-trip. Cleared whenever the alignment cache becomes stale.
   */
  private _cachedAlignedBlocks: string[] = [];
  private _cachedAlignedLabels: string[] = [];

  /** Inject the wizard's Mol* content panel as the Mol* viewer mount point.
   *  Must be called BEFORE any preview* method runs. */
  attachMolstarHost(el: HTMLElement): void {
    this._externalMolstarHost = el;
  }

  /** Inject the wizard's PDB-detail content panel. The orchestrator
   *  will render Stage 4 detail content directly here. */
  attachDetailHost(el: HTMLElement): void {
    this._externalDetailHost = el;
    this.pdbDetailHost = el;
  }

  /** Wire the per-interaction consensus-selection checkboxes (rendered in the
   *  Step 4 detail panel by `renderPdbDetail`) to the wizard's option state. */
  attachConsensusSelectionHooks(hooks: {
    isInteractionExcluded: (key: string) => boolean;
    toggleInteraction: (key: string) => void;
    isPdbExcluded: (pdbId: string) => boolean;
    isInputPdbExcluded: (pdbId: string) => boolean;
  }): void {
    this._consensusSelectionHooks = hooks;
  }

  /** How many interactions / PDBs the current options' exclusions leave for the
   *  consensus. Drives the Step 4 + Step 5 "consensus will use N / M" counts.
   *  Counts non-diagnostic feature rows after `filterFeaturesForConsensus`. */
  consensusUsedCount(opts: PipelineOptions): {interactions: number; pdbs: number} {
    const feats = this.previewCache?.ligandFeatures;
    if (!feats) return {interactions: 0, pdbs: 0};
    const filtered = filterFeaturesForConsensus(feats, opts);
    const skipCol = filtered.col('skip_reason');
    const pdbCol = filtered.col('pdb_id');
    const pdbs = new Set<string>();
    let n = 0;
    for (let i = 0; i < filtered.rowCount; i++) {
      if (skipCol && String(skipCol.get(i) ?? '').trim() !== '') continue;
      n += 1;
      if (pdbCol) pdbs.add(String(pdbCol.get(i) ?? '').toUpperCase());
    }
    return {interactions: n, pdbs: pdbs.size};
  }

  /** Drop the per-fingerprint cache. Called by the wizard when the user
   *  edits the PDB list or options (stale-cascade). */
  invalidatePreviewCache(): void {
    this.previewCache = null;
  }

  /** Full reset: drop the cache, close the BSV viewer instance, drop the
   *  cached summary DataFrame, unsubscribe from any per-PDB detail panel
   *  updates, and clear the detail panel content. Called by the wizard's
   *  Reset ribbon button. The wizard separately clears its own DOM regions
   *  (contentCenter, contentRight) and re-shows the placeholder. */
  resetWizardState(): void {
    this.previewCache = null;
    this._summaryDf = undefined;
    this._droppedPdbs = [];
    this._alignedStructureRefs.clear();
    this.clearOverlayRefs();
    this._cachedAlignedBlocks = [];
    this._cachedAlignedLabels = [];
    // Close the BSV viewer instance — without this, the Mol* plugin keeps
    // its WebGL context and old structures live in memory.
    if (this.viewer) {
      try { (this.viewer as any).close?.(); } catch (_e) { /* */ }
      this.viewer = undefined;
    }
    // Unsubscribe the per-PDB detail panel from its current summary df.
    if (this.pdbDetailSub) {
      try { this.pdbDetailSub.unsubscribe(); } catch (_e) { /* */ }
      this.pdbDetailSub = undefined;
    }
    // Clear the detail panel contents (the wizard host element is kept; we
    // just remove what was rendered into it).
    if (this.pdbDetailHost) ui.empty(this.pdbDetailHost);
  }

  /** Close the BSV viewer instance and empty the mounted host WITHOUT
   *  dropping the preview cache. Used by the wizard when navigating to a
   *  step that has no 3D content (Step 1 / Fetch) so a leftover render from
   *  a later step doesn't linger. Lighter than `resetWizardState`, which
   *  also wipes the cache and detail panel. The wizard re-adds its own
   *  placeholder element to the host after calling this. */
  clearViewer(): void {
    if (this.viewer) {
      try { (this.viewer as any).close?.(); } catch (_e) { /* */ }
      this.viewer = undefined;
    }
    this.viewerDockNode = undefined;
    this._alignedStructureRefs.clear();
    this.clearOverlayRefs();
    if (this._externalMolstarHost) ui.empty(this._externalMolstarHost);
  }

  /** Read-only view of the current preview cache (DataFrames per stage).
   *  Used by the wizard's per-step panels to render result tables without
   *  going through grok.shell.tableView. Returns null when no fingerprint
   *  has been computed yet. */
  previewCacheView(): PreviewCache | null {
    return this.previewCache;
  }

  /** Latest Stage 4 summary DataFrame (one row per PDB · ligand). Returns
   *  null until previewFeatures has run at least once for the current
   *  fingerprint. */
  previewSummaryDataFrame(): DG.DataFrame | null {
    return this._summaryDf ?? null;
  }

  /** Per-PDB drop records from the last Stage 1 run (resolution / method /
   *  ligand filters). Empty until Stage 1 has actually executed. */
  previewDroppedPdbs(): DroppedPdb[] {
    return this._droppedPdbs;
  }

  /** Assign the active DataFrame to the legacy TableView, if one exists.
   *  In the wizard path there is no TableView — the per-step panels in
   *  wizard-steps.ts pull their data from the preview cache directly via
   *  previewCacheView() and previewSummaryDataFrame(). */
  private setActiveDataFrame(df: DG.DataFrame): void {
    if (this.tableView) this.tableView.dataFrame = df;
  }


  async init(pdbIds?: string[]): Promise<DG.TableView> {
    // Start with an empty placeholder DataFrame; the docked viewer + dock layout
    // attach to the VIEW, not the DataFrame, so we can swap dataFrame between stages.
    const placeholder = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('pdb_id', []),
    ]);
    placeholder.name = 'Consensus Pharmacophore';

    this.tableView = grok.shell.addTableView(placeholder);
    this.tableView.name = 'Consensus Pharmacophore';

    // Mol* viewer is created lazily on the first Build — passing an empty `pdb`
    // option to `df.plot.fromType('Biostructure', ...)` is unverified at the time
    // of writing. We dock a placeholder div instead.

    // Stage-progress header has been retired — the same status is already
    // conveyed three other ways: (a) preview buttons gain a green ✓ class
    // after each stage finishes, (b) the task-bar progress indicator at the
    // bottom-left shows the active stage label, (c) success and error toasts
    // appear on completion. Two TOP docks were also crowding the legend on
    // typical (≤900px) browser windows. The buildStageHeader / setStageStatus
    // methods are retained for now — the `this.stageBadges` map is just empty,
    // which makes setStageStatus a no-op (its `if (!el) return;` guard
    // already handles that).

    // Family legend + Mol* caption strip — single top row, always visible so
    // users have a persistent reference for the 7 family CPK colors and a
    // one-line hint of what the right pane is currently showing.
    const legendBar = this.buildLegendBar();
    this.tableView.dockManager.dock(legendBar, DG.DOCK_TYPE.TOP, null, 'Legend & view', 0.09);

    // Input panel (left)
    const panel = this.buildInputPanel(pdbIds);
    this.tableView.dockManager.dock(panel, DG.DOCK_TYPE.LEFT, null, 'Input', 0.25);

    // Viewer placeholder (right) — replaced by the Mol* viewer after first Build.
    this.viewerPlaceholder = ui.divText(
      'Run the pipeline to render the consensus pharmacophore.',
      'cp-viewer-placeholder');
    this.viewerPlaceholderNode = this.tableView.dockManager.dock(
      this.viewerPlaceholder, DG.DOCK_TYPE.RIGHT, null, 'Mol*', 0.55);

    // Pre-dock the PDB-detail panel at root-level RIGHT (null neighbor) so
    // the dock manager splits the right area into Mol* on the left and the
    // detail panel on the far right. Reserving this space at init time
    // guarantees Mol*'s width is identical across Fetch / Align / Pocket /
    // Features / Consensus — previously Mol* shrank at Stage 4 when the
    // detail panel was docked late.
    //
    // Earlier attempts:
    //  - DG.DOCK_TYPE.RIGHT relative to viewerPlaceholderNode → detail was
    //    rendered as a TOP horizontal strip above Mol* (Datagrok's nested
    //    splits interpret RIGHT in the container's secondary axis).
    //  - Skipping init-time dock → Mol* width changed at Stage 4.
    this.pdbDetailHost = ui.div([], 'cp-pdb-detail-panel');
    this.pdbDetailHost.append(ui.divText(
      'Per-PDB interaction detail — click 4. Features to populate.',
      'cp-pdb-detail-empty'));
    this.pdbDetailNode = this.tableView.dockManager.dock(
      this.pdbDetailHost, DG.DOCK_TYPE.RIGHT, null, 'PDB detail', 0.25);

    // Family cell color — any `family` cell in the central grid (which is
    // re-bound to a different DataFrame at each pipeline stage) gets painted
    // with the canonical FAMILY_MAP color, so the grid mirrors the 3D legend.
    // Works across consensus_model, ligand_features, pdb_interaction_summary —
    // all of which carry a `family` column.
    this.installFamilyCellRenderer();

    return this.tableView;
  }

  /**
   * Hook grid.onCellRender to paint any cell in a column named `family`
   * with FAMILY_MAP[code].hexColor. Single registration on the persistent
   * tableView.grid — survives `view.dataFrame = newDf` swaps the
   * orchestrator does between stages, so the consensus-model grid, the
   * Stage-4 ligand_features detail grid, and the per-PDB summary grid all
   * get the same family coloring without per-stage wiring.
   */
  private installFamilyCellRenderer(): void {
    if (!this.tableView) return;
    this.tableView.grid.onCellRender.subscribe((args) => {
      const col = args.cell.gridColumn?.column;
      if (!col || col.name !== 'family') return;
      const value = args.cell.cell.value;
      if (value == null) return;
      const fam = resolveFamily(String(value));
      const g = args.g;
      const b = args.bounds;
      // Paint a subtle background tinted by family. Pure hexColor would be
      // too saturated to keep the cell text readable, so we mix in white
      // at the cell-renderer level by drawing the chip + then drawing the
      // text at the platform default color.
      g.fillStyle = fam.hexColor + '33'; // ~20% alpha hex suffix
      g.fillRect(b.x, b.y, b.width, b.height);
      // Family color dot, left-aligned in the cell.
      g.beginPath();
      g.fillStyle = fam.hexColor;
      const cy = b.y + b.height / 2;
      const r = Math.min(4, b.height / 3);
      g.arc(b.x + 8, cy, r, 0, Math.PI * 2);
      g.fill();
      // Family name, default text color, offset past the dot.
      g.fillStyle = '#222';
      g.font = '12px sans-serif';
      g.textBaseline = 'middle';
      g.fillText(String(value), b.x + 18, cy + 0.5);
      args.preventDefault();
    });
  }

  /**
   * Build the persistent legend + caption strip docked at the top. Two rows:
   *   row 1: caption ("Run the pipeline ..." until a preview rewrites it)
   *   row 2: 7 colored chips, one per pharmacophore family, generated from
   *          family-map.ts so renaming/recoloring a family in one place
   *          updates the legend automatically.
   */
  private buildLegendBar(): HTMLElement {
    const captionEl = ui.divText(
      'Run a preview or Build to render the consensus pharmacophore.',
      'cp-molstar-caption');
    captionEl.style.fontSize = '12px';
    captionEl.style.color = '#444';
    captionEl.style.padding = '4px 10px 2px 10px';
    captionEl.style.fontStyle = 'italic';
    this.molstarCaptionEl = captionEl;

    const chipRow = ui.divH([], 'cp-legend-row');
    chipRow.style.flexWrap = 'wrap';
    chipRow.style.gap = '8px';
    chipRow.style.padding = '0 10px 4px 10px';
    chipRow.style.alignItems = 'center';
    for (const code of FAMILY_CODES) {
      const fam = FAMILY_MAP[code];
      const item = ui.divH([], 'cp-legend-item');
      item.style.gap = '4px';
      item.style.alignItems = 'center';
      item.style.fontSize = '11px';
      const dot = ui.div();
      dot.style.width = '10px';
      dot.style.height = '10px';
      dot.style.borderRadius = '50%';
      dot.style.background = fam.hexColor;
      dot.style.border = '1px solid rgba(0,0,0,0.2)';
      const label = ui.divText(`${code} · ${fam.name}`);
      label.style.color = '#444';
      item.append(dot, label);
      item.title = `Family ${code} (${fam.name}) — rendered as element ${fam.element} in CPK colors`;
      chipRow.append(item);
    }

    const wrapper = ui.divV([captionEl, chipRow], 'cp-legend-bar');
    wrapper.style.padding = '0';
    return wrapper;
  }

  /** Rewrite the persistent Mol* caption. Each preview calls this with text
   *  appropriate to what's about to render (or just rendered). Used so a new
   *  user looking at the right pane understands what the spheres represent. */
  private setMolstarCaption(text: string): void {
    if (!this.molstarCaptionEl) return;
    this.molstarCaptionEl.textContent = text;
  }

  // -------------------------------------------------------------------------
  // UI building blocks
  // -------------------------------------------------------------------------

  private buildStageHeader(): HTMLElement {
    const stages = [
      {id: 'stage1', label: '1. Enrich'},
      {id: 'stage2a', label: '2a. Align'},
      {id: 'picker', label: 'Picker'},
      {id: 'stage3', label: '3. Pocket'},
      {id: 'stage2b', label: '2b. Re-align'},
      {id: 'stage4', label: '4. Features'},
      {id: 'stage5a', label: '5a. Consensus'},
      {id: 'stage5b', label: '5b. Render'},
    ];
    const badges = stages.map(({id, label}) => {
      const el = ui.divText(label, 'cp-stage-badge cp-stage-pending');
      el.title = 'Pending';
      this.stageBadges.set(id, el);
      return el;
    });
    return ui.divH(badges, 'cp-stage-header');
  }

  private setStageStatus(id: string, status: 'pending' | 'running' | 'done' | 'error'): void {
    const el = this.stageBadges.get(id);
    if (!el) return;
    el.classList.remove('cp-stage-pending', 'cp-stage-running', 'cp-stage-done', 'cp-stage-error');
    el.classList.add(`cp-stage-${status}`);
    el.title = status[0].toUpperCase() + status.slice(1);
  }

  private resetStageStatuses(): void {
    for (const id of this.stageBadges.keys()) this.setStageStatus(id, 'pending');
  }

  private buildInputPanel(seedPdbIds?: string[]): HTMLElement {
    const initial = (seedPdbIds ?? []).join('\n');

    // Empty label — the section header "PDB IDs" above the textarea already
    // serves as the visible label; ui.input's own inline label would
    // duplicate it.
    const pdbInput = ui.input.textArea('', {
      value: initial,
      size: {width: 240, height: 120},
    });
    (pdbInput.root.querySelector('textarea') as HTMLTextAreaElement | null)?.setAttribute(
      'placeholder', 'One per line. e.g.\n3W2S\n4WKQ\n5HG8');

    const resolutionInput = ui.input.float('Max resolution (A)', {value: DEFAULT_OPTIONS.maxResolution});
    const requireXrayInput = ui.input.bool('X-ray only', {value: DEFAULT_OPTIONS.requireXray});
    const mwInput = ui.input.float('Min ligand MW (Da)', {value: DEFAULT_OPTIONS.minLigandMw});
    const pocketMethodInput = ui.input.choice<PocketMethod>('Pocket method', {
      value: DEFAULT_OPTIONS.pocketMethod,
      items: ['cutoff', 'dbscan'] as PocketMethod[],
    });
    const pocketRadiusInput = ui.input.float('Pocket radius (A)', {value: DEFAULT_OPTIONS.pocketRadius});
    const pocketRepInput = ui.input.choice<PocketRep>('Pocket rep', {
      value: DEFAULT_OPTIONS.pocketRep,
      items: ['spacefill', 'gaussian-surface'] as PocketRep[],
    });
    const refPdbInput = ui.input.string('Reference PDB (optional)', {value: ''});

    // Stage 5a consensus knobs — exposed so the user can re-cluster without
    // editing the Python script. kq trades cluster count for cluster size;
    // minClusterSizeFraction is the support threshold (0.5 = "at least half
    // the ligands"); topClusterNumber caps clusters per family so a
    // pathological run with many small clusters doesn't flood the grid.
    const kqInput = ui.input.int('Avg features per cluster (kq)',
      {value: DEFAULT_OPTIONS.kq, min: 1, max: 30});
    kqInput.root.title = 'k-means n_clusters = ceil(n_features / kq). Lower = more, smaller clusters; ' +
      'higher = fewer, fatter clusters. Default 7 matches TeachOpenCADD T009.';
    const minFracInput = ui.input.float('Min ligand fraction',
      {value: DEFAULT_OPTIONS.minClusterSizeFraction, min: 0.0, max: 1.0});
    minFracInput.root.title = 'A cluster is kept only if at least this fraction of distinct ligands ' +
      'contributed a feature to it. Drop to 0.5 to surface less-conserved (but still meaningful) hotspots.';
    const topClusterInput = ui.input.int('Max clusters per family',
      {value: DEFAULT_OPTIONS.topClusterNumber, min: 1, max: 10});
    topClusterInput.root.title = 'After filtering by min ligand fraction, keep at most N clusters per ' +
      'family (sorted by number of contributing ligands, then by cluster size).';

    const demoLabel = ui.label('Demo: 5 active-state EGFR kinase structures',
      {classes: 'cp-demo-label'});
    demoLabel.style.display = seedPdbIds && seedPdbIds.length ? 'block' : 'none';

    // Status line under the textarea — shows accepted/ignored token counts and
    // the first few ignored tokens by name so the user is never surprised by a
    // silent filter (e.g. pasting "EGFR, kinase, 1XKK" used to keep only 1XKK).
    const pdbStatus = ui.divText('', 'cp-pdb-status');
    pdbStatus.style.fontSize = '11px';
    pdbStatus.style.color = '#777';
    pdbStatus.style.marginTop = '2px';
    pdbStatus.style.minHeight = '14px';

    const refreshPdbStatus = (): void => {
      const {accepted, rejected} = parsePdbIdsDetailed(pdbInput.value ?? '');
      if (accepted.length === 0 && rejected.length === 0) {
        pdbStatus.textContent = '';
        pdbStatus.style.color = '#777';
        return;
      }
      let msg = `${accepted.length} PDB ID${accepted.length === 1 ? '' : 's'} accepted`;
      if (rejected.length > 0) {
        const shown = rejected.slice(0, 3).join(', ');
        const more = rejected.length > 3 ? `, +${rejected.length - 3} more` : '';
        msg += ` · ${rejected.length} ignored (${shown}${more})`;
        pdbStatus.style.color = '#B26A00';
      } else {
        pdbStatus.style.color = '#4CAF50';
      }
      pdbStatus.textContent = msg;
    };

    // Clear the demo label and refresh the status line as soon as the user
    // edits the input area.
    pdbInput.onChanged.subscribe(() => {
      demoLabel.style.display = 'none';
      refreshPdbStatus();
    });

    // Hero demo button. Bigger, more prominent than the old text link so a
    // first-time visitor's path to a result is one obvious click. The small
    // text link is kept below the textarea as a discreet fallback for users
    // who already edited the textarea and want to revert.
    const demoHero = ui.bigButton('Load EGFR demo (5 PDBs)', async () => {
      const ids = await loadDemoPdbIds();
      pdbInput.value = ids.join('\n');
      demoLabel.style.display = 'block';
      this.previewCache = null;
      refreshPdbStatus();
    });
    demoHero.classList.add('cp-demo-hero-btn');
    demoHero.title = '1XKK / 3W2S / 4WKQ / 5HG5 / 5HG8 — five active-state EGFR ' +
      'kinase structures bound to known inhibitors. Use this for a fast first run.';

    const demoLink = ui.link('Use demo input', async () => {
      const ids = await loadDemoPdbIds();
      pdbInput.value = ids.join('\n');
      demoLabel.style.display = 'block';
      this.previewCache = null;
      refreshPdbStatus();
    });

    // Target lookup — UniProt accession OR gene/protein name. Resolved
    // PDB IDs replace the textarea contents (per the user's preference).
    // Name lookups that map to multiple species open a picker dialog so
    // the user disambiguates human / mouse / Drosophila / etc. before
    // we commit the resolved PDB list.
    const lookupStatus = ui.divText('', 'cp-target-lookup-status');
    lookupStatus.style.fontSize = '11px';
    lookupStatus.style.color = '#777';
    lookupStatus.style.marginTop = '4px';
    lookupStatus.style.minHeight = '14px';
    const lookupInput = ui.input.string('Target', {
      value: '',
      nullable: true,
    });
    (lookupInput.root.querySelector('input') as HTMLInputElement | null)?.setAttribute(
      'placeholder', 'UniProt ID or name (e.g. P00533, egfr)');
    const lookupBtn = ui.button('Find PDBs', async () => {
      const query = (lookupInput.value ?? '').trim();
      if (!query) {
        grok.shell.warning('Enter a UniProt accession (e.g. P00533) or a gene/protein name (e.g. egfr).');
        return;
      }
      lookupBtn.disabled = true;
      lookupStatus.textContent = 'Searching UniProt...';
      try {
        const hit = await pickUniProtHit(query, lookupStatus);
        if (!hit) {
          lookupStatus.textContent = 'No target selected.';
          return;
        }
        lookupStatus.textContent =
          `Fetching PDB cross-references for ${hit.accession}...`;
        const pdbIds = await getPdbIdsForAccession(hit.accession);
        if (pdbIds.length === 0) {
          lookupStatus.textContent =
            `No PDB structures cross-referenced from ${hit.accession}.`;
          grok.shell.warning(
            `UniProt ${hit.accession} (${hit.organism}) has no PDB cross-references.`);
          return;
        }
        pdbInput.value = pdbIds.join('\n');
        demoLabel.style.display = 'none';
        this.previewCache = null;
        const orgShort = hit.organism.split(' ').slice(0, 2).join(' ');
        lookupStatus.textContent =
          `${hit.accession} · ${orgShort} · ${pdbIds.length} PDBs — listed below.`;
        grok.shell.info(`${hit.accession} ${hit.geneName || hit.proteinName} ` +
          `(${hit.organism}): ${pdbIds.length} PDB ID${pdbIds.length === 1 ? '' : 's'} ` +
          `loaded. Edit the list below if you want a subset.`);
      } catch (e: any) {
        const msg = e?.message ?? String(e);
        lookupStatus.textContent = `Lookup failed: ${msg.slice(0, 120)}`;
        grok.shell.error(`Target lookup failed: ${msg}`);
      } finally {
        lookupBtn.disabled = false;
      }
    }, 'Resolve a UniProt accession or gene/protein name to the list of PDB ' +
       'structures cross-referenced from that entry. Multi-species hits open ' +
       'a picker so you choose human/mouse/etc.');
    lookupBtn.classList.add('cp-target-lookup-btn');
    // Enter in the Target input triggers Find PDBs — saves the round-trip
    // to the button for keyboard users. We bind on the raw <input> rather
    // than ui.input.string's onChanged because the latter only fires on
    // value commit (focus loss), which would miss the Enter intent.
    const lookupInputEl = lookupInput.root.querySelector('input') as HTMLInputElement | null;
    lookupInputEl?.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' && !e.shiftKey && !lookupBtn.disabled) {
        e.preventDefault();
        lookupBtn.click();
      }
    });
    const lookupRow = ui.divH([lookupInput.root, lookupBtn], 'cp-target-lookup-row');
    lookupRow.style.alignItems = 'flex-end';
    lookupRow.style.gap = '6px';

    const readIds = (): string[] => parsePdbIds(pdbInput.value ?? '');
    const readOpts = (): PipelineOptions => ({
      maxResolution: resolutionInput.value ?? DEFAULT_OPTIONS.maxResolution,
      requireXray:   requireXrayInput.value ?? DEFAULT_OPTIONS.requireXray,
      // Legacy UI exposes a single "X-ray only" checkbox; map it to the new
      // 4-flag scheme: if X-ray-only is true, accept only X-ray; otherwise
      // accept all 3 experimental methods. AlphaFold stays off.
      allowXray: true,
      allowNmr: !(requireXrayInput.value ?? DEFAULT_OPTIONS.requireXray),
      allowCryoEm: !(requireXrayInput.value ?? DEFAULT_OPTIONS.requireXray),
      allowAlphaFold: false,
      minLigandMw:   mwInput.value ?? DEFAULT_OPTIONS.minLigandMw,
      pocketMethod:  (pocketMethodInput.value as PocketMethod) ?? DEFAULT_OPTIONS.pocketMethod,
      pocketRadius:  pocketRadiusInput.value ?? DEFAULT_OPTIONS.pocketRadius,
      pocketRep:     (pocketRepInput.value as PocketRep) ?? DEFAULT_OPTIONS.pocketRep,
      refPdbId:      (refPdbInput.value ?? '').trim() || undefined,
      kq:                     kqInput.value ?? DEFAULT_OPTIONS.kq,
      minClusterSizeFraction: minFracInput.value ?? DEFAULT_OPTIONS.minClusterSizeFraction,
      topClusterNumber:       topClusterInput.value ?? DEFAULT_OPTIONS.topClusterNumber,
    });

    // Five incremental preview buttons — each shows progressively more of the
    // pipeline output without committing to the full Build. After a successful
    // run the button gets a green `.cp-preview-btn-done` class with a "✓ "
    // prefix so the user can see at a glance which stages they've already
    // exercised (mirrors the stage-badge done state at the top of the view).
    // The class is cleared on textarea edits because that invalidates the
    // previewCache — subsequent clicks recompute from scratch.
    const makePreviewBtn = (
      label: string, tooltip: string, run: (ids: string[], opts: PipelineOptions) => Promise<void>,
    ): HTMLButtonElement => {
      const btn = ui.button(label, async () => {
        const ids = readIds();
        if (ids.length === 0) {
          grok.shell.warning('Paste at least one PDB ID, or click "Load EGFR demo".');
          return;
        }
        previewButtons.forEach((b) => { b.disabled = true; });
        buildBtn.disabled = true;
        try {
          await run(ids, readOpts());
          btn.classList.add('cp-preview-btn-done');
        } catch (_e) {
          // Failed runs don't get the green check — the existing toast +
          // error stage badge already report the problem.
        } finally {
          previewButtons.forEach((b) => { b.disabled = false; });
          buildBtn.disabled = false;
        }
      }, tooltip);
      return btn;
    };

    const fetchBtn = makePreviewBtn('1. Fetch PDBs',
      'Download the input PDBs from RCSB and show them in their native crystallographic frames.',
      (ids, opts) => this.previewFetch(ids, opts));
    const alignBtn = makePreviewBtn('2. Align',
      'Superimpose the PDBs via Stage 2a (Calpha Kabsch). Reuses fetched data when inputs are unchanged.',
      (ids, opts) => this.previewAlign(ids, opts));
    const pocketBtn = makePreviewBtn('3. Pocket',
      'Run Stage 3 and engage Mol*\'s native binding-site renderer (highlights residues within the cutoff radius of each drug ligand).',
      (ids, opts) => this.previewPocket(ids, opts));
    const featuresBtn = makePreviewBtn('4. Features',
      'Run Stage 4 (ProLIF) and overlay productive interaction sites (chain F, colored by family). ~15-30s per PDB.',
      (ids, opts) => this.previewFeatures(ids, opts));
    const consensusBtn = makePreviewBtn('5. Consensus',
      'Run Stage 5a (k-means per family) and overlay the consensus pharmacophore (chain P, B-factor = frequency).',
      (ids, opts) => this.previewConsensus(ids, opts));
    const previewButtons: HTMLButtonElement[] = [fetchBtn, alignBtn, pocketBtn, featuresBtn, consensusBtn];

    const buildBtn = ui.bigButton('Run full pipeline', async () => {
      const ids = readIds();
      if (ids.length === 0) {
        grok.shell.warning('Paste at least one PDB ID, or click "Load EGFR demo".');
        return;
      }
      buildBtn.disabled = true;
      previewButtons.forEach((b) => { b.disabled = true; });
      try {
        await this.runPipeline(ids, readOpts());
      } finally {
        buildBtn.disabled = false;
        previewButtons.forEach((b) => { b.disabled = false; });
      }
    });
    buildBtn.title = 'Run all stages end-to-end. Unlike the preview buttons, ' +
      'this also opens the cluster picker (between Stage 2a and Stage 3) so ' +
      'you can drop conformationally mixed PDBs before pocket isolation.';

    // Editing the PDB list invalidates any cached intermediate results.
    // Also clear the `.cp-preview-btn-done` checkmarks so a stale green tick
    // doesn't suggest a finding that no longer reflects the (now-different)
    // PDB list — next click recomputes the whole stage from scratch.
    pdbInput.onChanged.subscribe(() => {
      this.previewCache = null;
      previewButtons.forEach((b) => b.classList.remove('cp-preview-btn-done'));
    });

    // Stage 5a knob changes invalidate ONLY the consensus_model cache —
    // Stage 1-4 outputs (pdbQc, aligned, pocketAtoms, ligandFeatures) are
    // unaffected so we don't want a full previewCache wipe (which would
    // force a re-run of the 30s ProLIF step). Same for the consensusBtn's
    // ✓ done marker so the user knows their last consensus is stale.
    const invalidateConsensus = () => {
      if (this.previewCache) this.previewCache.consensusModel = undefined;
      consensusBtn.classList.remove('cp-preview-btn-done');
    };
    kqInput.onChanged.subscribe(invalidateConsensus);
    minFracInput.onChanged.subscribe(invalidateConsensus);
    topClusterInput.onChanged.subscribe(invalidateConsensus);

    const resetBtn = this.buildResetBtn(
      pdbInput, lookupInput, pdbStatus, lookupStatus, previewButtons, demoLabel);

    // Initial paint of the status line for the seeded demo PDB list (if any).
    refreshPdbStatus();

    // Advanced options — collapsed by default. New users only see the
    // textarea + target lookup + preview/run buttons; power users open the
    // accordion to tune the QC / pocket / reference-PDB knobs.
    const advancedBody = ui.divV([
      ui.h3('QC'),
      resolutionInput.root,
      requireXrayInput.root,
      mwInput.root,
      ui.h3('Pocket'),
      pocketMethodInput.root,
      pocketRadiusInput.root,
      pocketRepInput.root,
      ui.h3('Consensus'),
      kqInput.root,
      minFracInput.root,
      topClusterInput.root,
      ui.h3('Output'),
      refPdbInput.root,
    ], 'cp-advanced-body');
    const advancedToggle = ui.divText('▸ Advanced options', 'cp-advanced-toggle');
    advancedToggle.style.cursor = 'pointer';
    advancedToggle.style.color = '#1976D2';
    advancedToggle.style.fontSize = '12px';
    advancedToggle.style.fontWeight = 'bold';
    advancedToggle.style.padding = '8px 0 4px 0';
    advancedToggle.style.userSelect = 'none';
    advancedBody.style.display = 'none';
    advancedToggle.addEventListener('click', () => {
      const open = advancedBody.style.display === 'none';
      advancedBody.style.display = open ? 'block' : 'none';
      advancedToggle.textContent = open ? '▾ Advanced options' : '▸ Advanced options';
    });

    // Caption above the preview-stage column — answers "what's the
    // difference between these 5 buttons and the big Run button below?".
    const previewHint = ui.divText(
      'Click any stage to preview just up to that point.', 'cp-preview-hint');
    previewHint.style.fontSize = '11px';
    previewHint.style.color = '#777';
    previewHint.style.marginBottom = '4px';

    return ui.divV([
      ui.h2('Consensus Pharmacophore'),
      demoLabel,
      demoHero,
      ui.h3('Target lookup'),
      lookupRow,
      lookupStatus,
      ui.h3('PDB IDs'),
      pdbInput.root,
      pdbStatus,
      ui.div([demoLink], 'cp-demo-link-row'),
      advancedToggle,
      advancedBody,
      ui.h3('Preview stages'),
      previewHint,
      ui.divV([fetchBtn, alignBtn, pocketBtn, featuresBtn, consensusBtn], 'cp-preview-btn-col'),
      ui.h3('Run'),
      buildBtn,
      resetBtn,
    ], 'cp-input-panel');
  }

  /** Build the reset button. Clears textarea + lookup field + status lines,
   *  resets stage badges, invalidates the preview cache, removes any
   *  preview "✓ done" markers, and closes the detail panel. Mol* state is
   *  left as-is — the next preview/run will overwrite it. */
  private buildResetBtn(
    pdbInput: DG.InputBase<string | null>,
    lookupInput: DG.InputBase<string | null>,
    pdbStatus: HTMLElement,
    lookupStatus: HTMLElement,
    previewButtons: HTMLButtonElement[],
    demoLabel: HTMLElement,
  ): HTMLButtonElement {
    const btn = ui.button('Reset', () => {
      pdbInput.value = '';
      lookupInput.value = '';
      pdbStatus.textContent = '';
      lookupStatus.textContent = '';
      demoLabel.style.display = 'none';
      this.previewCache = null;
      previewButtons.forEach((b) => b.classList.remove('cp-preview-btn-done'));
      this.resetStageStatuses();
      this.closePdbDetailPanel();
      this.setMolstarCaption(
        'Reset — paste PDB IDs (or use Load EGFR demo / target lookup) and ' +
        'click a preview button or Run full pipeline.');
    }, 'Clear the PDB list, target lookup, stage badges, and preview ' +
       'cache. Leaves the 3D viewer at its last state until you run again.');
    btn.classList.add('cp-reset-btn');
    return btn;
  }

  // -------------------------------------------------------------------------
  // Pipeline
  // -------------------------------------------------------------------------

  async runPipeline(pdbIds: string[], options: PipelineOptions): Promise<PipelineResult | null> {
    // tableView may be undefined in the wizard path — that's OK because
    // setActiveDataFrame() guards every dataFrame assignment.
    this.aborted = false;
    this.resetStageStatuses();
    this.closePdbDetailPanel();
    this.setMolstarCaption(
      'Running full pipeline — cluster picker will open after Stage 2a to let ' +
      'you drop conformationally mixed PDBs before pocket isolation.');

    // Wizard integration: populate previewCache as we go so the per-step
    // panels in wizard-steps.ts can render result tables for stages 1–5.
    const cache = this.ensurePreviewCache(this.previewFingerprint(pdbIds, options));

    const pi = DG.TaskBarProgressIndicator.create('Consensus Pharmacophore: running...');
    try {
      // Stage 1 — Enrich PDB metadata (RCSB-backed; replaces stubStage1)
      this.setStageStatus('stage1', 'running');
      pi.update(5, 'Stage 1: enriching PDB metadata...');
      const stage1 = await enrichPdbList(pdbIds, options);
      const pdbQc = stage1.accepted;
      if (pdbQc.rowCount === 0)
        throw new Error('No usable PDB entries after QC filtering. ' +
          'Check the PDB IDs, loosen the resolution cap, or disable the X-ray filter.');
      this.setActiveDataFrame(pdbQc);
      pdbQc.name = 'pdb_qc';
      cache.pdbQc = pdbQc;
      this._droppedPdbs = stage1.dropped;
      this.setStageStatus('stage1', 'done');
      if (this.aborted) return null;

      // Stage 2a — Kabsch pass 1
      this.setStageStatus('stage2a', 'running');
      pi.update(20, 'Stage 2a: aligning structures (pass 1)...');
      const aligned1 = await runStage2aAlign(pdbQc, options);
      this.setActiveDataFrame(aligned1);
      aligned1.name = 'aligned_structures (pass 1)';
      cache.aligned = aligned1;
      this.setStageStatus('stage2a', 'done');
      if (this.aborted) return null;

      // Cluster picker (TS dialog) — pass pdb_id -> ligand_smiles map from Stage 1
      // so the dialog can render ligand thumbnails on each cluster card.
      this.setStageStatus('picker', 'running');
      const smilesByPdbId = ligandSmilesMap(pdbQc);
      const pick = await pickClusterViaDialog(aligned1, smilesByPdbId);
      if (pick === null) {
        this.setStageStatus('picker', 'error');
        grok.shell.info('Pipeline cancelled at cluster picker.');
        return null;
      }
      this.setStageStatus('picker', 'done');
      if (this.aborted) return null;

      // Stage 3 — Pocket isolation
      this.setStageStatus('stage3', 'running');
      pi.update(40, 'Stage 3: isolating pocket...');
      const pocketAtoms = await runStage3IsolatePocket(aligned1, pick.clusterIds, options);
      this.setActiveDataFrame(pocketAtoms);
      pocketAtoms.name = 'pocket_atoms';
      cache.pocketAtoms = pocketAtoms;
      this.setStageStatus('stage3', 'done');
      if (this.aborted) return null;

      // Stage 2b — Kabsch pass 2 (pocket-only Calpha; filtered to picker's PDB list)
      this.setStageStatus('stage2b', 'running');
      pi.update(55, 'Stage 2b: re-aligning on pocket Calpha...');
      const aligned2 = await runStage2bAlignPocket(aligned1, pocketAtoms, pick.selectedPdbIds);
      this.setActiveDataFrame(aligned2);
      aligned2.name = 'aligned_structures (pass 2)';
      cache.alignedV2 = aligned2;
      this.setStageStatus('stage2b', 'done');
      if (this.aborted) return null;

      // Stage 4 — Protein-ligand interaction extraction (ProLIF)
      this.setStageStatus('stage4', 'running');
      pi.update(70, 'Stage 4: detecting protein-ligand interactions (ProLIF, ~15-30s/PDB)...');
      const ligandFeatures = await runStage4ExtractFeatures(aligned2, pocketAtoms, pick.selectedPdbIds);
      this.setActiveDataFrame(ligandFeatures);
      ligandFeatures.name = 'ligand_features';
      cache.ligandFeatures = ligandFeatures;
      // Cache the per-PDB summary so the wizard's Step 4 panel can render it.
      const smilesByPdb = cache.pdbQc ? ligandSmilesMap(cache.pdbQc) : undefined;
      this._summaryDf = aggregateFeaturesByPdb(ligandFeatures, smilesByPdb);
      this._summaryDf.name = 'pdb_interaction_summary';
      this.setStageStatus('stage4', 'done');
      if (this.aborted) return null;

      // Stage 5a — k-means consensus (honor any PDB/interaction exclusions)
      this.setStageStatus('stage5a', 'running');
      pi.update(85, 'Stage 5a: computing per-family consensus...');
      const consensus = await runStage5aConsensusKmeans(
        filterFeaturesForConsensus(ligandFeatures, options), options);
      this.setActiveDataFrame(consensus);
      consensus.name = 'consensus_model';
      cache.consensusModel = consensus;
      this.setStageStatus('stage5a', 'done');
      if (this.aborted) return null;

      // Stage 5b — TS renderer (REAL, not stub)
      this.setStageStatus('stage5b', 'running');
      pi.update(95, 'Stage 5b: rendering consensus pharmacophore...');
      const pdbBlock = await renderConsensusPharmacophore(consensus, options.refPdbId);
      await this.refreshDockedViewer(pdbBlock);
      this.setStageStatus('stage5b', 'done');

      pi.update(100, 'Done.');
      return {consensusModel: consensus, pdbBlock};
    } catch (e: any) {
      const failing = currentRunningStage(this.stageBadges) ?? 'pipeline';
      this.setStageStatus(failing, 'error');
      grok.shell.error(`${failing} failed: ${e?.message ?? e}`);
      throw e;
    } finally {
      pi.close();
    }
  }

  /**
   * Push a new PDB into the docked Mol* viewer; lazy-create the viewer on first call.
   * Defensive re-assertion per advisor Q3: `view.dataFrame =` swaps survive the dock binding
   * but we re-call setOptions to be safe.
   *
   * NOTE: B-factor coloring depends on Mol*'s default behaviour (Phase 0 V1 verification
   * pending). If V1 returns "PASS-with-toggle", surface the property here once the actual
   * knob is identified (likely a `representation` extension).
   */
  private async refreshDockedViewer(
    pdbBlock: string,
    extraOptions: Record<string, unknown> = {},
  ): Promise<void> {
    const allOptions = {pdb: pdbBlock, ...extraOptions};

    // The viewer is about to be torn down + rebuilt. Any cached Mol* state
    // refs (per-PDB aligned structures + per-PDB interaction overlays)
    // become invalid once the old viewer is closed — clear both maps
    // proactively so isolatePdb doesn't try to mutate dead refs.
    this._alignedStructureRefs.clear();
    this.clearOverlayRefs();

    // --- Wizard path: mount the Mol* viewer directly into the injected
    // host element, no dock manager. The host is the wizard's
    // .cp-wizard-content-center div, mounted ONCE for the lifetime of the
    // view — only the viewer child changes between stages.
    if (this._externalMolstarHost) {
      // Force-rebuild (same rationale as the legacy path — BSV's setOptions
      // viewSyncer reshuffle leaves cells in an inconsistent state).
      if (this.viewer) {
        try { (this.viewer as any).close?.(); } catch (_e) { /* */ }
        this.viewer = undefined;
      }
      ui.empty(this._externalMolstarHost);
      // A throw-away DataFrame so plot.fromType can be called. We could use
      // this.tableView.dataFrame here, but the wizard path has no tableView,
      // and plot.fromType is happy with a 1-column 1-row DataFrame.
      const tempDf = DG.DataFrame.fromColumns([
        DG.Column.fromStrings('pdb_id', ['']),
      ]);
      this.viewer = await tempDf.plot.fromType(
        'Biostructure', allOptions) as DG.Viewer;
      const vroot = (this.viewer as any).root ?? (this.viewer as DG.Viewer).root;
      if (vroot) {
        // Make the viewer fill its host. Without this, BSV's default sizing
        // produces a small canvas in the upper-left corner of the dark host.
        (vroot as HTMLElement).style.width = '100%';
        (vroot as HTMLElement).style.height = '100%';
        (vroot as HTMLElement).style.position = 'absolute';
        (vroot as HTMLElement).style.inset = '0';
        this._externalMolstarHost.append(vroot);
      }
      return;
    }

    if (!this.tableView) return;

    // Force-rebuild the Mol* viewer on every refresh. BSV's setOptions
    // performs an asynchronous viewSyncer reshuffle when `showBindingSite`
    // changes — empirically that leaves plugin.state.data.cells in a state
    // where iterating for cartoon/binding-site cells finds none (verified by
    // a [Mol* style] (0) console line that then retried 6× and gave up,
    // leaving the protein at full opacity and hiding the chain-F overlays
    // inside it). Tearing down and re-creating the viewer guarantees a
    // clean state tree and is only ~1–2 s of wall time on a 5-PDB demo.
    if (this.viewer && this.viewerDockNode) {
      try {
        const oldViewer = this.viewer;
        this.tableView.dockManager.close(this.viewerDockNode);
        try { (oldViewer as any).close?.(); } catch (_e) { /* */ }
      } catch (e) {
        console.warn('[refresh] failed to tear down old viewer cleanly:', e);
      }
      this.viewer = undefined;
      this.viewerDockNode = undefined;
    }

    // Close the empty-state placeholder BEFORE docking the new viewer.
    // Otherwise the dock manager has to fit three things in the right area
    // (placeholder + new viewer + detail), which leaves the input panel
    // squeezed to ~80 px wide and pushes detail past the viewport edge.
    if (this.viewerPlaceholderNode) {
      this.tableView.dockManager.close(this.viewerPlaceholderNode);
      this.viewerPlaceholderNode = undefined;
    }
    this.viewer = await this.tableView.dataFrame.plot.fromType(
      'Biostructure', allOptions) as DG.Viewer;
    // Dock the new viewer to the LEFT of the detail panel (which was
    // pre-docked at init time at root-level RIGHT). Using the detail node
    // as the neighbor keeps the layout stable: grid | Mol* | detail panel.
    // If we instead dock at root-level RIGHT with null neighbor, the dock
    // manager makes the new viewer the FAR-right sibling of the detail
    // panel, displacing detail into the middle (observed empirically).
    //
    // Ratio = 0.6 gives Mol* ~60% of the (Mol* + detail) horizontal area,
    // detail ~40%. The Stage-4 detail panel shows a 5-column grid (PDB ID,
    // ligand, residue, family, interaction) — at 40% of the available width
    // every row is visible without horizontal scrolling.
    if (this.pdbDetailNode) {
      this.viewerDockNode = this.tableView.dockManager.dock(
        this.viewer as unknown as DG.Viewer,
        DG.DOCK_TYPE.LEFT, this.pdbDetailNode, 'Mol*', 0.6);
    } else {
      // Fallback if detail node somehow wasn't created (defensive).
      this.viewerDockNode = this.tableView.dockManager.dock(
        this.viewer as unknown as DG.Viewer, DG.DOCK_TYPE.RIGHT, null, 'Mol*', 0.55);
    }
  }

  /**
   * Wait until BSV's MolstarViewer has built the binding-site components (i.e. the
   * `bindingSiteRefs` array is non-empty). BSV applies the binding site on a
   * `viewSyncer` queue after we toggle `showBindingSite: true`, so a brief poll
   * lets us patch the Mol* representation afterward without a race.
   */
  private async waitForBindingSiteApplied(
    previousRefsKey: string | null = null, maxMs = 4000,
  ): Promise<boolean> {
    if (!this.viewer) return false;
    const start = Date.now();
    while (Date.now() - start < maxMs) {
      const refs = (this.viewer as any).bindingSiteRefs;
      if (Array.isArray(refs) && refs.length > 0) {
        const curKey = JSON.stringify(refs);
        if (previousRefsKey === null || curKey !== previousRefsKey) {
          // Brief settling delay — viewSyncer may still be finalizing reps.
          await new Promise((r) => setTimeout(r, 300));
          return true;
        }
      }
      await new Promise((r) => setTimeout(r, 100));
    }
    return false;
  }

  private snapshotBindingSiteRefs(): string {
    return JSON.stringify((this.viewer as any)?.bindingSiteRefs ?? []);
  }

  /**
   * Poll until BSV's Mol* plugin has at least one Structure 3D cell loaded.
   * Needed when `showBindingSite: false` (so there's no binding-site refs to
   * wait on) and we still want to apply representation tweaks or load an
   * additional overlay structure. Without this poll, setMolstarStyle and
   * addInteractionOverlay can race against Mol*'s async PDB parser, hit a
   * not-yet-populated `plugin.state.data.cells`, and early-return — leaving
   * the viewer in its default style with no overlay.
   */
  private async waitForMolstarReady(maxMs = 6000): Promise<boolean> {
    if (!this.viewer) return false;
    const start = Date.now();
    while (Date.now() - start < maxMs) {
      const plugin = (this.viewer as any).viewer?.plugin;
      const cells = plugin?.state?.data?.cells;
      if (cells) {
        // Any Structure 3D cell means Mol* has at least loaded the PDB into
        // its trajectory pipeline. That's the earliest moment we can safely
        // mutate the state tree.
        for (const [, cell] of cells) {
          if (cell?.obj?.type?.name === 'Structure 3D') {
            await new Promise((r) => setTimeout(r, 150));
            return true;
          }
        }
      }
      await new Promise((r) => setTimeout(r, 100));
    }
    console.warn('[Mol* ready] timed out waiting for Structure 3D cells');
    return false;
  }

  /**
   * Style the docked Mol* for the Pocket view:
   *   - polymer cartoon alpha = `cartoonAlpha` (translucent for Pocket, opaque otherwise)
   *   - if binding-site components exist, switch their rep to `bindingSiteRep`
   *     (`spacefill` = dots, `gaussian-surface` = surface, `ball-and-stick` = Mol*'s default)
   *   - water / ion / non-standard alpha = `solventAlpha` (0 hides them entirely)
   *
   * Tags BSV uses for its built-in static components (see binding-site.ts):
   *   structure-component-static-{polymer,ligand,non-standard,water,ion,coarse,branched}
   */
  private async setMolstarStyle(opts: {
    cartoonAlpha?: number;
    bindingSiteRep?: 'spacefill' | 'gaussian-surface' | 'ball-and-stick' | null;
    bindingSiteAlpha?: number;
    solventAlpha?: number;
  } = {}): Promise<{nCartoon: number; nBinding: number; nSolvent: number}> {
    if (!this.viewer) return {nCartoon: 0, nBinding: 0, nSolvent: 0};
    const v = this.viewer as any;
    const plugin = v.viewer?.plugin;
    if (!plugin) return {nCartoon: 0, nBinding: 0, nSolvent: 0};
    const cartoonAlpha = opts.cartoonAlpha ?? 1.0;
    const bindingSiteRep = opts.bindingSiteRep ?? null;
    const bindingSiteAlpha = opts.bindingSiteAlpha ?? 1.0;
    const solventAlpha = opts.solventAlpha ?? 1.0;

    const bindingSiteRefs = new Set<string>(v.bindingSiteRefs ?? []);
    const SOLVENT_TAGS = new Set([
      'structure-component-static-water',
      'structure-component-static-ion',
      'structure-component-static-non-standard',
    ]);
    // Build parent-ref -> tag map so we know what each rep's component is.
    const parentTagByRef = new Map<string, string>();
    for (const [ref, cell] of plugin.state.data.cells) {
      const tag = cell?.transform?.tags?.[0];
      if (tag) parentTagByRef.set(ref, tag);
    }

    const update = plugin.build();
    let nCartoon = 0; let nBinding = 0; let nSolvent = 0;
    for (const [ref, cell] of plugin.state.data.cells) {
      if (cell?.obj?.type?.name !== 'Structure 3D') continue;
      const params = cell?.transform?.params;
      const parent = cell?.transform?.parent;
      const parentTag = parentTagByRef.get(parent);

      if (params?.type?.name === 'cartoon') {
        update.to(ref).update((old: any) => ({
          ...old,
          type: {...old.type, params: {...old.type.params, alpha: cartoonAlpha}},
        }));
        nCartoon += 1;
      }
      if (bindingSiteRep && bindingSiteRefs.has(parent)) {
        update.to(ref).update((old: any) => ({
          ...old,
          type: {name: bindingSiteRep, params: {...(old.type?.params ?? {}), alpha: bindingSiteAlpha}},
        }));
        nBinding += 1;
      }
      if (parentTag && SOLVENT_TAGS.has(parentTag)) {
        update.to(ref).update((old: any) => ({
          ...old,
          type: {...old.type, params: {...old.type.params, alpha: solventAlpha}},
        }));
        nSolvent += 1;
      }
    }
    await update.commit();
    console.log(`[Mol* style] cartoonAlpha=${cartoonAlpha} (${nCartoon}), ` +
      `bindingSiteRep=${bindingSiteRep ?? '(none)'} alpha=${bindingSiteAlpha} (${nBinding}), ` +
      `solventAlpha=${solventAlpha} (${nSolvent})`);
    return {nCartoon, nBinding, nSolvent};
  }

  /**
   * Retry setMolstarStyle with increasing delay until it finds the expected
   * number of cartoon cells. The plain setMolstarStyle race-loses against
   * Mol*'s async cell refresh especially after addPerRowOverlay adds many
   * structures back-to-back — the cartoon cell from the original protein
   * structure is momentarily absent while the state tree is being rebuilt.
   * This wrapper polls until we hit at least one cartoon cell or give up.
   */
  private async setMolstarStyleWithRetry(opts: {
    cartoonAlpha?: number;
    bindingSiteRep?: 'spacefill' | 'gaussian-surface' | 'ball-and-stick' | null;
    bindingSiteAlpha?: number;
    solventAlpha?: number;
  }, maxAttempts = 6): Promise<void> {
    for (let attempt = 0; attempt < maxAttempts; attempt++) {
      const result = await this.setMolstarStyle(opts);
      if (result.nCartoon > 0) return;
      // Cell tree wasn't ready yet — back off and try again. Each retry
      // doubles the wait (100, 200, 400, 800, ...) up to ~3s total.
      await new Promise((r) => setTimeout(r, 100 * Math.pow(2, attempt)));
    }
    console.warn(`[Mol* style] retried ${maxAttempts}× without finding a ` +
      `cartoon cell — protein may render at full opacity`);
  }

  /**
   * Load `features` as a SEPARATE Mol* structure (label `consensus-features`)
   * with a spacefill representation at `sizeFactor` (default ~1.0; we use 2.5
   * for Features/Consensus to make each interaction dot clearly bigger than
   * the binding-site spacefill spheres). Posting features as a second
   * structure rather than splicing them into the main combined PDB keeps Mol*
   * from absorbing them into the static-ligand / binding-site partitions of
   * the protein structure, which would force them to share a representation
   * with the drug ligand.
   *
   * Each call clears any previous `consensus-features` structure first so
   * repeated Features/Consensus clicks don't pile overlays up.
   */
  /**
   * Mount the chain-F interaction overlay as N separate Mol* structures
   * (one per PDB), tracking each ref in `_alignedInteractionRefs`. Used by
   * Step 4 (Features) so the row-click isolation flow can hide the overlay
   * of every PDB except the clicked one — without rebuilding the viewer.
   *
   * Cleans up any previous overlay structures (cells whose label starts
   * with `cp-overlay`) before mounting new ones. Each PDB gets its own
   * Data → Trajectory → Model → Structure → ligand component → spacefill
   * rep subtree, identical to the single-overlay path in
   * `addInteractionOverlay` but filtered to one PDB at a time.
   */
  private async addInteractionOverlayPerPdb(
    features: DG.DataFrame, sizeFactor = 2.5, chainLetter = 'F',
  ): Promise<void> {
    this.clearOverlayRefs();
    if (!this.viewer) return;
    const plugin = (this.viewer as any).viewer?.plugin;
    if (!plugin) {
      console.warn('[overlay] plugin not reachable; per-PDB overlay skipped');
      return;
    }

    // Clean up stale overlay cells (from a previous Step 4 / Step 5 run).
    const stale: string[] = [];
    for (const [ref, cell] of plugin.state.data.cells) {
      const label = cell?.transform?.params?.label ?? cell?.obj?.label;
      if (typeof label === 'string' && label.startsWith('cp-overlay'))
        stale.push(ref);
    }
    if (stale.length > 0) {
      const cleanup = plugin.build();
      for (const ref of stale) cleanup.delete(ref);
      try { await cleanup.commit(); } catch (e) { console.warn('[overlay] cleanup failed', e); }
    }

    const pdbIdCol = features.col('pdb_id');
    if (!pdbIdCol) {
      console.warn('[overlay] features has no pdb_id column; falling back to single-overlay path');
      await this.addInteractionOverlay(features, sizeFactor, {chain: chainLetter});
      return;
    }
    const famCol = features.col('family');
    const xCol = features.col('x');
    const yCol = features.col('y');
    const zCol = features.col('z');
    const skipCol = features.col('skip_reason');
    const f = (v: number): string => v.toFixed(3).padStart(8);

    // ONE overlay structure PER INTERACTION (not per family). Mol*'s
    // toggleVisibility works at structure granularity, so a single shared
    // per-family structure could not hide an individual contact — un-ticking
    // one interaction in the detail panel calls setInteractionHidden(key),
    // which hides exactly that structure. Each sphere is still uniform-colored
    // with its family's legend hex so the 3D view matches the legend chips.
    // Every interaction structure is tracked under its PDB id (for row-click
    // isolation) AND keyed by consensusFeatureKey (for per-interaction hide and
    // the isolate path's interaction-exclusion check).
    let totalAtoms = 0;
    for (let i = 0; i < features.rowCount; i++) {
      if (skipCol && String(skipCol.get(i) ?? '').trim() !== '') continue;
      const pdbId = String(pdbIdCol.get(i) ?? '').toUpperCase().trim();
      if (!pdbId) continue;
      const x = Number(xCol?.get(i));
      const y = Number(yCol?.get(i));
      const z = Number(zCol?.get(i));
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) continue;
      const fam = resolveFamily(String(famCol?.get(i) ?? ''));
      // Legend hex '#RRGGBB' → Mol* Color int (0xRRGGBB).
      const colorInt = parseInt(fam.hexColor.slice(1), 16);
      // Single-atom HETATM block (same column layout as the working consensus
      // overlay): element/resName from FAMILY_MAP drive the vdW sphere size.
      const atomNameField = fam.element.length === 2
        ? fam.element.padEnd(4)
        : ' ' + fam.element.padEnd(3);
      const block =
        `REMARK   Interaction ${pdbId} #${i} (${fam.name})\n` +
        `HETATM    1 ${atomNameField} ${fam.resName.padEnd(3)} ${chainLetter}   1    ` +
        f(x) + f(y) + f(z) +
        '  1.00  1.00          ' + fam.element.padStart(2) + '\nEND\n';
      const key = consensusFeatureKey(features, i);
      try {
        const data = await plugin.builders.data.rawData(
          {data: block, label: `cp-overlay-${pdbId}-${i}`});
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);
        const component = await plugin.builders.structure.tryCreateComponentStatic(
          structure, 'ligand', {label: `cp-overlay-component-${pdbId}-${i}`});
        if (component) {
          await plugin.builders.structure.representation.addRepresentation(component, {
            type: 'spacefill',
            // Slightly transparent so the spheres read as feature "clouds"
            // and the ligand/pocket behind them stays visible.
            typeParams: {sizeFactor, alpha: 0.7},
            color: 'uniform',
            colorParams: {value: colorInt},
          });
        }
        // Track the ROOT structure's ref directly from the builder return.
        const ref = (structure as any).ref;
        if (ref) {
          if (!this._alignedInteractionRefs.has(pdbId))
            this._alignedInteractionRefs.set(pdbId, []);
          this._alignedInteractionRefs.get(pdbId)!.push(ref);
          this._interactionRefByKey.set(key, ref);
          this._interactionKeyByRef.set(ref, key);
          totalAtoms += 1;
        }
      } catch (e) {
        console.error(`[overlay] failed for ${pdbId} interaction ${i}:`, e);
      }
    }
    console.log(`[overlay] tracked ${this._alignedInteractionRefs.size} PDB(s), ` +
      `${this._interactionRefByKey.size} per-interaction spheres (family-colored)`);
  }

  /**
   * Mount the Stage-3 pocket atoms as N separate Mol* structures (one per
   * PDB), tracking each ref in `_alignedInteractionRefs`. Used by Step 3
   * (Pocket) so the row-click isolation flow can hide the pocket overlay of
   * every PDB except the clicked one — exactly mirrors
   * `addInteractionOverlayPerPdb`, but draws Cα + ligand-seed pocket atoms
   * (chain P) instead of chain-F interactions. The rep type follows the
   * user's Pocket-rep choice (spacefill dots vs gaussian-surface).
   *
   * Reusing `_alignedInteractionRefs` is safe: only one step's overlays exist
   * at a time because navigation rebuilds the viewer (`syncViewerToStep`).
   */
  private async addPocketOverlayPerPdb(
    // `sizeFactor` multiplies each atom's vdW radius in Mol*'s spacefill rep.
    // ~0.5 renders the pocket Cα + ligand-seed atoms as fine "dots" (the look
    // the user asked to bring back); 1.6 made chunky, overlapping spheres that
    // read as a blob rather than dots. Bump back up for a denser spacefill.
    pocketAtoms: DG.DataFrame, rep: PocketRep, sizeFactor = 0.5,
  ): Promise<void> {
    this.clearOverlayRefs();
    if (!this.viewer) return;
    const plugin = (this.viewer as any).viewer?.plugin;
    if (!plugin) {
      console.warn('[pocket-overlay] plugin not reachable; per-PDB overlay skipped');
      return;
    }

    // Clean up stale overlay cells (cp-overlay*) from a previous step's render.
    const stale: string[] = [];
    for (const [ref, cell] of plugin.state.data.cells) {
      const label = cell?.transform?.params?.label ?? cell?.obj?.label;
      if (typeof label === 'string' && label.startsWith('cp-overlay')) stale.push(ref);
    }
    if (stale.length > 0) {
      const cleanup = plugin.build();
      for (const ref of stale) cleanup.delete(ref);
      try { await cleanup.commit(); } catch (e) { console.warn('[pocket-overlay] cleanup failed', e); }
    }

    const pdbIdCol = pocketAtoms.col('pdb_id');
    if (!pdbIdCol) {
      console.warn('[pocket-overlay] pocket_atoms has no pdb_id column; skipping per-PDB overlay');
      return;
    }
    const uniqueIds = new Set<string>();
    for (let i = 0; i < pocketAtoms.rowCount; i++) {
      const id = String(pdbIdCol.get(i) ?? '').toUpperCase().trim();
      if (id) uniqueIds.add(id);
    }

    const knownRefs = new Set<string>();
    for (const [ref, cell] of plugin.state.data.cells) {
      const tname = cell?.obj?.type?.name;
      if (tname === 'Structure' || tname === 'Molecular Structure') knownRefs.add(ref);
    }

    const isSurface = rep === 'gaussian-surface';
    let serial = 1;
    let totalAtoms = 0;
    for (const pdbId of uniqueIds) {
      const overlay = pocketAtomsToOverlayBlock(pocketAtoms, {
        chain: 'P', serialStart: serial, pdbIdFilter: pdbId,
      });
      if (!overlay) continue;
      const atomLines = overlay.split('\n').filter((l) => l.startsWith('HETATM')).length;
      serial += atomLines;
      totalAtoms += atomLines;
      const fullBlock = `REMARK   Pocket overlay for ${pdbId} (chain P)\n` + overlay + '\nEND\n';
      try {
        const data = await plugin.builders.data.rawData(
          {data: fullBlock, label: `cp-overlay-${pdbId}`});
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);
        const component = await plugin.builders.structure.tryCreateComponentStatic(
          structure, 'ligand', {label: `cp-overlay-component-${pdbId}`});
        if (component) {
          await plugin.builders.structure.representation.addRepresentation(component, {
            type: isSurface ? 'gaussian-surface' : 'spacefill',
            typeParams: isSurface ? {alpha: 0.55} : {sizeFactor, alpha: 1.0},
            color: 'element-symbol',
          });
        }
        for (const [ref, cell] of plugin.state.data.cells) {
          const tname = cell?.obj?.type?.name;
          if ((tname === 'Structure' || tname === 'Molecular Structure') &&
              !knownRefs.has(ref) && !this._alignedInteractionRefs.has(pdbId)) {
            // Pocket overlay = one structure per PDB → single-element array.
            this._alignedInteractionRefs.set(pdbId, [ref]);
            knownRefs.add(ref);
            break;
          }
        }
      } catch (e) {
        console.error(`[pocket-overlay] failed for ${pdbId}:`, e);
      }
    }
    console.log(`[pocket-overlay] tracked ${this._alignedInteractionRefs.size} per-PDB ` +
      `pocket overlays, ${totalAtoms} atoms total (rep: ${rep})`);
  }

  private async addInteractionOverlay(
    features: DG.DataFrame, sizeFactor = 2.5, opts: {chain?: string} = {},
  ): Promise<void> {
    if (!this.viewer) return;
    const plugin = (this.viewer as any).viewer?.plugin;
    if (!plugin) return;

    // Build the overlay PDB block. ligandFeaturesToOverlayBlock already maps
    // family→element (Donor→N, Acceptor→O, Aromatic→S, Hydrophobic→C,
    // Positive→Na, Negative→Cl, Halogen→Br); Mol*'s CPK element coloring
    // gives us a distinct color per family without any custom theme.
    const chain = opts.chain ?? 'F';
    const overlay = ligandFeaturesToOverlayBlock(features, {chain, serialStart: 1});
    if (!overlay) return;
    const fullBlock = `REMARK   Consensus interaction overlay (chain ${chain})\n` +
      overlay + '\nEND\n';

    // Delete any previous overlay structures (data cells with our marker).
    // plugin.build().delete(ref) is the idiomatic Mol* cleanup.
    const stale: string[] = [];
    for (const [ref, cell] of plugin.state.data.cells) {
      const label = cell?.transform?.params?.label ?? cell?.obj?.label;
      if (typeof label === 'string' && label.startsWith('cp-overlay'))
        stale.push(ref);
    }
    if (stale.length > 0) {
      const cleanup = plugin.build();
      for (const ref of stale) cleanup.delete(ref);
      try { await cleanup.commit(); } catch (e) { console.warn('[overlay] cleanup failed', e); }
    }

    // Per-row sizing path: when the features dataframe carries a
    // `cluster_radius_a` column (Stage 5a output), each row gets its own
    // Mol* structure + spacefill rep so the rendered sphere diameter can
    // approximate the cluster's actual spatial spread. This is the
    // "consensus pharmacophore" mode — the sphere visually covers the
    // positions of the features that voted for that cluster, instead of
    // being a constant blob bigger than the pocket.
    const radCol = features.col('cluster_radius_a');
    if (radCol != null) {
      await this.addPerRowOverlay(features, chain);
      return;
    }

    // Constant-size path (Stage 4 chain-F overlay): single Mol* structure
    // holding all atoms, one spacefill rep with one global sizeFactor. Lighter
    // on the state tree than the per-row path; sufficient when every point
    // represents the same kind of thing (a single interaction site).
    try {
      const data = await plugin.builders.data.rawData(
        {data: fullBlock, label: 'cp-overlay-data'});
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      // 'ligand' static component selects HETATMs — our overlay atoms are
      // all HETATM by construction, so this picks them up wholesale.
      const component = await plugin.builders.structure.tryCreateComponentStatic(
        structure, 'ligand', {label: 'cp-overlay-component'});
      if (component) {
        await plugin.builders.structure.representation.addRepresentation(component, {
          type: 'spacefill',
          typeParams: {sizeFactor, alpha: 1.0},
          color: 'element-symbol',
        });
      } else {
        console.warn('[overlay] tryCreateComponentStatic(ligand) returned null; ' +
          'no representation added — chain F atoms may be invisible.');
      }
      console.log(`[overlay] added ${features.rowCount} chain-${chain} atoms ` +
        `as separate structure with spacefill sizeFactor=${sizeFactor}`);
    } catch (e) {
      console.error('[overlay] failed to add Mol* overlay structure:', e);
    }
  }

  /**
   * Per-row chain-P consensus overlay. For each row of `features` (which
   * MUST carry `cluster_radius_a`), build a single-atom PDB, load it as its
   * own Mol* structure, and apply a spacefill rep whose sizeFactor is tuned
   * so the rendered sphere's WORLD-SPACE radius approximates the cluster's
   * max-distance-from-centroid value in Å. The rep is half-transparent
   * (alpha 0.5) so spheres overlap visibly instead of occluding each other
   * and the underlying drug ligands.
   *
   * Mol*'s spacefill rep multiplies `sizeFactor` by the element's vdW radius
   * (Å). For our family→element mapping (Donor→N 1.55, Acceptor→O 1.52,
   * Aromatic→S 1.8, Hydrophobic→C 1.7, Positive→Na 2.27, Negative→Cl 1.75,
   * Halogen→Br 1.85), we precompute sizeFactor = clusterRadiusA / vdw so the
   * rendered radius lands at clusterRadiusA Å. Falls back to vdW 1.7 (C) for
   * unknown elements.
   */
  private async addPerRowOverlay(
    features: DG.DataFrame, chain: string,
  ): Promise<void> {
    if (!this.viewer) return;
    const plugin = (this.viewer as any).viewer?.plugin;
    if (!plugin) return;
    const famCol = features.col('family');
    const xCol = features.col('x');
    const yCol = features.col('y');
    const zCol = features.col('z');
    const radCol = features.col('cluster_radius_a');
    if (!famCol || !xCol || !yCol || !zCol || !radCol) return;

    // Approximate vdW radii (Å) used by Mol*'s default spacefill size theme.
    const VDW: Record<string, number> = {
      N: 1.55, O: 1.52, S: 1.8, C: 1.7,
      NA: 2.27, CL: 1.75, BR: 1.85, I: 1.98,
    };
    let nDrawn = 0;
    for (let i = 0; i < features.rowCount; i++) {
      const famName = String(famCol.get(i) ?? '');
      const fam = resolveFamily(famName);
      const x = Number(xCol.get(i));
      const y = Number(yCol.get(i));
      const z = Number(zCol.get(i));
      const clusterRadius = Number(radCol.get(i));
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z))
        continue;
      const vdw = VDW[fam.element.toUpperCase()] ?? 1.7;
      // sizeFactor lower clamp 1.0 so the smallest sphere is at least vdw Å
      // radius (~1.5-2 Å, clearly visible on a 50 Å protein). Lower clamps
      // produced spheres too tiny to spot inside the kinase pocket when the
      // protein cartoon failed to dim (race condition with Mol*'s async cell
      // refresh). Upper clamp 2.5 keeps the largest spheres from filling the
      // pocket. Dynamic range now ~1.5–4 Å rendered radius across the
      // clamped clusterRadius input range of 0.8–3.5 Å.
      const sizeFactor = Math.max(1.0, Math.min(2.5,
        (Number.isFinite(clusterRadius) ? clusterRadius : 1.5) / vdw));

      // Single-atom PDB. resName/element come from FAMILY_MAP so the rep's
      // CPK coloring picks up the right family hue.
      const atomNameField = fam.element.length === 2
        ? fam.element.padEnd(4)
        : ' ' + fam.element.padEnd(3);
      const f = (v: number) => v.toFixed(3).padStart(8);
      const block =
        `REMARK   Consensus point ${i + 1} (${fam.name}, radius ${clusterRadius.toFixed(2)} Å)\n` +
        `HETATM    1 ${atomNameField} ${fam.resName.padEnd(3)} ${chain}   1    ` +
        f(x) + f(y) + f(z) +
        '  1.00  1.00          ' + fam.element.padStart(2) + '\nEND\n';

      try {
        const data = await plugin.builders.data.rawData(
          {data: block, label: `cp-overlay-row-${i}-data`});
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);
        const component = await plugin.builders.structure.tryCreateComponentStatic(
          structure, 'ligand', {label: `cp-overlay-row-${i}-component`});
        if (component) {
          await plugin.builders.structure.representation.addRepresentation(component, {
            type: 'spacefill',
            typeParams: {sizeFactor, alpha: 0.5},
            color: 'element-symbol',
          });
          nDrawn += 1;
        }
      } catch (e) {
        console.error(`[overlay] per-row overlay failed for row ${i}:`, e);
      }
    }
    console.log(`[overlay] added ${nDrawn} chain-${chain} consensus spheres ` +
      `(sized by cluster_radius_a, alpha=0.5)`);
  }

  /**
   * Build a per-PDB interaction-detail DG.Widget content host and dock it on
   * the RIGHT, between the main grid and the Mol* viewer. The host updates
   * when the user clicks a different row in the summary grid. Lazy-creates
   * the dock node on first call; subsequent calls reuse it. Callers must
   * pass the raw `features` (one row per interaction) and `summary` (one row
   * per PDB) — only `summary.onCurrentRowChanged` drives the panel updates.
   */
  private setupPdbDetailPanel(features: DG.DataFrame, summary: DG.DataFrame): void {
    // The dock node + host are pre-created in init() so Mol*'s width stays
    // constant across stages. If init() hasn't run yet, defensively create
    // them here. In the wizard path, _externalDetailHost was assigned to
    // pdbDetailHost via attachDetailHost(), and there is no dock node.
    if (!this.pdbDetailHost) {
      this.pdbDetailHost = ui.div([], 'cp-pdb-detail-panel');
    }
    if (!this.pdbDetailNode && this.tableView && !this._externalDetailHost) {
      this.pdbDetailNode = this.tableView.dockManager.dock(
        this.pdbDetailHost, DG.DOCK_TYPE.RIGHT, this.viewerDockNode, 'PDB detail', 0.4);
    }

    // Re-subscribe to onCurrentRowChanged on the NEW summary DataFrame.
    if (this.pdbDetailSub) {
      try { this.pdbDetailSub.unsubscribe(); } catch (_e) { /* */ }
      this.pdbDetailSub = undefined;
    }
    const update = () => this.renderPdbDetail(features, summary);
    this.pdbDetailSub = summary.onCurrentRowChanged.subscribe(update);
    // Initial paint — Datagrok keeps `currentRowIdx` at -1 by default, so
    // pre-select row 0 if any rows exist so the user immediately sees detail
    // for the first PDB rather than an empty panel.
    if (summary.rowCount > 0 && summary.currentRowIdx < 0)
      summary.currentRowIdx = 0;
    update();
  }

  /**
   * Tear down the per-PDB detail panel. Called before non-Features previews
   * (Pocket / Consensus) and on app abort so the panel doesn't linger with
   * stale content.
   */
  private closePdbDetailPanel(): void {
    // Unsubscribe from the previous summary's onCurrentRowChanged stream so
    // the old summary DataFrame doesn't pin live update callbacks. The dock
    // node itself stays (pre-created in init for layout stability) — we just
    // restore the empty-state hint.
    if (this.pdbDetailSub) {
      try { this.pdbDetailSub.unsubscribe(); } catch (_e) { /* */ }
      this.pdbDetailSub = undefined;
    }
    if (this.pdbDetailHost) {
      ui.empty(this.pdbDetailHost);
      this.pdbDetailHost.append(ui.divText(
        'Per-PDB interaction detail — click 4. Features to populate.',
        'cp-pdb-detail-empty'));
    }
  }

  /** Zoom the Mol* camera in on a single interaction's position — the same
   *  effect as clicking that atom in Mol*. Reachable purely through the
   *  runtime plugin (`plugin.managers.camera`), so no Mol* internals are
   *  bundled (hard rule #4). Pair with `resetCameraView()` to zoom back out
   *  when the user de-selects the interaction. */
  async focusInteraction(x: number, y: number, z: number): Promise<void> {
    if (!this.viewer) return;
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) return;
    const camera = (this.viewer as any).viewer?.plugin?.managers?.camera;
    if (!camera) return;
    try {
      // focusSphere frames a sphere of radius max(radius + extraRadius,
      // minRadius). radius 0 + a tight minRadius gives a close-up centred on
      // the feature (~9 Å of context around the contact point).
      camera.focusSphere({center: [x, y, z], radius: 0},
        {durationMs: 350, minRadius: 9, extraRadius: 2});
    } catch (e) {
      console.warn('[focus] focusSphere failed:', e);
    }
  }

  /** Zoom the Mol* camera back out to frame the whole scene — the "de-select"
   *  / standard-view counterpart to `focusInteraction()`. */
  async resetCameraView(): Promise<void> {
    if (!this.viewer) return;
    const camera = (this.viewer as any).viewer?.plugin?.managers?.camera;
    if (!camera) return;
    try {
      camera.reset(undefined, 350);
    } catch (e) {
      console.warn('[focus] camera reset failed:', e);
    }
  }

  /** Show or hide ALL of one PDB's structures (protein + its per-family
   *  interaction overlays) in Mol* without a viewer rebuild. Wired to the
   *  Step 4 "Use" checkbox: un-ticking a PDB removes it from the 3D view as
   *  well as from the Step 5 consensus. No-op when the fast-path refs aren't
   *  tracked yet. Uses the same single-transaction hierarchy.toggleVisibility
   *  path as `isolatePdb`. */
  async setPdbHidden(pdbId: string, hidden: boolean): Promise<void> {
    if (!this.viewer) return;
    const plugin = (this.viewer as any).viewer?.plugin;
    const hierarchy = plugin?.managers?.structure?.hierarchy;
    if (!plugin || !hierarchy?.current) return;
    const up = pdbId.toUpperCase();
    try {
      const structByRootRef = new Map<string, any>();
      for (const s of hierarchy.current.structures)
        structByRootRef.set(s.cell.transform.ref, s);
      const targets: any[] = [];
      const protRef = this._alignedStructureRefs.get(up);
      if (protRef && structByRootRef.has(protRef))
        targets.push(structByRootRef.get(protRef));
      for (const ref of (this._alignedInteractionRefs.get(up) ?? []))
        if (structByRootRef.has(ref)) targets.push(structByRootRef.get(ref));
      if (targets.length > 0)
        hierarchy.toggleVisibility(targets, hidden ? 'hide' : 'show');
    } catch (e) {
      console.warn('[setPdbHidden] failed:', e);
    }
  }

  /** Hide every PDB currently excluded via the Step 4 "Use" checkboxes. Called
   *  after the Step 4 viewer is (re)built so un-ticked PDBs don't reappear when
   *  the user navigates away and back (a rebuild starts every structure
   *  visible). */
  private async applyPdbExclusionsToViewer(): Promise<void> {
    const isCons = this._consensusSelectionHooks?.isPdbExcluded;
    const isInput = this._consensusSelectionHooks?.isInputPdbExcluded;
    if (!isCons && !isInput) return;
    for (const id of this._alignedStructureRefs.keys()) {
      const hide = (isCons && isCons(id)) || (isInput && isInput(id));
      if (hide) await this.setPdbHidden(id, true);
    }
  }

  /** Clear all overlay ref tracking (per-PDB + per-interaction). Called wherever
   *  the interaction overlay is torn down or rebuilt so stale refs never leak
   *  into isolatePdb / setInteractionHidden. Loop form so a global
   *  find-replace of the per-map .clear() calls never rewrites this body. */
  private clearOverlayRefs(): void {
    for (const m of [this._alignedInteractionRefs, this._interactionRefByKey,
      this._interactionKeyByRef]) m.clear();
  }

  /** Show or hide ONE interaction's sphere in Mol* (no rebuild). Wired to the
   *  per-interaction "use" checkbox in the Step 4 detail panel: un-ticking a
   *  contact removes its sphere from the 3D view as well as from the consensus.
   *  No-op when the interaction's structure ref isn't tracked. */
  async setInteractionHidden(key: string, hidden: boolean): Promise<void> {
    if (!this.viewer) return;
    const plugin = (this.viewer as any).viewer?.plugin;
    const hierarchy = plugin?.managers?.structure?.hierarchy;
    if (!plugin || !hierarchy?.current) return;
    const ref = this._interactionRefByKey.get(key);
    if (!ref) return;
    try {
      let target: any = null;
      for (const s of hierarchy.current.structures)
        if (s.cell.transform.ref === ref) { target = s; break; }
      if (target) hierarchy.toggleVisibility([target], hidden ? 'hide' : 'show');
    } catch (e) {
      console.warn('[setInteractionHidden] failed:', e);
    }
  }

  /** Hide every interaction the user has un-ticked. Called after the Step 4
   *  overlay is (re)built so excluded contacts don't reappear on rebuild. */
  private async applyInteractionExclusionsToViewer(): Promise<void> {
    const isExcluded = this._consensusSelectionHooks?.isInteractionExcluded;
    if (!isExcluded) return;
    for (const key of this._interactionRefByKey.keys())
      if (isExcluded(key)) await this.setInteractionHidden(key, true);
  }

  /** Render the per-PDB detail panel content for the currently-selected
   *  summary row. Pure DOM rewrite — replaces the existing host children. */
  private renderPdbDetail(features: DG.DataFrame, summary: DG.DataFrame): void {
    if (!this.pdbDetailHost) return;
    ui.empty(this.pdbDetailHost);
    const idx = summary.currentRowIdx;
    if (idx < 0 || idx >= summary.rowCount) {
      this.pdbDetailHost.append(ui.divText(
        'Select a PDB row to see per-interaction detail.', 'cp-pdb-detail-empty'));
      return;
    }
    const pdbId = String(summary.col('pdb_id')?.get(idx) ?? '');
    const compId = String(summary.col('ligand_comp_id')?.get(idx) ?? '');
    const nTotal = Number(summary.col('n_total')?.get(idx) ?? 0);
    const hasWarn = Boolean(summary.col('has_skip_warning')?.get(idx));

    // Header
    const header = ui.div([], 'cp-pdb-detail-header');
    header.style.fontWeight = 'bold';
    header.style.marginBottom = '8px';
    header.style.display = 'flex';
    header.style.alignItems = 'center';
    header.style.gap = '6px';
    header.append(ui.divText(pdbId));
    header.append(ui.divText('·'));
    header.append(ui.divText(compId));
    header.append(ui.divText(`— ${nTotal} interaction${nTotal === 1 ? '' : 's'}`));
    if (hasWarn) {
      const warn = ui.divText('⚠', 'cp-pdb-detail-warn');
      warn.title = 'One or more residues in this PDB were skipped (modified residue / parse error).';
      warn.style.color = '#FFA000';
      header.append(warn);
    }
    this.pdbDetailHost.append(header);

    // Family chips for this PDB
    const chipRow = ui.divH([], 'cp-pdb-detail-chips');
    chipRow.style.flexWrap = 'wrap';
    chipRow.style.gap = '4px';
    chipRow.style.marginBottom = '8px';
    const familyCounts: Array<[string, number]> = [];
    for (const code of FAMILY_CODES) {
      const colName = SUMMARY_FAMILY_COL[code];
      const n = Number(summary.col(colName)?.get(idx) ?? 0);
      if (n > 0) familyCounts.push([code, n]);
    }
    for (const [code, n] of familyCounts) {
      const fam = FAMILY_MAP[code];
      const chip = ui.divText(`${code} ×${n}`, 'cp-pdb-detail-chip');
      chip.style.background = fam.hexColor;
      chip.style.color = (code === 'H') ? '#222' : '#fff';
      chip.style.padding = '1px 6px';
      chip.style.borderRadius = '10px';
      chip.style.fontSize = '11px';
      chip.title = fam.name;
      chipRow.append(chip);
    }
    this.pdbDetailHost.append(chipRow);

    // Detail table — one mini-row per interaction.
    const list = ui.div([], 'cp-pdb-detail-list');
    list.style.display = 'flex';
    list.style.flexDirection = 'column';
    list.style.rowGap = '2px';

    const famCol = features.col('family');
    const intCol = features.col('interaction_type');
    const resCol = features.col('residue');
    const distCol = features.col('distance');
    const pdbCol = features.col('pdb_id');
    const skipCol = features.col('skip_reason');
    const xCol = features.col('x');
    const yCol = features.col('y');
    const zCol = features.col('z');
    // Highlighted interaction (clicked) — its row stays marked across renders.
    const highlightedKey = this._highlightedInteractionKey;

    let nShown = 0;
    for (let i = 0; i < features.rowCount && nShown < 200; i++) {
      if (pdbCol && String(pdbCol.get(i)) !== pdbId) continue;
      if (skipCol && String(skipCol.get(i) ?? '').trim() !== '') continue;
      const famName = String(famCol?.get(i) ?? '');
      const fam = resolveFamily(famName);
      const key = consensusFeatureKey(features, i);
      // One clickable row per interaction (same 5-col grid template as before:
      // checkbox · family dot · type · residue · distance), but a single
      // element so the WHOLE row highlights its Mol* sphere on click.
      const row = ui.div();
      row.classList.add('cp-detail-int-row');
      row.style.display = 'grid';
      // Fixed widths on the last two columns so independent per-row grids stay
      // aligned (residue / distance). checkbox=auto, dot=12px, type=1fr.
      row.style.gridTemplateColumns = 'auto 12px 1fr 88px 58px';
      row.style.columnGap = '6px';
      row.style.alignItems = 'center';
      row.style.cursor = 'pointer';
      row.style.borderRadius = '3px';
      row.style.padding = '1px 2px';
      row.title = 'Click to zoom to this interaction in the 3D viewer. ' +
        'Click again to zoom back out.';
      if (key === highlightedKey) row.style.background = 'rgba(33,150,243,0.18)';

      // Per-interaction "use" checkbox (tier 2) — exclude THIS contact from the
      // Step 5 consensus without dropping the whole PDB. stopPropagation so
      // ticking it doesn't also fire the row's 3D-highlight click.
      const selHooks = this._consensusSelectionHooks;
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.checked = selHooks ? !selHooks.isInteractionExcluded(key) : true;
      cb.disabled = !selHooks;
      cb.title = 'Include this interaction in the Step 5 consensus.';
      cb.style.cursor = 'pointer';
      cb.addEventListener('click', (e) => e.stopPropagation());
      cb.addEventListener('change', () => {
        selHooks?.toggleInteraction(key);
        // Un-ticking also hides this interaction's sphere in Mol* (re-ticking
        // shows it) — mirrors the per-PDB "Use" box.
        void this.setInteractionHidden(key, !cb.checked);
      });
      row.append(cb);
      const dot = ui.div();
      dot.style.width = '8px';
      dot.style.height = '8px';
      dot.style.borderRadius = '50%';
      dot.style.background = fam.hexColor;
      row.append(dot);
      row.append(ui.divText(String(intCol?.get(i) ?? '')));
      row.append(ui.divText(String(resCol?.get(i) ?? ''), 'cp-pdb-detail-res'));
      const d = Number(distCol?.get(i));
      row.append(ui.divText(Number.isFinite(d) ? d.toFixed(2) + ' Å' : ''));

      // Click anywhere on the row → zoom the 3D viewer in on this interaction
      // (like clicking the atom in Mol*). Click the SAME row again to toggle
      // off — de-select and zoom back out to the standard whole-scene view.
      const ix = Number(xCol?.get(i));
      const iy = Number(yCol?.get(i));
      const iz = Number(zCol?.get(i));
      row.addEventListener('click', () => {
        const wasFocused = this._highlightedInteractionKey === key;
        list.querySelectorAll<HTMLElement>('.cp-detail-int-row')
          .forEach((el) => { el.style.background = ''; });
        if (wasFocused) {
          this._highlightedInteractionKey = undefined;
          void this.resetCameraView();
        } else {
          this._highlightedInteractionKey = key;
          row.style.background = 'rgba(33,150,243,0.18)';
          void this.focusInteraction(ix, iy, iz);
        }
      });
      list.append(row);
      nShown += 1;
    }
    // Surface skip_reason diagnostics for THIS PDB so a 0-interaction result
    // explains itself (ProLIF conversion / bond-perception failure, ligand not
    // found, pdbfixer drop, etc.) instead of a bare "(no interactions)". This
    // is the user-facing "why did 5HG8 break?" answer.
    const skipReasons: string[] = [];
    if (skipCol && pdbCol) {
      const seen = new Set<string>();
      for (let i = 0; i < features.rowCount; i++) {
        if (String(pdbCol.get(i)) !== pdbId) continue;
        const r = String(skipCol.get(i) ?? '').trim();
        if (r && !seen.has(r)) { seen.add(r); skipReasons.push(r); }
      }
    }
    if (nShown === 0 && skipReasons.length === 0)
      list.append(ui.divText('(no interactions for this PDB)',
        'cp-pdb-detail-empty'));
    this.pdbDetailHost.append(list);
    if (skipReasons.length > 0) {
      const diag = ui.div([], 'cp-pdb-detail-skip');
      diag.style.marginTop = '8px';
      diag.style.fontSize = '11px';
      diag.style.color = '#B26A00';
      diag.append(ui.divText(nShown === 0
        ? '⚠ No interactions — ProLIF could not process this ligand:'
        : '⚠ Some ligand/residue rows were skipped:'));
      for (const r of skipReasons) {
        const line = ui.divText('• ' + r);
        line.style.fontFamily = 'monospace';
        line.style.whiteSpace = 'normal';
        diag.append(line);
      }
      this.pdbDetailHost.append(diag);
    }
  }

  // -------------------------------------------------------------------------
  // Preview cache + per-stage ensures
  // -------------------------------------------------------------------------

  /** Fingerprint of the inputs + every user-visible option, deterministic across
   *  reorderings. The wizard calls this to compare the live form state against
   *  `lastRunFingerprint()` — when they match, displayed stage results are still
   *  valid and the stale banner should clear. */
  previewFingerprint(pdbIds: string[], opts: PipelineOptions): string {
    return this._previewFingerprint(pdbIds, opts);
  }

  /** Fingerprint of the inputs + options that produced the currently-cached
   *  stage results. Returns null when the cache is empty (no run yet, or a
   *  full reset). */
  lastRunFingerprint(): string | null {
    return this.previewCache?.fingerprint ?? null;
  }

  private _previewFingerprint(pdbIds: string[], opts: PipelineOptions): string {
    return JSON.stringify({
      pdbIds: [...pdbIds].sort(),
      maxResolution: opts.maxResolution,
      allowXray: opts.allowXray,
      allowNmr: opts.allowNmr,
      allowCryoEm: opts.allowCryoEm,
      allowAlphaFold: opts.allowAlphaFold,
      minLigandMw: opts.minLigandMw,
      pocketMethod: opts.pocketMethod,
      pocketRadius: opts.pocketRadius,
      // refPdbId is the Stage 2a alignment template. Changing it produces
      // different RMSDs and a different aligned DataFrame, so it must be
      // part of the cache key — otherwise switching reference in the UI
      // would silently return the previous run's aligned result.
      refPdbId: (opts.refPdbId ?? '').trim().toUpperCase(),
      // Step 1 "Use" deselections — drop PDBs before Stage 2. Different
      // exclusion set = different aligned/pocket/features/consensus, so the
      // whole downstream cache must invalidate when the user toggles a box.
      excludedInputPdbs: [...(opts.excludedInputPdbs ?? [])]
        .map((s) => s.toUpperCase()).sort().join(','),
      // Dropped local files have no RCSB identity — key the cache on their
      // content so re-dropping an edited file invalidates downstream stages.
      localStore: localFingerprintFor(pdbIds),
    });
  }

  private ensurePreviewCache(fingerprint: string): PreviewCache {
    if (!this.previewCache || this.previewCache.fingerprint !== fingerprint)
      this.previewCache = {fingerprint};
    return this.previewCache;
  }

  /** Run Stage 1 (or use cached) and publish pdb_qc to the table view. */
  private async ensurePdbQc(
    cache: PreviewCache, pdbIds: string[], opts: PipelineOptions, pi: DG.ProgressIndicator,
  ): Promise<DG.DataFrame> {
    if (cache.pdbQc) return cache.pdbQc;
    this.setStageStatus('stage1', 'running');
    pi.update(10, 'Stage 1: enriching PDB metadata...');
    // enrichPdbList synthesizes QC rows for user-dropped local ids (no RCSB
    // round-trip) alongside the RCSB enrichment — one schema-stable build.
    const stage1 = await enrichPdbList(pdbIds, opts);
    const pdbQc = stage1.accepted;
    this._droppedPdbs = stage1.dropped;
    if (pdbQc.rowCount === 0)
      throw new Error('No usable PDB entries after QC filtering. ' +
        'Check the PDB IDs, loosen the resolution cap, or disable the X-ray filter.');
    pdbQc.name = 'pdb_qc';
    this.setStageStatus('stage1', 'done');
    cache.pdbQc = pdbQc;
    return pdbQc;
  }

  /** TS-side parallel fetch of raw PDB blocks. Caches per fingerprint. */
  private async ensureRawBlocks(
    cache: PreviewCache, pdbQc: DG.DataFrame, pi: DG.ProgressIndicator,
  ): Promise<Map<string, string>> {
    if (cache.rawBlocks) return cache.rawBlocks;
    const pdbIdCol = pdbQc.col('pdb_id');
    const seen = new Set<string>();
    const uniqueIds: string[] = [];
    if (pdbIdCol) {
      for (let i = 0; i < pdbQc.rowCount; i++) {
        const id = String(pdbIdCol.get(i));
        if (!seen.has(id)) { seen.add(id); uniqueIds.push(id); }
      }
    }
    pi.update(40, `Fetching ${uniqueIds.length} PDB file(s)...`);
    const blocks = await Promise.all(uniqueIds.map((id) => fetchPdbBlock(id)));
    cache.rawBlocks = new Map(uniqueIds.map((id, i) => [id, blocks[i]]));
    return cache.rawBlocks;
  }

  /** Run Stage 2a (or use cached) and publish aligned_structures. Filters
   *  cache.pdbQc by the Step 1 "Use" deselections before handing rows to
   *  Stage 2a — excluded PDBs never enter alignment / pocket / features /
   *  consensus. The exclusion set is in the fingerprint, so this filtered
   *  view is consistent with the cache slot. */
  private async ensureAligned(
    cache: PreviewCache, opts: PipelineOptions, pi: DG.ProgressIndicator,
  ): Promise<DG.DataFrame> {
    if (cache.aligned) return cache.aligned;
    if (!cache.pdbQc) throw new Error('ensureAligned: pdb_qc missing from cache');
    this.setStageStatus('stage2a', 'running');
    pi.update(50, 'Stage 2a: aligning structures (Kabsch)...');
    const stage2Input = filterPdbQcByInputExclusion(cache.pdbQc, opts);
    if (stage2Input.rowCount === 0)
      throw new Error('All accepted PDBs are un-ticked in Step 1’s "Use" column. ' +
        'Re-tick at least one PDB to continue.');
    const aligned = await runStage2aAlign(stage2Input, opts);
    aligned.name = 'aligned_structures (preview)';
    this.setStageStatus('stage2a', 'done');
    cache.aligned = aligned;
    return aligned;
  }

  /** Run Stage 3 (or use cached) and publish pocket_atoms. */
  private async ensurePocketAtoms(
    cache: PreviewCache, opts: PipelineOptions, pi: DG.ProgressIndicator,
  ): Promise<DG.DataFrame> {
    if (cache.pocketAtoms) return cache.pocketAtoms;
    if (!cache.aligned) throw new Error('ensurePocketAtoms: aligned missing from cache');
    this.setStageStatus('stage3', 'running');
    pi.update(70, 'Stage 3: isolating pocket...');
    const pocketAtoms = await runStage3IsolatePocket(cache.aligned, ['all'], opts);
    pocketAtoms.name = 'pocket_atoms (preview)';
    this.setStageStatus('stage3', 'done');
    cache.pocketAtoms = pocketAtoms;
    return pocketAtoms;
  }

  /** Run Stage 2b (or use cached). Refines transforms over the pocket Cα subset. */
  private async ensureAlignedV2(
    cache: PreviewCache, pi: DG.ProgressIndicator,
  ): Promise<DG.DataFrame> {
    if (cache.alignedV2) return cache.alignedV2;
    if (!cache.aligned) throw new Error('ensureAlignedV2: aligned missing from cache');
    if (!cache.pocketAtoms) throw new Error('ensureAlignedV2: pocketAtoms missing from cache');
    this.setStageStatus('stage2b', 'running');
    pi.update(78, 'Stage 2b: re-aligning on pocket Calpha...');
    // Preview path: don't filter by cluster (the cluster picker only runs in the
    // full pipeline; previews bypass it intentionally).
    const all: string[] = [];
    const v2 = await runStage2bAlignPocket(cache.aligned, cache.pocketAtoms, all);
    v2.name = 'aligned_structures (pass 2)';
    this.setStageStatus('stage2b', 'done');
    cache.alignedV2 = v2;
    return v2;
  }

  /** Run Stage 4 (or use cached) — ProLIF interaction extraction. */
  private async ensureLigandFeatures(
    cache: PreviewCache, pi: DG.ProgressIndicator,
  ): Promise<DG.DataFrame> {
    if (cache.ligandFeatures) return cache.ligandFeatures;
    if (!cache.alignedV2 || !cache.pocketAtoms)
      throw new Error('ensureLigandFeatures: alignedV2/pocketAtoms missing');
    this.setStageStatus('stage4', 'running');
    pi.update(85, 'Stage 4: detecting protein-ligand interactions (ProLIF, ~15-30s/PDB)...');
    const all: string[] = [];
    const features = await runStage4ExtractFeatures(cache.alignedV2, cache.pocketAtoms, all);
    features.name = 'ligand_features (preview)';
    this.setStageStatus('stage4', 'done');
    cache.ligandFeatures = features;
    return features;
  }

  /** Run Stage 5a (or use cached) — k-means consensus. */
  private async ensureConsensusModel(
    cache: PreviewCache, pi: DG.ProgressIndicator, opts: PipelineOptions,
  ): Promise<DG.DataFrame> {
    // Consensus-only signature = the k-means knobs + the user's PDB/interaction
    // exclusions. If unchanged and we have a cached model, reuse it. Otherwise
    // recompute from the FILTERED features. This makes a knob change OR a PDB
    // toggle re-run ONLY Stage 5a — Stages 1-4 stay cached (they're keyed by
    // the global fingerprint, which deliberately omits these consensus inputs).
    const sig = consensusSignature(opts);
    if (cache.consensusModel && cache.consensusSignature === sig)
      return cache.consensusModel;
    if (!cache.ligandFeatures) throw new Error('ensureConsensusModel: ligandFeatures missing');
    this.setStageStatus('stage5a', 'running');
    pi.update(92, 'Stage 5a: k-means consensus...');
    const features = filterFeaturesForConsensus(cache.ligandFeatures, opts);
    const consensus = await runStage5aConsensusKmeans(features, opts);
    consensus.name = 'consensus_model (preview)';
    this.setStageStatus('stage5a', 'done');
    cache.consensusModel = consensus;
    cache.consensusSignature = sig;
    return consensus;
  }

  /**
   * Preview 1/3 — fetch PDBs from RCSB (no alignment). Structures appear in their
   * native crystallographic frames; useful for sanity-checking the inputs.
   */
  async previewFetch(
    pdbIds: string[], options: PipelineOptions, fetchOpts: {render3d?: boolean} = {},
  ): Promise<void> {
    // tableView may be undefined in the wizard path — setActiveDataFrame()
    // guards the assignment.
    this.resetStageStatuses();
    this.closePdbDetailPanel();
    const cache = this.ensurePreviewCache(this.previewFingerprint(pdbIds, options));
    const pi = DG.TaskBarProgressIndicator.create('Fetch PDBs: starting...');
    try {
      const pdbQc = await this.ensurePdbQc(cache, pdbIds, options, pi);
      this.setActiveDataFrame(pdbQc);

      // Wizard path passes render3d:false — Step 1 (Fetch) is about QC, not
      // 3D. Showing the structures in their native crystallographic frames
      // here is actively confusing (each PDB sits at its own origin, so they
      // scatter across the canvas). The 3D story starts at Step 2 (Align).
      // Skipping the render also avoids fetching raw PDB blocks until they're
      // actually needed (Stage 2a fetches them anyway).
      if (fetchOpts.render3d === false) {
        pi.update(100, 'Done.');
        return;
      }

      // QC has succeeded and the Accepted/Dropped table is published. The 3D
      // preview of the (unaligned) structures is BEST-EFFORT: rendering many
      // structures in their native frames can fail inside Mol*/BSV, and such a
      // failure must NOT fail the whole Fetch step — QC is the point of Step 1,
      // and the real 3D story is at Align.
      try {
        const rawBlocks = await this.ensureRawBlocks(cache, pdbQc, pi);
        this.setMolstarCaption(
          'Fetched PDBs shown in their native crystallographic frames (no alignment yet).');

        // Unique pdb_ids preserved in pdb_qc first-occurrence order
        const labels: string[] = [];
        const blocks: string[] = [];
        const seen = new Set<string>();
        const pdbIdCol = pdbQc.col('pdb_id')!;
        for (let i = 0; i < pdbQc.rowCount; i++) {
          const id = String(pdbIdCol.get(i));
          if (seen.has(id)) continue;
          seen.add(id);
          labels.push(id);
          blocks.push(stripCofactorsFromPdb(rawBlocks.get(id) ?? ''));
        }

        pi.update(90, 'Rendering in Mol*...');
        // Each PDB is a SEPARATE Mol* structure so the Accepted-PDBs row-click
        // can isolate one and the "Use" checkbox can hide one. Native crystal
        // frames will visually scatter; per-PDB selection still works.
        await this.renderAlignedStructuresMulti(blocks, labels);
        this._cachedAlignedBlocks = blocks;
        this._cachedAlignedLabels = labels;
        await this.setMolstarStyle({cartoonAlpha: 1.0, solventAlpha: 1.0});
        grok.shell.info(`Fetch: ${labels.length} PDB(s) shown in their native frames. ` +
          'Click "Align" to superimpose them.');
      } catch (renderErr: any) {
        console.warn('[previewFetch] 3D preview failed (QC succeeded):', renderErr);
        grok.shell.warning('Structures fetched and QC complete — the 3D preview here ' +
          `couldn't render (${renderErr?.message ?? renderErr}); they'll show at Align.`);
      }
      pi.update(100, 'Done.');
    } catch (e: any) {
      const failing = currentRunningStage(this.stageBadges) ?? 'preview';
      this.setStageStatus(failing, 'error');
      grok.shell.error(`Fetch failed: ${e?.message ?? e}`);
      throw e;
    } finally {
      pi.close();
    }
  }

  /**
   * Preview 2/3 — Stage 1 + Stage 2a. Shows the input structures superimposed
   * by the Kabsch transform from the Python script.
   */
  async previewAlign(pdbIds: string[], options: PipelineOptions): Promise<void> {
    // tableView may be undefined in the wizard path.
    this.resetStageStatuses();
    this.closePdbDetailPanel();
    this.setMolstarCaption(
      'Stage 2a: structures superimposed on global Cα via Kabsch. Each PDB now ' +
      'shares a common frame; the first one is the reference.');
    const cache = this.ensurePreviewCache(this.previewFingerprint(pdbIds, options));
    const pi = DG.TaskBarProgressIndicator.create('Align: starting...');
    try {
      await this.ensurePdbQc(cache, pdbIds, options, pi);
      const aligned = await this.ensureAligned(cache, options, pi);
      this.setActiveDataFrame(aligned);

      const pdbCol = aligned.col('original_pdb')!;
      const tCol = aligned.col('transform_4x4_json')!;
      const idCol = aligned.col('pdb_id')!;
      const labels: string[] = [];
      const transformedBlocks: string[] = [];
      for (let i = 0; i < aligned.rowCount; i++) {
        const block = String(pdbCol.get(i) ?? '');
        const transformJson = String(tCol.get(i) ?? '');
        labels.push(String(idCol.get(i) ?? `row-${i}`));
        transformedBlocks.push(block ? stripCofactorsFromPdb(applyPdbTransform(block, transformJson)) : '');
      }

      pi.update(90, 'Rendering in Mol*...');
      // Render each PDB as a SEPARATE Mol* structure (not a concat). This
      // lets Step 2's row clicks toggle per-PDB visibility via Mol*'s
      // native `isHidden` flag — instant, no rebuild. Uses applyPreset on
      // each trajectory to match BSV's default cartoon+components setup.
      await this.renderAlignedStructuresMulti(transformedBlocks, labels);
      // Keep the cached blocks too — fallback path in case the multi-
      // structure setup doesn't populate `_alignedStructureRefs` (then
      // isolatePdb can still rebuild from blocks).
      this._cachedAlignedBlocks = transformedBlocks;
      this._cachedAlignedLabels = labels;
      await this.setMolstarStyle({cartoonAlpha: 1.0, solventAlpha: 1.0});
      pi.update(100, 'Done.');
      const refId = String(aligned.col('ref_pdb_id')?.get(0) ?? '?');
      grok.shell.info(`Align: ${labels.length} structure(s) Kabsch-aligned to template ${refId}.`);
    } catch (e: any) {
      const failing = currentRunningStage(this.stageBadges) ?? 'preview';
      this.setStageStatus(failing, 'error');
      grok.shell.error(`Align failed: ${e?.message ?? e}`);
      throw e;
    } finally {
      pi.close();
    }
  }

  /**
   * Preview 3/3 — Stage 1 + Stage 2a + Stage 3. Shows the aligned structures with
   * the pocket Cα atoms overlaid on chain P (rendered as a distinctly-colored
   * point cloud around each ligand).
   */
  async previewPocket(pdbIds: string[], options: PipelineOptions): Promise<void> {
    // tableView may be undefined in the wizard path.
    this.resetStageStatuses();
    this.closePdbDetailPanel();
    this.setMolstarCaption(
      `Stage 3: pocket residues highlighted (within ${options.pocketRadius} Å of any drug ligand). ` +
      'Cartoon translucent; drug ligand kept as the binding-site seed.');
    const cache = this.ensurePreviewCache(this.previewFingerprint(pdbIds, options));
    const pi = DG.TaskBarProgressIndicator.create('Pocket: starting...');
    try {
      await this.ensurePdbQc(cache, pdbIds, options, pi);
      const aligned = await this.ensureAligned(cache, options, pi);
      // pocket_atoms is computed and cached (Stage 2b/4 will consume it later)
      // but we no longer surface it in the grid — the user-facing output of
      // "Pocket" is the visual, not the atom list.
      const pocketAtoms = await this.ensurePocketAtoms(cache, options, pi);
      this.setActiveDataFrame(aligned);

      const pdbCol = aligned.col('original_pdb')!;
      const tCol = aligned.col('transform_4x4_json')!;
      const idCol = aligned.col('pdb_id')!;
      const labels: string[] = [];
      const transformedBlocks: string[] = [];
      for (let i = 0; i < aligned.rowCount; i++) {
        const block = String(pdbCol.get(i) ?? '');
        const transformJson = String(tCol.get(i) ?? '');
        labels.push(String(idCol.get(i) ?? `row-${i}`));
        transformedBlocks.push(block ? stripCofactorsFromPdb(applyPdbTransform(block, transformJson)) : '');
      }

      // Render proteins as N SEPARATE structures (same multi-structure path
      // as Step 2) so Step 3 rows can isolate one protein + its pocket via
      // hierarchy.toggleVisibility — no rebuild. The pocket "dots" are then
      // drawn as per-PDB overlays (chain P, Cα + ligand-seed atoms) instead
      // of BSV's whole-viewer binding-site renderer (which can't be isolated
      // per structure). See the `molstar-via-bsv` skill for the trade-off.
      pi.update(90, 'Rendering proteins in Mol*...');
      await this.renderAlignedStructuresMulti(transformedBlocks, labels);
      this._cachedAlignedBlocks = transformedBlocks;
      this._cachedAlignedLabels = labels;
      // Translucent cartoon so the pocket overlay reads clearly on top.
      await this.setMolstarStyle({cartoonAlpha: 0.25, solventAlpha: 0.0});
      pi.update(96, 'Adding pocket overlays...');
      await this.addPocketOverlayPerPdb(pocketAtoms, options.pocketRep);
      // Re-assert the translucent cartoon after the overlay (same race
      // mitigation as previewFeatures).
      await this.setMolstarStyleWithRetry({cartoonAlpha: 0.25, solventAlpha: 0.0});
      pi.update(100, 'Done.');
      const nResidues = countPocketResidues(pocketAtoms);
      grok.shell.info(`Pocket: ${nResidues} pocket residue(s) detected ` +
        `(${options.pocketMethod}, ${options.pocketRadius.toFixed(1)} A). ` +
        `Rep: ${options.pocketRep}; cartoon translucent; solvent hidden. ` +
        'Click a PDB row to isolate it.');
    } catch (e: any) {
      const failing = currentRunningStage(this.stageBadges) ?? 'preview';
      this.setStageStatus(failing, 'error');
      grok.shell.error(`Pocket failed: ${e?.message ?? e}`);
      throw e;
    } finally {
      pi.close();
    }
  }

  /**
   * Preview 4/5 — Stages 1 → 2a → 3 → 2b → 4. Layers raw ligand_features (one
   * HETATM per SMARTS match, colored by family via element symbol) on top of
   * the Pocket view.
   */
  async previewFeatures(pdbIds: string[], options: PipelineOptions): Promise<void> {
    // tableView may be undefined in the wizard path.
    this.resetStageStatuses();
    this.setMolstarCaption(
      'Stage 4: chain-F spheres = ProLIF protein-ligand interactions (one per ' +
      'detected contact). CPK colors map to pharmacophore family (see legend).');
    const cache = this.ensurePreviewCache(this.previewFingerprint(pdbIds, options));
    const pi = DG.TaskBarProgressIndicator.create('Features: starting...');
    try {
      await this.ensurePdbQc(cache, pdbIds, options, pi);
      await this.ensureAligned(cache, options, pi);
      await this.ensurePocketAtoms(cache, options, pi);
      // Stage 4 lifts each feature into the Stage 2b (pocket-Cα) frame, so the
      // rendered protein must use the same Stage 2b transforms — otherwise the
      // chain-F overlay drifts by the pass-2 refinement amount.
      const alignedV2 = await this.ensureAlignedV2(cache, pi);
      const features = await this.ensureLigandFeatures(cache, pi);

      // Option A: surface a per-PDB summary in the main grid (5 rows) and
      // show the per-interaction detail on a docked side panel that updates
      // when the user clicks a different summary row. The raw `features`
      // DataFrame stays in cache.ligandFeatures for downstream Stage 5a + the
      // Mol* chain-F overlay.
      const smilesByPdb = cache.pdbQc ? ligandSmilesMap(cache.pdbQc) : undefined;
      const summary = aggregateFeaturesByPdb(features, smilesByPdb);
      summary.name = 'pdb_interaction_summary';
      this._summaryDf = summary;
      this.setActiveDataFrame(summary);
      this.setupPdbDetailPanel(features, summary);

      // Render proteins as N SEPARATE Mol* structures (same pattern as
      // Step 2's `renderAlignedStructuresMulti`) so the Step 4 row click
      // can isolate one PDB instantly via the hierarchy.toggleVisibility
      // fast path — no viewer rebuild needed.
      const structureBlocks = this.buildAlignedBlocksFromCache(alignedV2);

      pi.update(95, 'Rendering proteins in Mol*...');
      await this.renderAlignedStructuresMulti(
        structureBlocks.blocks, structureBlocks.labels);
      // Keep the per-PDB blocks for the slow-path fallback in `isolatePdb`.
      this._cachedAlignedBlocks = structureBlocks.blocks;
      this._cachedAlignedLabels = structureBlocks.labels;
      await this.setMolstarStyle({
        cartoonAlpha: 0.25,
        solventAlpha: 0.0,
      });
      // Mount the chain-F interaction overlay as N separate structures
      // (one per PDB). This is what makes Step 4 isolation actually work —
      // each PDB's interaction spheres get their own ref so we can hide
      // them in lockstep with the protein.
      pi.update(98, 'Adding interaction overlays...');
      // sizeFactor 1.3 → distinct family-colored interaction "dots" rather than
      // chunky overlapping blobs; each sphere marks one ProLIF contact point.
      await this.addInteractionOverlayPerPdb(features, 1.3);
      // Re-apply any per-PDB "Use" exclusions to the freshly-built scene so a
      // PDB the user un-ticked stays hidden across viewer rebuilds (navigating
      // away and back to Step 4 rebuilds every structure visible).
      await this.applyPdbExclusionsToViewer();
      // Same for individually un-ticked interactions (per-interaction hide).
      await this.applyInteractionExclusionsToViewer();
      // Re-apply the cartoon dim AFTER the overlay (same race-mitigation as in
      // previewConsensus — see the comment there for the full reasoning).
      // Use the retry wrapper because a single setMolstarStyle pass can
      // still race-lose when many per-row overlay structures are being
      // added concurrently.
      await this.setMolstarStyleWithRetry({cartoonAlpha: 0.25, solventAlpha: 0.0});
      pi.update(100, 'Done.');
      grok.shell.info(`Features: ${features.rowCount} ProLIF interactions ` +
        `across ${summary.rowCount} PDB(s). ` +
        'Click a summary row to see the per-interaction detail in the side panel.');
    } catch (e: any) {
      const failing = currentRunningStage(this.stageBadges) ?? 'preview';
      this.setStageStatus(failing, 'error');
      grok.shell.error(`Features failed: ${e?.message ?? e}`);
      throw e;
    } finally {
      pi.close();
    }
  }

  /**
   * Preview 5/5 — Full pipeline up to consensus. Shows the consensus pharmacophore
   * in Mol* (chain P, B-factor = frequency) on top of the aligned + pocket view.
   * Uses the existing Stage 5b renderer.
   */
  async previewConsensus(pdbIds: string[], options: PipelineOptions): Promise<void> {
    // tableView may be undefined in the wizard path.
    this.resetStageStatuses();
    this.closePdbDetailPanel();
    this.setMolstarCaption(
      'Stage 5a: chain-P spheres = consensus pharmacophore (k-means centroids ' +
      'per family). Sphere radius = max spread of contributing features in Å ' +
      '(small = tight consensus across ligands; large = loose spatial cluster).');
    const cache = this.ensurePreviewCache(this.previewFingerprint(pdbIds, options));
    const pi = DG.TaskBarProgressIndicator.create('Consensus: starting...');
    try {
      await this.ensurePdbQc(cache, pdbIds, options, pi);
      await this.ensureAligned(cache, options, pi);
      await this.ensurePocketAtoms(cache, options, pi);
      // Stage 4 (feature extraction) and therefore Stage 5a (consensus) live in
      // the Stage 2b pocket-Cα frame; render the protein from the same frame so
      // the chain-P consensus HETATMs sit inside the pocket, not next to it.
      const alignedV2 = await this.ensureAlignedV2(cache, pi);
      await this.ensureLigandFeatures(cache, pi);
      const consensus = await this.ensureConsensusModel(cache, pi, options);
      this.setActiveDataFrame(consensus);

      // Same separate-structure overlay path as previewFeatures: load proteins
      // only into the main Mol* structure, then add chain P consensus as a
      // separate structure with bigger spheres. Same reasoning — keep Mol*
      // from absorbing the consensus HETATMs into binding-site / static-ligand
      // partitions where they'd be indistinguishable from pocket residues.
      const structureBlocks = this.buildAlignedBlocksFromCache(alignedV2);
      const concat = concatPdbStructures(structureBlocks.blocks, structureBlocks.labels,
        {terminate: true});

      pi.update(95, 'Rendering proteins in Mol*...');
      await this.refreshDockedViewer(concat.body, {showBindingSite: false});
      await this.waitForMolstarReady();
      await this.setMolstarStyle({cartoonAlpha: 0.25, solventAlpha: 0.0});
      pi.update(98, 'Adding consensus overlay...');
      // sizeFactor=3.0 for consensus (slightly bigger than the per-PDB feature
      // overlay's 2.5) so the consensus centres pop above the underlying
      // crowd in Build-mode views. Chain P, with frequency baked into B-factor
      // (B-factor coloring kicks in if a Mol* color scheme other than
      // element-symbol is selected by the user).
      await this.addInteractionOverlay(consensus, 3.0, {chain: 'P'});
      // Re-apply style with polling retry. Per-row overlays add many
      // structures back-to-back, and a single setMolstarStyle pass
      // race-loses (consistently observed: nCartoon=0 in console). The
      // wrapper retries with exponential backoff (100..3200ms) until at
      // least one cartoon cell is found and dimmed. Tier 2 fix: render
      // queue.
      await this.setMolstarStyleWithRetry({cartoonAlpha: 0.25, solventAlpha: 0.0});
      pi.update(100, 'Done.');
      grok.shell.info(`Consensus: ${consensus.rowCount} consensus pharmacophore points ` +
        `across ${countDistinct(consensus, 'family')} families. ` +
        'Chain P holds the consensus centroids (large spheres, CPK-colored).');
    } catch (e: any) {
      const failing = currentRunningStage(this.stageBadges) ?? 'preview';
      this.setStageStatus(failing, 'error');
      grok.shell.error(`Consensus failed: ${e?.message ?? e}`);
      throw e;
    } finally {
      pi.close();
    }
  }

  /** Helper: rebuild the (labels, transformed+stripped PDB blocks) tuple from a Stage 2a/2b DataFrame. */
  private buildAlignedBlocksFromCache(aligned: DG.DataFrame): {labels: string[]; blocks: string[]} {
    const pdbCol = aligned.col('original_pdb')!;
    const tCol = aligned.col('transform_4x4_json')!;
    const idCol = aligned.col('pdb_id')!;
    const labels: string[] = [];
    const blocks: string[] = [];
    for (let i = 0; i < aligned.rowCount; i++) {
      const block = String(pdbCol.get(i) ?? '');
      const transformJson = String(tCol.get(i) ?? '');
      labels.push(String(idCol.get(i) ?? `row-${i}`));
      blocks.push(block ? stripCofactorsFromPdb(applyPdbTransform(block, transformJson)) : '');
    }
    return {labels, blocks};
  }

  /**
   * Render the aligned structures as N SEPARATE Mol* structures (one per
   * PDB) instead of a single concat'd structure. This is what `previewAlign`
   * uses so that Step 2 can toggle per-PDB visibility via Mol*'s native
   * `isHidden` flag — no viewer rebuild needed when the user clicks rows
   * in the RMSD table.
   *
   * Strategy: load the first PDB via BSV (sets up the viewer chrome), then
   * add each remaining PDB via `plugin.builders.structure.hierarchy.applyPreset`
   * which is the SAME high-level path BSV uses internally — it creates a
   * Trajectory → Model → Structure subtree plus the default polymer +
   * representations automatically. This is more reliable than manually
   * composing `tryCreateComponentStatic('polymer')` + `addRepresentation`,
   * which was attempted previously and produced structures that didn't
   * actually render in the canvas.
   *
   * Side effect: populates `this._alignedStructureRefs` with PDB id (upper)
   * → state ref so `isolatePdb()` can find structures to hide/show.
   */
  private async renderAlignedStructuresMulti(blocks: string[], labels: string[]): Promise<void> {
    // Drop structures with no PDB-format block (e.g. RCSB serves only mmCIF for
    // very large entries, so fetchPdbBlock returned ''). Feeding an empty block
    // to BSV/Mol* throws "Bad state: No element", which previously killed the
    // whole render — and blocks[0] being empty is the common trigger.
    const kept = blocks
      .map((b, i) => ({b, l: labels[i]}))
      .filter((p) => typeof p.b === 'string' && p.b.trim().length > 0);
    if (kept.length < blocks.length) {
      const dropped = labels.filter((_, i) => !(blocks[i] && blocks[i].trim().length > 0));
      console.warn('[render] skipping structures with no PDB-format block:', dropped);
    }
    blocks = kept.map((p) => p.b);
    labels = kept.map((p) => p.l);
    if (blocks.length === 0) return;
    // The viewer is being rebuilt — any interaction-highlight marker is gone,
    // so clear the tracked key (else the detail row would show as highlighted
    // with no matching sphere).
    this._highlightedInteractionKey = undefined;
    // 1. Mount the viewer chrome with the first structure BSV can load. A
    //    structure Mol*/BSV chokes on (very large, unusual records) must NOT
    //    kill the whole multi-structure render — fall through to the next one.
    let baseIdx = -1;
    for (let i = 0; i < blocks.length; i++) {
      try {
        await this.refreshDockedViewer(blocks[i], {showBindingSite: false});
        if (this.viewer) { baseIdx = i; break; }
      } catch (e) {
        console.error(`[renderAlignedStructuresMulti] base load failed for ${labels[i]}:`, e);
      }
    }
    if (baseIdx < 0 || !this.viewer) return;
    await this.waitForMolstarReady();

    // 2. Snapshot the structure refs that exist AFTER BSV's first load.
    //    Anything new we add after this point will produce a new structure
    //    ref that we can pair with the corresponding PDB id.
    //    NOTE: BSV exposes the Mol* PluginContext at `.viewer.plugin`, not
    //    `.plugin` directly. Getting this wrong returns undefined and
    //    silently breaks tracking — keep this path consistent with
    //    `waitForMolstarReady` / `setMolstarStyle`.
    const plugin = (this.viewer as any).viewer?.plugin;
    if (!plugin) {
      console.warn('[renderAlignedStructuresMulti] plugin not reachable; ' +
        'isolatePdb will fall back to rebuild path');
      return;
    }
    // Pick out ROOT Molecular Structure cells only — i.e. those whose parent
    // cell is a Model, not another Structure. BSV's default preset, like the
    // 'default' applyPreset path, creates a fan-out under each root Structure
    // (polymer / ligand / water / branched component "structures") that also
    // claim type name 'Molecular Structure'. Tracking the root parent gives
    // us a single ref per PDB whose `isHidden` toggle propagates to every
    // descendant — exactly what isolatePdb wants.
    const collectRootStructRefs = (): Set<string> => {
      const roots = new Set<string>();
      for (const [ref, cell] of plugin.state.data.cells) {
        const tname = cell?.obj?.type?.name;
        if (tname !== 'Structure' && tname !== 'Molecular Structure') continue;
        const parentRef = cell?.transform?.parent;
        const parentCell = parentRef ? plugin.state.data.cells.get(parentRef) : null;
        const parentType = parentCell?.obj?.type?.name;
        // A root structure's parent is the Model it was derived from.
        // Component structures have another Structure as their parent.
        if (parentType === 'Model' || parentType === 'Molecule.Model') roots.add(ref);
      }
      return roots;
    };

    const preExistingRoots = collectRootStructRefs();
    this._alignedStructureRefs.clear();
    // Pair labels[baseIdx] with the single root structure BSV produced from it.
    if (preExistingRoots.size === 1) {
      const [firstRef] = preExistingRoots;
      this._alignedStructureRefs.set(labels[baseIdx].toUpperCase(), firstRef);
    } else {
      console.warn(`[renderAlignedStructuresMulti] expected exactly 1 root structure ` +
        `after BSV load, found ${preExistingRoots.size}`);
    }

    // 3. Add the remaining PDBs via applyPreset — produces the same kind
    //    of subtree (components + cartoon reps) that BSV created for the
    //    first PDB.
    const trackedRoots = new Set(preExistingRoots);
    for (let i = 0; i < blocks.length; i++) {
      if (i === baseIdx) continue; // already mounted as the base structure
      const block = blocks[i];
      const label = labels[i];
      if (!block) continue;
      try {
        const data = await plugin.builders.data.rawData(
          {data: block, label: `cp-aligned-${label}`});
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        // applyPreset auto-creates Model → Structure → Components → Reps
        // following the default preset (cartoon for polymer, etc.).
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        // Find the freshly-created ROOT Structure cell — the one whose
        // parent is the new Model cell.
        const rootsNow = collectRootStructRefs();
        for (const ref of rootsNow) {
          if (!trackedRoots.has(ref)) {
            this._alignedStructureRefs.set(label.toUpperCase(), ref);
            trackedRoots.add(ref);
            break;
          }
        }
      } catch (e) {
        console.error(`[renderAlignedStructuresMulti] failed for ${label}:`, e);
      }
    }
    console.log(`[renderAlignedStructuresMulti] tracked ${this._alignedStructureRefs.size} / ` +
      `${blocks.length} structures: ${Array.from(this._alignedStructureRefs.keys()).join(', ')}`);
  }

  /**
   * Show only the named PDB in the Mol* viewer (hide all others), or
   * release the isolation (`pdbId === null` → show all).
   *
   * Fast path: if `_alignedStructureRefs` has all the expected PDBs
   * tracked, toggle their `isHidden` flag via `plugin.state.data.updateCellState`.
   * This is a state-tree update with no parser invocation — ~10 ms.
   *
   * Slow path: if the refs map is incomplete (e.g. `renderAlignedStructuresMulti`
   * didn't populate every structure), fall back to a viewer rebuild with
   * the filtered concat block from the cached aligned data (~1 s).
   */
  async isolatePdb(pdbId: string | null): Promise<void> {
    const blocks = this._cachedAlignedBlocks;
    const labels = this._cachedAlignedLabels;
    if (blocks.length === 0) return;

    // Fast path — only valid when EVERY label has a tracked structure ref.
    const fastPathValid = this.viewer && labels.length > 0 &&
      labels.every((l) => this._alignedStructureRefs.has(l.toUpperCase()));
    if (fastPathValid) {
      // Same `.viewer.plugin` path as renderAlignedStructuresMulti — see
      // the note there.
      const plugin = (this.viewer as any).viewer?.plugin;
      const hierarchy = plugin?.managers?.structure?.hierarchy;
      if (plugin && hierarchy?.current) {
        try {
          // Map our tracked root refs (proteins + per-PDB interaction
          // overlays) to the hierarchy manager's StructureRef objects.
          // The hierarchy manager has a toggleVisibility([refs], 'show'|'hide')
          // API that walks the structure subtree internally and flips
          // visibility on every component + representation in a SINGLE
          // state-tree transaction — exactly what we need to avoid per-cell
          // re-render storms.
          const structByRootRef = new Map<string, any>();
          for (const s of hierarchy.current.structures)
            structByRootRef.set(s.cell.transform.ref, s);

          // Collect (pdbId, rootRef) tuples from BOTH maps. A single PDB
          // typically contributes 2 structures in Step 4 (protein + overlay);
          // 1 in Steps 2/3 (protein only). Isolating must hide/show both
          // in lockstep so the interaction spheres of hidden PDBs don't
          // hang in space.
          const allTracked: Array<[string, string]> = [
            ...this._alignedStructureRefs.entries(),
          ];
          // _alignedInteractionRefs maps a PDB → MANY overlay refs (Step 4 has
          // one per family); flatten so every overlay toggles with its protein.
          for (const [id, refs] of this._alignedInteractionRefs.entries())
            for (const ref of refs) allTracked.push([id, ref]);
          const toHide: any[] = [];
          const toShow: any[] = [];
          const isPdbExcl = this._consensusSelectionHooks?.isPdbExcluded;
          const isInputExcl = this._consensusSelectionHooks?.isInputPdbExcluded;
          const isIntExcl = this._consensusSelectionHooks?.isInteractionExcluded;
          for (const [id, rootRef] of allTracked) {
            const s = structByRootRef.get(rootRef);
            if (!s) continue;
            // PDB-level: isolation (pdbId != null) shows only the targeted PDB;
            // show-all (pdbId == null) keeps PDBs the user un-ticked hidden so
            // releasing isolation doesn't resurrect an excluded PDB. BOTH the
            // Step 1 "Use" (input) and Step 4 "Use" (consensus) checkboxes
            // contribute to hide-on-show-all.
            let hidden = pdbId != null
              ? id !== pdbId.toUpperCase()
              : !!((isPdbExcl && isPdbExcl(id)) ||
                   (isInputExcl && isInputExcl(id)));
            // Interaction-level: an un-ticked interaction stays hidden even
            // when its PDB is shown or isolated.
            if (!hidden) {
              const k = this._interactionKeyByRef.get(rootRef);
              if (k && isIntExcl && isIntExcl(k)) hidden = true;
            }
            (hidden ? toHide : toShow).push(s);
          }
          if (toHide.length > 0) hierarchy.toggleVisibility(toHide, 'hide');
          if (toShow.length > 0) hierarchy.toggleVisibility(toShow, 'show');
          // Auto-reframe the camera on what's NOW visible. Without this the
          // camera keeps its previous framing (the whole-scene bounding box),
          // so an isolated PDB sitting at a distant native-frame origin
          // (Step 1, no alignment yet) ends up off-screen — the user sees an
          // empty viewport even though the structure is technically rendered.
          // requestCameraReset() with no snapshot = auto-frame visible content.
          try { plugin.managers.camera.reset(undefined, 350); } catch (_e) { /* */ }
          return;
        } catch (e) {
          console.warn('[isolatePdb] hierarchy.toggleVisibility failed, falling back:', e);
        }
      }
    }

    // Slow path — rebuild the viewer with all or one PDB.
    console.log('[isolatePdb] fast path unavailable; falling back to rebuild');
    if (pdbId == null) {
      const concat = concatPdbStructures(blocks, labels, {terminate: true});
      await this.refreshDockedViewer(concat.body, {showBindingSite: false});
    } else {
      const idx = labels.findIndex((l) => l.toUpperCase() === pdbId.toUpperCase());
      if (idx < 0) return;
      const filtered = concatPdbStructures([blocks[idx]], [labels[idx]], {terminate: true});
      await this.refreshDockedViewer(filtered.body, {showBindingSite: false});
    }
    this._cachedAlignedBlocks = blocks;
    this._cachedAlignedLabels = labels;
  }

  /** External hook — Phase 6 has no cancel button yet but the orchestrator honours an aborted flag. */
  abort(): void { this.aborted = true; }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/**
 * Family-code → summary count column name. Drives the per-PDB summary
 * DataFrame produced by `aggregateFeaturesByPdb` AND the chip rendering
 * inside the per-PDB detail panel. Aligned 1:1 with FAMILY_CODES (D/A/a/H/P/N/X)
 * so adding a family elsewhere doesn't drift this mapping silently.
 */
export const SUMMARY_FAMILY_COL: Readonly<Record<string, string>> = {
  D: 'n_donor',
  A: 'n_acceptor',
  a: 'n_aromatic',
  H: 'n_hydrophobic',
  P: 'n_positive',
  N: 'n_negative',
  X: 'n_halogen',
};

/**
 * Resolve the user's target-lookup query to a single UniProt hit. If the
 * query is a syntactically valid accession, fetch it directly. Otherwise
 * search by gene/protein name and either (a) auto-pick the only hit, or
 * (b) open a small ui.dialog with a radio list so the user disambiguates
 * across species (human / mouse / fly / ...). Returns `null` only if the
 * user cancels the dialog; throws on network / API errors.
 *
 * The `status` element is updated in place with brief progress text so the
 * user sees what's happening on long-running searches.
 */
export async function pickUniProtHit(
  query: string, status: HTMLElement,
): Promise<UniProtHit | null> {
  // Fast path — query is already a UniProt accession.
  if (isUniProtAccession(query)) {
    status.textContent = `Fetching ${query.toUpperCase()} from UniProt...`;
    const hit = await fetchUniProtById(query);
    if (!hit) throw new Error(`UniProt accession ${query} not found.`);
    return hit;
  }

  status.textContent = `Searching UniProt for "${query}"...`;
  const hits = await searchUniProt(query);
  if (hits.length === 0)
    throw new Error(`No SwissProt entries matched "${query}". Try a UniProt accession instead.`);
  if (hits.length === 1) return hits[0];

  // Multi-hit picker. Build a label-keyed map so ui.input.choice (which
  // works in string space) can round-trip back to the chosen `UniProtHit`.
  const labelOf = (h: UniProtHit): string =>
    `${h.organism || '?'} · ${h.geneName || '(no gene)'} · ` +
    `${h.accession}${h.pdbCount > 0 ? ` (${h.pdbCount} PDB)` : ''}`;
  const labels = hits.map(labelOf);
  const byLabel = new Map<string, UniProtHit>(hits.map((h, i) => [labels[i], h]));
  const choice = ui.input.choice<string>('Target', {
    value: labels[0],
    items: labels,
  });
  // Tooltip carries the full protein name so the user can disambiguate
  // ambiguous gene names without leaving the dialog.
  const tooltipRow = ui.divText('', 'cp-target-picker-tooltip');
  tooltipRow.style.fontSize = '11px';
  tooltipRow.style.color = '#555';
  tooltipRow.style.marginTop = '6px';
  tooltipRow.style.maxWidth = '380px';
  tooltipRow.style.whiteSpace = 'normal';
  const updateTooltip = () => {
    const h = byLabel.get((choice.value ?? '') as string);
    tooltipRow.textContent = h?.proteinName ?? '';
  };
  updateTooltip();
  choice.onChanged.subscribe(updateTooltip);
  const body = ui.divV([
    ui.divText(`Found ${hits.length} SwissProt entries matching "${query}". ` +
      'Pick one — the species matters because PDB structures differ across ' +
      'orthologs.', 'cp-target-picker-hint'),
    choice.root,
    tooltipRow,
  ]);
  const dlg = ui.dialog({title: 'Choose target organism'});
  dlg.add(body);
  status.textContent = 'Waiting for organism choice...';
  return new Promise<UniProtHit | null>((resolve) => {
    let picked: UniProtHit | null = null;
    dlg.onOK(() => {
      picked = byLabel.get((choice.value ?? '') as string) ?? null;
    });
    dlg.onCancel(() => { picked = null; });
    dlg.onClose.subscribe(() => resolve(picked));
    dlg.show();
  });
}


/**
 * Aggregate Stage 4 `ligand_features` into a per-PDB summary. One row per
 * (pdb_id, ligand_comp_id), with per-family counts mapped via SUMMARY_FAMILY_COL,
 * an `n_total` of non-diagnostic rows, a `top_residues` comma-joined list of
 * the three most-frequent residues, and `has_skip_warning` flagging PDBs that
 * had any skip_reason row (so the detail-panel header can surface a warning
 * icon without scanning the raw features dataframe again).
 *
 * Diagnostic rows (skip_reason non-empty) are EXCLUDED from the counts but
 * still set `has_skip_warning`. PDBs that produced only diagnostic rows still
 * appear in the summary so the user can see the warning rather than wonder
 * where their PDB went.
 */
function aggregateFeaturesByPdb(
  features: DG.DataFrame, smilesByPdb?: Map<string, string>,
): DG.DataFrame {
  const pdbCol = features.col('pdb_id');
  const compCol = features.col('ligand_comp_id');
  const famCol = features.col('family');
  const resCol = features.col('residue');
  const skipCol = features.col('skip_reason');
  if (!pdbCol || !famCol)
    return DG.DataFrame.fromColumns([
      DG.Column.fromStrings('pdb_id', []),
    ]);

  interface Acc {
    pdb_id: string;
    ligand_comp_id: string;
    counts: Record<string, number>;
    nTotal: number;
    residueCounts: Map<string, number>;
    hasSkip: boolean;
  }
  const accByKey = new Map<string, Acc>();

  for (let i = 0; i < features.rowCount; i++) {
    const pid = String(pdbCol.get(i) ?? '');
    if (!pid) continue;
    const comp = String(compCol?.get(i) ?? '');
    const key = `${pid}::${comp}`;
    let acc = accByKey.get(key);
    if (!acc) {
      acc = {
        pdb_id: pid, ligand_comp_id: comp,
        counts: {D: 0, A: 0, a: 0, H: 0, P: 0, N: 0, X: 0},
        nTotal: 0,
        residueCounts: new Map(),
        hasSkip: false,
      };
      accByKey.set(key, acc);
    }
    if (skipCol && String(skipCol.get(i) ?? '').trim() !== '') {
      acc.hasSkip = true;
      continue;
    }
    const famName = String(famCol.get(i) ?? '');
    if (!famName) continue;
    const fam = resolveFamily(famName);
    acc.counts[fam.code] = (acc.counts[fam.code] ?? 0) + 1;
    acc.nTotal += 1;
    const residue = String(resCol?.get(i) ?? '').trim();
    if (residue)
      acc.residueCounts.set(residue, (acc.residueCounts.get(residue) ?? 0) + 1);
  }

  // Roll up CSX/MSE-style skip-only ligand rows into the real drug-ligand row
  // for the SAME pdb_id. ProLIF was called once per (pdb_id, comp_id) seeded
  // by Stage 3's ligand_seed flag; modified residues that survived the
  // STANDARD_AA filter (CSX = S-OXY-CYSTEINE, MSE = selenomethionine, etc.)
  // produce a single diagnostic row and no real interactions. We want their
  // warning to surface on the PDB's real ligand row, not as a ghost row of
  // its own.
  const skipOnlyByPdb = new Map<string, boolean>();
  for (const acc of accByKey.values()) {
    if (acc.nTotal === 0 && acc.hasSkip)
      skipOnlyByPdb.set(acc.pdb_id, true);
  }
  for (const key of [...accByKey.keys()]) {
    const acc = accByKey.get(key)!;
    if (acc.nTotal === 0 && acc.hasSkip && skipOnlyByPdb.has(acc.pdb_id)) {
      // Find a sibling (same pdb_id, has interactions) to absorb the warning into.
      const sibling = [...accByKey.values()].find(
        (a) => a.pdb_id === acc.pdb_id && a.nTotal > 0);
      if (sibling) {
        sibling.hasSkip = true;
        accByKey.delete(key);
      }
      // If no sibling exists (PDB produced only skip rows), keep the entry so
      // the user sees the PDB and its warning in the summary.
    }
  }

  // Emit rows in stable, alphabetical pdb_id order so the table matches the
  // input ordering the user sees in the textarea.
  const sortedAccs = [...accByKey.values()].sort((a, b) => {
    if (a.pdb_id !== b.pdb_id) return a.pdb_id < b.pdb_id ? -1 : 1;
    return a.ligand_comp_id < b.ligand_comp_id ? -1 : 1;
  });

  const rows = sortedAccs.map((acc) => {
    const topRes = [...acc.residueCounts.entries()]
      .sort((a, b) => b[1] - a[1])
      .slice(0, 3)
      .map(([r]) => r)
      .join(', ');
    const row: Record<string, unknown> = {
      pdb_id: acc.pdb_id,
      ligand_comp_id: acc.ligand_comp_id,
      ligand_smiles: smilesByPdb?.get(acc.pdb_id) ?? '',
      n_total: acc.nTotal,
      top_residues: topRes,
      has_skip_warning: acc.hasSkip,
    };
    for (const code of FAMILY_CODES)
      row[SUMMARY_FAMILY_COL[code]] = acc.counts[code] ?? 0;
    return row;
  });

  return DG.DataFrame.fromObjects(rows) ?? DG.DataFrame.fromColumns([
    DG.Column.fromStrings('pdb_id', []),
  ]);
}


/**
 * Build a pdb_id -> ligand_smiles map from a pdb_qc DataFrame (first-occurrence
 * wins; subsequent multi-ligand rows are ignored for the lookup).
 */
function ligandSmilesMap(pdbQc: DG.DataFrame): Map<string, string> {
  const out = new Map<string, string>();
  const idCol = pdbQc.col('pdb_id');
  const smiCol = pdbQc.col('ligand_smiles');
  if (!idCol || !smiCol) return out;
  for (let i = 0; i < pdbQc.rowCount; i++) {
    const id = String(idCol.get(i));
    if (out.has(id)) continue;
    const smiles = String(smiCol.get(i) ?? '').trim();
    if (smiles) out.set(id, smiles);
  }
  return out;
}

/** Count distinct values in a string column of a DataFrame. */
function countDistinct(df: DG.DataFrame, colName: string): number {
  const c = df.col(colName);
  if (!c) return 0;
  const set = new Set<string>();
  for (let i = 0; i < df.rowCount; i++) set.add(String(c.get(i) ?? ''));
  return set.size;
}

/** Count unique (pdb_id, chain, res_seq) tuples in pocket_atoms (i.e. pocket residues). */
function countPocketResidues(pocketAtoms: DG.DataFrame): number {
  const pdbCol = pocketAtoms.col('pdb_id');
  const chainCol = pocketAtoms.col('chain');
  const resSeqCol = pocketAtoms.col('res_seq');
  const seedCol = pocketAtoms.col('ligand_seed');
  if (!pdbCol || !chainCol || !resSeqCol) return pocketAtoms.rowCount;
  const seen = new Set<string>();
  for (let i = 0; i < pocketAtoms.rowCount; i++) {
    if (seedCol && Boolean(seedCol.get(i))) continue; // exclude ligand atoms from residue count
    seen.add(`${pdbCol.get(i)}|${chainCol.get(i)}|${resSeqCol.get(i)}`);
  }
  return seen.size;
}

export function parsePdbIds(raw: string): string[] {
  return parsePdbIdsDetailed(raw).accepted;
}

/**
 * Parse the PDB-ID textarea into `accepted` (valid 4-char tokens, uppercased
 * and deduped) and `rejected` (non-empty tokens that didn't match the PDB
 * format — e.g. "EGFR", "kinase", "1AB" too-short). The orchestrator surfaces
 * the rejected list under the textarea so silent token drops don't confuse
 * the user (e.g. pasting "EGFR, kinase, 1XKK" used to leave only "1XKK" on
 * screen without any explanation).
 */
export function parsePdbIdsDetailed(raw: string): {accepted: string[]; rejected: string[]} {
  const accepted: string[] = [];
  const rejected: string[] = [];
  const seen = new Set<string>();
  for (const rawToken of raw.split(/[\s,;]+/)) {
    const t = rawToken.trim();
    if (!t) continue;
    const upper = t.toUpperCase();
    // Accept RCSB 4-char codes and registered local (dropped-file) ids.
    if (/^[0-9][A-Z0-9]{3}$/.test(upper) || isLocalId(upper)) {
      if (!seen.has(upper)) {
        seen.add(upper);
        accepted.push(upper);
      }
    } else {
      rejected.push(t);
    }
  }
  return {accepted, rejected};
}

function currentRunningStage(badges: Map<string, HTMLElement>): string | null {
  for (const [id, el] of badges)
    if (el.classList.contains('cp-stage-running')) return id;
  return null;
}

/** Load demo PDB IDs from the bundled demo-kinase-pdb-ids.csv. */
export async function loadDemoPdbIds(): Promise<string[]> {
  const df = await _package.files.readCsv('demo-kinase-pdb-ids.csv');
  const ids: string[] = [];
  const col = df.col('pdb_id');
  if (!col) return ids;
  for (let i = 0; i < df.rowCount; i++) ids.push(String(col.get(i)).toUpperCase());
  return ids;
}

// ---------------------------------------------------------------------------
// STUB-PHASE-6 — each function below is replaced by a grok.functions.call to a
// Python script in the corresponding phase. Schemas match Appendix B.
// ---------------------------------------------------------------------------

/**
 * Stage 2a — Kabsch global Cα alignment. Calls the Python script via
 * `grok.functions.call(FN.STAGE2A, ...)`. On failure (no Python runtime, RCSB
 * down inside the container, MDAnalysis missing, etc.), logs a warning and
 * falls back to the TS-side stub so the UI stays usable.
 */
async function runStage2aAlign(
  pdbQc: DG.DataFrame, opts: PipelineOptions,
): Promise<DG.DataFrame> {
  const refPdbIdArg = (opts.refPdbId ?? '').trim();
  console.log(`[Stage 2a] calling ${FN.STAGE2A} with ${pdbQc.rowCount} input rows, ` +
    `reference_pdb_id=${JSON.stringify(refPdbIdArg) || '"" (auto)'}`);
  try {
    const result = await grok.functions.call(FN.STAGE2A, {
      pdb_qc: pdbQc,
      chain_selection: 'auto',
      alt_loc_filter: true,
      min_occupancy: 0.5,
      // User-overridable reference PDB; empty string = auto-pick (Python
      // falls back to highest-resolution entry).
      reference_pdb_id: refPdbIdArg,
    });
    console.log('[Stage 2a] raw result type:', result?.constructor?.name, result);
    const df = (result instanceof DG.DataFrame) ? result : (result as any)?.aligned_structures;
    console.log(`[Stage 2a] resolved DataFrame? ${df instanceof DG.DataFrame}, rows=${df?.rowCount ?? 'n/a'}`);
    if (df && df.rowCount > 0) {
      const tCol = df.col('transform_4x4_json');
      if (tCol) {
        for (let i = 0; i < Math.min(df.rowCount, 5); i++)
          console.log(`[Stage 2a] row ${i} transform = ${String(tCol.get(i)).slice(0, 80)}`);
      } else
        console.warn('[Stage 2a] result DataFrame has no transform_4x4_json column!');
      df.col('original_pdb')!.semType = DG.SEMTYPE.MOLECULE3D;
      df.col('original_pdb')!.setTag(DG.TAGS.UNITS, 'pdb');
      return df;
    }
    throw new Error('Stage 2a Python returned an empty DataFrame.');
  } catch (e: any) {
    console.error('[Stage 2a] Python call failed:', e);
    grok.shell.warning(`Stage 2a Python unavailable (${e?.message ?? e}). ` +
      'Falling back to TS stub (no real alignment — structures stay in native frames).');
    return stubStage2aAlign(pdbQc, opts);
  }
}

async function stubStage2aAlign(
  pdbQc: DG.DataFrame, _opts: PipelineOptions,
): Promise<DG.DataFrame> {
  // STUB-PHASE-6: real impl is `grok.functions.call(FN.STAGE2A, {pdb_qc: pdbQc, ...})`.
  // We fetch real PDB blocks (so the Preview button has something to show, and so
  // downstream stubs / Phase 3 work against realistic data). Transforms remain
  // identity — only the real Python Stage 2a produces non-identity transforms.
  const pdbIdCol = pdbQc.col('pdb_id');
  if (!pdbIdCol) throw new Error('stage2a stub: missing pdb_id');
  const identity = JSON.stringify(IDENTITY_4X4);

  // Dedup: Stage 1 emits one row per (pdb_id, ligand). For alignment we want one
  // row per pdb_id — preserve the order of first occurrence.
  const seen = new Set<string>();
  const uniqueIds: string[] = [];
  for (let i = 0; i < pdbQc.rowCount; i++) {
    const id = String(pdbIdCol.get(i));
    if (seen.has(id)) continue;
    seen.add(id);
    uniqueIds.push(id);
  }
  if (uniqueIds.length === 0) throw new Error('stage2a stub: no pdb ids to align');
  const refPdbId = uniqueIds[0];

  const pdbBlocks = await Promise.all(uniqueIds.map((id) => fetchPdbBlock(id)));
  const missing = uniqueIds.filter((_, i) => !pdbBlocks[i]);
  if (missing.length)
    grok.shell.warning(`Could not download PDB(s): ${missing.join(', ')}. Continuing with empty blocks.`);

  const cols = [
    DG.Column.fromStrings('pdb_id', uniqueIds),
    DG.Column.fromStrings('original_pdb', pdbBlocks),
    DG.Column.fromStrings('transform_4x4_json', uniqueIds.map(() => identity)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'rmsd_to_ref', uniqueIds.map((_, i) => (i === 0 ? 0.0 : 1.2))),
    DG.Column.fromStrings('ref_pdb_id', uniqueIds.map(() => refPdbId)),
    DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'bootstrap_failed', uniqueIds.map(() => false)),
    DG.Column.fromStrings('cluster_id', uniqueIds.map(() => 'C1')),
  ];
  const df = DG.DataFrame.fromColumns(cols);
  df.col('original_pdb')!.semType = DG.SEMTYPE.MOLECULE3D;
  df.col('original_pdb')!.setTag(DG.TAGS.UNITS, 'pdb');
  return df;
}

/**
 * Stage 3 — Pocket isolation. Calls FN.STAGE3 (Python) and falls back to the
 * TS stub if the call fails. Selects between cutoff and DBSCAN paths via
 * `options.pocketMethod`.
 */
async function runStage3IsolatePocket(
  aligned: DG.DataFrame, _clusterIds: string[], opts: PipelineOptions,
): Promise<DG.DataFrame> {
  console.log(`[Stage 3] calling ${FN.STAGE3} with ${aligned.rowCount} aligned rows ` +
    `(method=${opts.pocketMethod}, radius=${opts.pocketRadius}A)`);
  try {
    const result = await grok.functions.call(FN.STAGE3, {
      aligned_structures: aligned,
      pocket_method: opts.pocketMethod,
      pocket_radius: opts.pocketRadius,
      dbscan_eps: 4.0,
      dbscan_min_samples: 6,
    });
    const df = (result instanceof DG.DataFrame) ? result : (result as any)?.pocket_atoms;
    console.log(`[Stage 3] resolved DataFrame? ${df instanceof DG.DataFrame}, rows=${df?.rowCount ?? 'n/a'}`);
    if (df && df.rowCount > 0) return df;
    throw new Error('Stage 3 Python returned an empty DataFrame.');
  } catch (e: any) {
    console.error('[Stage 3] Python call failed:', e);
    grok.shell.warning(`Stage 3 Python unavailable (${e?.message ?? e}). ` +
      'Falling back to TS stub (mock pocket atoms).');
    return stubStage3IsolatePocket(aligned, _clusterIds, opts);
  }
}

async function stubStage3IsolatePocket(
  aligned: DG.DataFrame, _clusterIds: string[], _opts: PipelineOptions,
): Promise<DG.DataFrame> {
  // STUB-PHASE-6: real impl is `grok.functions.call(FN.STAGE3, {aligned, ...})`.
  const pdbIdCol = aligned.col('pdb_id');
  const ids: string[] = [];
  if (pdbIdCol) for (let i = 0; i < aligned.rowCount; i++) ids.push(String(pdbIdCol.get(i)));

  // Emit 6 mock pocket residues per structure.
  const rows: any[] = [];
  for (const pdbId of ids) {
    for (let k = 0; k < 6; k++) {
      rows.push({
        pdb_id: pdbId, chain: 'A', res_name: 'ASP', res_seq: 720 + k, atom_name: 'CA',
        element: 'C', x_consensus: k * 1.5, y_consensus: 0, z_consensus: 0,
        ligand_seed: false, cluster_label: 0,
      });
    }
  }
  return DG.DataFrame.fromObjects(rows) ?? DG.DataFrame.fromColumns([
    DG.Column.fromStrings('pdb_id', []),
  ]);
}

/**
 * Stage 2b — pocket-only Cα Kabsch refinement. Calls FN.STAGE2B (Python). Falls
 * back to the TS stub if the call fails. `selectedPdbIds` is the flattened PDB
 * list from the cluster picker — an empty list (or ['all']) means "use all".
 */
async function runStage2bAlignPocket(
  aligned1: DG.DataFrame, pocketAtoms: DG.DataFrame, selectedPdbIds: string[],
): Promise<DG.DataFrame> {
  console.log(`[Stage 2b] calling ${FN.STAGE2B} with ${aligned1.rowCount} aligned rows, ` +
    `${pocketAtoms.rowCount} pocket atoms, ${selectedPdbIds.length} selected PDB(s)`);
  try {
    const result = await grok.functions.call(FN.STAGE2B, {
      aligned_structures: aligned1,
      pocket_atoms: pocketAtoms,
      selected_pdb_ids: JSON.stringify(selectedPdbIds),
    });
    const df = (result instanceof DG.DataFrame) ? result : (result as any)?.aligned_structures_v2;
    console.log(`[Stage 2b] resolved DataFrame? ${df instanceof DG.DataFrame}, rows=${df?.rowCount ?? 'n/a'}`);
    if (df && df.rowCount > 0) {
      df.col('original_pdb')!.semType = DG.SEMTYPE.MOLECULE3D;
      df.col('original_pdb')!.setTag(DG.TAGS.UNITS, 'pdb');
      return df;
    }
    throw new Error('Stage 2b Python returned an empty DataFrame.');
  } catch (e: any) {
    console.error('[Stage 2b] Python call failed:', e);
    grok.shell.warning(`Stage 2b Python unavailable (${e?.message ?? e}). ` +
      'Falling back to TS stub (pass-2 transforms = scaled pass-1).');
    return stubStage2bAlignPocket(aligned1, pocketAtoms, selectedPdbIds);
  }
}

async function stubStage2bAlignPocket(
  aligned1: DG.DataFrame, _pocketAtoms: DG.DataFrame, _clusterIds: string[],
): Promise<DG.DataFrame> {
  // STUB-PHASE-6: pass-2 produces the same shape as 2a; we mark pass-2 RMSDs slightly tighter.
  const out = aligned1.clone();
  const rmsdCol = out.col('rmsd_to_ref');
  if (rmsdCol) {
    for (let i = 0; i < out.rowCount; i++) {
      const v = Number(rmsdCol.get(i)) || 0;
      rmsdCol.set(i, Math.max(0, v * 0.7));
    }
  }
  return out;
}

/**
 * Stage 4 — protein-ligand interaction extraction via ProLIF. Calls FN.STAGE4
 * (Python in the BSV-shared conda env: rdkit + mdanalysis + prolif + pdbfixer +
 * openmm + openbabel). Replaces the earlier SMARTS-only extractor that ran
 * against ligand atoms in isolation; the SMARTS-only version is archived at
 * `files/_reference/04_extract_features_smarts_only.py`.
 *
 * ProLIF runs once per PDB. For our 5-PDB EGFR demo that's roughly 75-150s
 * end-to-end (pdbfixer + ProLIF is ~15-30s per PDB, single-threaded).
 * The bundled SMARTS CSV is no longer passed to Stage 4 — ProLIF uses the
 * Datagrok-aligned SMARTS baked into the script itself.
 */
async function runStage4ExtractFeatures(
  aligned2: DG.DataFrame, pocketAtoms: DG.DataFrame, selectedPdbIds: string[],
): Promise<DG.DataFrame> {
  console.log(`[Stage 4 ProLIF] calling ${FN.STAGE4} with ${aligned2.rowCount} aligned rows, ` +
    `${pocketAtoms.rowCount} pocket atoms, ${selectedPdbIds.length} selected (~15-30s/PDB)`);
  try {
    const result = await grok.functions.call(FN.STAGE4, {
      aligned_structures: aligned2,
      pocket_atoms: pocketAtoms,
      selected_pdb_ids: JSON.stringify(selectedPdbIds),
    });
    const df = (result instanceof DG.DataFrame) ? result : (result as any)?.ligand_features;
    console.log(`[Stage 4] resolved DataFrame? ${df instanceof DG.DataFrame}, rows=${df?.rowCount ?? 'n/a'}`);
    if (df && df.rowCount >= 0) return df; // can be 0 if all ligands failed; let pipeline decide
    throw new Error('Stage 4 Python returned an empty DataFrame.');
  } catch (e: any) {
    console.error('[Stage 4] Python call failed:', e);
    grok.shell.warning(`Stage 4 Python unavailable (${e?.message ?? e}). ` +
      'Falling back to TS stub (mock 4-feature ligands).');
    return stubStage4ExtractFeatures(aligned2, pocketAtoms, selectedPdbIds);
  }
}

async function stubStage4ExtractFeatures(
  aligned2: DG.DataFrame, _pocketAtoms: DG.DataFrame, _clusterIds: string[],
): Promise<DG.DataFrame> {
  // STUB-PHASE-6: emit 4 features per ligand spread across 4 families so Stage 5a has
  // enough to cluster. Coordinates lie in a small cube around the origin.
  const pdbIdCol = aligned2.col('pdb_id');
  const ids: string[] = [];
  if (pdbIdCol) for (let i = 0; i < aligned2.rowCount; i++) ids.push(String(pdbIdCol.get(i)));

  const familySeq = ['Donor', 'Acceptor', 'Aromatic', 'Hydrophobic'];
  const rows: any[] = [];
  for (const pdbId of ids) {
    familySeq.forEach((fam, idx) => {
      rows.push({
        pdb_id: pdbId,
        ligand_comp_id: 'STI',
        family: fam,
        feature_idx: idx,
        atom_indices_json: '[]',
        x: idx * 1.2,
        y: 0,
        z: 0,
        bond_order_corrected: true,
        skip_reason: '',
      });
    });
  }
  return DG.DataFrame.fromObjects(rows) ?? DG.DataFrame.fromColumns([
    DG.Column.fromStrings('pdb_id', []),
  ]);
}

/**
 * Stable identifier for one ProLIF interaction row. Used by BOTH the Step 4
 * per-interaction selection checkboxes (tier 2) and `filterFeaturesForConsensus`
 * so the wizard and the orchestrator agree on what "this interaction" means.
 * Composite of the identifying fields: PDB, protein residue, interaction type,
 * and the ligand-centroid coords (rounded) to disambiguate multiple contacts
 * to the same residue.
 */
export function consensusFeatureKey(features: DG.DataFrame, i: number): string {
  const g = (c: string): string => String(features.col(c)?.get(i) ?? '');
  const n = (c: string): string => {
    const v = Number(features.col(c)?.get(i));
    return Number.isFinite(v) ? v.toFixed(2) : '';
  };
  return `${g('pdb_id').toUpperCase()}|${g('residue')}|${g('interaction_type')}|` +
    `${n('x')},${n('y')},${n('z')}`;
}

/**
 * Consensus-only cache signature: the Stage 5a k-means knobs plus the user's
 * PDB and per-interaction exclusions. `ensureConsensusModel` recomputes the
 * consensus whenever this changes — without invalidating Stages 1-4.
 */
function consensusSignature(opts: PipelineOptions): string {
  return JSON.stringify({
    kq: opts.kq,
    minFrac: opts.minClusterSizeFraction,
    top: opts.topClusterNumber,
    exclPdbs: [...(opts.consensusExcludedPdbs ?? [])].map((s) => s.toUpperCase()).sort(),
    exclInts: [...(opts.consensusExcludedInteractions ?? [])].sort(),
  });
}

/**
 * Drop the rows the user excluded from the consensus: whole PDBs
 * (`consensusExcludedPdbs`) and individual interactions
 * (`consensusExcludedInteractions`, keyed by `consensusFeatureKey`). Returns
 * the input unchanged when nothing is excluded. Diagnostic (skip_reason) rows
 * are left as-is — Stage 5a already ignores them.
 */
function filterFeaturesForConsensus(
  features: DG.DataFrame, opts: PipelineOptions,
): DG.DataFrame {
  const exclPdbs = new Set((opts.consensusExcludedPdbs ?? []).map((s) => s.toUpperCase()));
  const exclInts = new Set(opts.consensusExcludedInteractions ?? []);
  if (exclPdbs.size === 0 && exclInts.size === 0) return features;
  const pdbCol = features.col('pdb_id');
  if (!pdbCol) return features;
  const keep = new Array<boolean>(features.rowCount);
  let nKept = 0;
  for (let i = 0; i < features.rowCount; i++) {
    const pid = String(pdbCol.get(i) ?? '').toUpperCase();
    let drop = exclPdbs.has(pid);
    if (!drop && exclInts.size > 0) drop = exclInts.has(consensusFeatureKey(features, i));
    keep[i] = !drop;
    if (keep[i]) nKept++;
  }
  if (nKept === features.rowCount) return features;
  const mask = DG.BitSet.create(features.rowCount, (i) => keep[i]);
  return features.clone(mask);
}

/** Drop rows from pdb_qc whose pdb_id appears in opts.excludedInputPdbs (the
 *  Step 1 "Use" deselections). Returns the original df when nothing is
 *  excluded so the common case stays zero-cost. Used right before Stage 2a so
 *  excluded PDBs never enter alignment / pocket / features / consensus. */
function filterPdbQcByInputExclusion(
  pdbQc: DG.DataFrame, opts: PipelineOptions,
): DG.DataFrame {
  const excl = new Set((opts.excludedInputPdbs ?? []).map((s) => s.toUpperCase()));
  if (excl.size === 0) return pdbQc;
  const pdbCol = pdbQc.col('pdb_id');
  if (!pdbCol) return pdbQc;
  const keep = new Array<boolean>(pdbQc.rowCount);
  let nKept = 0;
  for (let i = 0; i < pdbQc.rowCount; i++) {
    keep[i] = !excl.has(String(pdbCol.get(i) ?? '').toUpperCase());
    if (keep[i]) nKept++;
  }
  if (nKept === pdbQc.rowCount) return pdbQc;
  const mask = DG.BitSet.create(pdbQc.rowCount, (i) => keep[i]);
  return pdbQc.clone(mask);
}

/**
 * Stage 5a — per-family k-means consensus. Calls FN.STAGE5A. Defaults come
 * from TeachOpenCADD T009 (kq=7, top_cluster_number=4, min_cluster_size_fraction=0.75)
 * but the orchestrator now passes the live PipelineOptions values so the user
 * can tune them from the Advanced options accordion. Falls back to the TS
 * stub on failure.
 */
async function runStage5aConsensusKmeans(
  ligandFeatures: DG.DataFrame, opts: PipelineOptions,
): Promise<DG.DataFrame> {
  console.log(`[Stage 5a] calling ${FN.STAGE5A} with ${ligandFeatures.rowCount} feature rows ` +
    `(kq=${opts.kq}, min_frac=${opts.minClusterSizeFraction}, top=${opts.topClusterNumber})`);
  try {
    const result = await grok.functions.call(FN.STAGE5A, {
      ligand_features: ligandFeatures,
      kq: opts.kq,
      top_cluster_number: opts.topClusterNumber,
      min_cluster_size_fraction: opts.minClusterSizeFraction,
      kmeans_random_state: 42,
    });
    const df = (result instanceof DG.DataFrame) ? result : (result as any)?.consensus_model;
    console.log(`[Stage 5a] resolved DataFrame? ${df instanceof DG.DataFrame}, rows=${df?.rowCount ?? 'n/a'}`);
    if (df && df.rowCount >= 0) return df;
    throw new Error('Stage 5a Python returned an empty DataFrame.');
  } catch (e: any) {
    console.error('[Stage 5a] Python call failed:', e);
    grok.shell.warning(`Stage 5a Python unavailable (${e?.message ?? e}). ` +
      'Falling back to TS stub (one mock row per family).');
    return stubStage5aConsensusKmeans(ligandFeatures);
  }
}

async function stubStage5aConsensusKmeans(_ligandFeatures: DG.DataFrame): Promise<DG.DataFrame> {
  // STUB-PHASE-6: emit one row per family across the consensus frame so the Mol*
  // viewer renders 7 colored HETATMs at varied B-factors.
  const rows = FAMILY_CODES.map((code, i) => ({
    family: FAMILY_MAP[code].name,
    x: i * 2.0,
    y: 0,
    z: 0,
    frequency: 0.4 + 0.1 * i, // 0.4..1.0 — varied so B-factor coloring (if it fires) is visible
    n_ligands: 5,
  }));
  return DG.DataFrame.fromObjects(rows)!;
}

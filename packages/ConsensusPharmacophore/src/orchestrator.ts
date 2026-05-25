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
import {enrichPdbList} from './rcsb-client';
import {applyPdbTransform, concatPdbStructures, fetchPdbBlock, IDENTITY_4X4,
  ligandFeaturesToOverlayBlock, stripCofactorsFromPdb} from './pdb-utils';
import {fetchUniProtById, getPdbIdsForAccession, isUniProtAccession,
  searchUniProt, UniProtHit} from './uniprot-client';

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
}

export class ConsensusPharmacophoreApp {
  tableView?: DG.TableView;
  private viewer?: DG.Viewer;
  private viewerPlaceholder?: HTMLElement;
  private viewerPlaceholderNode?: DG.DockNode;
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
      this.viewerPlaceholder, DG.DOCK_TYPE.RIGHT, null, 'Mol*', 0.4);

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

    const pdbInput = ui.input.textArea('PDB IDs', {
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
    if (!this.tableView) throw new Error('runPipeline called before init()');
    this.aborted = false;
    this.resetStageStatuses();
    this.closePdbDetailPanel();
    this.setMolstarCaption(
      'Running full pipeline — cluster picker will open after Stage 2a to let ' +
      'you drop conformationally mixed PDBs before pocket isolation.');

    const pi = DG.TaskBarProgressIndicator.create('Consensus Pharmacophore: running...');
    try {
      // Stage 1 — Enrich PDB metadata (RCSB-backed; replaces stubStage1)
      this.setStageStatus('stage1', 'running');
      pi.update(5, 'Stage 1: enriching PDB metadata...');
      const pdbQc = await enrichPdbList(pdbIds, options);
      if (pdbQc.rowCount === 0)
        throw new Error('No usable PDB entries after QC filtering. ' +
          'Check the PDB IDs, loosen the resolution cap, or disable the X-ray filter.');
      this.tableView.dataFrame = pdbQc;
      pdbQc.name = 'pdb_qc';
      this.setStageStatus('stage1', 'done');
      if (this.aborted) return null;

      // Stage 2a — Kabsch pass 1
      this.setStageStatus('stage2a', 'running');
      pi.update(20, 'Stage 2a: aligning structures (pass 1)...');
      const aligned1 = await runStage2aAlign(pdbQc, options);
      this.tableView.dataFrame = aligned1;
      aligned1.name = 'aligned_structures (pass 1)';
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
      this.tableView.dataFrame = pocketAtoms;
      pocketAtoms.name = 'pocket_atoms';
      this.setStageStatus('stage3', 'done');
      if (this.aborted) return null;

      // Stage 2b — Kabsch pass 2 (pocket-only Calpha; filtered to picker's PDB list)
      this.setStageStatus('stage2b', 'running');
      pi.update(55, 'Stage 2b: re-aligning on pocket Calpha...');
      const aligned2 = await runStage2bAlignPocket(aligned1, pocketAtoms, pick.selectedPdbIds);
      this.tableView.dataFrame = aligned2;
      aligned2.name = 'aligned_structures (pass 2)';
      this.setStageStatus('stage2b', 'done');
      if (this.aborted) return null;

      // Stage 4 — Protein-ligand interaction extraction (ProLIF)
      this.setStageStatus('stage4', 'running');
      pi.update(70, 'Stage 4: detecting protein-ligand interactions (ProLIF, ~15-30s/PDB)...');
      const ligandFeatures = await runStage4ExtractFeatures(aligned2, pocketAtoms, pick.selectedPdbIds);
      this.tableView.dataFrame = ligandFeatures;
      ligandFeatures.name = 'ligand_features';
      this.setStageStatus('stage4', 'done');
      if (this.aborted) return null;

      // Stage 5a — k-means consensus
      this.setStageStatus('stage5a', 'running');
      pi.update(85, 'Stage 5a: computing per-family consensus...');
      const consensus = await runStage5aConsensusKmeans(ligandFeatures, options);
      this.tableView.dataFrame = consensus;
      consensus.name = 'consensus_model';
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
    if (!this.tableView) return;
    const allOptions = {pdb: pdbBlock, ...extraOptions};
    if (!this.viewer) {
      // First render — create the Mol* viewer with the consensus PDB, dock it on the right,
      // remove the placeholder.
      this.viewer = await this.tableView.dataFrame.plot.fromType(
        'Biostructure', allOptions) as DG.Viewer;
      this.tableView.dockManager.dock(
        this.viewer as unknown as DG.Viewer, DG.DOCK_TYPE.RIGHT, null, 'Mol*', 0.4);
      if (this.viewerPlaceholderNode) {
        this.tableView.dockManager.close(this.viewerPlaceholderNode);
        this.viewerPlaceholderNode = undefined;
      }
      return;
    }
    (this.viewer as any).setOptions(allOptions);
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
    if (!this.tableView) return;
    if (!this.pdbDetailHost) {
      this.pdbDetailHost = ui.div([], 'cp-pdb-detail-panel');
      this.pdbDetailHost.style.padding = '8px 10px';
      this.pdbDetailHost.style.minWidth = '240px';
      this.pdbDetailHost.style.fontSize = '12px';
      this.pdbDetailHost.style.overflowY = 'auto';
    }
    // Dock on first call. ratio 0.2 = ~20% of remaining width — leaves Mol* room.
    if (!this.pdbDetailNode) {
      this.pdbDetailNode = this.tableView.dockManager.dock(
        this.pdbDetailHost, DG.DOCK_TYPE.RIGHT, null, 'PDB detail', 0.2);
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
    if (this.pdbDetailSub) {
      try { this.pdbDetailSub.unsubscribe(); } catch (_e) { /* */ }
      this.pdbDetailSub = undefined;
    }
    if (this.pdbDetailNode && this.tableView) {
      try { this.tableView.dockManager.close(this.pdbDetailNode); } catch (_e) { /* */ }
      this.pdbDetailNode = undefined;
    }
    this.pdbDetailHost = undefined;
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
    list.style.display = 'grid';
    list.style.gridTemplateColumns = '12px 1fr auto auto';
    list.style.columnGap = '6px';
    list.style.rowGap = '3px';
    list.style.alignItems = 'center';

    const famCol = features.col('family');
    const intCol = features.col('interaction_type');
    const resCol = features.col('residue');
    const distCol = features.col('distance');
    const pdbCol = features.col('pdb_id');
    const skipCol = features.col('skip_reason');

    let nShown = 0;
    for (let i = 0; i < features.rowCount && nShown < 200; i++) {
      if (pdbCol && String(pdbCol.get(i)) !== pdbId) continue;
      if (skipCol && String(skipCol.get(i) ?? '').trim() !== '') continue;
      const famName = String(famCol?.get(i) ?? '');
      const fam = resolveFamily(famName);
      const dot = ui.div();
      dot.style.width = '8px';
      dot.style.height = '8px';
      dot.style.borderRadius = '50%';
      dot.style.background = fam.hexColor;
      list.append(dot);
      list.append(ui.divText(String(intCol?.get(i) ?? '')));
      list.append(ui.divText(String(resCol?.get(i) ?? ''), 'cp-pdb-detail-res'));
      const d = Number(distCol?.get(i));
      list.append(ui.divText(Number.isFinite(d) ? d.toFixed(2) + ' Å' : ''));
      nShown += 1;
    }
    if (nShown === 0)
      list.append(ui.divText('(no interactions for this PDB)',
        'cp-pdb-detail-empty'));
    this.pdbDetailHost.append(list);
  }

  // -------------------------------------------------------------------------
  // Preview cache + per-stage ensures
  // -------------------------------------------------------------------------

  private previewFingerprint(pdbIds: string[], opts: PipelineOptions): string {
    return JSON.stringify({
      pdbIds: [...pdbIds].sort(),
      maxResolution: opts.maxResolution,
      requireXray: opts.requireXray,
      minLigandMw: opts.minLigandMw,
      pocketMethod: opts.pocketMethod,
      pocketRadius: opts.pocketRadius,
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
    const pdbQc = await enrichPdbList(pdbIds, opts);
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

  /** Run Stage 2a (or use cached) and publish aligned_structures. */
  private async ensureAligned(
    cache: PreviewCache, opts: PipelineOptions, pi: DG.ProgressIndicator,
  ): Promise<DG.DataFrame> {
    if (cache.aligned) return cache.aligned;
    if (!cache.pdbQc) throw new Error('ensureAligned: pdb_qc missing from cache');
    this.setStageStatus('stage2a', 'running');
    pi.update(50, 'Stage 2a: aligning structures (Kabsch)...');
    const aligned = await runStage2aAlign(cache.pdbQc, opts);
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
    if (cache.consensusModel) return cache.consensusModel;
    if (!cache.ligandFeatures) throw new Error('ensureConsensusModel: ligandFeatures missing');
    this.setStageStatus('stage5a', 'running');
    pi.update(92, 'Stage 5a: k-means consensus...');
    const consensus = await runStage5aConsensusKmeans(cache.ligandFeatures, opts);
    consensus.name = 'consensus_model (preview)';
    this.setStageStatus('stage5a', 'done');
    cache.consensusModel = consensus;
    return consensus;
  }

  /**
   * Preview 1/3 — fetch PDBs from RCSB (no alignment). Structures appear in their
   * native crystallographic frames; useful for sanity-checking the inputs.
   */
  async previewFetch(pdbIds: string[], options: PipelineOptions): Promise<void> {
    if (!this.tableView) throw new Error('previewFetch called before init()');
    this.resetStageStatuses();
    this.closePdbDetailPanel();
    this.setMolstarCaption(
      'Fetched PDBs shown in their native crystallographic frames (no alignment yet).');
    const cache = this.ensurePreviewCache(this.previewFingerprint(pdbIds, options));
    const pi = DG.TaskBarProgressIndicator.create('Fetch PDBs: starting...');
    try {
      const pdbQc = await this.ensurePdbQc(cache, pdbIds, options, pi);
      this.tableView.dataFrame = pdbQc;
      const rawBlocks = await this.ensureRawBlocks(cache, pdbQc, pi);

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
      await this.refreshDockedViewer(
        concatPdbStructures(blocks, labels).body,
        {showBindingSite: false},
      );
      // Reset to opaque cartoon + show solvent (in case Pocket was clicked previously).
      await this.setMolstarStyle({cartoonAlpha: 1.0, solventAlpha: 1.0});
      pi.update(100, 'Done.');
      grok.shell.info(`Fetch: ${labels.length} PDB(s) shown in their native frames. ` +
        'Click "Align" to superimpose them.');
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
    if (!this.tableView) throw new Error('previewAlign called before init()');
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
      this.tableView.dataFrame = aligned;

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
      await this.refreshDockedViewer(
        concatPdbStructures(transformedBlocks, labels).body,
        {showBindingSite: false},
      );
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
    if (!this.tableView) throw new Error('previewPocket called before init()');
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
      this.tableView.dataFrame = aligned;

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

      // Hand the aligned structures to Mol* and turn on its native binding-site
      // renderer (BSV's `showBindingSite` flag wires up `applyBindingSiteView()`,
      // which picks the largest non-cofactor ligand and highlights residues
      // within `bindingSiteRadius`). No custom HETATM overlay needed.
      pi.update(90, 'Rendering in Mol*...');
      const oldRefsKey = this.snapshotBindingSiteRefs();
      await this.refreshDockedViewer(
        concatPdbStructures(transformedBlocks, labels).body,
        {
          showBindingSite: true,
          bindingSiteRadius: options.pocketRadius,
          bindingSiteWholeResidues: true,
        },
      );
      // Wait for BSV to FINISH applying the binding site (clear-then-build runs on a
      // viewSyncer queue, so we wait until bindingSiteRefs CHANGE from their old
      // value — otherwise our style edits land on cells BSV is about to destroy).
      await this.waitForBindingSiteApplied(oldRefsKey);
      const isSurface = options.pocketRep === 'gaussian-surface';
      await this.setMolstarStyle({
        cartoonAlpha: 0.25,
        bindingSiteRep: options.pocketRep,
        bindingSiteAlpha: isSurface ? 0.55 : 1.0,
        solventAlpha: 0.0,
      });
      pi.update(100, 'Done.');
      const nResidues = countPocketResidues(pocketAtoms);
      grok.shell.info(`Pocket: ${nResidues} pocket residue(s) detected ` +
        `(${options.pocketMethod}, ${options.pocketRadius.toFixed(1)} A). ` +
        `Rep: ${options.pocketRep}; cartoon translucent; solvent hidden.`);
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
    if (!this.tableView) throw new Error('previewFeatures called before init()');
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
      this.tableView.dataFrame = summary;
      this.setupPdbDetailPanel(features, summary);

      // Render proteins ONLY (no chain F splice). Chain F gets added afterward as
      // a SEPARATE Mol* structure so we can style it with bigger spheres without
      // colliding with BSV's static-ligand / binding-site partitions.
      const structureBlocks = this.buildAlignedBlocksFromCache(alignedV2);
      const concat = concatPdbStructures(structureBlocks.blocks, structureBlocks.labels,
        {terminate: true});

      pi.update(95, 'Rendering proteins in Mol*...');
      await this.refreshDockedViewer(concat.body, {
        // Binding site stays OFF — the user only wants the interaction dots,
        // not the pocket-residue sphere mesh on top of them.
        showBindingSite: false,
      });
      // With showBindingSite:false there's no binding-site refs to wait on,
      // so we poll for plugin readiness ourselves — otherwise setMolstarStyle
      // and addInteractionOverlay early-return on a not-yet-initialized plugin.
      await this.waitForMolstarReady();
      await this.setMolstarStyle({
        cartoonAlpha: 0.25,
        solventAlpha: 0.0,
      });
      // Overlay chain F as a separate Mol* structure (big spacefill, no
      // binding-site interference). sizeFactor 2.5 vs Mol*'s default 1.0 makes
      // each interaction sphere about 2.5× the radius of a regular spacefill
      // atom — clearly bigger than what the binding site used to render at.
      pi.update(98, 'Adding interaction overlay...');
      await this.addInteractionOverlay(features, 2.5);
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
    if (!this.tableView) throw new Error('previewConsensus called before init()');
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
      this.tableView.dataFrame = consensus;

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
const SUMMARY_FAMILY_COL: Readonly<Record<string, string>> = {
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
async function pickUniProtHit(
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

function parsePdbIds(raw: string): string[] {
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
    if (/^[0-9][A-Z0-9]{3}$/.test(upper)) {
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
  console.log(`[Stage 2a] calling ${FN.STAGE2A} with ${pdbQc.rowCount} input rows`);
  try {
    const result = await grok.functions.call(FN.STAGE2A, {
      pdb_qc: pdbQc,
      chain_selection: 'auto',
      alt_loc_filter: true,
      min_occupancy: 0.5,
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

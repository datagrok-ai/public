import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as chemSearches from '../chem-searches';
import {similarityMetric} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import $ from 'cash-dom';
import {Fingerprint} from '../utils/chem-common';
import {renderMolecule} from '../rendering/render-molecule';
import {ChemSearchBaseViewer, SIMILARITY, RowSourceTypes} from './chem-search-base-viewer';
import {malformedDataWarning} from '../utils/malformed-data-utils';
import '../../css/chem.css';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {Subscription} from 'rxjs';

/** Similarity can display up to this many results; the shared base default (and the diversity
 * viewer, whose selection is O(n^2)) keep a much smaller cap. */
const MAX_LIMIT_SIMILARITY = 10000;
/** Hard cap on how many cards are built into the DOM at once, independent of `limit`. Select &
 * filter act on the full above-cutoff set, so capping the display loses nothing functionally. */
const DISPLAY_CAP = 200;

export class ChemSimilarityViewer extends ChemSearchBaseViewer {
  followCurrentRow: boolean;
  sketchButton: HTMLElement;
  sketchedMolecule: string = '';
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  cutoff: number;
  searchAsYouSketch: boolean;
  sketcherDebounceMs: number;
  targetMoleculeIdx: number = 0;
  // --- select / filter controls ---
  similarSetBitset: DG.BitSet | null = null;
  filterActive: boolean = false;
  selectBtn: HTMLElement | null = null;
  filterBtn: HTMLElement | null = null;
  private _filterSub: Subscription | null = null;
  private _molObserver: IntersectionObserver | null = null;
  private _pushedSummary: string | null = null;
  private _skipContribution: boolean = false;
  private _suppressTimer: ReturnType<typeof setTimeout> | null = null;
  private _searchGen: number = 0;

  get targetMolecule(): string {
    return this.isEditedFromSketcher ?
      this.sketchedMolecule :
      this.moleculeColumn?.get(this.targetMoleculeIdx);
  }

  constructor() {
    super(SIMILARITY);
    // Search/show within the filtered rows by default so the result cards always respect
    // the active filter. With no filter active, 'Filtered' covers all rows (same as 'All').
    this.rowSource = RowSourceTypes.Filtered;
    // Similarity may show up to MAX_LIMIT_SIMILARITY results (the shared base/diversity cap stays small).
    // Set the property's max here in the constructor so a restored limit is never clamped against a stale cap.
    const limitProp = this.getProperties().find((p) => p.name === 'limit');
    if (limitProp)
      limitProp.max = MAX_LIMIT_SIMILARITY;
    this.cutoff = this.float('cutoff', 0.01, {min: 0, max: 1, category: 'Similarity search'});
    this.followCurrentRow = this.bool('followCurrentRow', true,
      {description: 'Re-compute similarity search when changing current row'});
    this.searchAsYouSketch = this.bool('searchAsYouSketch', true, {
      description: 'Re-run similarity search dynamically while drawing in the sketcher (debounced)',
    });
    this.sketcherDebounceMs = this.int('sketcherDebounceMs', 100, {
      min: 50,
      max: 2000,
      description: 'Debounce delay (ms) for search-as-you-sketch updates',
    });
    this.sketchButton = ui.icons.edit(() => {
      const sketcher = new grok.chem.Sketcher();
      const savedMolecule = this.targetMolecule;
      const savedIsEditedFromSketcher = this.isEditedFromSketcher;
      const savedSketchedMolecule = this.sketchedMolecule;
      sketcher.setMolecule(this.targetMolecule);

      const applySketcherState = (isEdited: boolean, mol: string) => {
        this.isEditedFromSketcher = isEdited;
        this.sketchedMolecule = mol;
        this.gridSelect = false;
        this.render();
      };

      let liveSub: Subscription | null = null;
      if (this.searchAsYouSketch) {
        liveSub = DG.debounce(sketcher.onChanged, this.sketcherDebounceMs).subscribe(() => {
          const mol = sketcher.getMolFile();
          if (DG.chem.Sketcher.isEmptyMolfile(mol))
            return;
          applySketcherState(true, mol);
        });
      }

      const dialog = ui.dialog()
        .add(sketcher.root)
        .onOK(() => {
          const editedMolecule = sketcher.getMolFile();
          if (DG.chem.Sketcher.isEmptyMolfile(editedMolecule)) {
            grok.shell.error(`Empty molecule cannot be used for similarity search`);
            this.isEditedFromSketcher = true;
            this.sketchedMolecule = savedMolecule;
          } else
            applySketcherState(true, editedMolecule);
        })
        .onCancel(() => applySketcherState(savedIsEditedFromSketcher, savedSketchedMolecule));
      dialog.show();
      dialog.onClose.subscribe(() => liveSub?.unsubscribe());
    }, 'Edit');
    this.sketchButton.classList.add('chem-similarity-search-edit');
    this.sketchButton.classList.add('chem-mol-view-icon');

    // Select all shown (the similar set, i.e. rows with score >= cutoff). Deselect = Escape (platform built-in).
    this.selectBtn = ui.iconFA('th-large', () => {
      if (!this.similarSetBitset || !this.dataFrame ||
        this.similarSetBitset.trueCount === 0 ||
        this.similarSetBitset.length !== this.dataFrame.filter.length)
        return; // nothing similar → don't wipe the user's existing selection with an all-false set
      // Respect the active filter: select only similar molecules that are currently visible
      // (similar set AND the table's filter), not similar molecules hidden by other filters.
      const sel = DG.BitSet.create(this.dataFrame.rowCount).copyFrom(this.similarSetBitset);
      sel.and(this.dataFrame.filter, false);
      this.dataFrame.selection.copyFrom(sel);
    }, 'Select all similar');
    this.selectBtn.classList.add('chem-similarity-action-btn');

    // Collaborative filter toggle (off by default).
    this.filterBtn = ui.iconFA('filter', () => this._toggleFilter(),
      'Filter table to similar set (toggle)');
    this.filterBtn.classList.add('chem-similarity-action-btn');

    // Accessibility: make the icon <i> elements behave as buttons (focusable + keyboard-activatable).
    for (const b of [this.selectBtn, this.filterBtn]) {
      b.setAttribute('role', 'button');
      b.setAttribute('tabindex', '0');
      b.addEventListener('keydown', (e: KeyboardEvent) => {
        if (e.key === 'Enter' || e.key === ' ') {
          e.preventDefault();
          b.click();
        }
      });
    }
    this.filterBtn.setAttribute('aria-pressed', 'false');

    this.metricsDiv!.appendChild(ui.divH([this.selectBtn, this.filterBtn], 'chem-similarity-action-group'));
    this.updateMetricsLink(this, {});
  }

  init(): void {
    this.isEditedFromSketcher = false;
    this.followCurrentRow = true;
    this.initialized = true;
  }

  isReferenceMolecule(idx: number): boolean {
    return idx === this.targetMoleculeIdx && !this.isEditedFromSketcher;
  }

  private _filterSummary(): string {
    return `Most similar structures: ${this.distanceMetric} ≥ ${this.cutoff}`;
  }

  /** Builds the "similar set": rows within the current row source whose similarity
   * score meets the cutoff. Malformed rows carry a 100.0 sentinel and are excluded. */
  private _buildSimilarSetBitset(allDistances: number[], rowSourceIdxs: DG.BitSet): DG.BitSet {
    const cutoff = this.cutoff;
    // A sketched query is not a table row, so it excludes nothing; otherwise exclude the reference row.
    const refIdx = this.isEditedFromSketcher ? -1 : this.targetMoleculeIdx;
    // Read the row source as a packed buffer once (one interop call) and test bits in JS, rather than
    // one rowSourceIdxs.get(i) interop per row — matters on 50k–100k tables. Single round-trip build.
    const srcBuf = rowSourceIdxs.getBuffer();
    return DG.BitSet.create(this.dataFrame!.rowCount, (i) =>
      i !== refIdx && ((srcBuf[i >>> 5] & (1 << (i & 31))) !== 0) &&
      allDistances[i] >= cutoff && allDistances[i] <= 1);
  }

  /** Single shared, cancellable suppress window. Both the probe and the re-apply call this, so the
   * one timer is rescheduled by the last requestFilter — the flag is never cleared mid-search, and
   * the two callers never race to clear each other's window. 300 ms > the 50 ms onFilterChanged debounce. */
  private _holdSuppress(): void {
    this._suppressFilterRender = true;
    if (this._suppressTimer)
      clearTimeout(this._suppressTimer);
    this._suppressTimer = setTimeout(() => {
      this._suppressFilterRender = false;
      this._suppressTimer = null;
    }, 300);
  }

  /** Re-runs the filter cycle while suppressing the viewer's own onFilterChanged re-render,
   * so applying our own mask does not trigger a recompute loop. */
  private _reapplyFilter(): void {
    if (!this.dataFrame)
      return;
    this._holdSuppress();
    this.dataFrame.rows.requestFilter();
  }

  private _removePushedSummary(): void {
    if (this._pushedSummary && this.dataFrame && this.dataFrame.rows.filters.includes(this._pushedSummary))
      this.dataFrame.rows.filters.remove(this._pushedSummary);
    this._pushedSummary = null;
  }

  private _toggleFilter(): void {
    if (!this.dataFrame || !this.similarSetBitset)
      return;
    const df = this.dataFrame;
    this.filterActive = !this.filterActive;
    this.filterBtn?.classList.toggle('active', this.filterActive);
    this.filterBtn?.setAttribute('aria-pressed', String(this.filterActive));
    if (this.filterActive) {
      this._filterSub = df.onRowsFiltering.subscribe(() => {
        // Skip our contribution during a "probe" cycle used to read the other filters' result
        // (no AND, no label push — so the probe emits no intermediate "others-only" filter state).
        if (this._skipContribution)
          return;
        // Don't hide everything if there are no similar molecules — only AND a non-empty mask.
        if (this.similarSetBitset && this.similarSetBitset.length === df.filter.length &&
          this.similarSetBitset.trueCount > 0) {
          df.filter.and(this.similarSetBitset, false);
          // Surface the filter in the filter panel (replace a stale summary when the cutoff changes).
          const s = this._filterSummary();
          if (this._pushedSummary !== s)
            this._removePushedSummary();
          if (!df.rows.filters.includes(s))
            df.rows.filters.push(s);
          this._pushedSummary = s;
        }
      });
      this._reapplyFilter();
    } else {
      this._filterSub?.unsubscribe();
      this._filterSub = null;
      this._removePushedSummary();
      // Re-filter WITHOUT suppression so the viewer re-renders against the now-wider set.
      df.rows.requestFilter();
    }
  }

  /** Reads the result of all OTHER filters (excluding this viewer's own contribution) by running
   * one probe filter cycle with our contribution skipped, then restoring it. requestFilter runs the
   * onRowsFiltering handlers synchronously, so df.filter reflects each pass before the call returns,
   * and the restore happens in the same tick — so the grid never repaints the intermediate state.
   * This is order-independent (unlike snapshotting before our own AND in the handler). */
  private _computeOtherFiltersMask(): DG.BitSet {
    const df = this.dataFrame;
    this._holdSuppress();
    let others: DG.BitSet;
    this._skipContribution = true;
    try {
      df.rows.requestFilter(); // df.filter = result of all other filters (our AND skipped)
      others = DG.BitSet.create(df.rowCount).copyFrom(df.filter);
    } finally {
      // Always restore the flag, even if a third-party filter handler throws, so our mask is
      // never permanently dropped (which would leave the toggle on but filtering nothing).
      this._skipContribution = false;
    }
    df.rows.requestFilter(); // restore df.filter = others AND ours
    this._holdSuppress(); // reschedule the single window to also cover the restore cycle
    return others;
  }

  /** Nearest scrollable ancestor, used as the IntersectionObserver root for lazy molecule rendering
   * (the viewer scrolls inside the context panel, not the document viewport). */
  private _findScrollParent(el: HTMLElement): HTMLElement | null {
    let p: HTMLElement | null = el.parentElement;
    while (p) {
      const oy = getComputedStyle(p).overflowY;
      if (oy === 'auto' || oy === 'scroll')
        return p;
      p = p.parentElement;
    }
    return null;
  }

  override getRowSourceIndexes(): DG.BitSet {
    // While our own filter is active, search against the OTHER filters' result so changing the
    // reference row does not progressively shrink the set (order-independent probe, not a snapshot).
    if (this.filterActive && this.rowSource === RowSourceTypes.Filtered)
      return this._computeOtherFiltersMask();
    return super.getRowSourceIndexes();
  }

  override async onTableAttached(): Promise<void> {
    // Keep the modest default limit — do NOT auto-inflate it to the whole table. The render loop
    // builds DOM nodes per result up front, so a huge limit (with a low cutoff) would freeze the UI.
    // Select & filter operate on the full above-cutoff set regardless of how many cards are shown,
    // so a small display limit loses nothing functionally. Only clamp a too-large (e.g. restored
    // from a saved project) value down to the table size.
    if (this.dataFrame && this.limit > this.dataFrame.rowCount)
      this.limit = Math.max(1, this.dataFrame.rowCount);
    await super.onTableAttached();
  }

  override detach(): void {
    this._molObserver?.disconnect();
    this._filterSub?.unsubscribe();
    if (this._suppressTimer)
      clearTimeout(this._suppressTimer);
    this._suppressTimer = null;
    this._suppressFilterRender = false; // don't carry a stale suppress into a later re-attach
    this._removePushedSummary();
    // If we were filtering, re-run the filter cycle (our listener is gone now) so the table is
    // restored instead of being left filtered with no viewer or label to explain why.
    const wasFiltering = this.filterActive;
    this.filterActive = false;
    // Reset the toggle's visual state in case the platform reuses this instance for another table.
    this.filterBtn?.classList.remove('active');
    this.filterBtn?.setAttribute('aria-pressed', 'false');
    if (wasFiltering && this.dataFrame)
      this.dataFrame.rows.requestFilter();
    super.detach();
  }

  async renderInternal(computeData: boolean): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.moleculeColumn && this.dataFrame) {
      if (this.moleculeColumn.type !== DG.TYPE.STRING) {
        this.closeWithError('Incorrect target column type');
        return;
      }
      let progressBar: DG.TaskBarProgressIndicator | null = null;
      this.curIdx = this.dataFrame.currentRowIdx === -1 ? 0 : this.dataFrame.currentRowIdx;
      if (computeData && (!this.gridSelect && this.followCurrentRow || this.isEditedFromSketcher)) {
        progressBar = DG.TaskBarProgressIndicator.create(`Similarity search running...`);
        this.isComputing = true;
        this.error = '';
        this.root.classList.remove(`chem-malformed-molecule-error`);
        this.targetMoleculeIdx = this.dataFrame.currentRowIdx === -1 ? 0 : this.dataFrame.currentRowIdx;
        if (DG.chem.Sketcher.isEmptyMolfile(this.targetMolecule)) {
          this.closeWithError(`Empty molecule cannot be used for similarity search`, progressBar);
          return;
        }
        const searchGen = ++this._searchGen;
        try {
          const rowSourceIdxs = this.getRowSourceIndexes();
          rowSourceIdxs.set(this.targetMoleculeIdx, true);
          const result = await chemSimilaritySearchEx(this.dataFrame!, this.moleculeColumn!,
            this.targetMolecule, this.distanceMetric as BitArrayMetrics, this.limit, this.cutoff,
            this.fingerprint as Fingerprint, rowSourceIdxs);
          if (searchGen !== this._searchGen)
            return; // a newer search started while this one ran — let it own the result/render
          if (!result) {
            this.closeWithError(`Malformed molecule cannot be used for similarity search`); // finally closes the progress bar
            return;
          }
          this.molCol = result.df.getCol('smiles');
          this.idxs = result.df.getCol('indexes');
          this.scores = result.df.getCol('score');
          this.similarSetBitset = this._buildSimilarSetBitset(result.allDistances, rowSourceIdxs);
          if (this.filterActive)
            this._reapplyFilter();
        } catch (e: any) {
          grok.shell.error(e.message);
          return;
        } finally {
          progressBar?.close();
        }
      } else if (this.gridSelect)
        this.gridSelect = false;
      if (this.error) {
        this.closeWithError(this.error, progressBar);
        return;
      }
      this.clearResults();
      // Grey out the action icons until a search has produced a non-empty result set.
      const noSimilar = !this.similarSetBitset || this.similarSetBitset.trueCount === 0;
      this.selectBtn?.classList.toggle('chem-similarity-action-disabled', noSimilar);
      this.filterBtn?.classList.toggle('chem-similarity-action-disabled', noSimilar);
      const panel = [];
      const grids = [];
      const lazyHosts: {host: HTMLElement; molStr: string}[] = [];
      let cnt = 0; let cnt2 = 0;
      panel[cnt++] = this.metricsDiv;
      if (this.molCol && this.idxs && this.scores) {
        if (this.isEditedFromSketcher) {
          const label = this.sketchButton;
          const grid = ui.div([
            renderMolecule(
              this.targetMolecule, { width: this.sizesMap[this.size].width, height: this.sizesMap[this.size].height }),
            label], { style: { position: 'relative' } });
          let divClass = 'd4-flex-col';
          divClass += ' d4-current';
          grid.style.boxShadow = '0px 0px 1px var(--grey-6)';
          $(grid).addClass(divClass);
          grids[cnt2++] = grid;
        }
        const displayCount = Math.min(this.molCol.length, DISPLAY_CAP);
        for (let i = 0; i < displayCount; ++i) {
          const idx = this.idxs.get(i);
          const similarity = this.scores.get(i).toPrecision(2);
          const refMolecule = this.isReferenceMolecule(idx);
          const label = refMolecule ? this.sketchButton : ui.div();
          const molProps = this.createMoleculePropertiesDiv(idx, refMolecule, similarity);
          const molHost = ui.div([]);
          molHost.style.width = `${this.sizesMap[this.size].width}px`;
          molHost.style.height = `${this.sizesMap[this.size].height}px`;
          lazyHosts.push({host: molHost, molStr: this.molCol!.get(i)});
          const grid = ui.div([
            molHost,
            label,
            molProps], { style: { position: 'relative' } });
          let divClass = 'd4-flex-col';
          if (idx === this.curIdx) {
            divClass += ' d4-current';
            grid.style.backgroundColor = '#ddffd9';
          }
          if (idx === this.targetMoleculeIdx && !this.isEditedFromSketcher) {
            divClass += ' d4-current';
            grid.style.boxShadow = '0px 0px 1px var(--grey-6)';
          }
          if (this.dataFrame?.selection.get(idx)) {
            divClass += ' d4-selected';
            if (divClass === 'd4-flex-col d4-selected')
              grid.style.backgroundColor = '#f8f8df';
            else
              grid.style.backgroundColor = '#d3f8bd';
          }
          $(grid).addClass(divClass);
          grid.addEventListener('click', (event: MouseEvent) => {
            if (this.dataFrame && this.idxs) {
              if (event.shiftKey || event.altKey)
                this.dataFrame.selection.set(idx, true);
              else if (event.metaKey) {
                const selected = this.dataFrame.selection;
                this.dataFrame.selection.set(idx, !selected.get(idx));
              } else {
                this.dataFrame.currentRowIdx = idx;
                this.gridSelect = true;
              }
            }
          });
          grids[cnt2++] = grid;
        }
        if (this.molCol.length > DISPLAY_CAP) {
          grids[cnt2++] = ui.divText(`+ ${this.molCol.length - DISPLAY_CAP} more — use Select / Filter to act on all`,
            {style: {alignSelf: 'center', padding: '0 8px', color: 'var(--grey-5)', fontSize: '11px'}});
        }
      }
      const gridDiv = ui.divH(grids, 'chem-viewer-grid');
      panel[cnt++] = gridDiv;
      this.root.appendChild(ui.panel([ui.divV(panel)]));
      // Lazy-render molecule structures: only draw a card when it scrolls into view, rooted at the
      // viewer's actual scroll container (the context panel), so the deferral works on panel-scroll.
      // If no scrollable ancestor is found, fall back to the viewport with a 0px margin (a 300px
      // margin against the viewport would mark every panel card as intersecting and defeat lazy load).
      const scrollRoot = this._findScrollParent(gridDiv);
      this._molObserver?.disconnect();
      this._molObserver = new IntersectionObserver((entries) => {
        for (const e of entries) {
          if (!e.isIntersecting)
            continue;
          const host = e.target as HTMLElement & {_molStr?: string; _rendered?: boolean};
          this._molObserver!.unobserve(host);
          if (host._rendered)
            continue;
          host._rendered = true;
          host.appendChild(renderMolecule(host._molStr!,
            {width: this.sizesMap[this.size].width, height: this.sizesMap[this.size].height}));
        }
      }, {root: scrollRoot, rootMargin: scrollRoot ? '300px' : '0px'});
      for (const lh of lazyHosts) {
        (lh.host as HTMLElement & {_molStr?: string})._molStr = lh.molStr;
        this._molObserver.observe(lh.host);
      }
    }
  }
}

export interface SimilaritySearchResult {
  /** 3-column dataframe (smiles, score, indexes) of the displayed top-N. */
  df: DG.DataFrame;
  /** Similarity score for every source row (length = table.rowCount); malformed rows = 100.0. */
  allDistances: number[];
}

/** Backward-compatible wrapper returning only the display dataframe. */
export async function chemSimilaritySearch(
  table: DG.DataFrame,
  smiles: DG.Column,
  molecule: string,
  metricName: BitArrayMetrics,
  limit: number,
  minScore: number,
  fingerprint: Fingerprint,
  rowSourceIndexes: DG.BitSet,
) : Promise<DG.DataFrame | null> {
  const res = await chemSimilaritySearchEx(table, smiles, molecule, metricName, limit, minScore,
    fingerprint, rowSourceIndexes);
  return res ? res.df : null;
}

export async function chemSimilaritySearchEx(
  table: DG.DataFrame,
  smiles: DG.Column,
  molecule: string,
  metricName: BitArrayMetrics,
  limit: number,
  minScore: number,
  fingerprint: Fingerprint,
  rowSourceIndexes: DG.BitSet,
) : Promise<SimilaritySearchResult | null> {
  const targetFingerprint = chemSearches.chemGetFingerprint(molecule, fingerprint, () => {return null;});
  if (!targetFingerprint)
    return null; //returning null in case target molecule is malformed
  const fingerprintCol = await chemSearches.chemGetFingerprints(smiles, fingerprint, false);
  malformedDataWarning(fingerprintCol, smiles);
  const distances: number[] = [];

  const fpSim = similarityMetric[metricName];
  for (let row = 0; row < fingerprintCol.length; row++) {
    const fp = fingerprintCol[row];
    distances[row] = (!fp || fp!.allFalse) ? 100.0 : fpSim(targetFingerprint, fp!);
  }

  function range(end: number) {
    return Array(end).fill(0).map((_, idx) => idx);
  }

  function compare(i1: number, i2: number) {
    if (distances[i1] > distances[i2])
      return -1;

    if (distances[i1] < distances[i2])
      return 1;

    return 0;
  }

  const indexes = range(table.rowCount)
    .filter((idx) => fingerprintCol[idx] && !fingerprintCol[idx]!.allFalse && rowSourceIndexes.get(idx))
    .sort(compare);
  const molsList = [];
  const scoresList = [];
  const molsIdxs = [];
  limit = Math.min(indexes.length, limit);
  for (let n = 0; n < limit; n++) {
    const idx = indexes[n];
    const score = distances[idx];
    if (score < minScore)
      break;

    molsIdxs[n] = idx;
    molsList[n] = smiles.get(idx);
    scoresList[n] = score;
  }
  const mols = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'smiles', molsList);
  mols.semType = DG.SEMTYPE.MOLECULE;
  const scores = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'score', scoresList);
  const newIndexes = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'indexes', molsIdxs);
  return {df: DG.DataFrame.fromColumns([mols, scores, newIndexes]), allDistances: distances};
}

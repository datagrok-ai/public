import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as chemSearches from '../chem-searches';
import {similarityMetric} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import $ from 'cash-dom';
import {Fingerprint} from '../utils/chem-common';
import {renderMolecule} from '../rendering/render-molecule';
import {ChemSearchBaseViewer, SIMILARITY, RowSourceTypes} from './chem-search-base-viewer';
import {getSearchProgressEventName} from '../constants';
import {malformedDataWarning} from '../utils/malformed-data-utils';
import '../../css/chem.css';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {Subscription} from 'rxjs';

/** Display cap; the shared base default and diversity viewer (O(n^2) selection) keep a smaller cap. */
const MAX_LIMIT_SIMILARITY = 10000;
/** Hard cap on cards built into the DOM, independent of `limit` (Select/filter still act on the full set). */
const DISPLAY_CAP = 200;
/** Default similar molecules shown (base's 12 is too small — this viewer is cutoff-driven). */
const DEFAULT_LIMIT_SIMILARITY = 200;
/** Prefix of the filter-panel label this viewer pushes; used to find/remove its own label. */
const FILTER_LABEL_PREFIX = 'Most similar structures:';
/** Suppress window for our own requestFilter cycles; must exceed the base's 50 ms debounce. */
const SUPPRESS_WINDOW_MS = 300;
/** Debounce coalescing a burst of filter cycles into a single re-search once activity settles. */
const FILTER_RECHECK_DEBOUNCE_MS = 300;
/** Watchdog: un-block the recheck after this quiet period if a substructure search dies before 100%. */
const EXTERNAL_SEARCH_WATCHDOG_MS = 2000;
/** Dataframes a viewer is currently probing, so a SECOND viewer doesn't mistake a sibling's probe for an
 * external change. Guards only the PROBE: with two active similarity filters a _reapplyFilter() still
 * fires the sibling's onRowsFiltering, which is correct — _sameMask dedupe converges in ~1–2 rounds. */
const _probingDataframes = new WeakSet<DG.DataFrame>();

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
  /** Cutoff/metric the current similarSetBitset was built with — the label reads these (not live values)
   * so a fast slider drag can't show a threshold the mask doesn't match. */
  private _bitsetCutoff: number = 0;
  private _bitsetMetric: string = '';
  filterActive: boolean = false;
  selectBtn!: HTMLElement; // assigned in the constructor
  filterBtn!: HTMLElement;
  private _filterSub: Subscription | null = null;
  private _molObserver: IntersectionObserver | null = null;
  // Snapshot of the last full card render, so a non-recompute render (selection/current-row/metadata)
  // refreshes per-card highlights IN PLACE instead of rebuilding the panel (which flickers every card).
  private _renderedCards: {grid: HTMLElement; idx: number}[] = [];
  private _renderedLazyHosts: {host: HTMLElement; molStr: string}[] = [];
  private _renderedGridDiv: HTMLElement | null = null;
  private _renderedMolCol: DG.Column | null = null;
  private _renderedSize: string = '';
  private _pushedSummary: string | null = null;
  /** Subscription to the chem substructure-search progress event for our molecule column. */
  private _searchProgressSub: Subscription | null = null;
  /** Live while a substructure search streams on our column (it calls requestFilter() per batch). We
   * defer and run a single search at 100% instead of one per batch; also a watchdog if it dies early. */
  private _externalSearchTimer: ReturnType<typeof setTimeout> | null = null;
  private _skipContribution: boolean = false;
  /** True only during the viewer's OWN synchronous requestFilter calls (probe/reapply), so the handler
   * tells its cycles from an external change precisely (a clear mid-processing isn't dropped as ours). */
  private _ownFilterCycle: boolean = false;
  /** Open while our OWN requestFilter cycles run so onExternalFilterChanged skips them (no feedback
   * loop). Time-windowed via _suppressTimer; see _holdSuppress. */
  private _suppressFilterRender: boolean = false;
  private _suppressTimer: ReturnType<typeof setTimeout> | null = null;
  private _searchGen: number = 0;
  /** Row source the current similarSetBitset was last searched against — lets a filter change be detected
   * even when net df.filter is unchanged (e.g. clearing a filter). */
  private _builtAgainstFilter: DG.BitSet | null = null;
  private _filterRecheckTimer: ReturnType<typeof setTimeout> | null = null;
  /** Hand-off of the mask the recheck just probed, so the render it triggers reuses it instead of
   * probing again. Valid only for the immediately following render (see _scheduleFilterRecheck). */
  private _recheckRowSource: DG.BitSet | null = null;

  get targetMolecule(): string {
    return this.isEditedFromSketcher ?
      this.sketchedMolecule :
      this.moleculeColumn?.get(this.targetMoleculeIdx);
  }

  constructor() {
    super(SIMILARITY);
    // Per-viewer property overrides: raise the limit cap, default row source to 'Filtered'. Done here (not
    // onTableAttached) so a restored value isn't clamped against a stale max.
    const limitProp = this.getProperty('limit');
    if (limitProp) {
      limitProp.max = MAX_LIMIT_SIMILARITY;
      limitProp.defaultValue = DEFAULT_LIMIT_SIMILARITY; // 'Reset to default' / project default → 200
    }
    this.limit = DEFAULT_LIMIT_SIMILARITY;
    const rowSourceProp = this.getProperty('rowSource');
    if (rowSourceProp)
      rowSourceProp.defaultValue = RowSourceTypes.Filtered;
    this.rowSource = RowSourceTypes.Filtered;
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

    // Select all shown (score >= cutoff); deselect = Escape. Reuse the grid "select all" icon to match.
    this.selectBtn = ui.iconSvg('select-all', () => {
      // Guard empty set + length mismatch (BitSet.and requires equal length).
      if (!this.dataFrame || this.similarCount === 0 || this.similarSetBitset!.length !== this.dataFrame.filter.length)
        return;
      // Select only similar rows passing the active filter (similar set ∩ table filter).
      const sel = this.similarSetBitset!.clone().and(this.dataFrame.filter, false);
      if (sel.trueCount === 0)
        return; // every similar row is hidden by another filter — don't wipe the existing selection
      const changed = !this._sameMask(sel, this.dataFrame.selection); // skip the toast on a repeat click
      this.dataFrame.selection.copyFrom(sel);
      // Confirm the off-panel action — the selection happens in the grid, away from the click.
      if (changed)
        grok.shell.info(`${sel.trueCount} similar ${sel.trueCount === 1 ? 'row' : 'rows'} selected`);
    }, 'Select all similar');
    this.selectBtn.classList.add('chem-similarity-action-btn');

    // Collaborative filter toggle (off by default).
    this.filterBtn = ui.icons.filter(() => this._toggleFilter(),
      'Filter table to similar set (toggle)');
    this.filterBtn.classList.add('chem-similarity-action-btn');

    // Accessibility: explicit aria-label up front (renderInternal later adds the count, _toggleFilter the
    // on/off state) so a screen reader never sees an unlabelled button; then make the <i> behave as buttons.
    this.selectBtn.setAttribute('aria-label', 'Select all similar');
    this.filterBtn.setAttribute('aria-label', 'Filter table to similar set');
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

  override onPropertyChanged(property: DG.Property): void {
    // setOptions bypasses the slider's [0, 1] bounds; clamp so a stray cutoff can't mark every (0) or
    // nothing (1) similar. (limit is clamped in the base.)
    if (property?.name === 'cutoff')
      this.cutoff = Math.max(0, Math.min(1, this.cutoff));
    super.onPropertyChanged(property);
    // The progress event is column-scoped — re-subscribe when the molecule column changes.
    if (property?.name === 'moleculeColumnName')
      this._subscribeSearchProgress();
  }

  isReferenceMolecule(idx: number): boolean {
    return idx === this.targetMoleculeIdx && !this.isEditedFromSketcher;
  }

  /** Rows in the current similar set (0 when none meet the cutoff or no search has run). */
  private get similarCount(): number {
    return this.similarSetBitset?.trueCount ?? 0;
  }

  private _filterSummary(): string {
    // Read the cutoff/metric the mask was BUILT with (not live) so a mid-drag label matches the applied
    // mask. Round so 0.3 doesn't render as 0.30000000000000004.
    const cutoff = Number(this._bitsetCutoff.toFixed(3));
    return `${FILTER_LABEL_PREFIX} ${this._bitsetMetric} ≥ ${cutoff}`;
  }

  /** Builds the "similar set": rows in the row source whose score meets the cutoff (malformed rows carry
   * a 100.0 sentinel and are excluded). */
  private _buildSimilarSetBitset(allDistances: number[], rowSourceIdxs: DG.BitSet, cutoff: number): DG.BitSet {
    // A sketched query isn't a table row (excludes nothing); otherwise exclude the reference row.
    const refIdx = this.isEditedFromSketcher ? -1 : this.targetMoleculeIdx;
    // Read the row source as a packed buffer once and test bits in JS — avoids a get(i) interop per row.
    const srcBuf = rowSourceIdxs.getBuffer();
    // Bound by scored rows: a row added during a fast (cached) search has no score yet — exclude it
    // (rowCount can exceed allDistances.length; don't read OOB).
    const scored = allDistances.length;
    return DG.BitSet.create(this.dataFrame!.rowCount, (i) =>
      i < scored && i !== refIdx && ((srcBuf[i >>> 5] & (1 << (i & 31))) !== 0) &&
      allDistances[i] >= cutoff && allDistances[i] <= 1);
  }

  /** Open/extend a single suppress window so our own requestFilter cycles don't trigger a re-render.
   * One shared cancellable timer means probe and re-apply never race to clear each other's window. */
  private _holdSuppress(): void {
    this._suppressFilterRender = true;
    if (this._suppressTimer)
      clearTimeout(this._suppressTimer);
    this._suppressTimer = setTimeout(() => {
      this._suppressFilterRender = false;
      this._suppressTimer = null;
    }, SUPPRESS_WINDOW_MS);
  }

  /** Re-run the filter cycle while suppressing our own re-render, so applying our mask doesn't loop. */
  private _reapplyFilter(): void {
    if (!this.dataFrame)
      return;
    this._holdSuppress();
    this._ownFilterCycle = true;
    try {
      this.dataFrame.rows.requestFilter();
    } finally {
      this._ownFilterCycle = false;
    }
  }

  private _removePushedSummary(): void {
    if (this._pushedSummary && this.dataFrame && this.dataFrame.rows.filters.includes(this._pushedSummary))
      this.dataFrame.rows.filters.remove(this._pushedSummary);
    this._pushedSummary = null;
  }

  /** Drop the similar set and its label together — a failed/empty search must leave neither a stale mask
   * (which onRowsFiltering would keep ANDing) nor a label contradicting the error state. */
  private _clearSimilarSet(): void {
    this.similarSetBitset = null;
    this._removePushedSummary();
  }

  private _toggleFilter(): void {
    // Block activation on an empty set, but always allow deactivation so the toggle can't get stuck on.
    if (!this.dataFrame || (!this.filterActive && this.similarCount === 0))
      return;
    const df = this.dataFrame;
    this.filterActive = !this.filterActive;
    this.filterBtn.classList.toggle('active', this.filterActive);
    this.filterBtn.setAttribute('aria-pressed', String(this.filterActive));
    // Distinct labels per state so a screen reader doesn't announce the same thing on and off.
    this.filterBtn.setAttribute('aria-label',
      this.filterActive ? 'Deactivate similarity filter' : 'Filter table to similar set');
    if (this.filterActive) {
      this._filterSub = df.onRowsFiltering.subscribe(() => {
        // Skip our contribution during a probe cycle (reads the other filters' result without ours).
        if (this._skipContribution)
          return;
        // An external cycle (not our probe/reapply) may have changed the population — schedule one recheck.
        // Essential for CLEARING a filter (our subset mask leaves net df.filter unchanged, so onFilterChanged
        // never fires). Gate on _ownFilterCycle (not the suppress window) so a clear isn't dropped; skip
        // while any viewer is probing.
        if (!this._ownFilterCycle && !_probingDataframes.has(df))
          this._scheduleFilterRecheck();
        // Only AND a non-empty mask, so an empty similar set never hides every row.
        if (this.similarCount > 0 && this.similarSetBitset!.length === df.filter.length) {
          df.filter.and(this.similarSetBitset!, false);
          // Surface the filter in the panel (replace a stale summary when the cutoff changes).
          const s = this._filterSummary();
          if (this._pushedSummary !== s)
            this._removePushedSummary();
          if (!df.rows.filters.includes(s))
            df.rows.filters.push(s);
          this._pushedSummary = s;
        } else {
          // Not contributing (empty set, or length mismatch right after a row add/remove before rebuild)
          // — drop the label so the panel never shows a filter that isn't applied.
          this._removePushedSummary();
        }
      });
      this._reapplyFilter();
    } else {
      this._filterSub?.unsubscribe();
      this._filterSub = null;
      // Drop any pending recheck + watchdog: they self-guard on `filterActive`, but cancelling here
      // avoids a spurious search if the user rapidly toggles back on mid-stream.
      if (this._filterRecheckTimer) {
        clearTimeout(this._filterRecheckTimer);
        this._filterRecheckTimer = null;
      }
      if (this._externalSearchTimer) {
        clearTimeout(this._externalSearchTimer);
        this._externalSearchTimer = null;
      }
      this._removePushedSummary();
      // Un-filter the TABLE but suppress the re-render: the similar set is built against the OTHER filters
      // (probe excludes our contribution), so removing our toggle doesn't change which cards show —
      // rebuilding them would just flicker the panel.
      this._holdSuppress();
      df.rows.requestFilter();
    }
  }

  /** Result of all OTHER filters, excluding this viewer's contribution. requestFilter runs onRowsFiltering
   * synchronously, so we probe with our contribution skipped, snapshot df.filter, then restore — all in
   * one tick (no intermediate repaint). Order-independent, unlike snapshotting before our own AND. */
  private _computeOtherFiltersMask(): DG.BitSet {
    const df = this.dataFrame;
    // ASSUMES requestFilter() fires ALL onRowsFiltering handlers synchronously (true for every DG.Filter
    // today). An async future filter would let the probe snapshot the all-true reset state — wrong result.
    this._holdSuppress();
    this._ownFilterCycle = true; // both requestFilter cycles below are ours, not external
    _probingDataframes.add(df!); // tell any sibling viewer these cycles are a probe, not a change
    let others: DG.BitSet;
    this._skipContribution = true;
    try {
      df.rows.requestFilter(); // df.filter = result of all other filters (our AND skipped)
      others = df.filter.clone();
    } finally {
      this._skipContribution = false;
      try {
        df.rows.requestFilter(); // restore df.filter = others AND ours
        this._holdSuppress(); // reschedule the window to also cover the restore cycle
      } finally {
        // Reset the cycle flags FIRST so a throwing re-AND below can't leave the viewer deaf to external
        // changes (the AND uses notify=false — nothing observes the flags in between).
        this._ownFilterCycle = false;
        _probingDataframes.delete(df!);
        // Defensively re-apply our mask: idempotent on the happy path, but recovers our contribution if a
        // third-party handler THREW during the restore (else the table strands on the others-only mask).
        if (this.filterActive && this.similarSetBitset && this.similarSetBitset.length === df!.filter.length)
          df!.filter.and(this.similarSetBitset, false);
      }
    }
    return others;
  }

  /** Nearest scrollable ancestor — the IntersectionObserver root for lazy rendering (the viewer scrolls
   * inside the context panel, not the document viewport). */
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

  /** Set a card's current-row/selected highlight (class + background) from the live selection and
   * reference row. Used both when building cards and refreshing them in place (see renderInternal). */
  private _applyCardStyle(grid: HTMLElement, idx: number): void {
    const isCurRow = idx === this.curIdx;
    const isRef = idx === this.targetMoleculeIdx && !this.isEditedFromSketcher;
    const isSelected = !!this.dataFrame?.selection.get(idx);
    grid.classList.toggle('d4-current', isCurRow || isRef);
    grid.classList.toggle('d4-selected', isSelected);
    grid.style.boxShadow = isRef ? '0px 0px 1px var(--grey-6)' : '';
    grid.style.backgroundColor = isSelected ?
      ((isCurRow || isRef) ? '#d3f8bd' : '#f8f8df') : (isCurRow ? '#ddffd9' : '');
  }

  /** (Re)arm the lazy-render observer over `lazyHosts`, rooted at the scroll container with a 300px
   * prefetch margin (100px against the viewport with no scrollable ancestor). Already-rendered hosts are
   * skipped, so safe to call again on an in-place refresh. */
  private _observeLazyHosts(gridDiv: HTMLElement, lazyHosts: {host: HTMLElement; molStr: string}[]): void {
    const {width, height} = this.sizesMap[this.size];
    const scrollRoot = this._findScrollParent(gridDiv);
    this._molObserver?.disconnect();
    // Close over the local `observer`, not this._molObserver (which the next render reassigns) — so a late
    // callback unobserves on its OWN observer and never touches a newer instance.
    const observer = new IntersectionObserver((entries) => {
      for (const e of entries) {
        if (!e.isIntersecting)
          continue;
        const host = e.target as HTMLElement & {_molStr?: string; _rendered?: boolean};
        observer.unobserve(host);
        if (host._rendered || !host.isConnected)
          continue; // skip detached/already-drawn hosts — no wasted RDKit work after a fast scroll-then-close
        host._rendered = true;
        host.appendChild(renderMolecule(host._molStr!, {width, height}));
      }
    }, {root: scrollRoot, rootMargin: scrollRoot ? '300px' : '100px'});
    this._molObserver = observer;
    for (const lh of lazyHosts) {
      (lh.host as HTMLElement & {_molStr?: string})._molStr = lh.molStr;
      observer.observe(lh.host);
    }
  }

  /** True when two masks are equal — compares packed word buffers directly (no clone), bailing early. */
  private _sameMask(a: DG.BitSet | null, b: DG.BitSet | null): boolean {
    if (!a || !b || a.length !== b.length)
      return false;
    const bufA = a.getBuffer(); const bufB = b.getBuffer();
    for (let i = 0; i < bufA.length; i++) {
      if (bufA[i] !== bufB[i])
        return false;
    }
    return true;
  }

  /** Coalesce filter-driven re-searches: after activity settles, re-run the search at most once, and only
   * if the other-filters result actually changed. Keeps any live filter from re-searching every cycle. */
  private _scheduleFilterRecheck(): void {
    // A substructure search is streaming (requestFilter per batch) — its progress handler triggers a
    // single recheck at 100%. The watchdog timer is live only mid-stream, so its presence = "in progress".
    if (this._externalSearchTimer)
      return;
    if (this._filterRecheckTimer)
      clearTimeout(this._filterRecheckTimer);
    this._filterRecheckTimer = setTimeout(() => {
      this._filterRecheckTimer = null;
      this._recheckRowSource = null; // discard any stale hand-off
      if (!this.filterActive || !this.dataFrame)
        return;
      const filterBased = this.rowSource === RowSourceTypes.Filtered ||
        this.rowSource === RowSourceTypes.FilteredSelected;
      // Filter-based modes depend on the other filters — skip the expensive search when unchanged. Other
      // modes can't be altered by a filter change, so just re-render.
      if (filterBased) {
        const current = this.getRowSourceIndexes(); // runs the probe ONCE
        if (this._sameMask(current, this._builtAgainstFilter))
          return;
        // Hand the just-probed mask to the render so it doesn't probe again (consumed before render's
        // first await; the post-clear drops it on a no-recompute path so it can't leak).
        this._recheckRowSource = current;
        this.render(true);
        this._recheckRowSource = null;
      } else {
        // All/Selected: filter-independent set — re-render (refresh highlights) WITHOUT an O(N) search.
        this.render(false);
      }
    }, FILTER_RECHECK_DEBOUNCE_MS);
  }

  /** Subscribe to the substructure-search progress event for the current molecule column. Such a filter
   * calls requestFilter() per batch; while streaming we hold off the recheck and fire exactly one search
   * at 100% (or after a quiet period if it was terminated early). */
  private _subscribeSearchProgress(): void {
    this._searchProgressSub?.unsubscribe();
    this._searchProgressSub = null;
    if (this._externalSearchTimer) {
      clearTimeout(this._externalSearchTimer);
      this._externalSearchTimer = null;
    }
    if (!this.dataFrame || !this.moleculeColumn)
      return;
    const eventName = getSearchProgressEventName(this.dataFrame.name, this.moleculeColumn.name);
    this._searchProgressSub = grok.events.onCustomEvent(eventName).subscribe((progress: number) => {
      if (this._externalSearchTimer) {
        clearTimeout(this._externalSearchTimer);
        this._externalSearchTimer = null;
      }
      if (progress < 100) {
        // Search in progress: suppress the recheck and drop any batch-triggered pending one.
        if (this._filterRecheckTimer) {
          clearTimeout(this._filterRecheckTimer);
          this._filterRecheckTimer = null;
        }
        // Watchdog: if terminated before 100%, recover after a quiet period so the recheck isn't blocked
        // forever. While set, _scheduleFilterRecheck stays gated.
        this._externalSearchTimer = setTimeout(() => {
          this._externalSearchTimer = null;
          if (this.filterActive)
            this._scheduleFilterRecheck();
        }, EXTERNAL_SEARCH_WATCHDOG_MS);
      } else {
        // Search complete (timer cleared above): run one similarity search against the settled set.
        if (this.filterActive)
          this._scheduleFilterRecheck();
      }
    });
  }

  /** Route filter-change re-renders through the debounced, deduped recheck. One of TWO paths that schedule
   * it — onRowsFiltering is the other (and the only one that fires when net df.filter is unchanged, e.g.
   * clearing). Load-bearing redundancy: the suppress window can swallow this path for a genuine external
   * change, and onRowsFiltering covers that. Don't remove either without re-checking clear-and-re-expand. */
  override async onExternalFilterChanged(): Promise<void> {
    // Our own requestFilter cycles open a short suppress window so they don't loop back as a re-render.
    if (this._suppressFilterRender)
      return;
    if (this.filterActive)
      this._scheduleFilterRecheck();
    else
      await super.onExternalFilterChanged();
  }

  override getRowSourceIndexes(): DG.BitSet {
    // PLATFORM CONTRACT: _computeOtherFiltersMask assumes requestFilter() runs every onRowsFiltering
    // subscriber SYNCHRONOUSLY (true for all current DG filters). An async future filter would snapshot
    // the reset state → wrong result. Reuse the mask the pending recheck just probed instead of re-probing.
    if (this._recheckRowSource) {
      const m = this._recheckRowSource;
      this._recheckRowSource = null;
      return m;
    }
    // Search against the OTHER filters' result so changing the reference row doesn't progressively shrink
    // the set. Covers both filter-based sources; FilteredSelected also intersects the selection.
    if (this.filterActive &&
      (this.rowSource === RowSourceTypes.Filtered || this.rowSource === RowSourceTypes.FilteredSelected)) {
      const others = this._computeOtherFiltersMask();
      return this.rowSource === RowSourceTypes.FilteredSelected ?
        others.and(this.dataFrame!.selection, false) : others;
    }
    return super.getRowSourceIndexes();
  }

  /** Reset all per-table filter-runtime state (incl. the three recheck-gating timers). Shared by detach()
   * and onTableAttached()'s re-attach guard so a field added to one path isn't forgotten in the other.
   * Callers handle their own extras (progress/observer subs, requestFilter restore, limit clamp, etc.). */
  private _resetFilterState(): void {
    this._searchGen++; // abandon any in-flight search so it can't write state after the reset
    this._filterSub?.unsubscribe();
    this._filterSub = null;
    this.filterActive = false;
    // Drop the previous table's similar set: a stale trueCount could slip past the toggle's activation
    // guard during the async attach-render, briefly showing an "active" filter that ANDs nothing.
    this.similarSetBitset = null;
    // Drop the cached card snapshot so the next render full-rebuilds (and we don't hold detached DOM refs).
    this._renderedCards = [];
    this._renderedLazyHosts = [];
    this._renderedGridDiv = null;
    this._renderedMolCol = null;
    this._builtAgainstFilter = null;
    this._recheckRowSource = null;
    this._bitsetCutoff = 0;
    this._bitsetMetric = '';
    this._ownFilterCycle = false;
    this._suppressFilterRender = false;
    if (this._suppressTimer) {
      clearTimeout(this._suppressTimer);
      this._suppressTimer = null;
    }
    if (this._filterRecheckTimer) {
      clearTimeout(this._filterRecheckTimer);
      this._filterRecheckTimer = null;
    }
    // The watchdog also gates _scheduleFilterRecheck, so clear it too (its progress SUBSCRIPTION is torn
    // down separately, by detach / _subscribeSearchProgress).
    if (this._externalSearchTimer) {
      clearTimeout(this._externalSearchTimer);
      this._externalSearchTimer = null;
    }
    this._removePushedSummary();
    this.filterBtn.classList.remove('active');
    this.filterBtn.setAttribute('aria-pressed', 'false');
    // Reset the label too, so a screen reader doesn't announce "Deactivate…" on a now-off button.
    this.filterBtn.setAttribute('aria-label', 'Filter table to similar set');
  }

  override async onTableAttached(): Promise<void> {
    // Idempotent re-attach guard: if the platform reused this instance WITHOUT a detach, drop any stale
    // filter sub/state from the previous dataframe (no double-subscribe, no carried-over toggle). A no-op
    // in the normal detach-then-attach path.
    this._resetFilterState();
    // Keep the modest default limit (select/filter act on the full above-cutoff set); only clamp a
    // too-large restored value down to the table size.
    if (this.dataFrame && this.limit > this.dataFrame.rowCount)
      this.limit = Math.max(1, this.dataFrame.rowCount);
    // The toggle is session-only — strip any of our labels left in the filter panel by a restored project
    // so it never lingers without a live filter behind it.
    if (this.dataFrame) {
      for (const s of [...this.dataFrame.rows.filters]) {
        if (s.startsWith(FILTER_LABEL_PREFIX))
          this.dataFrame.rows.filters.remove(s);
      }
    }
    // Match the other Chem filters: "Reset all filters" must clear ours too. Register AFTER super (its
    // first act is `this.subs = []`, which would drop a sub pushed before it), and in `finally` so a
    // throwing first render can't skip it. The base's subs-reset prevents accumulation on re-attach.
    try {
      await super.onTableAttached();
    } finally {
      this.subs.push(grok.events.onResetFilterRequest.subscribe(() => {
        // onResetFilterRequest is GLOBAL (no df payload), so it also fires for other tables' "Reset all
        // filters". Only honour it when our table is in view, so resetting table B doesn't switch off A.
        if (this.filterActive && this.dataFrame === grok.shell.tv?.dataFrame)
          this._toggleFilter();
      }));
      // Re-render on row add (similarity only — out of the base to avoid an O(n^2) diversity recompute):
      // a mask sized to the old row count would silently turn the filter into a no-op. Toggle on → route
      // through the recheck (rebuilds the mask to the new length); else a plain re-render.
      if (this.dataFrame) {
        this.subs.push(DG.debounce(this.dataFrame.onRowsAdded, 50).subscribe(async () => {
          if (this.filterActive)
            this._scheduleFilterRecheck();
          else
            await this.render();
        }));
      }
      // Watch substructure-search progress so a streaming search doesn't trigger one search per batch
      // (super has detected this.moleculeColumn by now).
      this._subscribeSearchProgress();
    }
  }

  override detach(): void {
    const wasFiltering = this.filterActive; // capture before _resetFilterState clears it
    // Teardown beyond the shared filter state: lazy-render observer + progress sub (watchdog timer is
    // cleared by _resetFilterState).
    this._molObserver?.disconnect();
    this._molObserver = null;
    this._searchProgressSub?.unsubscribe();
    this._searchProgressSub = null;
    this._resetFilterState();
    // If we were filtering, re-run the filter cycle (our listener is gone) so the table is restored
    // instead of left filtered with no viewer or label to explain why.
    if (wasFiltering && this.dataFrame)
      this.dataFrame.rows.requestFilter();
    super.detach();
  }

  async renderInternal(computeData: boolean): Promise<void> {
    // Tear down the previous render's lazy observer up front: the early-exit paths below never reach the
    // success-path disconnect, so a prior observer would otherwise keep its detached card hosts alive.
    this._molObserver?.disconnect();
    this._molObserver = null;
    if (!this.beforeRender())
      return;
    if (this.moleculeColumn && this.dataFrame) {
      if (this.moleculeColumn.type !== DG.TYPE.STRING) {
        this.closeWithError('Incorrect target column type');
        return;
      }
      if (this.dataFrame.rowCount === 0) {
        this.closeWithError('No rows to search'); // clearer than "Empty molecule…"
        return;
      }
      let progressBar: DG.TaskBarProgressIndicator | null = null;
      this.curIdx = this.dataFrame.currentRowIdx === -1 ? 0 : this.dataFrame.currentRowIdx;
      // _recheckRowSource is set ONLY by a filter-driven recheck. Force the recompute even when
      // followCurrentRow is off, else the set stays built against the old population (stale mask ANDed).
      const filterRecheck = this._recheckRowSource != null;
      if (computeData && (filterRecheck || (!this.gridSelect && this.followCurrentRow) || this.isEditedFromSketcher)) {
        progressBar = DG.TaskBarProgressIndicator.create(`Similarity search running...`);
        this.isComputing = true;
        this.error = '';
        this.root.classList.remove(`chem-malformed-molecule-error`);
        // Move the reference to the current row only when tracking it; a filter-driven recheck must keep
        // the pinned reference, not jump to the current row when followCurrentRow is off.
        if (!filterRecheck || this.followCurrentRow)
          this.targetMoleculeIdx = this.curIdx;
        if (DG.chem.Sketcher.isEmptyMolfile(this.targetMolecule)) {
          this.closeWithError(`Empty molecule cannot be used for similarity search`, progressBar);
          return;
        }
        const searchGen = ++this._searchGen;
        // Capture cutoff AND metric once: search and mask must use the SAME values even if a slider drag
        // or metric switch changes them during the await (else set, mask, and label could disagree).
        const cutoff = this.cutoff;
        const distanceMetric = this.distanceMetric;
        try {
          const rowSourceIdxs = this.getRowSourceIndexes();
          // Snapshot the row source this search runs against, so a later filter change is detected even
          // when net df.filter is unchanged — only while the toggle is on (the recheck dedupe is its sole
          // reader; cloning otherwise is waste).
          this._builtAgainstFilter = this.filterActive ? rowSourceIdxs.clone() : null;
          rowSourceIdxs.set(this.targetMoleculeIdx, true);
          const result = await chemSimilaritySearchEx(this.dataFrame!, this.moleculeColumn!,
            this.targetMolecule, distanceMetric as BitArrayMetrics, this.limit, cutoff,
            this.fingerprint as Fingerprint, rowSourceIdxs);
          if (searchGen !== this._searchGen)
            return; // a newer search started while this one ran — let it own the result/render
          if (!result) {
            // Clear the stale set + label so table/panel don't contradict the error state (no requestFilter
            // here — it could re-trigger the failing search; a later filter cycle un-narrows).
            this._clearSimilarSet();
            this.closeWithError(`Malformed molecule cannot be used for similarity search`);
            return;
          }
          this.molCol = result.df.getCol('smiles');
          this.idxs = result.df.getCol('indexes');
          this.scores = result.df.getCol('score');
          this.similarSetBitset = this._buildSimilarSetBitset(result.allDistances, rowSourceIdxs, cutoff);
          // Record the captured cutoff/metric this mask was built with, so the label stays in sync with
          // the applied mask rather than a mid-drag live value.
          this._bitsetCutoff = cutoff;
          this._bitsetMetric = distanceMetric;
          if (this.filterActive)
            this._reapplyFilter();
        } catch (e: unknown) {
          this._clearSimilarSet(); // clear stale set + label so neither contradicts the error state
          this.clearResults(); // don't strand the previous render's cards under the error toast
          grok.shell.error(e instanceof Error ? e.message : String(e)); // non-Error throw → readable text
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
      // Displayed set unchanged (non-recompute render → SAME molCol): refresh highlights in place and
      // re-arm the observer instead of rebuilding (which flashes every card). A recompute yields a NEW
      // molCol → full rebuild below.
      if (this.molCol === this._renderedMolCol && this.size === this._renderedSize && this._renderedCards.length) {
        for (const {grid, idx} of this._renderedCards)
          this._applyCardStyle(grid, idx);
        if (this._renderedGridDiv)
          this._observeLazyHosts(this._renderedGridDiv, this._renderedLazyHosts);
        progressBar?.close();
        return;
      }
      this.clearResults();
      // Grey out the action icons when there's nothing to act on. The Filter toggle stays clickable while
      // ON even on an empty set (the disabled CSS sets pointer-events:none, which would trap it on).
      const noSimilar = this.similarCount === 0;
      const setDisabled = (b: HTMLElement, disabled: boolean) => {
        b.classList.toggle('chem-similarity-action-disabled', disabled);
        b.setAttribute('aria-disabled', String(disabled));
        b.setAttribute('tabindex', disabled ? '-1' : '0');
      };
      setDisabled(this.selectBtn, noSimilar);
      setDisabled(this.filterBtn, noSimilar && !this.filterActive);
      // Set empties while toggle on → handler stops ANDing; drop the now-meaningless label.
      if (noSimilar && this.filterActive)
        this._removePushedSummary();
      // Surface the full similar-set size: Select/Filter act on ALL of it, which can far exceed the shown
      // cards (the "see more than 50" ask).
      const selectMsg = noSimilar ? 'Select all similar' : `Select all ${this.similarCount} similar`;
      ui.tooltip.bind(this.selectBtn, selectMsg);
      this.selectBtn.setAttribute('aria-label', selectMsg);
      const {width, height} = this.sizesMap[this.size];
      const panel = [];
      const grids = [];
      const lazyHosts: {host: HTMLElement; molStr: string}[] = [];
      const renderedCards: {grid: HTMLElement; idx: number}[] = []; // for later in-place restyle
      let cnt = 0; let cnt2 = 0;
      panel[cnt++] = this.metricsDiv;
      if (this.molCol && this.idxs && this.scores) {
        if (this.isEditedFromSketcher) {
          const label = this.sketchButton;
          const grid = ui.div([
            renderMolecule(this.targetMolecule, {width, height}),
            label], {style: {position: 'relative'}});
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
          molHost.style.width = `${width}px`;
          molHost.style.height = `${height}px`;
          lazyHosts.push({host: molHost, molStr: this.molCol!.get(i)});
          const grid = ui.div([
            molHost,
            label,
            molProps], { style: { position: 'relative' } });
          grid.classList.add('d4-flex-col');
          this._applyCardStyle(grid, idx); // current-row/selected highlight (also reused by in-place path)
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
          renderedCards.push({grid, idx});
        }
        // The full set can exceed the shown cards (capped at `limit`/DISPLAY_CAP). Count displayed cards
        // actually in the set (the reference card occupies a slot but isn't in similarCount) so the "more"
        // delta is exact either way.
        let shownSimilar = 0;
        for (let i = 0; i < displayCount; i++) {
          if (this.similarSetBitset?.get(this.idxs.get(i)))
            shownSimilar++;
        }
        const more = this.similarCount - shownSimilar;
        if (more > 0) {
          grids[cnt2++] = ui.divText(
            `+ ${more} more — Select / Filter act on all ${this.similarCount}`,
            {style: {alignSelf: 'center', padding: '0 8px', color: 'var(--grey-5)', fontSize: '11px'}});
        }
      }
      const gridDiv = ui.divH(grids, 'chem-viewer-grid');
      panel[cnt++] = gridDiv;
      this.root.appendChild(ui.panel([ui.divV(panel)]));
      this._observeLazyHosts(gridDiv, lazyHosts);
      // Snapshot this render so a later non-recompute render can refresh highlights in place.
      this._renderedCards = renderedCards;
      this._renderedLazyHosts = lazyHosts;
      this._renderedGridDiv = gridDiv;
      this._renderedMolCol = this.molCol;
      this._renderedSize = this.size;
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
  // An unknown distanceMetric (e.g. set via setOptions) → clear error rather than a cryptic TypeError later.
  if (typeof fpSim !== 'function')
    throw new Error(`Unknown similarity metric: ${metricName}`);
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

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

/** Similarity can display up to this many results; the shared base default (and the diversity
 * viewer, whose selection is O(n^2)) keep a much smaller cap. */
const MAX_LIMIT_SIMILARITY = 10000;
/** Hard cap on how many cards are built into the DOM at once, independent of `limit`. Select &
 * filter act on the full above-cutoff set, so capping the display loses nothing functionally. */
const DISPLAY_CAP = 200;
/** Default number of similar molecules shown (the base viewer's default of 12 is too small here —
 * the similarity viewer is cutoff-driven and shows up to the display cap). */
const DEFAULT_LIMIT_SIMILARITY = 200;
/** Prefix of the filter-panel label this viewer pushes; used to find/remove its own label. */
const FILTER_LABEL_PREFIX = 'Most similar structures:';
/** Suppress-window duration: our own requestFilter cycles shouldn't trigger a re-render. Must
 * comfortably exceed the base's 50 ms onFilterChanged debounce. */
const SUPPRESS_WINDOW_MS = 300;
/** Debounce before re-running the similarity search once filter activity settles — coalesces a burst
 * of filter cycles into a single search. */
const FILTER_RECHECK_DEBOUNCE_MS = 300;
/** Watchdog: if a substructure search is terminated before emitting 100%, un-block the recheck after
 * this quiet period so it isn't suppressed forever. */
const EXTERNAL_SEARCH_WATCHDOG_MS = 2000;
/** Dataframes a similarity viewer is currently probing (the probe's requestFilter fires every viewer's
 * handler). Shared so a SECOND similarity viewer on the same dataframe doesn't mistake a sibling's probe
 * for an external filter change and trigger a mutual recheck cascade.
 * KNOWN LIMITATION (two active similarity filters on one dataframe): this guards only the PROBE; a
 * viewer's _reapplyFilter() still fires the sibling's onRowsFiltering, so each re-search triggers the
 * other's. That is correct, not a bug — when viewer A's contribution changes the filtered population,
 * viewer B (searching "similar within the filtered set") genuinely must re-search. The _sameMask dedupe
 * makes it converge in ~1–2 debounced rounds to the mutually-consistent fixed point; the only visible
 * cost is a brief card reflash. Suppressing the cross-viewer recheck would be WRONG (B would show stale
 * results), so this rare config is left to converge rather than guarded. */
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
  /** The cutoff/metric the current similarSetBitset was built with — the filter-panel label reads these
   * (not the live this.cutoff/this.distanceMetric) so a fast slider drag can't show a threshold that
   * disagrees with the mask actually being applied until the next search settles. */
  private _bitsetCutoff: number = 0;
  private _bitsetMetric: string = '';
  filterActive: boolean = false;
  selectBtn!: HTMLElement; // assigned in the constructor
  filterBtn!: HTMLElement;
  private _filterSub: Subscription | null = null;
  private _molObserver: IntersectionObserver | null = null;
  // Snapshot of the last full card render, so a non-recompute render (selection / current-row / metadata
  // change) can refresh the per-card highlights IN PLACE instead of clearing + rebuilding the whole panel
  // (which flashes every card and re-pops its lazily-rendered molecule — the "flicker" on Select/Filter).
  private _renderedCards: {grid: HTMLElement; idx: number}[] = [];
  private _renderedLazyHosts: {host: HTMLElement; molStr: string}[] = [];
  private _renderedGridDiv: HTMLElement | null = null;
  private _renderedMolCol: DG.Column | null = null;
  private _renderedSize: string = '';
  private _pushedSummary: string | null = null;
  /** Subscription to the chem substructure-search progress event for our molecule column. */
  private _searchProgressSub: Subscription | null = null;
  /** Live (non-null) while a substructure search is streaming results on our molecule column. The
   * substructure filter calls requestFilter() per progress batch, so without this we'd launch one
   * similarity search per batch; instead we defer and run a single search when it reaches 100%. Doubles
   * as a watchdog so a search terminated before 100% can't block the recheck forever. */
  private _externalSearchTimer: ReturnType<typeof setTimeout> | null = null;
  private _skipContribution: boolean = false;
  /** True only during the viewer's OWN synchronous requestFilter calls (probe/reapply). Lets the
   * onRowsFiltering handler tell its own cycles apart from an external filter change precisely, so a
   * user clearing a filter mid-processing is never mistaken for our own cycle and dropped. */
  private _ownFilterCycle: boolean = false;
  /** Open while our OWN requestFilter cycles run, so onExternalFilterChanged skips them (avoids a feedback
   * loop when we apply our own mask). Time-windowed via _suppressTimer; see _holdSuppress. */
  private _suppressFilterRender: boolean = false;
  private _suppressTimer: ReturnType<typeof setTimeout> | null = null;
  private _searchGen: number = 0;
  /** The row-source (other-filters result) the current similarSetBitset was last searched against.
   * Lets a filter change be detected even when the net df.filter is unchanged (e.g. clearing a filter). */
  private _builtAgainstFilter: DG.BitSet | null = null;
  private _filterRecheckTimer: ReturnType<typeof setTimeout> | null = null;
  /** Hand-off of the row source the recheck just probed, so the render it triggers reuses it instead
   * of probing again. Valid only for the immediately following render (see _scheduleFilterRecheck). */
  private _recheckRowSource: DG.BitSet | null = null;

  get targetMolecule(): string {
    return this.isEditedFromSketcher ?
      this.sketchedMolecule :
      this.moleculeColumn?.get(this.targetMoleculeIdx);
  }

  constructor() {
    super(SIMILARITY);
    // Per-viewer property overrides (kept in the subclass, not the base): raise the limit cap, and default
    // the row source to 'Filtered' so cards respect the active filter. Setting them here (not in
    // onTableAttached) means a restored value is never clamped against a stale max, and 'Reset to default'
    // honours the similarity-specific defaults.
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

    // Select all shown (the similar set, i.e. rows with score >= cutoff). Deselect = Escape (platform built-in).
    // Reuse the platform's grid "select all" icon (svg sprite) so it matches the table's own select-all control.
    this.selectBtn = ui.iconSvg('select-all', () => {
      // Guard against an empty set and a length mismatch (BitSet.and requires equal length).
      if (!this.dataFrame || this.similarCount === 0 || this.similarSetBitset!.length !== this.dataFrame.filter.length)
        return;
      // Select only similar rows that pass the active filter (similar set ∩ table filter).
      const sel = this.similarSetBitset!.clone().and(this.dataFrame.filter, false);
      if (sel.trueCount === 0)
        return; // every similar row is hidden by another filter — don't wipe the existing selection
      const changed = !this._sameMask(sel, this.dataFrame.selection); // skip the toast on a repeat click
      this.dataFrame.selection.copyFrom(sel);
      // Confirm the off-panel action — the selection happens in the grid, away from where the user clicked.
      if (changed)
        grok.shell.info(`${sel.trueCount} similar ${sel.trueCount === 1 ? 'row' : 'rows'} selected`);
    }, 'Select all similar');
    this.selectBtn.classList.add('chem-similarity-action-btn');

    // Collaborative filter toggle (off by default).
    this.filterBtn = ui.icons.filter(() => this._toggleFilter(),
      'Filter table to similar set (toggle)');
    this.filterBtn.classList.add('chem-similarity-action-btn');

    // Accessibility: give each icon an explicit aria-label up front (renderInternal later refines the
    // select label with a count, _toggleFilter the filter label with on/off state) — don't rely on the
    // ui.icons.* helpers to set it, or a screen reader sees an unlabelled button before the first
    // render/toggle. Then make the <i> elements behave as buttons (focusable + keyboard-activatable).
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
    // setOptions bypasses the slider's [0, 1] bounds; clamp so a stray cutoff can't mark every molecule
    // similar (cutoff 0) or nothing (cutoff 1) from out-of-range input. (limit is clamped in the base.)
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
    // Read the cutoff/metric the mask was BUILT with (not the live values) so a mid-drag label can't show
    // a threshold the applied mask doesn't match. Round so a slider value like 0.3 doesn't render as
    // 0.30000000000000004.
    const cutoff = Number(this._bitsetCutoff.toFixed(3));
    return `${FILTER_LABEL_PREFIX} ${this._bitsetMetric} ≥ ${cutoff}`;
  }

  /** Builds the "similar set": rows within the current row source whose similarity
   * score meets the cutoff. Malformed rows carry a 100.0 sentinel and are excluded. */
  private _buildSimilarSetBitset(allDistances: number[], rowSourceIdxs: DG.BitSet, cutoff: number): DG.BitSet {
    // A sketched query is not a table row, so it excludes nothing; otherwise exclude the reference row.
    const refIdx = this.isEditedFromSketcher ? -1 : this.targetMoleculeIdx;
    // Read the row source as a packed buffer once and test bits in JS — avoids one get(i) interop per row.
    const srcBuf = rowSourceIdxs.getBuffer();
    // Bound by the number of scored rows: if rows were added during a fast (cached) search, rowCount can
    // exceed allDistances.length; those new rows have no score yet and must be excluded (don't read OOB).
    const scored = allDistances.length;
    return DG.BitSet.create(this.dataFrame!.rowCount, (i) =>
      i < scored && i !== refIdx && ((srcBuf[i >>> 5] & (1 << (i & 31))) !== 0) &&
      allDistances[i] >= cutoff && allDistances[i] <= 1);
  }

  /** Open/extend a single suppress window so our own requestFilter cycles don't trigger a re-render.
   * One shared, cancellable timer means the probe and re-apply never race to clear each other's window.
   * 300 ms comfortably covers the 50 ms onFilterChanged debounce. */
  private _holdSuppress(): void {
    this._suppressFilterRender = true;
    if (this._suppressTimer)
      clearTimeout(this._suppressTimer);
    this._suppressTimer = setTimeout(() => {
      this._suppressFilterRender = false;
      this._suppressTimer = null;
    }, SUPPRESS_WINDOW_MS);
  }

  /** Re-runs the filter cycle while suppressing the viewer's own onFilterChanged re-render,
   * so applying our own mask does not trigger a recompute loop. */
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

  /** Drop the similar set and its filter-panel label together — a failed/empty search must leave neither a
   * stale mask (which onRowsFiltering would keep ANDing) nor a label that contradicts the error state. */
  private _clearSimilarSet(): void {
    this.similarSetBitset = null;
    this._removePushedSummary();
  }

  private _toggleFilter(): void {
    // Block activation on an empty set, but always allow deactivation (e.g. Reset all filters, or
    // turning it off after the cutoff emptied the set) so the toggle can never get stuck on.
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
        // An external filter cycle (NOT one of our own probe/reapply cycles) may have changed the
        // filtered population — schedule a single debounced recheck rather than searching inline.
        // Essential for CLEARING a filter (our mask is a subset of the other filters, so removing one
        // leaves the net df.filter unchanged and onFilterChanged never fires); the debounce + dedupe also
        // collapses a live substructure search into ONE similarity search after it settles. Gate on the
        // precise _ownFilterCycle flag (not the time-windowed suppress) so a user clearing a filter while
        // we're still processing a prior change is honoured instead of dropped. Also skip while ANY
        // similarity viewer is probing this dataframe, so two viewers don't trigger a mutual recheck cascade.
        if (!this._ownFilterCycle && !_probingDataframes.has(df))
          this._scheduleFilterRecheck();
        // Only AND a non-empty mask, so an empty similar set never hides every row.
        if (this.similarCount > 0 && this.similarSetBitset!.length === df.filter.length) {
          df.filter.and(this.similarSetBitset!, false);
          // Surface the filter in the filter panel (replace a stale summary when the cutoff changes).
          const s = this._filterSummary();
          if (this._pushedSummary !== s)
            this._removePushedSummary();
          if (!df.rows.filters.includes(s))
            df.rows.filters.push(s);
          this._pushedSummary = s;
        } else {
          // Not contributing this cycle (empty set, or a length mismatch right after a row add/remove
          // before the mask is rebuilt) — drop the label so the panel never shows a filter that isn't applied.
          this._removePushedSummary();
        }
      });
      this._reapplyFilter();
    } else {
      this._filterSub?.unsubscribe();
      this._filterSub = null;
      // Drop any pending filter-driven recheck (and the substructure-search watchdog): they're self-guarding
      // on `filterActive` when they fire, but cancelling here avoids a spurious search if the user rapidly
      // toggles back on while a substructure search is still streaming.
      if (this._filterRecheckTimer) {
        clearTimeout(this._filterRecheckTimer);
        this._filterRecheckTimer = null;
      }
      if (this._externalSearchTimer) {
        clearTimeout(this._externalSearchTimer);
        this._externalSearchTimer = null;
      }
      this._removePushedSummary();
      // Un-filter the TABLE, but suppress the viewer's re-render: the similar set is built against the OTHER
      // filters (the probe excludes our own contribution), so removing our toggle doesn't change which cards
      // are shown — recomputing + rebuilding them would just flicker the panel for no visible change.
      this._holdSuppress();
      df.rows.requestFilter();
    }
  }

  /** Result of all OTHER filters, excluding this viewer's own contribution. requestFilter runs the
   * onRowsFiltering handlers synchronously, so we probe with our contribution skipped, snapshot
   * df.filter, then restore — all in one tick, so the grid never repaints the intermediate state.
   * Order-independent, unlike snapshotting before our own AND inside the handler. */
  private _computeOtherFiltersMask(): DG.BitSet {
    const df = this.dataFrame;
    // ASSUMES requestFilter() fires ALL onRowsFiltering handlers synchronously (true for every
    // DG.Filter today). If a future filter does async work in its handler, the probe would snapshot the
    // all-true reset state and we'd search the full table while believing we respect the other filters.
    this._holdSuppress();
    this._ownFilterCycle = true; // both requestFilter cycles below are ours, not external
    _probingDataframes.add(df!); // tell any sibling similarity viewer these cycles are a probe, not a change
    let others: DG.BitSet;
    this._skipContribution = true;
    try {
      df.rows.requestFilter(); // df.filter = result of all other filters (our AND skipped)
      others = df.filter.clone();
    } finally {
      this._skipContribution = false;
      try {
        df.rows.requestFilter(); // restore df.filter = others AND ours
        this._holdSuppress(); // reschedule the single window to also cover the restore cycle
      } finally {
        // Reset the cycle flags FIRST so that even if the defensive re-AND below throws, the viewer never
        // stays deaf to external changes (the AND uses notify=false, so nothing observes the flags between
        // here and it — a zero-risk reorder that lets the unreachable "AND throws" case self-recover).
        this._ownFilterCycle = false;
        _probingDataframes.delete(df!);
        // Defensively re-apply our mask. On the happy path the restore cycle already ANDed it (this is a
        // harmless idempotent re-AND); but if a third-party onRowsFiltering handler THREW during the
        // restore, df.filter would otherwise be left without our contribution — this recovers it instead
        // of stranding the table on the others-only mask until the next filter cycle. (If the restore
        // throws after the probe also threw, the restore's error supersedes the probe's — a finally-throw
        // replaces the original; acceptable for that narrow double-throw case.)
        if (this.filterActive && this.similarSetBitset && this.similarSetBitset.length === df!.filter.length)
          df!.filter.and(this.similarSetBitset, false);
      }
    }
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

  /** Set a result card's current-row / selected highlight (class + background) from the live selection and
   * reference row. Used both when building cards and when refreshing them in place (see renderInternal). */
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

  /** (Re)arm the lazy-render observer over `lazyHosts`, rooted at the viewer's scroll container with a 300px
   * prefetch margin (100px against the viewport when there's no scrollable ancestor — enough to prefetch
   * just-off-screen cards without marking every card as intersecting). Already-rendered hosts are skipped,
   * so this is safe to call again on an in-place refresh to cover cards not yet scrolled into view. */
  private _observeLazyHosts(gridDiv: HTMLElement, lazyHosts: {host: HTMLElement; molStr: string}[]): void {
    const {width, height} = this.sizesMap[this.size];
    const scrollRoot = this._findScrollParent(gridDiv);
    this._molObserver?.disconnect();
    // Close over the local `observer`, not this._molObserver (which the next render reassigns) — so a late
    // callback from this render unobserves on its OWN observer and never touches a newer instance.
    const observer = new IntersectionObserver((entries) => {
      for (const e of entries) {
        if (!e.isIntersecting)
          continue;
        const host = e.target as HTMLElement & {_molStr?: string; _rendered?: boolean};
        observer.unobserve(host);
        if (host._rendered || !host.isConnected)
          continue; // skip detached/already-drawn hosts — don't waste RDKit work after a fast scroll-then-close
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

  /** True when two row-source masks are equal — compares the packed word buffers directly (no clone),
   * bailing on the first differing word. */
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

  /** Coalesce filter-driven re-searches: after filter activity settles, re-run the similarity search at
   * most once, and only if the other-filters result actually changed since the current set was built.
   * Keeps a live substructure (or any) filter from launching a similarity search on every cycle. */
  private _scheduleFilterRecheck(): void {
    // A substructure search is streaming results (calling requestFilter per batch) — don't re-search on
    // every batch; the progress handler will trigger a single recheck once the search reaches 100%. The
    // watchdog timer is live exactly while such a search is mid-stream (set on progress < 100, cleared on
    // completion/quiesce), so its presence is the "in progress" flag.
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
      // In filter-based modes the similar set depends on the other filters — skip the (expensive) search
      // when they are unchanged. In other modes a filter change can't alter the set, so just re-render.
      if (filterBased) {
        const current = this.getRowSourceIndexes(); // runs the probe ONCE
        if (this._sameMask(current, this._builtAgainstFilter))
          return;
        // Hand the just-probed mask to the render so it doesn't probe again (halves the filter cycles).
        // render()'s synchronous prefix consumes it before its first await; the post-clear below drops
        // it if the render took a no-recompute path, so it can never leak to an unrelated caller.
        this._recheckRowSource = current;
        this.render(true);
        this._recheckRowSource = null;
      } else {
        // All/Selected: the similar set is filter-independent, so re-render (refresh highlights) WITHOUT
        // an O(N) similarity search.
        this.render(false);
      }
    }, FILTER_RECHECK_DEBOUNCE_MS);
  }

  /** Subscribe to the chem substructure-search progress event for the current molecule column. A
   * substructure (or "is similar") filter on that column calls requestFilter() on every progress batch;
   * while that stream is running we hold off the recheck (see _scheduleFilterRecheck) and fire exactly
   * one similarity search when the search reaches 100% (or after a quiet period, if it was terminated
   * before completing). */
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
        // Watchdog: if the search is terminated before 100% (no further progress), recover after a quiet
        // period so the recheck isn't blocked forever. While it is set, _scheduleFilterRecheck stays gated.
        this._externalSearchTimer = setTimeout(() => {
          this._externalSearchTimer = null;
          if (this.filterActive)
            this._scheduleFilterRecheck();
        }, EXTERNAL_SEARCH_WATCHDOG_MS);
      } else {
        // Search complete (timer already cleared above): run one similarity search against the settled set.
        if (this.filterActive)
          this._scheduleFilterRecheck();
      }
    });
  }

  /** While our filter is active, route filter-change re-renders through the debounced, deduped recheck
   * so a burst of filter cycles produces at most one similarity search after the filter settles.
   * NOTE: this is one of TWO paths that schedule the recheck — the onRowsFiltering handler is the other
   * (and the only one that fires when net df.filter is unchanged, e.g. clearing a filter). The
   * redundancy is intentional and load-bearing: the 300ms _suppressFilterRender window can swallow this
   * onFilterChanged path for a genuine external change, and the onRowsFiltering path covers that case.
   * Don't remove either without re-checking the clear-and-re-expand behaviour. */
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
    // PLATFORM CONTRACT: the filter-active branch below probes via _computeOtherFiltersMask, which assumes
    // requestFilter() runs every onRowsFiltering subscriber SYNCHRONOUSLY (true for all current DG filters,
    // but undocumented). If a future filter ever contributes asynchronously, the probe would snapshot the
    // reset (all-true) state and this viewer would silently search the whole table believing it respects the
    // other filters — a wrong-result, not a crash. Revisit here if DG filters ever go async.
    // Reuse the mask the pending recheck just probed (see _scheduleFilterRecheck) instead of probing again.
    if (this._recheckRowSource) {
      const m = this._recheckRowSource;
      this._recheckRowSource = null;
      return m;
    }
    // While our own filter is active, search against the OTHER filters' result so changing the
    // reference row does not progressively shrink the set (order-independent probe, not a snapshot).
    // Covers both filter-based row sources; FilteredSelected also intersects the selection. In All/
    // Selected modes the similar set is filter-independent, so no probe — the handler still ANDs our
    // mask into df.filter, but an external filter change can't change which molecules ARE similar.
    if (this.filterActive &&
      (this.rowSource === RowSourceTypes.Filtered || this.rowSource === RowSourceTypes.FilteredSelected)) {
      const others = this._computeOtherFiltersMask();
      return this.rowSource === RowSourceTypes.FilteredSelected ?
        others.and(this.dataFrame!.selection, false) : others;
    }
    return super.getRowSourceIndexes();
  }

  /** Reset all per-table filter-runtime state, including the three timers that gate the recheck. Shared by
   * detach() and onTableAttached()'s re-attach guard so a field added to one path can't be silently
   * forgotten in the other. Callers handle their own extras (detach: progress/observer subs + requestFilter
   * restore; onTableAttached: limit clamp / label strip). */
  private _resetFilterState(): void {
    this._searchGen++; // abandon any in-flight search so it can't write state after the reset
    this._filterSub?.unsubscribe();
    this._filterSub = null;
    this.filterActive = false;
    // Drop the previous table's similar set: a stale non-zero trueCount would otherwise slip past the
    // toggle's activation guard during the async attach-render, briefly showing an "active" filter that
    // ANDs nothing (the length guard skips it) until the render rebuilds the mask.
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
    // The substructure-search watchdog also gates _scheduleFilterRecheck, so clear it here too (the progress
    // SUBSCRIPTION that drives it is torn down separately, by detach / _subscribeSearchProgress).
    if (this._externalSearchTimer) {
      clearTimeout(this._externalSearchTimer);
      this._externalSearchTimer = null;
    }
    this._removePushedSummary();
    this.filterBtn.classList.remove('active');
    this.filterBtn.setAttribute('aria-pressed', 'false');
    // Reset the label too, so a screen reader doesn't announce "Deactivate…" on a now-off button after a
    // detach-while-active.
    this.filterBtn.setAttribute('aria-label', 'Filter table to similar set');
  }

  override async onTableAttached(): Promise<void> {
    // Idempotent re-attach guard: if the platform reused this instance WITHOUT a detach, drop any stale
    // filter subscription/state bound to the previous dataframe so we neither double-subscribe nor carry
    // the toggle over. (In the normal detach-then-attach path this is all already reset — a no-op.)
    this._resetFilterState();
    // Keep the modest default limit (select & filter act on the full above-cutoff set regardless of
    // how many cards show); only clamp a too-large restored value down to the table size.
    if (this.dataFrame && this.limit > this.dataFrame.rowCount)
      this.limit = Math.max(1, this.dataFrame.rowCount);
    // The toggle is session-only (off on a fresh attach) — strip any of our labels left in the filter
    // panel by a restored project so it never lingers without a live filter behind it.
    if (this.dataFrame) {
      for (const s of [...this.dataFrame.rows.filters]) {
        if (s.startsWith(FILTER_LABEL_PREFIX))
          this.dataFrame.rows.filters.remove(s);
      }
    }
    // Match the other Chem filters: the filter panel's "Reset all filters" must clear ours too. Register
    // AFTER super (super's first act is `this.subs = []`, which would otherwise drop a sub pushed before
    // it), and in a `finally` so a throwing first render still can't skip it. The base's subs-reset already
    // prevents accumulation on re-attach, so this single registration is safe.
    try {
      await super.onTableAttached();
    } finally {
      this.subs.push(grok.events.onResetFilterRequest.subscribe(() => {
        // onResetFilterRequest is GLOBAL (no dataframe payload), so with several tables open it also fires
        // for other tables' "Reset all filters". Only honour it when our table is the one in view, so
        // resetting filters on table B doesn't silently switch our filter off on table A.
        if (this.filterActive && this.dataFrame === grok.shell.tv?.dataFrame)
          this._toggleFilter();
      }));
      // Re-render on row add (similarity only — kept out of the base so it doesn't force an O(n^2) diversity
      // recompute on every addition): a mask sized to the old row count would otherwise silently turn the
      // filter into a no-op while its label lingers in the filter panel. While the toggle is on, route
      // through the recheck (which rebuilds the mask to the new length even when followCurrentRow is off,
      // avoiding a ~300ms zombie-filter flicker); otherwise a plain re-render suffices.
      if (this.dataFrame) {
        this.subs.push(DG.debounce(this.dataFrame.onRowsAdded, 50).subscribe(async () => {
          if (this.filterActive)
            this._scheduleFilterRecheck();
          else
            await this.render();
        }));
      }
      // Watch substructure-search progress on our molecule column so a streaming search doesn't trigger
      // one similarity search per batch (super has detected this.moleculeColumn by now).
      this._subscribeSearchProgress();
    }
  }

  override detach(): void {
    const wasFiltering = this.filterActive; // capture before _resetFilterState clears it
    // Teardown beyond the shared filter state: the lazy-render observer and the progress subscription
    // (its watchdog timer is cleared by _resetFilterState).
    this._molObserver?.disconnect();
    this._molObserver = null;
    this._searchProgressSub?.unsubscribe();
    this._searchProgressSub = null;
    this._resetFilterState();
    // If we were filtering, re-run the filter cycle (our listener is gone now) so the table is restored
    // instead of being left filtered with no viewer or label to explain why.
    if (wasFiltering && this.dataFrame)
      this.dataFrame.rows.requestFilter();
    super.detach();
  }

  async renderInternal(computeData: boolean): Promise<void> {
    // Tear down the previous render's lazy-render observer up front: the early-exit paths below (no/invalid
    // column, empty table, malformed query) never reach the success-path disconnect, so without this an
    // observer from a prior populated render would keep its (now-detached) card hosts alive until the next
    // successful render. The success path creates a fresh observer.
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
        this.closeWithError('No rows to search'); // clearer than "Empty molecule cannot be used…"
        return;
      }
      let progressBar: DG.TaskBarProgressIndicator | null = null;
      this.curIdx = this.dataFrame.currentRowIdx === -1 ? 0 : this.dataFrame.currentRowIdx;
      // _recheckRowSource is set ONLY by a filter-driven recheck (see _scheduleFilterRecheck). Force the
      // recompute in that case even when followCurrentRow is off — otherwise the similar set stays built
      // against the old filtered population and the contribution handler keeps ANDing a stale mask.
      const filterRecheck = this._recheckRowSource != null;
      if (computeData && (filterRecheck || (!this.gridSelect && this.followCurrentRow) || this.isEditedFromSketcher)) {
        progressBar = DG.TaskBarProgressIndicator.create(`Similarity search running...`);
        this.isComputing = true;
        this.error = '';
        this.root.classList.remove(`chem-malformed-molecule-error`);
        // Move the reference to the current row only when actually tracking it; a filter-driven recheck must
        // keep the existing (pinned) reference, not jump it to the current row when followCurrentRow is off.
        if (!filterRecheck || this.followCurrentRow)
          this.targetMoleculeIdx = this.curIdx;
        if (DG.chem.Sketcher.isEmptyMolfile(this.targetMolecule)) {
          this.closeWithError(`Empty molecule cannot be used for similarity search`, progressBar);
          return;
        }
        const searchGen = ++this._searchGen;
        // Capture cutoff AND metric once: the search and the mask must use the SAME values even if a slider
        // drag or metric switch changes this.cutoff / this.distanceMetric during the await (otherwise the
        // displayed set, the mask, and the filter-panel label could disagree).
        const cutoff = this.cutoff;
        const distanceMetric = this.distanceMetric;
        try {
          const rowSourceIdxs = this.getRowSourceIndexes();
          // Snapshot the row source (other-filters result) this search runs against, so a later filter
          // change is detected even when net df.filter is unchanged — but only while the toggle is on, since
          // _builtAgainstFilter is read solely by the (filter-active) recheck dedupe. Cloning a full-width
          // bitset on every search when the toggle is off is pure waste.
          this._builtAgainstFilter = this.filterActive ? rowSourceIdxs.clone() : null;
          rowSourceIdxs.set(this.targetMoleculeIdx, true);
          const result = await chemSimilaritySearchEx(this.dataFrame!, this.moleculeColumn!,
            this.targetMolecule, distanceMetric as BitArrayMetrics, this.limit, cutoff,
            this.fingerprint as Fingerprint, rowSourceIdxs);
          if (searchGen !== this._searchGen)
            return; // a newer search started while this one ran — let it own the result/render
          if (!result) {
            // Clear the stale set + label so the table/panel don't contradict the error state (a later
            // filter cycle then un-narrows; no requestFilter here, which could re-trigger the failing search).
            this._clearSimilarSet();
            this.closeWithError(`Malformed molecule cannot be used for similarity search`);
            return;
          }
          this.molCol = result.df.getCol('smiles');
          this.idxs = result.df.getCol('indexes');
          this.scores = result.df.getCol('score');
          this.similarSetBitset = this._buildSimilarSetBitset(result.allDistances, rowSourceIdxs, cutoff);
          // Record the captured cutoff/metric this mask was built with, so the filter-panel label stays in
          // sync with the applied mask rather than a mid-drag this.cutoff / this.distanceMetric.
          this._bitsetCutoff = cutoff;
          this._bitsetMetric = distanceMetric;
          if (this.filterActive)
            this._reapplyFilter();
        } catch (e: unknown) {
          this._clearSimilarSet(); // clear the stale set + label so neither contradicts the error state
          this.clearResults(); // don't leave the previous render's cards stranded under the error toast
          grok.shell.error(e instanceof Error ? e.message : String(e)); // a non-Error throw → readable text
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
      // Displayed set unchanged (a non-recompute render: selection / current-row / metadata change produces
      // the SAME molCol). Refresh the per-card highlights in place and re-arm the lazy observer instead of
      // clearing + rebuilding the whole panel — which would flash every card and re-pop its molecule. A
      // recompute yields a NEW molCol, so it falls through to the full rebuild below.
      if (this.molCol === this._renderedMolCol && this.size === this._renderedSize && this._renderedCards.length) {
        for (const {grid, idx} of this._renderedCards)
          this._applyCardStyle(grid, idx);
        if (this._renderedGridDiv)
          this._observeLazyHosts(this._renderedGridDiv, this._renderedLazyHosts);
        progressBar?.close();
        return;
      }
      this.clearResults();
      // Grey out the action icons when there is nothing to act on. Select is disabled on an empty set.
      // The Filter toggle is disabled only when it's OFF — if it's already ON and the set just emptied,
      // it MUST stay clickable so the user can turn it back off (the CSS sets pointer-events:none on the
      // disabled class, which would otherwise trap a stuck-on toggle that _toggleFilter is happy to clear).
      const noSimilar = this.similarCount === 0;
      const setDisabled = (b: HTMLElement, disabled: boolean) => {
        b.classList.toggle('chem-similarity-action-disabled', disabled);
        b.setAttribute('aria-disabled', String(disabled));
        b.setAttribute('tabindex', disabled ? '-1' : '0');
      };
      setDisabled(this.selectBtn, noSimilar);
      setDisabled(this.filterBtn, noSimilar && !this.filterActive);
      // If the set empties while the toggle is on, the handler stops ANDing (table un-filters) — drop
      // the now-meaningless label too, so it doesn't linger until the user toggles off.
      if (noSimilar && this.filterActive)
        this._removePushedSummary();
      // Surface the full similar-set size: Select/Filter act on ALL of it, which can far exceed the
      // shown cards (this is the "see more than 50" ask — make the real count visible).
      const selectMsg = noSimilar ? 'Select all similar' : `Select all ${this.similarCount} similar`;
      ui.tooltip.bind(this.selectBtn, selectMsg);
      this.selectBtn.setAttribute('aria-label', selectMsg);
      const {width, height} = this.sizesMap[this.size];
      const panel = [];
      const grids = [];
      const lazyHosts: {host: HTMLElement; molStr: string}[] = [];
      const renderedCards: {grid: HTMLElement; idx: number}[] = []; // idx cards, for later in-place restyle
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
          this._applyCardStyle(grid, idx); // current-row / selected highlight (also reused by the in-place path)
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
        // The full similar set can exceed the shown cards (which are capped at `limit` and DISPLAY_CAP) —
        // tell the user, and that Select / Filter act on the whole set, not just what is rendered. Count
        // displayed cards that are actually in the similar set (the reference card occupies a slot but is
        // excluded from similarCount) so the "more" delta is exact whether or not the reference is shown.
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
      // Snapshot this render so a later non-recompute render can refresh highlights in place (see above).
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

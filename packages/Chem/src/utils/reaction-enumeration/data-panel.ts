import {Subscription} from 'rxjs';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EnumeratorConfig} from './config';
import {PerRoundOverride} from './enumerate';
import {MountedViewerRegistry} from './viewer-mount';
import {CHANGED_DOT_STYLE, clampRounds, DataKey, detectChemSemTypes, Mode} from './shared';

export type TabBadge = {el: HTMLSpanElement; refresh: (n: number | null) => void};
export const makeTabBadge = (): TabBadge => {
  const el = document.createElement('span');
  el.className = 'chem-enum-tab-badge';
  return {el, refresh: (n: number | null) => {
    el.textContent = n != null ? String(n) : '';
    el.style.display = n != null ? '' : 'none';
  }};
};

// Table names aren't unique in a workspace, and a ChoiceInput's `.value` getter re-wraps the table
// on every read — `.id` (assigned on registration) is the only reliable identity.
function isSameTrackedTable(a: DG.DataFrame | null | undefined, b: DG.DataFrame | null | undefined): boolean {
  return a?.id != null && a.id === b?.id;
}

// ui.input.table is a ChoiceInput over the workspace tables list, which silently rejects a value
// that isn't in it. Register first, then dip through null to drop cached pointer-equality.
function assignTableInput(input: DG.InputBase<DG.DataFrame | null>, df: DG.DataFrame): void {
  grok.shell.addTable(df);
  try {input.value = null;} catch {/* nullable: false rejects */}
  input.value = df;
}

const SELECT_ROWS_OR_FILTER = 'Select rows (Ctrl/Shift+click) or apply a filter';

function closeTableSafe(t: DG.DataFrame | null | undefined, label: string): void {
  if (!t) return;
  try {grok.shell.closeTable(t);} catch (e) {console.warn(`Could not close ${label}: ${e}`);}
}

/** Clones `df` down to its selection, falling back to the active filter when nothing is selected
 * (so filtering alone is enough to subset). Returns null, after an info toast, when there's
 * nothing to narrow. */
function cloneSubsetByRows(df: DG.DataFrame, emptyMsg: string): DG.DataFrame | null {
  const sel = df.selection;
  const mask = sel.anyTrue ? sel : df.filter;
  if (!mask.anyTrue) {
    grok.shell.info(emptyMsg);
    return null;
  }
  if (!mask.anyFalse) {
    grok.shell.info('All rows are selected — nothing to subset.');
    return null;
  }
  const subset = df.clone(mask);
  subset.name = `${df.name} (subset, ${subset.rowCount}/${df.rowCount} rows)`;
  detectChemSemTypes(subset);
  // clone(mask) leaves the source's own .selection set to `mask`, and the clone's .filter all-false.
  df.selection.setAll(false, false);
  subset.filter.setAll(true, false);
  return subset;
}

// Breadth-first draws from every earlier product, so it never reads the per-round BB pool.
function bbOverrideSuppressedInBreadth(isBuildingBlocks: boolean, mode: Mode): boolean {
  return isBuildingBlocks && mode === 'breadth';
}

interface SubsetState {
  prev: DG.DataFrame | null; original: DG.DataFrame | null; freshClone: DG.DataFrame | null;
}
interface StepState { df: DG.DataFrame | null; sub: Subscription | null; committed: boolean; }

export interface DataPanelOpts {
  idx: number; noun: string;
  input: DG.InputBase<DG.DataFrame | null>;
  badge?: TabBadge;
  noTableMsg?: string;
  emptyMsg?: string;
  /** How to fold this component's per-round work df into a PerRoundOverride fragment. */
  apply: (o: PerRoundOverride, work: DG.DataFrame, cfg: EnumeratorConfig) => void;
}

export interface DataPanelDeps {
  view: DG.View;
  viewerHost: MountedViewerRegistry;
  getConfig: () => EnumeratorConfig;
  currentMode: () => Mode;
  currentRounds: () => number;
  refreshValidation: () => void;
  refreshCfgRibbon: () => void;
}

/** One grid plus a round-tab strip: "All rounds" shows the full library (with the global "Subset by
 * selection"); "Round k" shows a display-only clone whose rows are that round's subset. Switching
 * tabs swaps what the single grid displays — there is never a second grid. */
export class DataPanel {
  readonly panel: HTMLElement;

  private selStep = 0; // 0 = All rounds; 1..rounds = that round's subset
  private filtersOn = false;
  // host -> what was last mounted there, so the second of the two renders a table swap triggers
  // (input.onChanged plus the caller's own explicit render) is a no-op.
  private readonly lastMounted = new Map<HTMLElement, {identity: string | DG.DataFrame; filtersOn: boolean}>();
  private readonly state: SubsetState = {prev: null, original: null, freshClone: null};
  // Index k-1 = round k. `committed` is explicit, never inferred from row count: a clone can drift
  // from the global table without any onChanged firing, and inferring from that drift let deleted
  // rows resurrect in a round's reactant pool.
  private readonly stepState: StepState[] = [];
  private readonly stepTabsHost: HTMLElement;
  private stepTabsSub: Subscription | null = null;
  private readonly stepDots: (HTMLElement | null)[] = []; // index 0 ("All rounds") unused
  // Each pane owns its own barHost/gridHost — relocating a live grid between panes corrupts it.
  private paneHosts: {barHost: HTMLElement; gridHost: HTMLElement}[] = [];

  constructor(private readonly opts: DataPanelOpts, private readonly deps: DataPanelDeps) {
    this.stepTabsHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', display: 'flex',
      flexDirection: 'column', overflow: 'hidden'}});
    this.panel = ui.div([this.stepTabsHost], {style: {
      height: '100%', display: 'flex', flexDirection: 'column', minHeight: '0',
      background: 'var(--white)', overflow: 'hidden'}});
    this.buildStepTabs(0);
  }

  render(): void {this.renderGrid(); this.renderBar(); this.updateDots();}
  refreshDisplay(): void {this.renderBar(); this.updateDots();}

  onTableChanged(): void {
    const df = this.opts.input.value;
    if (df) detectChemSemTypes(df);
    // Matching neither prev nor freshClone means the input changed from outside both — e.g. a
    // different file — so the old subset bookkeeping would restore the wrong table's `original`.
    if (df && !isSameTrackedTable(df, this.state.prev) && !isSameTrackedTable(df, this.state.freshClone)) {
      this.state.original = null;
      this.state.prev = null;
      this.state.freshClone = null;
    }
    for (const s of this.stepState) s?.sub?.unsubscribe();
    this.stepState.length = 0;
    this.buildStepTabs(0);
    this.deps.refreshValidation();
  }

  // Deliberately keeps stepState entries beyond the new round count: typing "10" over "5"
  // transiently passes through "1", which would otherwise destroy overrides on rounds 2-5.
  onRoundsChanged(): void {this.buildStepTabs(this.selStep);}

  /** Applies this component's round-r override (if any) onto `out`; returns whether it did. */
  applyOverrideForRound(r: number, out: PerRoundOverride, cfg: EnumeratorConfig): boolean {
    const entry = this.stepState[r - 1];
    if (!entry?.committed || !entry.df) return false;
    this.opts.apply(out, entry.df, cfg);
    return true;
  }

  hasAnyOverride(): boolean {
    for (let k = 1; k <= this.roundCount(); k++) if (this.hasOverride(k)) return true;
    return false;
  }

  private hasOverride(k: number): boolean {
    if (!this.stepState[k - 1]?.committed) return false;
    if (bbOverrideSuppressedInBreadth(this.opts.idx === 1, this.deps.currentMode())) return false;
    return true;
  }

  private roundCount(): number {
    return clampRounds(this.deps.currentRounds());
  }

  private updateDots(): void {
    for (let k = 1; k <= this.stepDots.length - 1; k++) {
      const dot = this.stepDots[k];
      if (dot) dot.style.display = this.hasOverride(k) ? '' : 'none';
    }
  }

  // Preserves the existing `committed` value; callers that need to change it do so right after.
  private setStepWork(k: number, work: DG.DataFrame | null): void {
    this.stepState[k - 1]?.sub?.unsubscribe();
    const sub = work ? work.onSelectionChanged.subscribe(() => {this.updateDots(); this.renderBar();}) : null;
    if (sub) this.deps.view.subs.push(sub);
    const entry = this.stepState[k - 1];
    if (entry) {entry.df = work; entry.sub = sub;} else this.stepState[k - 1] = {df: work, sub, committed: false};
  }

  private stepClone(k: number): DG.DataFrame | null {
    const existing = this.stepState[k - 1]?.df;
    if (existing) return existing;
    const global = this.opts.input.value;
    if (!global) return null;
    // Display-only clone, never registered in the workspace.
    const work = global.clone(null);
    work.name = `${global.name} · round ${k}`;
    detectChemSemTypes(work);
    this.setStepWork(k, work);
    return work;
  }

  // Selecting rows in a round's grid is only staging; the round narrows once this commits it.
  private subsetStepBySelection(k: number): void {
    const w = this.stepState[k - 1]?.df;
    if (!w) return;
    const subset = cloneSubsetByRows(w,
      `${SELECT_ROWS_OR_FILTER} to use only those ${this.opts.noun} in round ${k}.`);
    if (!subset) return;
    this.setStepWork(k, subset);
    this.stepState[k - 1].committed = true;
    this.render();
    this.deps.refreshCfgRibbon();
    this.deps.viewerHost.deferredFilterReset(subset);
  }

  // Drops the clone entirely so the round re-derives from "All rounds" lazily.
  private useAllForStep(k: number): void {
    this.setStepWork(k, null);
    this.stepState[k - 1].committed = false;
    this.render();
    this.deps.refreshCfgRibbon();
    const w = this.stepState[k - 1]?.df;
    if (w) this.deps.viewerHost.deferredFilterReset(w);
  }

  /** Commits `newDf` into the tracked table input and closes the tables it replaces. `keepIfSameAs`
   * skips closing a previous table that is the same tracked table as it. renderGrid() is called
   * explicitly because assignTableInput's null-then-real assignment doesn't reliably reach here
   * through the input's own onChanged. */
  private commitTableSwap(
    newDf: DG.DataFrame, prev: DG.DataFrame | null, prevFresh: DG.DataFrame | null,
    keepIfSameAs?: DG.DataFrame | null,
  ): void {
    assignTableInput(this.opts.input, newDf);
    this.renderGrid();
    this.deps.viewerHost.deferredFilterReset(newDf);
    this.deps.refreshValidation();
    const shouldClose = (t: DG.DataFrame | null): boolean =>
      keepIfSameAs === undefined || !isSameTrackedTable(t, keepIfSameAs);
    if (shouldClose(prev)) closeTableSafe(prev, 'prev subset');
    if (shouldClose(prevFresh)) closeTableSafe(prevFresh, 'prev fresh clone');
  }

  private doSubsetBySelection(): void {
    const df = this.opts.input.value;
    if (!df) {
      if (this.opts.noTableMsg) grok.shell.info(this.opts.noTableMsg);
      return;
    }
    const subset = cloneSubsetByRows(df,
      `Select rows in the ${this.opts.noun} grid (Ctrl/Shift+click) or apply a filter first.`);
    if (!subset) return;
    if (!isSameTrackedTable(df, this.state.prev)) this.state.original = df;
    const prev = this.state.prev;
    const prevFresh = this.state.freshClone;
    this.state.prev = subset;
    this.state.freshClone = null;
    this.commitTableSwap(subset, prev, prevFresh, df);
  }

  private doRestoreFullTable(): void {
    const clearOrInform = (current: DG.DataFrame | null): void => {
      if (current && current.filter.trueCount < current.rowCount) current.filter.setAll(true, true);
      else grok.shell.info(`The full ${this.opts.noun} library is already in use.`);
    };
    const orig = this.state.original;
    if (!orig) return clearOrInform(this.opts.input.value);
    if (this.state.freshClone && isSameTrackedTable(this.opts.input.value, this.state.freshClone))
      return clearOrInform(this.state.freshClone);
    const prev = this.state.prev;
    const prevFresh = this.state.freshClone;
    this.state.prev = null;
    // A distinctly-named clone, not `orig` itself: reusing it carries forward per-column tags (e.g.
    // CHEM_APPLY_FILTER_SYNC) that can hang a re-run substructure search, and even a tag-free clone
    // sharing its name gets the platform's remembered filter/sketch state reapplied.
    const fresh = orig.clone(null);
    // `orig` can itself be an earlier "Use all" clone, so strip before appending.
    fresh.name = `${orig.name.replace(/ \(full\)$/, '')} (full)`;
    this.state.freshClone = fresh;
    this.commitTableSwap(fresh, prev, prevFresh);
  }

  // TabControl has no add/remove-pane API, only clear()+re-addPane().
  private buildStepTabs(initialStep = 0): void {
    this.stepTabsSub?.unsubscribe();
    // Close viewers BEFORE wiping stepTabsHost, or their hosts leave the DOM while still registered
    // and the Viewer instances are orphaned instead of released.
    for (const ph of this.paneHosts) if (ph) this.deps.viewerHost.close(ph.gridHost);
    this.lastMounted.clear();
    this.stepTabsHost.innerHTML = '';
    this.stepDots.length = 0;
    this.paneHosts = [];
    const tc = ui.tabControl(null, false);
    tc.root.style.width = '100%';
    tc.root.style.flex = '1 1 0';
    tc.root.style.minHeight = '0';
    tc.root.style.overflow = 'hidden';
    // Recorded at its own fixed index k, so it works whether addPane's factory runs eagerly or lazily.
    const makePaneContent = (k: number): () => HTMLElement => () => {
      const barHost = ui.div([], {style: {display: 'flex', alignItems: 'center', gap: '8px', flex: '0 0 auto',
        padding: '4px 8px 5px', borderBottom: '1px solid var(--grey-2)'}});
      const gridHost = ui.div([], {style: {display: 'flex', flexDirection: 'column', flex: '1 1 0',
        minHeight: '0', overflow: 'hidden'}});
      this.paneHosts[k] = {barHost, gridHost};
      return ui.div([barHost, gridHost], {style: {
        height: '100%', display: 'flex', flexDirection: 'column', overflow: 'hidden'}});
    };
    const allPane = tc.addPane('All rounds', makePaneContent(0));
    ui.tooltip.bind(allPane.header, 'Narrow this component for one round only. "All rounds" is the full ' +
      'library used by every round; pick a round to restrict just that round. Per-round building-block ' +
      'subsets apply in depth-first / reagents mode; in breadth-first mode a round draws from all earlier ' +
      'products, so a BB subset has no effect. Resets when you swap the input file (in-range overrides ' +
      'survive a round-count change).');
    this.stepDots.push(null);
    for (let k = 1; k <= this.roundCount(); k++) {
      const pane = tc.addPane(`Round ${k}`, makePaneContent(k));
      // Absolute, with a positive left offset (a negative one bleeds into the flush adjacent tab)
      // and a fixed top (the header box shrinks ~7px when selected, so a percentage would drift).
      const dot = ui.div([], {style: {...CHANGED_DOT_STYLE, position: 'absolute', left: '5px', top: '12px',
        display: 'none'}});
      pane.header.style.position = 'relative';
      pane.header.appendChild(dot);
      this.stepDots.push(dot);
    }
    this.stepTabsSub = tc.onTabChanged.subscribe(() => {
      // Insertion order above makes the pane index equal selStep, so the label never needs parsing.
      const idx = tc.currentPane ? tc.panes.indexOf(tc.currentPane) : -1;
      this.selStep = idx < 0 ? 0 : idx;
      this.render();
    });
    // Each rebuild retires the prior instance, but the last one needs view close to reach it.
    this.deps.view.subs.push(this.stepTabsSub);
    this.stepTabsHost.appendChild(tc.root);
    // Select explicitly — onTabChanged may not fire if the target is already the control's default.
    const clamped = Math.min(Math.max(0, initialStep), this.roundCount());
    this.selStep = clamped;
    const target = tc.panes[clamped] ?? allPane;
    if (target !== tc.currentPane) tc.currentPane = target;
    this.render();
  }

  private renderBar(): void {
    const barHost = this.paneHosts[this.selStep]?.barHost;
    if (!barHost) return;
    barHost.innerHTML = '';
    const hintEl = (t: string): HTMLElement =>
      ui.divText(t, {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '1 1 auto', marginRight: '4px'}});
    const filterIcon = ui.iconFA('filter',
      () => {this.filtersOn = !this.filtersOn; this.renderGrid(); this.renderBar();},
      this.filtersOn ? 'Hide filters' : 'Show filters');
    filterIcon.style.cursor = 'pointer';
    filterIcon.style.padding = '2px 5px';
    filterIcon.style.flex = '0 0 auto';
    filterIcon.style.color = this.filtersOn ? 'var(--blue-2)' : 'var(--grey-5)';
    if (this.selStep === 0) {
      // Only warn when the click actually swapped the table AND there were overrides to lose — a
      // no-op click already gets its own toast from doSubsetBySelection/doRestoreFullTable.
      const doGlobalAction = (action: () => void, clearedSuffix: string): void => {
        const hadOverride = this.hasAnyOverride();
        const prevValue = this.opts.input.value;
        action();
        if (hadOverride && !isSameTrackedTable(prevValue, this.opts.input.value)) {
          grok.shell.info(`Per-round ${this.opts.noun} overrides were cleared — every round now uses the ` +
            `${clearedSuffix}.`);
        }
      };
      const btn = ui.link('Subset by selection', () => doGlobalAction(
        () => this.doSubsetBySelection(), `new ${this.opts.noun} subset`));
      ui.tooltip.bind(btn, `Replace the ${this.opts.noun} library with only the selected rows, or — if nothing ` +
        `is selected — the currently filtered rows (applies to every round). Click "Use all" to restore the ` +
        `full set.`);
      const useAll = ui.link('Use all', () => doGlobalAction(
        () => this.doRestoreFullTable(), `full ${this.opts.noun} library again`));
      ui.tooltip.bind(useAll, `Restore the full ${this.opts.noun} library (undo "Subset by selection").`);
      barHost.append(hintEl(`Full ${this.opts.noun} library — used by every round unless a round overrides it.`),
        filterIcon, btn, useAll);
    } else {
      const entry = this.stepState[this.selStep - 1];
      const w = entry?.df;
      // While suppressed, `w` is still the committed subset clone — reporting its row count as
      // "all" would be both the wrong number and the wrong word.
      const suppressed = !!entry?.committed &&
        bbOverrideSuppressedInBreadth(this.opts.idx === 1, this.deps.currentMode());
      let statusText = '';
      if (w) {
        if (this.hasOverride(this.selStep))
          statusText = `using ${w.rowCount} / ${this.opts.input.value?.rowCount ?? w.rowCount}`;
        else {
          const fullCount = this.opts.input.value?.rowCount ?? w.rowCount;
          statusText = `all ${fullCount}` + (suppressed ? ' (subset ignored in breadth-first)' : '');
        }
      }
      const status = ui.divText(statusText,
        {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto'}});
      const btn = ui.link('Subset by selection', () => this.subsetStepBySelection(this.selStep));
      ui.tooltip.bind(btn, `Narrow round ${this.selStep} to only the selected rows (Ctrl/Shift+click), or — if ` +
        `nothing is selected — the currently filtered rows. Click "Use all" to go back to the full ` +
        `${this.opts.noun} library.`);
      const useAll = ui.link('Use all', () => this.useAllForStep(this.selStep));
      ui.tooltip.bind(useAll, `Undo "Subset by selection" so round ${this.selStep} uses the full ${this.opts.noun} ` +
        `library (same as "All rounds").`);
      barHost.append(
        hintEl(`${SELECT_ROWS_OR_FILTER}, then "Subset by selection" to use only those ${this.opts.noun} in ` +
          `round ${this.selStep}.`),
        status, filterIcon, btn, useAll);
    }
  }

  private renderGrid(): void {
    const gridHost = this.paneHosts[this.selStep]?.gridHost;
    if (!gridHost) return;
    const currentDf = this.selStep === 0 ? this.opts.input.value : this.stepClone(this.selStep);
    if (!currentDf) {
      this.lastMounted.delete(gridHost);
      // Close before wiping the DOM: closing after innerHTML='' hands the viewer a detached
      // container, which throws deep in the Dart-side close path.
      this.deps.viewerHost.close(gridHost);
      gridHost.innerHTML = '';
      gridHost.appendChild(ui.div(
        [ui.divText(this.opts.emptyMsg ?? `No ${this.opts.noun} table selected.`,
          {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}})],
        {style: {overflowY: 'auto', overflowX: 'hidden', height: '100%'}}));
      if (this.selStep === 0) this.opts.badge?.refresh(null);
      return;
    }
    // Round clones are never registered, so they have no `.id`, and stepClone() names every fresh
    // clone identically — but it does guarantee reference stability, so use the reference itself.
    const identity = this.selStep === 0 ? currentDf.id : currentDf;
    const prevMounted = this.lastMounted.get(gridHost);
    if (prevMounted && prevMounted.identity === identity && prevMounted.filtersOn === this.filtersOn) {
      if (this.selStep === 0) this.opts.badge?.refresh(currentDf.rowCount);
      return;
    }
    this.deps.viewerHost.mountDf(gridHost, currentDf, this.filtersOn);
    this.lastMounted.set(gridHost, {identity, filtersOn: this.filtersOn});
    if (this.selStep === 0) this.opts.badge?.refresh(currentDf.rowCount);
  }
}

/** A round with no narrowed component falls back to the global set; the whole result is undefined
 * when nothing is overridden. */
export function buildPerRoundOverrides(panels: DataPanel[], cfg: EnumeratorConfig): PerRoundOverride[] | undefined {
  const overrides: PerRoundOverride[] = [];
  let any = false;
  for (let r = 0; r < clampRounds(cfg.enumeration.num_rounds); r++) {
    const o: PerRoundOverride = {};
    for (const panel of panels)
      if (panel.applyOverrideForRound(r + 1, o, cfg)) any = true;

    overrides.push(o);
  }
  return any ? overrides : undefined;
}

/** Shared by the Strategy summary and the Preview recap. */
export function overrideCountFor(
  overrides: PerRoundOverride[] | undefined, mode: Mode, r: number, key: DataKey,
): number | null {
  if (bbOverrideSuppressedInBreadth(key === 'buildingBlocks', mode)) return null;
  const list = overrides?.[r - 1]?.[key];
  return list ? list.length : null;
}

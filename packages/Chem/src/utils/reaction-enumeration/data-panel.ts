import {Subscription} from 'rxjs';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EnumeratorConfig} from './config';
import {PerRoundOverride} from './enumerate';
import {MountedViewerRegistry} from './viewer-mount';
import {detectChemSemTypes, MAX_ROUNDS, OVERRIDE_DOT_COLOR} from './enumerator-app';

type Mode = 'depth' | 'breadth' | 'reagents';
type DataKey = 'templates' | 'buildingBlocks' | 'reagents';

export type TabBadge = {el: HTMLSpanElement; refresh: (n: number | null) => void};
export const makeTabBadge = (): TabBadge => {
  const el = document.createElement('span');
  el.className = 'chem-enum-tab-badge';
  return {el, refresh: (n: number | null) => {
    el.textContent = n != null ? String(n) : '';
    el.style.display = n != null ? '' : 'none';
  }};
};

// Table names aren't unique in a workspace, so `.name` can't tell "same table" from "same-named
// table" — compare `.id` instead, which the platform assigns once a table is registered.
function isSameTrackedTable(a: DG.DataFrame | null | undefined, b: DG.DataFrame | null | undefined): boolean {
  return a?.id != null && a.id === b?.id;
}

// ui.input.table is a ChoiceInput over the workspace tables list; setting .value to a DataFrame
// not in that list is silently rejected. Register via grok.shell.addTable first, and dip through
// null to drop any cached pointer-equality with the previous selection.
function assignTableInput(input: DG.InputBase<DG.DataFrame | null>, df: DG.DataFrame): void {
  grok.shell.addTable(df);
  try {input.value = null;} catch {/* nullable: false rejects */}
  input.value = df;
}

// Shared wording for the two identically-phrased "how to subset" prompts (per-step's empty-
// selection message and its hint text) — kept as one constant so the two stay in sync.
const SELECT_ROWS_OR_FILTER = 'Select rows (Ctrl/Shift+click) or apply a filter';

// Shared clone/rename/detect core for both the global ("All steps") and per-step "Subset by
// selection" actions — returns null (after an info toast) when there's nothing to subset.
// Falls back to the active filter when nothing is explicitly selected, so filtering alone is
// enough to subset without an extra manual "select all" step — explicit selection still wins
// when both are present, since it's the more deliberate action.
function cloneSubsetByRows(df: DG.DataFrame, emptyMsg: string): DG.DataFrame | null {
  const sel = df.selection;
  const mask = sel.anyTrue ? sel : df.filter;
  if (!mask.anyTrue) {
    grok.shell.info(emptyMsg);
    return null;
  }
  // Also covers the trivial "filter matches every row" case — nothing left to narrow either way.
  if (!mask.anyFalse) {
    grok.shell.info('All rows are selected — nothing to subset.');
    return null;
  }
  const subset = df.clone(mask);
  subset.name = `${df.name} (subset, ${subset.rowCount}/${df.rowCount} rows)`;
  detectChemSemTypes(subset);
  // df.clone(mask) leaves the SOURCE df's own .selection set to `mask` as a side effect — left
  // uncleared, a later filter-only "Subset by selection" on the restored original table picks up
  // this stale selection instead of the new filter, since selection is checked first above.
  df.selection.setAll(false, false);
  // Same clone-carries-BitSet-state issue on the other side: the new subset's OWN .filter comes
  // back all-false when cloning off a filtered source — reset it so the subset isn't born with
  // every row hidden.
  subset.filter.setAll(true, false);
  return subset;
}

// Breadth-first ignores per-round BB overrides entirely (see `eligibleSmiles` in enumerate.ts) —
// one shared rule instead of every caller re-encoding the same check.
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
  // How to fold this component's per-round work df into a PerRoundOverride fragment.
  apply: (o: PerRoundOverride, work: DG.DataFrame, cfg: EnumeratorConfig) => void;
}

// Shared across every DataPanel instance (all three read the same config/mode/rounds and mediate
// through the same orchestrator-level refresh functions).
export interface DataPanelDeps {
  view: DG.View;
  viewerHost: MountedViewerRegistry;
  getConfig: () => EnumeratorConfig;
  currentMode: () => Mode;
  currentRounds: () => number;
  refreshValidation: () => void;
  refreshCfgRibbon: () => void;
}

/** Single-grid per-component panel with a per-step strip. Each data tab shows ONE grid plus a
 * horizontal step strip: "All steps" shows the full library (with the global "Subset by
 * selection"); "Step k" shows a display-only clone whose row selection is that round's subset.
 * Switching chips swaps what the single grid displays — no second grid. */
export class DataPanel {
  readonly panel: HTMLElement;

  private selStep = 0; // 0 = All steps (full library); 1..rounds = that round's subset
  private filtersOn = false; // funnel toggle: show the standard Datagrok filters panel next to the grid
  private currentDf: DG.DataFrame | null = null;
  // host -> last (table identity, filtersOn) successfully mounted there. "Subset by selection"/"Use
  // all" both trigger a render twice (see the comment above assignTableInput) — renderGrid no-ops
  // the second call once it sees nothing actually changed, comparing by `.id`, not reference.
  private readonly lastMounted = new Map<HTMLElement, {id: string; filtersOn: boolean}>();
  private readonly state: SubsetState = {prev: null, original: null, freshClone: null};
  // Per-round state (index k-1 = round k) — one array instead of parallel ones, so df/sub/committed
  // can't drift out of sync. `committed` is an explicit flag, NOT inferred from row-count: a step's
  // clone can silently drift from the global table if rows are added/removed in place (no onChanged
  // fires), and inferring "override" from that drift was confirmed to let deleted rows resurrect in
  // a round's reactant pool. Set only by subsetStepBySelection, cleared only by useAllForStep.
  private readonly stepState: StepState[] = [];
  // Step selector: a real ui.tabControl. Each pane owns its own persistent barHost/gridHost (built
  // once via its addPane factory) — relocating a live grid between panes corrupts it, so nothing is
  // shared. `paneHosts[selStep]` gives renderBar/renderGrid the current pane's own hosts.
  private readonly stepTabsHost: HTMLElement;
  private stepTabsSub: Subscription | null = null;
  private readonly stepDots: (HTMLElement | null)[] = []; // index k = 1..rounds; index 0 unused
  private paneHosts: {barHost: HTMLElement; gridHost: HTMLElement}[] = [];

  constructor(private readonly opts: DataPanelOpts, private readonly deps: DataPanelDeps) {
    this.stepTabsHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', display: 'flex',
      flexDirection: 'column', overflow: 'hidden'}});
    // min-height:0 — same reason as tabPanel()'s wrapper: this is a flex child of the platform's
    // own .d4-tab-content, which without it refuses to shrink below content height.
    this.panel = ui.div([this.stepTabsHost], {style: {
      height: '100%', display: 'flex', flexDirection: 'column', minHeight: '0',
      background: 'var(--white)', overflow: 'hidden'}});
    this.buildStepTabs(0); // also mounts the grid once hosts are in the DOM
  }

  render(): void { this.renderGrid(); this.renderBar(); this.updateDots(); }
  refreshDisplay(): void { this.renderBar(); this.updateDots(); }

  onTableChanged(): void {
    const df = this.opts.input.value;
    if (df) detectChemSemTypes(df);
    // A match against neither prev nor freshClone means the input changed from outside either of
    // them — e.g. a different file — so the old subset bookkeeping no longer applies and must be
    // cleared, or a later "Use all" would restore the wrong table's `original`.
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

  // Deliberately does NOT drop stepState entries beyond the new round count — Dart int inputs fire
  // onChanged per keystroke, so typing "10" over "5" transiently passes through "1", which would
  // eagerly destroy committed overrides on steps 2-5 before the user finishes typing.
  onRoundsChanged(): void { this.buildStepTabs(this.selStep); }

  // Applies this component's round-r override (if any) onto `out`; returns whether it did. Lets
  // buildPerRoundOverrides reach each panel's own stepState without a shared indexed array.
  applyOverrideForRound(r: number, out: PerRoundOverride, cfg: EnumeratorConfig): boolean {
    const entry = this.stepState[r - 1];
    if (!entry?.committed || !entry.df) return false; // !entry.df shouldn't happen if committed
    this.opts.apply(out, entry.df, cfg);
    return true;
  }

  // Whether any round has a custom subset for this component — drives its ribbon chip's dot.
  hasAnyOverride(): boolean {
    for (let k = 1; k <= this.roundCount(); k++) if (this.hasOverride(k)) return true;
    return false;
  }

  // A step "has an override" only once stepState[k-1].committed is explicitly set — see the flag's
  // own comment above for why row-count inference was wrong.
  private hasOverride(k: number): boolean {
    if (!this.stepState[k - 1]?.committed) return false;
    // BB overrides don't apply in breadth-first (see bbOverrideSuppressedInBreadth) — don't show the dot.
    if (bbOverrideSuppressedInBreadth(this.opts.idx === 1, this.deps.currentMode())) return false;
    return true;
  }

  // Capped defensively regardless of validation state — an invalid (too-high) input value still
  // blocks Run via `validate()`, but must not make buildStepTabs try to build hundreds of tabs.
  private roundCount(): number {
    return Math.min(MAX_ROUNDS, Math.max(1, this.deps.currentRounds()));
  }

  private updateDots(): void {
    for (let k = 1; k <= this.stepDots.length - 1; k++) {
      const dot = this.stepDots[k];
      if (dot) dot.style.display = this.hasOverride(k) ? '' : 'none';
    }
  }

  // Shared bookkeeping behind every place a step's clone changes: unsubscribe the old selection
  // listener, subscribe the new one (if any), and replace the round's state record. Preserves the
  // existing `committed` value — callers that need to change it do so explicitly right after.
  private setStepWork(k: number, work: DG.DataFrame | null): void {
    this.stepState[k - 1]?.sub?.unsubscribe();
    const sub = work ? work.onSelectionChanged.subscribe(() => { this.updateDots(); this.renderBar(); }) : null;
    if (sub) this.deps.view.subs.push(sub);
    const entry = this.stepState[k - 1];
    if (entry) { entry.df = work; entry.sub = sub; }
    else this.stepState[k - 1] = {df: work, sub, committed: false};
  }

  private stepClone(k: number): DG.DataFrame | null {
    const existing = this.stepState[k - 1]?.df;
    if (existing) return existing;
    const global = this.opts.input.value;
    if (!global) return null;
    // Display-only clone, never registered in the workspace; its selection carries the subset.
    const work = global.clone(null);
    work.name = `${global.name} · round ${k}`;
    detectChemSemTypes(work);
    this.setStepWork(k, work);
    return work;
  }

  // Per-step mirror of the global "Subset by selection": selecting rows in a step's grid is just
  // staging, the round only narrows once this commits it by swapping in a new clone.
  private subsetStepBySelection(k: number): void {
    const w = this.stepState[k - 1]?.df;
    if (!w) return; // grid must be mounted (via stepClone) before this button is reachable
    const subset = cloneSubsetByRows(w,
      `${SELECT_ROWS_OR_FILTER} to use only those ${this.opts.noun} in round ${k}.`);
    if (!subset) return;
    this.setStepWork(k, subset);
    this.stepState[k - 1].committed = true;
    this.renderGrid(); this.renderBar(); this.updateDots();
    this.deps.refreshCfgRibbon(); // this component's ribbon chip dot depends on hasAnyOverride()
    this.deps.viewerHost.deferredFilterReset(subset);
  }

  // Undo: drop the clone entirely so the step falls back to (re-derives from) "All steps" lazily.
  private useAllForStep(k: number): void {
    this.setStepWork(k, null);
    this.stepState[k - 1].committed = false;
    this.renderGrid(); this.renderBar(); this.updateDots();
    this.deps.refreshCfgRibbon(); // this component's ribbon chip dot depends on hasAnyOverride()
    // stepClone(k) can inherit the global table's active filter — reset it.
    const w = this.stepState[k - 1]?.df;
    if (w) this.deps.viewerHost.deferredFilterReset(w);
  }

  // assignTableInput's null-then-real two-step assignment does NOT reliably re-render via the
  // input's own onChanged -> onTableChanged reaction alone — relying on that path alone left the
  // grid empty after "Subset by selection". The explicit calls below are intentionally redundant
  // with onTableChanged's own re-render — correctness over the extra render.
  private doSubsetBySelection(): void {
    const df = this.opts.input.value;
    if (!df) {
      if (this.opts.noTableMsg) grok.shell.info(this.opts.noTableMsg);
      return;
    }
    const subset = cloneSubsetByRows(df,
      `Select rows in the ${this.opts.noun} grid (Ctrl/Shift+click) or apply a filter first.`);
    if (!subset) return;
    if (!isSameTrackedTable(df, this.state.prev)) this.state.original = df; // remember the user's own table for "Use all"
    const prev = this.state.prev;
    const prevFresh = this.state.freshClone;
    this.state.prev = subset;
    this.state.freshClone = null; // no longer showing a "Use all"-produced copy
    assignTableInput(this.opts.input, subset);
    this.renderGrid();
    // The just-mounted Filters viewer can clobber subset.filter to all-false (e.g. zero-variance
    // column after subsetting) — reset it once that settles.
    this.deps.viewerHost.deferredFilterReset(subset);
    this.deps.refreshValidation();
    // Close the previous subset only after the input has switched away from it.
    if (prev && !isSameTrackedTable(prev, df))
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
    if (prevFresh && !isSameTrackedTable(prevFresh, df))
      try {grok.shell.closeTable(prevFresh);} catch (e) {console.warn(`Could not close prev fresh clone: ${e}`);}
  }

  // Undo "Subset by selection": put the remembered full table back into the input and close the
  // subset clone. No-op (with a hint) when the full library is already in use.
  private doRestoreFullTable(): void {
    // Shared "no swap needed" case: just clear an active filter directly, or say so if there isn't one.
    const clearOrInform = (current: DG.DataFrame | null): void => {
      if (current && current.filter.trueCount < current.rowCount) current.filter.setAll(true, true);
      else grok.shell.info(`The full ${this.opts.noun} library is already in use.`);
    };
    const orig = this.state.original;
    if (!orig) return clearOrInform(this.opts.input.value); // no subset was ever taken from this table
    if (this.state.freshClone && isSameTrackedTable(this.opts.input.value, this.state.freshClone))
      return clearOrInform(this.state.freshClone); // already showing a "Use all" copy — no need to re-clone
    const prev = this.state.prev;
    const prevFresh = this.state.freshClone;
    this.state.prev = null;
    // "Use all" swaps in a fresh, distinctly-named clone of `orig` rather than reusing it directly.
    // Reusing `orig` carries forward per-column tags (e.g. chem's CHEM_APPLY_FILTER_SYNC) that can hang
    // a re-run substructure search — and even a tag-free clone NAMED the same as `orig` still gets the
    // platform's own remembered filter/sketch state (keyed by table+column name) reapplied to it.
    // clone(null) plus a distinct name sidesteps both. `orig` stays untouched so repeated clicks keep
    // deriving from the same known-good source.
    const fresh = orig.clone(null);
    // Strip a prior "(full)" suffix instead of appending blindly — `orig` can itself be an earlier
    // "Use all" clone (subsetting from it re-points state.original at it), and re-suffixing every cycle
    // would otherwise grow "X (full) (full) (full)...".
    fresh.name = `${orig.name.replace(/ \(full\)$/, '')} (full)`;
    this.state.freshClone = fresh;
    assignTableInput(this.opts.input, fresh);
    this.renderGrid();
    // `fresh` can still inherit a stale filter bitset from `orig` — reset it once the mount settles.
    this.deps.viewerHost.deferredFilterReset(fresh);
    this.deps.refreshValidation();
    if (prev)
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
    if (prevFresh)
      try {grok.shell.closeTable(prevFresh);} catch (e) {console.warn(`Could not close prev fresh clone: ${e}`);}
  }

  // (Re)build the step tab strip, landing on `initialStep` (clamped). TabControl has no
  // add/remove-pane API, only clear()+re-addPane(); `tc.panes` insertion order lines up with selStep.
  private buildStepTabs(initialStep = 0): void {
    this.stepTabsSub?.unsubscribe();
    // Close each pane's mounted viewer(s) BEFORE wiping stepTabsHost — otherwise their gridHost
    // divs are dropped from the DOM while still registered in `mountedViewers`, orphaning the
    // Viewer instances (never closed) instead of releasing their Dart-side resources.
    for (const ph of this.paneHosts) if (ph) this.deps.viewerHost.close(ph.gridHost);
    // Every gridHost below is about to be discarded and rebuilt from scratch — drop renderGrid's
    // per-host dedup entries too, or they'd hold the old (now-detached) elements alive forever.
    this.lastMounted.clear();
    this.stepTabsHost.innerHTML = '';
    this.stepDots.length = 0;
    this.paneHosts = [];
    const tc = ui.tabControl(null, false);
    // tc.root must fill available height, not size to its header-strip content — each pane's own
    // content div (built below) needs the real space to lay its grid out in.
    tc.root.style.width = '100%';
    tc.root.style.flex = '1 1 0';
    tc.root.style.minHeight = '0';
    tc.root.style.overflow = 'hidden';
    // Builds one pane's persistent content (barHost+gridHost), recorded in paneHosts at its OWN
    // fixed index k (not push-order) — works whether TabControl's addPane factory runs eagerly or
    // lazily, since position k always lands at paneHosts[k] regardless of firing order.
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
      // Position dot absolutely (header stacks children vertically, so marginLeft lands below the
      // label). Positive left offset only — a negative one bleeds into the adjacent tab, since tabs
      // sit flush. Fixed top offset, not top:50% — the header's box shrinks ~7px when selected
      // (underline indicator), so a percentage-based center drifts on selection.
      const dot = ui.div([], {style: {position: 'absolute', left: '5px', top: '12px',
        width: '6px', height: '6px', borderRadius: '50%',
        background: OVERRIDE_DOT_COLOR, display: 'none'}});
      pane.header.style.position = 'relative';
      pane.header.appendChild(dot);
      this.stepDots.push(dot);
    }
    const renderCurrent = (): void => { this.renderGrid(); this.renderBar(); this.updateDots(); };
    this.stepTabsSub = tc.onTabChanged.subscribe(() => {
      // Index into tc.panes IS selStep (0 = "All steps", k = "Step k", by insertion order above)
      // — no need to parse the pane's label, which would break if it were ever reworded.
      const idx = tc.currentPane ? tc.panes.indexOf(tc.currentPane) : -1;
      this.selStep = idx < 0 ? 0 : idx;
      renderCurrent();
    });
    // Each rebuild unsubscribes the prior instance, but the LAST one has no later rebuild to retire
    // it — track it in view.subs too so view close still reaches it.
    this.deps.view.subs.push(this.stepTabsSub);
    this.stepTabsHost.appendChild(tc.root);
    // Select explicitly — onTabChanged may not fire if the target is already the control's default.
    const clamped = Math.min(Math.max(0, initialStep), this.roundCount());
    this.selStep = clamped;
    const target = tc.panes[clamped] ?? allPane;
    if (target !== tc.currentPane) tc.currentPane = target;
    renderCurrent();
  }

  private renderBar(): void {
    const barHost = this.paneHosts[this.selStep]?.barHost;
    if (!barHost) return; // pane not built yet (shouldn't happen once buildStepTabs has run once)
    barHost.innerHTML = '';
    const hintEl = (t: string): HTMLElement =>
      ui.divText(t, {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '1 1 auto', marginRight: '4px'}});
    // Funnel toggle — shows the standard Datagrok filters panel for the visible grid. Off by default.
    const filterIcon = ui.iconFA('filter',
      () => { this.filtersOn = !this.filtersOn; this.renderGrid(); this.renderBar(); },
      this.filtersOn ? 'Hide filters' : 'Show filters');
    filterIcon.style.cursor = 'pointer';
    filterIcon.style.padding = '2px 5px';
    filterIcon.style.flex = '0 0 auto';
    filterIcon.style.color = this.filtersOn ? 'var(--blue-2)' : 'var(--grey-5)';
    if (this.selStep === 0) {
      // Warn only when the click actually swaps the table AND overrides existed to lose — a no-op
      // click already gets its own info toast from doSubsetBySelection/doRestoreFullTable.
      const doGlobalAction = (action: () => void, clearedSuffix: string): void => {
        const hadOverride = this.hasAnyOverride();
        const prevValue = this.opts.input.value;
        action();
        if (hadOverride && this.opts.input.value !== prevValue) {
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
      const w = this.stepState[this.selStep - 1]?.df;
      const status = ui.divText(
        w ? (this.hasOverride(this.selStep) ? `using ${w.rowCount} / ${this.opts.input.value?.rowCount ?? w.rowCount}` :
          `all ${w.rowCount}`) : '',
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
    if (!gridHost) return; // pane not built yet (shouldn't happen once buildStepTabs has run once)
    this.currentDf = this.selStep === 0 ? this.opts.input.value : this.stepClone(this.selStep);
    if (!this.currentDf) {
      this.lastMounted.delete(gridHost);
      // Close mounted viewers BEFORE wiping the DOM — closing after innerHTML='' hands the viewer
      // a detached container, which throws ("Cannot read properties of null") deep in the Dart-side
      // close path. Under rapid re-triggering (e.g. the filter icon clicked several times in quick
      // succession) that cascades badly enough to crash the tab's renderer.
      this.deps.viewerHost.close(gridHost);
      gridHost.innerHTML = '';
      // gridHost itself stays overflow:hidden (correct once a real grid — which scrolls itself —
      // is mounted there); this empty-state text gets its own scrollable wrapper instead, since it
      // can outgrow a short window on its own.
      gridHost.appendChild(ui.div(
        [ui.divText(this.opts.emptyMsg ?? `No ${this.opts.noun} table selected.`,
          {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}})],
        {style: {overflowY: 'auto', overflowX: 'hidden', height: '100%'}}));
      if (this.selStep === 0) this.opts.badge?.refresh(null);
      return;
    }
    // Step override clones (selStep != 0) never get a registered `.id` — they're tracked by
    // stepState, not SubsetState — so skip that lookup and fall back to name for those.
    const key = {id: (this.selStep === 0 ? this.currentDf.id : null) ?? this.currentDf.name, filtersOn: this.filtersOn};
    const prevMounted = this.lastMounted.get(gridHost);
    if (prevMounted && prevMounted.id === key.id && prevMounted.filtersOn === key.filtersOn) {
      if (this.selStep === 0) this.opts.badge?.refresh(this.currentDf.rowCount);
      return;
    }
    this.deps.viewerHost.mountDf(gridHost, this.currentDf, this.filtersOn); // mountDf itself closes-then-clears the host
    this.lastMounted.set(gridHost, key);
    if (this.selStep === 0) this.opts.badge?.refresh(this.currentDf.rowCount);
  }
}

// A step's clone IS the subset once committed via "Subset by selection" — no deriving from a live
// .selection bitset at run time. A round with no narrowed component falls back to the global set;
// the whole result is undefined when nothing is overridden.
export function buildPerRoundOverrides(panels: DataPanel[], cfg: EnumeratorConfig): PerRoundOverride[] | undefined {
  const overrides: PerRoundOverride[] = [];
  let any = false;
  for (let r = 0; r < cfg.enumeration.num_rounds; r++) {
    const o: PerRoundOverride = {};
    for (const panel of panels) {
      if (panel.applyOverrideForRound(r + 1, o, cfg)) any = true;
    }
    overrides.push(o);
  }
  return any ? overrides : undefined;
}

// Shared by the Strategy summary and the Preview recap, so the "does round r have a custom
// subset" logic lives in exactly one place.
export function overrideCountFor(
  overrides: PerRoundOverride[] | undefined, mode: Mode, r: number, key: DataKey,
): number | null {
  if (bbOverrideSuppressedInBreadth(key === 'buildingBlocks', mode)) return null;
  const list = overrides?.[r - 1]?.[key];
  return list ? list.length : null;
}

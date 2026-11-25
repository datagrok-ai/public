import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {Engine, isServerEngine} from '../engines/engine';
import {ModelTask} from '../generated/db';
import {METRIC_LABELS, MetricId} from '../metrics/metrics';
import {isRecord} from '../preparation/preparation-options';
import {DEFAULT_BUDGET_SECONDS, defaultMetric, defaultUncheckReason, exclusionOf, isQuickLimit, isQuickMethod,
  MethodProbabilities, probabilitiesOf, SEARCH_METRICS, SearchChoice} from '../training/search';
import {LimitInputs, StoredLimit, TRIALS} from './limit-inputs';
import {NOT_SEARCHED, STATUS_TOOLTIPS, TuneSettings} from './tune-settings';

/** The key of Find best model's settings in the `grok.userSettings` storage of {@link TuneSettings}. */
export const FIND_KEY = 'Find best';
export const DEFAULT_SEARCH = 'Default';
export const CUSTOM_SEARCH = 'Custom';
type SearchKind = typeof DEFAULT_SEARCH | typeof CUSTOM_SEARCH;
const DEFAULT_LIMIT: StoredLimit = {limit: TRIALS, minutes: DEFAULT_BUDGET_SECONDS / 60, trials: 100};
const NO_METHOD_CHECKED = 'Check at least one method';
const NO_PROBABILITY_METHOD = 'No checked method gives class probabilities';
const SERVER_TOOLTIP = 'Runs on the server; not searched in the browser';
const METHODS_TOOLTIP = 'The methods to try, each with its defaults first';
const SEARCH_TOOLTIP = 'Search settings for the best model';
const SEARCH_TOOLTIPS: Record<SearchKind, string> = {[DEFAULT_SEARCH]: 'Search with the recommended settings',
  [CUSTOM_SEARCH]: 'Choose the metric, the limit and the methods yourself'};

/** A method of the Candidates list: whether it retrains live on the data and runs on the server, and its
 * probabilities. */
export interface MethodEntry { name: string; isLive: boolean; isServer: boolean; probabilities: MethodProbabilities }

/** What is remembered per user: Search (Default or Custom), the metric, the limit and the user's toggles of the
 * methods. */
interface StoredFind extends StoredLimit {
  search: SearchKind;
  metric?: string;
  methods: {[name: string]: boolean};
}

/** The Find best model block under **Method** while Auto is chosen: the status line (Not searched), **Search**
 * (Default / Custom) and the Find bar. Default searches with the task's default metric, the default limit (100 trials)
 * and the methods checked by their defaults. Custom shows the grid of Metric, Limit, Max time or Max trials, and
 * Candidates (the browser methods that can learn from the data, cheapest first; a server method, one without class
 * probabilities while Predict probability is on, or one without them under a probability metric is not listed, its
 * reason in the label's tooltip; a slow one, or a method whose probabilities are a train input (SVM) on a classifier,
 * unchecked unless the time limit is above a minute; SVM unchecked on a table of more than 1000 rows at any limit).
 * Both share the limit equally between the methods. Search, Metric, the limit and the user's method toggles are
 * remembered per user; Default keeps Custom's settings as they are. */
export class FindSettings {
  readonly root = ui.div([], 'forge-find');
  readonly statusText: HTMLSpanElement;
  readonly search: DG.InputBase<string | null>;
  readonly grid = ui.div([], 'forge-tune-grid');
  readonly metric: HTMLSelectElement;
  /** Limit and Max time or Max trials. */
  readonly limit: LimitInputs;
  readonly candidates: DG.MultiChoiceInput<string>;
  readonly candidatesLabel: HTMLElement;
  /** The subscriptions of the block; the view drops them with the form. */
  readonly subs: rxjs.Subscription[] = [];
  private readonly host: {onChanged: () => void};
  private entries: MethodEntry[] = [];
  private readonly overrides: {[name: string]: boolean};
  private task: ModelTask | undefined;
  private readonly storedMetric: string | undefined;
  private readonly saver: {save: () => void; flush: () => void};
  /** Whether the limit was quick at the last listing ({@link isQuickLimit}): the methods' default checks follow it. */
  private isQuick = true;
  private predictProbability = false;
  private rowCount = 0;
  private isLocked = false;
  private isLoaded = false;
  private isRebuilding = false;

  /** [findButton] goes into the Find bar; the grids' label column follows [formLabel]; [host] is told of every change
   * of Search, the methods, the metric or the budget. */
  constructor(findButton: HTMLElement, formLabel: HTMLElement, host: {onChanged: () => void}) {
    this.host = host;
    const stored = FindSettings.load();
    this.storedMetric = stored.metric;
    this.overrides = stored.methods;
    this.saver = TuneSettings.saver(FIND_KEY, () => this.record());
    this.subs.push(new rxjs.Subscription(() => this.flushSave()));
    const status = TuneSettings.statusLine();
    this.statusText = status.text;
    this.statusText.textContent = NOT_SEARCHED;
    ui.tooltip.bind(this.statusText, STATUS_TOOLTIPS[NOT_SEARCHED]);
    this.search = ui.input.radio('Search', {items: [DEFAULT_SEARCH, CUSTOM_SEARCH], value: stored.search,
      nullable: false, onValueChanged: () => {
        ui.setDisplay(this.grid, this.isCustom);
        this.rebuild();
        this.onSearchChanged();
      }});
    this.search.root.classList.add('forge-inline-radio');
    const search = LimitInputs.row(this.search, SEARCH_TOOLTIP);
    LimitInputs.bindOptionTooltips(this.search,
      (option) => SEARCH_TOOLTIPS[option.textContent === CUSTOM_SEARCH ? CUSTOM_SEARCH : DEFAULT_SEARCH]);
    // The metric decides which methods are listed.
    const metric = TuneSettings.metricRow(() => {
      this.rebuild();
      this.onSearchChanged();
    });
    this.metric = metric.select;
    TuneSettings.setMetrics(this.metric, SEARCH_METRICS.classification, defaultMetric('classification'));
    this.candidates = ui.input.multiChoice<string>('Candidates', {items: [], value: [],
      onValueChanged: () => this.onMethodsChanged()});
    this.limit = new LimitInputs(stored, DEFAULT_LIMIT, () => this.onLimitChanged(),
      (trials) => this.trialsProblem(trials));
    this.candidates.root.classList.add('forge-tune-input', 'forge-find-candidates');
    this.candidatesLabel = ui.div(['Candidates'], 'forge-tune-label');
    ui.tooltip.bind(this.candidatesLabel, () => this.candidatesTooltip());
    this.grid.append(metric.label, metric.cell, ...this.limit.elements, this.candidatesLabel,
      ui.div([this.candidates.root], 'forge-tune-cell'));
    ui.setDisplay(this.grid, this.isCustom);
    this.root.append(status.row, ui.div([ui.div([search.label, search.cell], 'forge-tune-grid'), this.grid],
      'forge-tune-auto'), ui.div([findButton], 'forge-tune-bar'));
    this.subs.push(TuneSettings.alignLabels(this.root, formLabel));
    this.rebuild();
    this.updateErrors();
    this.isLoaded = true;
  }

  /** Search is Custom: the grid's settings are used; Default uses the defaults. */
  get isCustom(): boolean {
    return this.search.value === CUSTOM_SEARCH;
  }

  /** The platform inputs of the block: Search, Limit, Max time, Max trials and Candidates. */
  get inputs(): DG.InputBase[] {
    return [this.search, ...this.limit.inputs, this.candidates];
  }

  /** The inputs whose values the search uses now: the limit only in Custom. */
  get usedInputs(): DG.InputBase[] {
    return [this.search, ...this.isCustom ? this.limit.inputs : [], this.candidates];
  }

  /** The checked methods, in the list's order. */
  get checkedMethods(): string[] {
    const checked = this.candidates.value ?? [];
    return this.candidates.items.filter((name) => checked.includes(name));
  }

  /** Why Predict probability cannot be used: no checked method gives class probabilities; null while one does, while
   * none is checked (Find is blocked by that) or while the list is not known yet. */
  get probabilityReason(): string | null {
    const checked = this.checkedMethods;
    return this.entries.length === 0 || checked.length === 0 || this.entries.some((e) => e.probabilities !== 'none' &&
      checked.includes(e.name)) ? null : NO_PROBABILITY_METHOD;
  }

  /** Lists [entries]: the browser methods, cheapest first, then the server ones. */
  setMethods(entries: MethodEntry[]): void {
    this.entries = entries;
    this.rebuild();
    this.updateErrors();
  }

  /** Predict probability is [on]: the methods it keeps out, and the defaults, follow. */
  setPredictProbability(on: boolean): void {
    if (on === this.predictProbability)
      return;
    this.predictProbability = on;
    this.rebuild();
    this.updateErrors();
  }

  /** The table has [rowCount] rows: the default checks follow it ({@link defaultUncheckReason}). */
  setRowCount(rowCount: number): void {
    if (rowCount === this.rowCount)
      return;
    this.rowCount = rowCount;
    this.rebuild();
    this.updateErrors();
  }

  /** Lists the metrics of [task], keeping the chosen one when it applies, else the remembered one, else the default. */
  setTask(task: ModelTask | undefined): void {
    if (task === undefined || task === this.task)
      return;
    this.task = TuneSettings.listMetrics(this.metric, this.task, task, this.storedMetric);
    this.rebuild();
    this.onSearchChanged();
  }

  /** Writes a pending save at once. */
  flushSave(): void {
    this.saver.flush();
  }

  /** Disables Search, the inputs and Metric while the form is locked. */
  setLocked(isLocked: boolean): void {
    this.isLocked = isLocked;
    for (const input of [this.search, this.candidates])
      input.enabled = !isLocked;
    this.limit.setLocked(isLocked);
    this.metric.disabled = isLocked;
  }

  /** Why Find cannot search: no method checked, or the limit used; null when it can. */
  problem(): string | null {
    return this.checkedMethods.length === 0 ? NO_METHOD_CHECKED :
      this.isCustom ? this.limit.problem() : this.trialsProblem(DEFAULT_LIMIT.trials);
  }

  /** What the search is asked: the metric, an equal split, the limit and the checked methods; Default's defaults or
   * Custom's settings. */
  choice(): SearchChoice {
    return {metric: this.metricId, strategy: 'equal', ...this.usedLimit, methods: this.checkedMethods};
  }

  /** The limit the search uses: Custom's, or Default's trials. */
  private get usedLimit(): {budgetSeconds: number} | {maxTrials: number} {
    return this.isCustom ? this.limit.limit() : {maxTrials: DEFAULT_LIMIT.trials};
  }

  /** The list entry of [engine], which retrains live on the data when [isLive]. */
  static entryOf(engine: Engine, isLive: boolean): MethodEntry {
    return {name: engine.name, isLive, isServer: isServerEngine(engine), probabilities: probabilitiesOf(engine)};
  }

  /** Candidates again: the entries that can take part, each checked by its default, in Custom by the user's toggle
   * first. */
  private rebuild(): void {
    this.isQuick = this.isQuickNow;
    const listed = this.entries.filter((e) => this.reasonOf(e) === null);
    this.isRebuilding = true;
    try {
      // Cleared first: the platform keeps a checked item that is no longer listed.
      this.candidates.value = [];
      this.candidates.items = listed.map((e) => e.name);
      const overrides = this.isCustom ? this.overrides : {};
      this.candidates.value = listed.filter((e) => overrides[e.name] ?? this.defaultCheck(e)).map((e) => e.name);
    } finally {
      this.isRebuilding = false;
    }
    // The new checkboxes do not know of the lock.
    if (this.isLocked) {
      for (const check of Array.from(this.candidates.root.querySelectorAll('input')))
        check.disabled = true;
    }
  }

  /** The tooltip of the Candidates label: what the list is, then why each unlisted method is not listed. */
  private candidatesTooltip(): HTMLElement {
    const reasons = this.entries.map((e) => this.reasonOf(e)).filter((r): r is string => r !== null);
    const unchecked = this.entries.filter((e) => this.reasonOf(e) === null)
      .map((e) => defaultUncheckReason(e.name, this.rowCount)).filter((r): r is string => r !== null);
    return ui.divV([METHODS_TOOLTIP, ...reasons.map((r) => `Not listed: ${r}`), ...unchecked]
      .map((text) => ui.divText(text)));
  }

  /** Why [entry] cannot take part: it runs on the server, or it gives no class probabilities while Predict
   * probability is on or the metric needs them; null when it can. */
  private reasonOf(entry: MethodEntry): string | null {
    if (entry.isServer)
      return `${entry.name}: ${SERVER_TOOLTIP}`;
    const metric = this.metricId;
    const exclusion = exclusionOf(entry.probabilities, this.predictProbability, metric);
    return exclusion === 'predictProbability' ? `${entry.name} gives no class probabilities` :
      exclusion === 'metric' ? `${entry.name} gives no class probabilities, which ${METRIC_LABELS[metric]} needs` :
        null;
  }

  /** Checked above the quick limit; at it, only a quick method ({@link isQuickMethod}); never one
   * {@link defaultUncheckReason} names. */
  private defaultCheck(entry: MethodEntry): boolean {
    return defaultUncheckReason(entry.name, this.rowCount) === null &&
      (!this.isQuick || isQuickMethod(entry.isLive, entry.probabilities, this.task === 'classification'));
  }

  /** Whether the limit used is quick now: a trial limit (Default's), or a time limit of at most
   * {@link QUICK_BUDGET_SECONDS}. */
  private get isQuickNow(): boolean {
    const limit = this.usedLimit;
    return isQuickLimit('budgetSeconds' in limit ? limit.budgetSeconds : undefined);
  }

  /** The limit changed: Candidates again only when it became quick or stopped being quick. */
  private onLimitChanged(): void {
    if (this.isQuickNow !== this.isQuick)
      this.rebuild();
    this.onSearchChanged();
  }

  /** A method was toggled: a toggle against its default is kept as the user's, one back to it is forgotten. */
  private onMethodsChanged(): void {
    if (this.isRebuilding)
      return;
    const checked = this.checkedMethods;
    for (const entry of this.entries) {
      if (this.reasonOf(entry) !== null)
        continue;
      const isChecked = checked.includes(entry.name);
      if (isChecked !== this.defaultCheck(entry))
        this.overrides[entry.name] = isChecked;
      else
        delete this.overrides[entry.name];
    }
    this.onSearchChanged();
  }

  /** Saved, the error line shown, the host told; nothing while the block is being built. */
  private onSearchChanged(): void {
    if (!this.isLoaded)
      return;
    this.updateErrors();
    this.save();
    this.host.onChanged();
  }

  /** Custom's limit validated; Default does not use it. */
  private updateErrors(): void {
    if (this.isCustom)
      this.limit.validate();
  }

  /** An equal split needs a trial per checked method. */
  private trialsProblem(trials: number): string | null {
    const checked = this.checkedMethods.length;
    return trials < checked ? `At least one trial per checked method: ${checked}` : null;
  }

  private get metricId(): MetricId {
    return this.isCustom ? TuneSettings.metricOf(this.metric, this.task) : defaultMetric(this.task ?? 'classification');
  }

  private save(): void {
    this.saver.save();
  }

  private record(): StoredFind {
    return {search: this.isCustom ? CUSTOM_SEARCH : DEFAULT_SEARCH,
      metric: TuneSettings.metricOf(this.metric, this.task), ...this.limit.stored(), methods: this.overrides};
  }

  /** The remembered settings, each item checked: a malformed or stale one falls back to its default (Search to
   * Default); the limit as {@link LimitInputs.load} reads it. */
  static load(): StoredFind {
    const value = TuneSettings.storedRecord(FIND_KEY);
    const {metric, methods} = value;
    const toggles: {[name: string]: boolean} = {};
    for (const [name, checked] of Object.entries(isRecord(methods) ? methods : {})) {
      if (typeof checked === 'boolean')
        toggles[name] = checked;
    }
    return {search: value.search === CUSTOM_SEARCH ? CUSTOM_SEARCH : DEFAULT_SEARCH,
      metric: typeof metric === 'string' ? metric : undefined, ...LimitInputs.load(value, DEFAULT_LIMIT),
      methods: toggles};
  }
}

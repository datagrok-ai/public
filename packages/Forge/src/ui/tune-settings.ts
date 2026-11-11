import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {defaultHyperparameters, defaultValuesOf, Engine, Hyperparameters, visibleHyperparametersOf}
  from '../engines/engine';
import {errorMessage} from '../forge-error';
import {ModelTask} from '../generated/db';
import {METRIC_LABELS, MetricId} from '../metrics/metrics';
import {isRecord} from '../preparation/preparation-options';
import {chosenOf, DEFAULT_BUDGET_SECONDS, defaultMetric, rangeProblem, SEARCH_METRICS, SearchChoice, SearchSpace,
  searchSpaceOf, SpaceDimension, TunedRange} from '../training/search';
import {LimitInputs, StoredLimit, TRIALS} from './limit-inputs';

export const MANUAL = 'Manual';
export const AUTO_TUNE = 'Auto-tune';
const DEFAULT_LIMIT: StoredLimit = {limit: TRIALS, minutes: DEFAULT_BUDGET_SECONDS / 60, trials: 50};
export const FIX_RANGE = 'Fix the highlighted range';
export const ALL_FIXED = 'All settings are fixed';
const CHOOSE_ONE = 'Choose at least one';
export type TuneStatus = 'Default' | 'Auto-tuned' | 'Customized';
export const CUSTOMIZED = 'Customized';
/** The status under Method while Auto is chosen. */
export const NOT_SEARCHED = 'Not searched';
export const STATUS_TOOLTIPS: Record<TuneStatus | typeof NOT_SEARCHED, string> = {
  'Default': 'The method\'s recommended settings', 'Auto-tuned': 'Found by the search', [CUSTOMIZED]: 'Changed by you',
  [NOT_SEARCHED]: 'Click Find to search'};
const PIN_TOOLTIP = 'Pin to exclude from tuning';
const SETTINGS_TOOLTIP = 'How the settings are chosen: by hand or by a search on this data';
const METRIC_TOOLTIP = 'Metric the search optimizes';
const FIXED_TOOLTIP = 'Fixed, can’t be tuned';
/** How long a settings save waits for the next change. */
const SAVE_DELAY_MS = 300;

type Control = HTMLInputElement | HTMLSelectElement | HTMLButtonElement;

/** What is remembered per user and method. */
interface StoredTune extends StoredLimit {
  mode: typeof MANUAL | typeof AUTO_TUNE;
  expanded: boolean;
  pinned: string[];
  ranges: {[name: string]: TunedRange};
  metric?: string;
}

/** A setting's row in **Search options**: its label, then on one line its value (a range, the choices, or the
 * setting's own value) and its pin, and the error line under them. */
export interface TuneRow {
  readonly property: DG.Property;
  /** The setting's input of Manual; [value] follows it both ways. */
  readonly own: DG.InputBase;
  readonly dim?: SpaceDimension;
  readonly label: HTMLElement;
  readonly line: HTMLElement;
  readonly from?: HTMLInputElement;
  readonly to?: HTMLInputElement;
  readonly choicesButton?: HTMLButtonElement;
  readonly choices?: HTMLInputElement[];
  readonly value: HTMLInputElement | HTMLSelectElement;
  readonly pin: HTMLButtonElement;
  readonly error: HTMLElement;
  isPinned: boolean;
}

/** A row with the parts a pin swaps (the range or choices, the own value with its box) and its problem as last
 * shown. */
interface Row extends TuneRow {
  readonly slot: HTMLElement;
  readonly tuned?: HTMLElement;
  readonly valueBox: HTMLElement;
  problem: string | null;
}

/** The data a tune ran on: the open table (compared by its `dart`), the target's and the features' names. */
export interface TuneContext {
  table: DG.DataFrame;
  target: string;
  features: string[];
}

/** What the settings ask of the Train view: the setting values, the form's data (undefined without a target), the
 * last tune of the method with its data ([isFound] when Find best model chose them), and the call after a change of
 * the mode or the search (pins, ranges, metric, budget); [affectsRules] when the mode or a pin changed, which decide
 * the settings the method ignores. */
interface TuneHost {
  values: () => Hyperparameters;
  context: () => TuneContext | undefined;
  lastTuned: () => {values: Hyperparameters; context: TuneContext; isFound: boolean} | undefined;
  onChanged: (affectsRules: boolean) => void;
}

let gridCount = 0;

/** The method's **Settings** (Manual / Auto-tune) with the reset button and the status under them; in Auto-tune the
 * **Search options** toggle and grid (Metric, Limit and Max time or Max trials, then each setting as a range, choices
 * or its own value with a pin) and the Tune bar. Manual shows the settings' own inputs, which are the view's. Both
 * layouts are built once; the mode only shows one of them. The Search options block is laid out by its own grid (the
 * user's decision: the platform form laid platform inputs out wrongly there), the limit's platform inputs in its cells
 * without their own labels; its label column follows the form's. Everything but the values is remembered per user and
 * method. */
export class TuneSettings {
  /** The `grok.userSettings` storage of the settings, one key per method. */
  static storage = 'Forge Tune';
  readonly root = ui.div([], 'forge-tune');
  readonly mode: DG.InputBase<string | null>;
  readonly resetButton: HTMLButtonElement;
  readonly statusText: HTMLSpanElement;
  readonly toggle: HTMLButtonElement;
  readonly grid = ui.div([], 'forge-tune-grid');
  readonly metric: HTMLSelectElement;
  /** Limit and Max time or Max trials. */
  readonly limit: LimitInputs;
  /** The subscriptions of the settings; the view drops them with the method. */
  readonly subs: rxjs.Subscription[] = [];
  /** The search space the method declares; undefined without one (Auto-tune is unavailable). */
  readonly space: SearchSpace | undefined;
  /** Why Auto-tune is unavailable: the method declares no search space, or a malformed one. */
  readonly noSpaceReason: string | undefined;
  private readonly engine: Engine;
  private readonly host: TuneHost;
  private readonly tuneBar: HTMLElement;
  private readonly autoBlock: HTMLElement;
  private readonly autoItem: HTMLElement | undefined;
  private readonly showExpanded: (isExpanded: boolean) => void;
  private readonly rows = new Map<string, Row>();
  /** The default values as {@link status} compares them. */
  private readonly defaultsText: string | undefined;
  private inactive: ReadonlySet<string> = new Set();
  private task: ModelTask | undefined;
  private storedMetric: string | undefined;
  private readonly saver: {save: () => void; flush: () => void};
  private isLocked = false;
  private isLoaded = false;
  private isResetting = false;

  constructor(engine: Engine, own: ReadonlyMap<string, DG.InputBase>, tuneButton: HTMLElement, host: TuneHost) {
    this.engine = engine;
    this.host = host;
    const properties = visibleHyperparametersOf(engine);
    let space: SearchSpace | undefined;
    try {
      space = searchSpaceOf(engine, properties);
      this.noSpaceReason = space === undefined ? `${engine.name} declares no search space.` : undefined;
    } catch (e) {
      this.noSpaceReason = errorMessage(e);
    }
    this.space = space;
    const stored = TuneSettings.load(engine, space, properties);
    this.storedMetric = stored.metric;
    this.saver = TuneSettings.saver(engine.name, () => this.record());
    this.subs.push(new rxjs.Subscription(() => this.flushSave()));

    this.resetButton = TuneSettings.iconButton('undo');
    this.resetButton.addEventListener('click', () => this.reset());
    ui.tooltip.bind(this.resetButton, () => this.resetButton.getAttribute('aria-label'));
    this.mode = ui.input.radio('Settings', {items: [MANUAL, AUTO_TUNE], value: stored.mode, nullable: false,
      onValueChanged: () => {
        this.showMode();
        this.save();
        this.host.onChanged(true);
      }});
    this.mode.root.classList.add('forge-inline-radio', 'forge-tune-mode');
    this.mode.addOptions(this.resetButton);
    const isUnavailable = (option: HTMLElement) => option.textContent === AUTO_TUNE && this.noSpaceReason !== undefined;
    const options = LimitInputs.bindOptionTooltips(this.mode,
      (option) => isUnavailable(option) ? this.noSpaceReason! : SETTINGS_TOOLTIP);
    this.autoItem = options.find((item) => item.textContent === AUTO_TUNE);
    if (this.autoItem !== undefined && isUnavailable(this.autoItem))
      this.autoItem.classList.add('forge-tune-unavailable');
    ui.tooltip.bind(this.mode.captionLabel, SETTINGS_TOOLTIP);
    const statusRow = TuneSettings.statusLine();
    this.statusText = statusRow.text;
    ui.tooltip.bind(this.statusText, () => this.statusTooltip());

    const metric = TuneSettings.metricRow(() => this.onSearchChanged());
    this.metric = metric.select;
    TuneSettings.setMetrics(this.metric, SEARCH_METRICS.classification, defaultMetric('classification'));
    this.limit = new LimitInputs(stored, DEFAULT_LIMIT, () => this.onSearchChanged());
    this.grid.append(metric.label, metric.cell, ...this.limit.elements);

    for (const property of properties) {
      const input = own.get(property.name);
      if (input !== undefined)
        this.rows.set(property.name, this.rowOf(property, input, stored));
    }
    this.defaultsText = this.textOf(defaultValuesOf(properties));

    ({toggle: this.toggle, show: this.showExpanded} = TuneSettings.optionsToggle(this.grid));
    this.toggle.addEventListener('click', () => this.setExpanded(!this.isExpanded));
    this.tuneBar = ui.div([tuneButton], 'forge-tune-bar');
    this.autoBlock = ui.div([this.toggle, this.grid], 'forge-tune-auto');
    // A method without settings shows no Settings row at all.
    if (this.rows.size > 0) {
      this.root.append(this.mode.root, statusRow.row, ...[...this.rows.values()].map((row) => row.own.root),
        this.autoBlock);
    }
    this.root.append(this.tuneBar);

    this.subs.push(TuneSettings.alignLabels(this.root, this.mode.captionLabel));
    this.setExpanded(stored.expanded);
    this.showMode();
    this.updateMode();
    this.updateStatus();
    this.isLoaded = true;
  }

  get isAuto(): boolean {
    return this.mode.value === AUTO_TUNE;
  }

  get isExpanded(): boolean {
    return this.toggle.getAttribute('aria-expanded') === 'true';
  }

  /** The platform inputs of the settings: the mode, Limit, Max time and Max trials. */
  get inputs(): DG.InputBase[] {
    return [this.mode, ...this.limit.inputs];
  }

  /** Every control of the settings but the platform inputs and the Search options toggle: the reset button, each
   * row's fields, choices and pin, and Metric. */
  get controls(): Control[] {
    return [this.resetButton, ...[...this.rows.values()].flatMap((row) => [...TuneSettings.fieldsOf(row), row.pin]),
      this.metric];
  }

  /** The settings the search chooses: those in the space and not pinned. */
  get tunedNames(): Set<string> {
    return new Set([...this.rows].filter(([, row]) => !row.isPinned).map(([name]) => name));
  }

  /** Auto-tuned when the values are the last tune's and the form's data is the tune's ({@link isSameContext}), unless
   * they are the method's defaults and a Tune (not Find best model) chose them; else Default for the defaults, else
   * Customized; the values compared as text over every setting. */
  get status(): TuneStatus {
    const tuned = this.tunedContext();
    const isTuned = tuned !== undefined && TuneSettings.isSameContext(tuned, this.host.context());
    return isTuned && (this.host.lastTuned()!.isFound || !this.isDefault) ? 'Auto-tuned' :
      this.isDefault ? 'Default' : CUSTOMIZED;
  }

  /** The values are the method's defaults. */
  private get isDefault(): boolean {
    return this.textOf(this.host.values()) === this.defaultsText;
  }

  /** The row of the setting [name] in Search options. */
  row(name: string): TuneRow | undefined {
    return this.rows.get(name);
  }

  /** Whether the method uses the setting [name] with the form's values. */
  isActive(name: string): boolean {
    return !this.inactive.has(name);
  }

  /** The settings the method ignores with the form's values: disabled, their ranges not checked. */
  setInactive(inactive: ReadonlySet<string>): void {
    if (inactive.size === this.inactive.size && [...inactive].every((name) => this.inactive.has(name)))
      return;
    this.inactive = inactive;
    this.updateControls();
    this.updateErrors();
  }

  /** Disables every control and the mode while the form is locked; the unlock gives each its own state back. */
  setLocked(isLocked: boolean): void {
    this.isLocked = isLocked;
    this.updateMode();
    this.updateControls();
  }

  /** Lists the metrics of [task], keeping the chosen one when it applies, else the remembered one, else the default. */
  setTask(task: ModelTask | undefined): void {
    if (task === undefined || task === this.task)
      return;
    this.task = TuneSettings.listMetrics(this.metric, this.task, task, this.storedMetric);
  }

  /** Switches the mode as the user would. */
  setMode(mode: typeof MANUAL | typeof AUTO_TUNE): void {
    if (this.mode.value !== mode)
      this.mode.value = mode;
  }

  /** Shows or hides the Search options grid. */
  setExpanded(isExpanded: boolean): void {
    this.showExpanded(isExpanded);
    this.save();
  }

  /** Re-reads the status: its text, and the reset button's look. */
  updateStatus(): void {
    this.statusText.textContent = this.status;
    this.resetButton.setAttribute('aria-disabled', String(this.isDefault));
    this.resetButton.setAttribute('aria-label', this.resetText());
  }

  /** Why Tune cannot search: a tuned setting's range or choices (as the error lines show them), nothing to tune, or the
   * limit; null when it can. */
  problem(): string | null {
    const tuned = [...this.rows].filter(([name, row]) => !row.isPinned && !this.inactive.has(name));
    if (tuned.some(([, row]) => row.problem !== null))
      return FIX_RANGE;
    return tuned.length === 0 ? ALL_FIXED : this.limit.problem();
  }

  /** What the search is asked: the metric, the limit, the pinned settings and the tuned ones' ranges. */
  choice(): SearchChoice {
    const pinned: string[] = [];
    const ranges: {[name: string]: TunedRange} = {};
    for (const [name, row] of this.rows) {
      if (row.dim === undefined)
        continue;
      if (row.isPinned)
        pinned.push(name);
      else
        ranges[name] = TuneSettings.rangeOf(row);
    }
    return {metric: this.metricId, ...this.limit.limit(), space: {pinned, ranges}};
  }

  /** Pins or unpins [name], as a click on its pin does: its row shows its own value, or its range or choices. */
  togglePin(name: string): void {
    const row = this.rows.get(name);
    if (row === undefined || row.dim === undefined || this.isLocked)
      return;
    // The pin's tooltip changes with the click, and a hovered control leaving the document gets no mouseleave.
    ui.tooltip.hide();
    row.isPinned = !row.isPinned;
    TuneSettings.showPin(row);
    this.onSearchChanged(true);
  }

  /** Why the row of [name] cannot be searched as shown: its range or choices when tuned, its own value's validity
   * when pinned or fixed; null for a setting the method ignores. */
  rowProblem(name: string): string | null {
    const row = this.rows.get(name);
    if (row === undefined || this.inactive.has(name))
      return null;
    if (row.isPinned)
      return row.own.validity;
    if (row.choices !== undefined)
      return row.choices.some((c) => c.checked) ? null : CHOOSE_ONE;
    return rangeProblem(row.dim!, row.property, TuneSettings.numberOf(row.from!.value),
      TuneSettings.numberOf(row.to!.value));
  }

  /** The reset button's tooltip and label, naming the method. */
  private resetText(): string {
    return this.isDefault ? `${this.engine.name} settings are already default` :
      `Reset ${this.engine.name} settings to defaults`;
  }

  private get metricId(): MetricId {
    return TuneSettings.metricOf(this.metric, this.task);
  }

  /** [values] as {@link status} compares them: every setting's value as text. */
  private textOf(values: Hyperparameters | undefined): string | undefined {
    return values === undefined ? undefined : JSON.stringify([...this.rows.keys()].map((name) => String(values[name])));
  }

  /** The data of the method's last tune when the values are that tune's; undefined otherwise. */
  private tunedContext(): TuneContext | undefined {
    const tuned = this.host.lastTuned();
    return tuned !== undefined && this.textOf(this.host.values()) === this.textOf(tuned.values) ? tuned.context :
      undefined;
  }

  /** The status's tooltip; tuned values on other data name the data they were tuned for. */
  private statusTooltip(): string {
    const status = this.status;
    const tuned = this.tunedContext();
    return status === CUSTOMIZED && tuned !== undefined ?
      `Tuned for ${tuned.target} on ${tuned.table.name}; tune again for this data` : STATUS_TOOLTIPS[status];
  }

  /** A pin, a range, the choices or the limit changed: saved, the errors shown, the view told. */
  private onSearchChanged(affectsRules = false): void {
    this.updateErrors();
    this.save();
    this.host.onChanged(affectsRules);
  }

  /** Manual shows the settings' own inputs; Auto-tune the Search options block and the Tune bar, which stays in the
   * document in Manual, so the tooltip overlay of a disabled Tune lives on. */
  private showMode(): void {
    const isAuto = this.isAuto;
    for (const row of this.rows.values())
      ui.setDisplay(row.own.root, !isAuto);
    ui.setDisplay(this.autoBlock, isAuto);
    ui.setDisplay(this.tuneBar, isAuto);
    this.updateErrors();
  }

  /** The mode follows the lock; Auto-tune stays unavailable for a method without a space, whatever enabled the
   * radio. */
  private updateMode(): void {
    this.mode.enabled = !this.isLocked;
    const radio = this.autoItem?.querySelector('input');
    if (radio !== null && radio !== undefined && this.noSpaceReason !== undefined)
      radio.disabled = true;
  }

  private rowOf(property: DG.Property, own: DG.InputBase, stored: StoredTune): Row {
    const name = property.name;
    const dim = this.space?.find((d) => d.name === name);
    const caption = own.caption;
    const describe = (element: HTMLElement) => property.description ? ui.tooltip.bind(element, property.description) :
      element;
    const label = describe(ui.div([caption], 'forge-tune-label'));
    const pin = TuneSettings.iconButton('thumbtack');
    pin.classList.add('forge-tune-pin');
    const error = ui.div([], 'forge-tune-err');
    error.setAttribute('role', 'alert');
    const {value, valueBox} = this.valueOf(own, property);
    describe(value);
    const slot = ui.div([], 'forge-tune-slot');
    const line = ui.div([slot, pin], 'forge-tune-row');
    this.grid.append(label, ui.div([line, error], 'forge-tune-cell'));
    const base = {property, own, label, line, value, pin, error, slot, valueBox, problem: null};
    if (dim === undefined) {
      slot.append(valueBox);
      pin.querySelector('i')?.classList.replace('fal', 'fas');
      pin.setAttribute('aria-disabled', 'true');
      pin.setAttribute('aria-pressed', 'true');
      pin.setAttribute('aria-label', `${caption}: ${FIXED_TOOLTIP}`);
      ui.tooltip.bind(pin, `${caption}: ${FIXED_TOOLTIP}`);
      return {...base, isPinned: true};
    }
    pin.addEventListener('click', () => this.togglePin(name));
    ui.tooltip.bind(pin, () => `${PIN_TOOLTIP}. ${pin.getAttribute('aria-label')}`);
    const range = stored.ranges[name];
    const isPinned = stored.pinned.includes(name);
    let row: Row;
    if (dim.type === 'categorical') {
      const chosen = (range?.choices ?? dim.choices!).map(String);
      const {box, button, checks} = this.choicesOf(caption, dim.choices!.map(String), chosen);
      for (const element of [button, ...checks.map((check) => check.parentElement!)])
        describe(element);
      row = {...base, dim, tuned: box, choicesButton: button, choices: checks, isPinned};
    } else {
      const from = TuneSettings.field(`${caption} from`, `${range?.low ?? dim.low}`);
      const to = TuneSettings.field(`${caption} to`, `${range?.high ?? dim.high}`);
      for (const field of [from, to]) {
        describe(field);
        field.addEventListener('input', () => this.onSearchChanged());
      }
      row = {...base, dim, tuned: ui.div([from, ui.span(['–'], 'forge-tune-hint'), to], 'forge-tune-range'), from, to,
        isPinned};
    }
    TuneSettings.showPin(row);
    return row;
  }

  /** The control of the setting's own value: a select for a choice, a checkbox for a yes/no, a field otherwise; it
   * sets [own] as typed and follows it, whoever changes it. */
  private valueOf(own: DG.InputBase, property: DG.Property):
    {value: HTMLInputElement | HTMLSelectElement; valueBox: HTMLElement} {
    const choices = property.choices ?? [];
    const isCheck = choices.length === 0 && property.propertyType === DG.TYPE.BOOL;
    let value: HTMLInputElement | HTMLSelectElement;
    let valueBox: HTMLElement;
    if (choices.length > 0) {
      const select: HTMLSelectElement = ui.element('select');
      select.append(...choices.map((c) => TuneSettings.option(c, c)));
      select.setAttribute('aria-label', own.caption);
      value = select;
      valueBox = TuneSettings.selectBox(select);
    } else if (isCheck) {
      const check: HTMLInputElement = ui.element('input');
      check.type = 'checkbox';
      check.setAttribute('aria-label', own.caption);
      value = check;
      valueBox = ui.div([check], 'forge-tune-check');
    } else {
      value = TuneSettings.field(own.caption, '');
      valueBox = value;
    }
    const show = () => {
      if (isCheck)
        (value as HTMLInputElement).checked = own.stringValue === 'true';
      else if (value.value !== (own.stringValue ?? ''))
        value.value = own.stringValue ?? '';
    };
    show();
    value.addEventListener(isCheck || choices.length > 0 ? 'change' : 'input', () => {
      own.stringValue = isCheck ? String((value as HTMLInputElement).checked) : value.value;
      // The input validates after it tells its listeners: the error line reads the validity once it is current.
      own.validate();
      this.updateErrors();
    });
    this.subs.push(own.onChanged.subscribe(() => {
      show();
      if (this.isResetting)
        return;
      this.updateErrors();
      this.updateStatus();
    }));
    return {value, valueBox};
  }

  /** The choices of a tuned categorical setting: a field listing the checked ones that opens a list of checkboxes. */
  private choicesOf(caption: string, choices: string[], chosen: string[]):
    {box: HTMLElement; button: HTMLButtonElement; checks: HTMLInputElement[]} {
    const text = ui.span([], 'forge-tune-choices-text');
    const button: HTMLButtonElement = ui.element('button', 'forge-tune-choices-button');
    button.type = 'button';
    button.classList.add('forge-tune-control');
    button.setAttribute('aria-label', `${caption} choices`);
    button.setAttribute('aria-haspopup', 'true');
    button.append(text);
    const checks = choices.map((choice) => {
      const check: HTMLInputElement = ui.element('input');
      check.type = 'checkbox';
      check.checked = chosen.includes(choice);
      return check;
    });
    const list = ui.div(checks.map((check, i) => {
      const item = ui.element('label', 'forge-tune-choice');
      item.append(check, ui.span([choices[i]]));
      return item;
    }), 'forge-tune-choice-list');
    const box = ui.div([button, list], 'forge-tune-choices');
    const showText = () => {
      const checked = choices.filter((_, i) => checks[i].checked);
      text.textContent = checked.length > 0 ? checked.join(', ') : 'None';
    };
    const onOutside = (e: MouseEvent) => {
      if (!box.contains(e.target as Node))
        setOpen(false);
    };
    const setOpen = (isOpen: boolean) => {
      ui.setDisplay(list, isOpen);
      button.setAttribute('aria-expanded', String(isOpen));
      if (isOpen)
        document.addEventListener('mousedown', onOutside, true);
      else
        document.removeEventListener('mousedown', onOutside, true);
    };
    button.addEventListener('click', () => setOpen(list.style.display === 'none'));
    box.addEventListener('keydown', (e) => {
      if (e.key === 'Escape')
        setOpen(false);
    });
    for (const check of checks) {
      check.addEventListener('change', () => {
        showText();
        this.onSearchChanged();
      });
    }
    showText();
    setOpen(false);
    return {box, button, checks};
  }

  /** Each control's state: all disabled while locked; a setting the method ignores keeps only its pin. */
  private updateControls(): void {
    for (const [name, row] of this.rows) {
      const isOff = this.isLocked || this.inactive.has(name);
      for (const control of TuneSettings.fieldsOf(row))
        control.disabled = isOff;
      row.pin.disabled = this.isLocked;
    }
    for (const control of [this.resetButton, this.metric])
      control.disabled = this.isLocked;
    this.limit.setLocked(this.isLocked);
  }

  /** In Auto-tune, the error line of each row and of the limit; Search options opens on one. */
  private updateErrors(): void {
    if (!this.isAuto)
      return;
    let hasError = false;
    for (const [name, row] of this.rows) {
      row.problem = this.rowProblem(name);
      hasError ||= row.problem !== null;
      row.error.textContent = row.problem ?? '';
      ui.setDisplay(row.error, row.problem !== null);
      for (const field of [row.from, row.to, row.choicesButton, row.value])
        field?.setAttribute('aria-invalid', String(row.problem !== null));
    }
    hasError = !this.limit.validate() || hasError;
    if (hasError && !this.isExpanded)
      this.setExpanded(true);
  }

  /** Sets the settings to the method's defaults at once; the mode, pins, ranges, metric and limit stay. */
  private reset(): void {
    if (this.isDefault || this.isLocked)
      return;
    const defaults = defaultHyperparameters(this.engine);
    this.isResetting = true;
    try {
      for (const [name, row] of this.rows) {
        const value = defaults[name];
        if (value !== undefined && String(row.own.value) !== String(value))
          row.own.value = value;
      }
    } finally {
      this.isResetting = false;
    }
    this.updateErrors();
    this.updateStatus();
  }

  private save(): void {
    if (this.isLoaded)
      this.saver.save();
  }

  /** Writes a pending save at once. */
  flushSave(): void {
    this.saver.flush();
  }

  /** What is remembered: everything but the values. */
  private record(): StoredTune {
    const pinned: string[] = [];
    const ranges: {[name: string]: TunedRange} = {};
    for (const [name, row] of this.rows) {
      if (row.dim === undefined)
        continue;
      ranges[name] = TuneSettings.rangeOf(row);
      if (row.isPinned)
        pinned.push(name);
    }
    return {mode: this.isAuto ? AUTO_TUNE : MANUAL, expanded: this.isExpanded, pinned, ranges,
      metric: this.metricId, ...this.limit.stored()};
  }

  /** A save of [record] under [key] of the storage, written {@link SAVE_DELAY_MS} after the last request, into the
   * storage of that request; [flush] writes a pending one at once. */
  static saver(key: string, record: () => object): {save: () => void; flush: () => void} {
    let storage: string | undefined;
    let timer: ReturnType<typeof setTimeout> | undefined;
    const flush = () => {
      clearTimeout(timer);
      timer = undefined;
      if (storage !== undefined)
        grok.userSettings.add(storage, key, JSON.stringify(record()));
      storage = undefined;
    };
    const save = () => {
      storage = TuneSettings.storage;
      clearTimeout(timer);
      timer = setTimeout(flush, SAVE_DELAY_MS);
    };
    return {save, flush};
  }

  /** The `grok.userSettings` record of [key] in the storage; empty when it is missing or malformed. */
  static storedRecord(key: string): {[key: string]: unknown} {
    try {
      const parsed: unknown = JSON.parse(grok.userSettings.getValue(TuneSettings.storage, key) ?? '');
      return isRecord(parsed) ? parsed : {};
    } catch {
      return {};
    }
  }

  /** The metric of [task] chosen in [select]; the task's default while none of its metrics is. */
  static metricOf(select: HTMLSelectElement, task: ModelTask | undefined): MetricId {
    const of = task ?? 'classification';
    return SEARCH_METRICS[of].find((id) => id === select.value) ?? defaultMetric(of);
  }

  /** Lists the metrics of [task] in [select] (after [previous]), keeping the chosen one when it applies, else
   * [stored], else the default; [task]. */
  static listMetrics(select: HTMLSelectElement, previous: ModelTask | undefined, task: ModelTask,
    stored: string | undefined): ModelTask {
    const chosen = previous === undefined ? undefined : TuneSettings.metricOf(select, previous);
    const ids = SEARCH_METRICS[task];
    TuneSettings.setMetrics(select, ids, ids.find((x) => x === chosen) ?? ids.find((x) => x === stored) ??
      defaultMetric(task));
    return task;
  }

  /** The italic status line under an input (Settings' and the classification threshold's): an empty label in the
   * form's label column, then the text, announced as it changes. */
  static statusLine(): {row: HTMLDivElement; text: HTMLSpanElement} {
    const text = ui.span([], 'forge-tune-status-text');
    text.setAttribute('aria-live', 'polite');
    return {row: ui.div([ui.label(' ', 'ui-input-label'), text], 'ui-input-root forge-tune-status'), text};
  }

  /** The **Search options** toggle of [grid] (Tune's and Find best model's): an empty chevron box as the groups'
   * headers have, so the chevron starts where the Method header's text does, the chevron and the caption; [show] shows
   * or hides the grid. The caller handles the click. */
  static optionsToggle(grid: HTMLElement): {toggle: HTMLButtonElement; show: (isExpanded: boolean) => void} {
    const expandIcon = ui.iconFA('chevron-right');
    const collapseIcon = ui.iconFA('chevron-down');
    grid.id = `forge-tune-grid-${gridCount++}`;
    const toggle: HTMLButtonElement = ui.element('button', 'forge-tune-toggle');
    toggle.type = 'button';
    toggle.append(ui.div([], 'forge-group-chevron'), expandIcon, collapseIcon, ui.span(['Search options']));
    toggle.setAttribute('aria-controls', grid.id);
    ui.tooltip.bind(toggle, () => toggle.getAttribute('aria-expanded') === 'true' ? 'Hide search options' :
      'Show search options');
    const show = (isExpanded: boolean) => {
      ui.setDisplay(grid, isExpanded);
      toggle.setAttribute('aria-expanded', String(isExpanded));
      ui.setDisplay(expandIcon, !isExpanded);
      ui.setDisplay(collapseIcon, isExpanded);
    };
    return {toggle, show};
  }

  /** The Metric row of a Search options grid: its label and its cell with the select, which calls [onChange]. */
  static metricRow(onChange: () => void): {label: HTMLElement; cell: HTMLElement; select: HTMLSelectElement} {
    const select: HTMLSelectElement = ui.element('select');
    select.setAttribute('aria-label', 'Metric');
    select.addEventListener('change', onChange);
    const label = ui.div(['Metric'], 'forge-tune-label');
    for (const element of [label, select])
      ui.tooltip.bind(element, METRIC_TOOLTIP);
    return {label, cell: TuneSettings.selectBox(select), select};
  }

  /** Lists the metrics [ids] in [select], [id] chosen. */
  static setMetrics(select: HTMLSelectElement, ids: MetricId[], id: MetricId): void {
    select.replaceChildren(...ids.map((x) => TuneSettings.option(x, METRIC_LABELS[x])));
    select.value = id;
  }

  /** The label column of [root]'s grid is as wide as [label] of the form (the platform sets the labels' widths). */
  static alignLabels(root: HTMLElement, label: HTMLElement): rxjs.Subscription {
    return ui.onSizeChanged(label).subscribe(() => {
      const width = label.getBoundingClientRect().width;
      if (width > 0)
        root.style.setProperty('--forge-tune-label', `${width}px`);
    });
  }

  /** [row]'s slot shows its own value when pinned, else its range or choices; its pin shows the state. */
  private static showPin(row: Row): void {
    row.slot.replaceChildren(row.isPinned ? row.valueBox : row.tuned!);
    row.pin.classList.toggle('forge-tune-pinned', row.isPinned);
    row.pin.querySelector('i')?.classList.replace(row.isPinned ? 'fal' : 'fas', row.isPinned ? 'fas' : 'fal');
    row.pin.setAttribute('aria-pressed', String(row.isPinned));
    row.pin.setAttribute('aria-label', `${row.own.caption}: ` +
      (row.isPinned ? 'pinned, excluded from tuning' : 'will be tuned'));
  }

  /** The row's fields that set its value: the range, the choices and the own value. */
  private static fieldsOf(row: TuneRow): Control[] {
    return [row.from, row.to, row.choicesButton, ...row.choices ?? [], row.value]
      .filter((c): c is Control => c !== undefined);
  }

  /** Whether [a] and [b] are the same data: the same open table, target and set of features. */
  private static isSameContext(a: TuneContext, b: TuneContext | undefined): boolean {
    return b !== undefined && a.table.dart === b.table.dart && a.target === b.target &&
      a.features.length === b.features.length && a.features.every((name) => b.features.includes(name));
  }

  /** The number typed in a field; null when it is empty. */
  private static numberOf(text: string): number | null {
    return text.trim() === '' ? null : Number(text);
  }

  private static rangeOf(row: TuneRow): TunedRange {
    if (row.choices !== undefined)
      return {choices: row.dim!.choices!.filter((_, i) => row.choices![i].checked)};
    return {low: TuneSettings.numberOf(row.from!.value) ?? undefined,
      high: TuneSettings.numberOf(row.to!.value) ?? undefined};
  }

  private static field(ariaLabel: string, value: string): HTMLInputElement {
    const field: HTMLInputElement = ui.element('input', 'forge-tune-field');
    field.type = 'text';
    field.classList.add('forge-tune-control');
    field.value = value;
    field.setAttribute('aria-label', ariaLabel);
    return field;
  }

  /** [select] as the platform's choice input draws it (`ui-input-editor`), in a box that fills its cell. */
  private static selectBox(select: HTMLSelectElement): HTMLElement {
    select.classList.add('ui-input-editor', 'forge-tune-control');
    return ui.div([select], 'forge-tune-select');
  }

  private static option(value: string, text: string): HTMLOptionElement {
    const option: HTMLOptionElement = ui.element('option');
    option.value = value;
    option.textContent = text;
    return option;
  }

  private static iconButton(icon: string): HTMLButtonElement {
    const button: HTMLButtonElement = ui.element('button', 'forge-tune-icon');
    button.type = 'button';
    button.append(ui.iconFA(icon));
    return button;
  }

  /** The remembered settings of [engine], each item checked against its space and settings ([properties]): a
   * malformed or stale item (a setting no longer in the space, a range outside its min and max, a mode the method
   * cannot have) falls back to its default; the limit as {@link LimitInputs.load} reads it. */
  static load(engine: Engine, space: SearchSpace | undefined, properties: DG.Property[]): StoredTune {
    const value = TuneSettings.storedRecord(engine.name);
    const names = (space ?? []).map((d) => d.name);
    const ranges: {[name: string]: TunedRange} = {};
    const storedRanges = isRecord(value.ranges) ? value.ranges : {};
    for (const dim of space ?? []) {
      const range = storedRanges[dim.name];
      const property = properties.find((p) => p.name === dim.name);
      if (!isRecord(range) || property === undefined)
        continue;
      if (dim.type === 'categorical') {
        const chosen = Array.isArray(range.choices) ? chosenOf(dim, range.choices) : [];
        if (chosen.length > 0)
          ranges[dim.name] = {choices: chosen};
      } else if (typeof range.low === 'number' && typeof range.high === 'number' &&
        rangeProblem(dim, property, range.low, range.high) === null)
        ranges[dim.name] = {low: range.low, high: range.high};
    }
    return {mode: value.mode === AUTO_TUNE && space !== undefined && properties.length > 0 ? AUTO_TUNE : MANUAL,
      expanded: value.expanded === true,
      pinned: Array.isArray(value.pinned) ? names.filter((name) => (value.pinned as unknown[]).includes(name)) : [],
      ranges, metric: typeof value.metric === 'string' ? value.metric : undefined,
      ...LimitInputs.load(value, DEFAULT_LIMIT)};
  }
}

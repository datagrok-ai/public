import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {visibleHyperparametersOf} from '../engines/engine';
import {fixed, METRIC_DESCRIPTIONS, METRIC_IDS, METRIC_LABELS, MetricId, metricOfLabel, Overfitting, overfittingOf,
  roundTo3} from '../metrics/metrics';
import {ValidationSettings} from '../training/k-fold';
import {defaultMetric, directionOf} from '../training/search';
import {CollapsibleGroup} from './collapsible-group';
import {metricIdsOf, valueColumns} from './metrics-table';
import {isSameInputs, SessionRun} from './runs-panel';
import {TuneSettings, TuneStatus} from './tune-settings';

const GROUP = 'Trained models';
const PIN_TOOLTIP = 'Pin to compare with other models';
const UNPIN_TOOLTIP = 'Remove from the comparison';
const COMPARE_TOOLTIP = 'The metric the cards show, sort and compare by.';
const SETTINGS_KINDS: Record<TuneStatus, string> = {'Default': 'Default settings',
  'Customized': 'Customized settings', 'Auto-tuned': 'Auto-tuned settings'};
const AUTO_TUNED_MODEL = 'Auto-tuned model';
const NO_OVERFITTING = 'No overfitting';
const OVERFITTING = 'Overfitting';
const SAME_ON_NEW_DATA = 'Should perform about as well on new data.';
const WORSE_ON_NEW_DATA = 'May perform worse on new data.';

/** A setting as a card lists it. */
export interface CardSetting {
  caption: string;
  value: string;
  /** The setting's Manual input tooltip; '' without one. */
  description: string;
}

interface Card {
  run: SessionRun;
  key: string;
  root: HTMLElement;
  /** The current card's "+", hidden while a pinned card has the run's inputs. */
  pin?: HTMLElement;
  best: HTMLElement;
  metric: HTMLElement;
  value: HTMLElement;
  status: HTMLElement;
  metricIds: MetricId[];
  /** The metric the card was last filled with. */
  filled?: MetricId;
}

/** The **Trained models** group: Compare by and the cards of the form's data, sorted by the chosen metric. */
export class ModelCards {
  readonly root = ui.divV([], 'forge-cards');
  readonly compareBy: DG.ChoiceInput<string | null>;
  readonly group: CollapsibleGroup;
  private currentCard: Card | undefined;
  private pinnedCards: Card[] = [];
  private key: string | undefined;
  private metric: MetricId | undefined;
  private isMetricChosen = false;

  constructor(private readonly onRestore: (run: SessionRun) => void) {
    this.compareBy = ui.input.choice<string | null>('Compare by', {items: [], value: null, nullable: false,
      tooltipText: COMPARE_TOOLTIP, onValueChanged: (label) => {
        this.metric = metricOfLabel(label);
        this.isMetricChosen = true;
        this.render();
      }});
    this.group = new CollapsibleGroup(GROUP, [this.compareBy.root, this.root]);
    ui.setDisplay(this.group.root, false);
  }

  get current(): SessionRun | undefined {
    return this.currentCard?.run;
  }

  /** The runs of the pinned cards, newest first. */
  get pinned(): SessionRun[] {
    return this.pinnedCards.map((c) => c.run);
  }

  /** The metric **Compare by** shows; undefined without a shown card. */
  get metricId(): MetricId | undefined {
    return this.cards().some((c) => c.key === this.key) ? this.metric : undefined;
  }

  /** The data a run or the form trains on: the target, the validation scheme with its folds or test fraction,
   * stratified, shuffle, and the seed when shuffled. */
  static keyOf(targetName: string, {scheme, folds, testFraction, stratified, shuffle}: ValidationSettings,
    seed: number): string {
    return JSON.stringify([targetName, scheme, scheme === 'holdout' ? testFraction : folds, stratified, shuffle,
      shuffle ? seed : null]);
  }

  /** Shows [run] as the current card, in place of the previous one. */
  setCurrent(run: SessionRun): void {
    const card = this.build(run, false);
    if (this.currentCard !== undefined)
      ui.remove(this.currentCard.root);
    this.root.prepend(card.root);
    this.currentCard = card;
    this.refreshPin();
    this.render();
  }

  /** Pins [run]; nothing while a pinned card has its inputs. */
  pin(run: SessionRun): void {
    if (this.pinnedCards.some((c) => isSameInputs(c.run, run)))
      return;
    this.pinnedCards.unshift(this.build(run, true));
    this.refreshPin();
    this.render();
  }

  /** Removes the pinned card of [run]. */
  unpin(run: SessionRun): void {
    const index = this.pinnedCards.findIndex((c) => c.run === run);
    if (index < 0)
      return;
    ui.remove(this.pinnedCards[index].root);
    this.pinnedCards.splice(index, 1);
    this.refreshPin();
    this.render();
  }

  /** Shows the cards trained on [key], the form's data ({@link keyOf}), and hides the others; undefined hides all. */
  update(key: string | undefined): void {
    if (key === this.key)
      return;
    this.key = key;
    this.render();
  }

  /** Removes every card: a Table change starts a new list. */
  clear(): void {
    for (const card of this.cards())
      ui.remove(card.root);
    this.currentCard = undefined;
    this.pinnedCards = [];
    this.render();
  }

  /** For each of [values] (undefined for a card without the value), whether it is at the best value as shown
   * ([direction] 1 the highest, -1 the lowest), ties included; none for fewer than two values. */
  static bestMarks(values: (number | undefined)[], direction: 1 | -1): boolean[] {
    const rounded = values.map((v) => v === undefined ? undefined : roundTo3(v));
    const defined = rounded.filter((v): v is number => v !== undefined);
    if (values.length < 2 || defined.length === 0)
      return values.map(() => false);
    const best = direction === 1 ? Math.max(...defined) : Math.min(...defined);
    return rounded.map((v) => v === best);
  }

  /** How [run]'s settings were chosen: `Default settings`, `Customized settings`, `Auto-tuned settings`, or
   * `Auto-tuned model` for the winner of a Find best model search over two methods or more. */
  static kindLine({settingsStatus, search}: Pick<SessionRun, 'settingsStatus' | 'search'>): string {
    return settingsStatus === 'Auto-tuned' && search?.mode === 'best' && (search.methods?.length ?? 0) > 1 ?
      AUTO_TUNED_MODEL : SETTINGS_KINDS[settingsStatus];
  }

  /** Every visible setting of [run] with its caption, value and description (its Manual input's tooltip); none for a
   * method without settings. */
  static settingsOf({engine, hyperparameters}: Pick<SessionRun, 'engine' | 'hyperparameters'>): CardSetting[] {
    return visibleHyperparametersOf(engine).map((p) => ({caption: ModelCards.captionOf(p),
      value: String(hyperparameters[p.name] ?? ''), description: p.description ?? ''}));
  }

  /** [kind] on its own line, then [settings] as a name-value table, two pairs a row; nothing without settings. */
  static settingsBlock(kind: string, settings: CardSetting[]): HTMLElement[] {
    if (settings.length === 0)
      return [];
    const described = (setting: CardSetting, text: string) => {
      const element = ui.span([text]);
      if (setting.description !== '')
        ui.tooltip.bind(element, setting.description);
      return element;
    };
    const rows = Array.from({length: Math.ceil(settings.length / 2)}, (_, i) => [settings[2 * i], settings[2 * i + 1]]);
    const cells = (setting: CardSetting | undefined) =>
      setting === undefined ? [] : [described(setting, setting.caption), described(setting, setting.value)];
    const table = ui.table(rows, ([first, second]) => [...cells(first), ...cells(second)]);
    return [ui.div([ui.span([kind], 'forge-card-kind')], 'forge-card-settings'),
      ui.div([table], 'forge-card-settings-table')];
  }

  /** Every metric of [run] with its Train and Validation values (never the Test of a search with a hold out), the row
   * of [metric] highlighted. */
  static metricsTooltip(run: SessionRun, metric: MetricId): HTMLElement {
    if (run.metrics === undefined)
      return ui.divText('No metrics');
    const columns = valueColumns({metrics: run.metrics}).filter(([column]) => column !== 'Test');
    const ids = metricIdsOf({metrics: run.metrics});
    const table = ui.table(ids, (id) => [METRIC_LABELS[id], ...columns.map(([, v]) => v[id]?.toFixed(3) ?? '')],
      ['Metric', ...columns.map(([column]) => column)]);
    table.classList.add('forge-card-tip-table');
    if (ids.includes(metric))
      table.rows[ids.indexOf(metric) + 1].classList.add('forge-card-tip-selected');
    return table;
  }

  /** The overfitting status of [metric] of [run]; null without both values. */
  static statusTooltip(run: SessionRun, metric: MetricId): HTMLElement | null {
    const overfitting = ModelCards.overfittingOfRun(run, metric);
    if (overfitting === undefined)
      return null;
    const {kind, gap, ratio, isOverfitting} = overfitting;
    const table = ui.table<[string, string]>([['Train', fixed(run.metrics?.train[metric])],
      ['Validation', fixed(run.metrics?.validation[metric])],
      kind === 'gap' ? ['Gap', fixed(gap)] : ['Ratio', ratio === undefined ? '-' : ratio.toFixed(2)]],
    (row) => row);
    table.classList.add('forge-card-tip-table', 'forge-card-status-tip-table');
    return ui.divV([ui.divText(`${METRIC_LABELS[metric]}:`), table,
      ui.divText(isOverfitting ? WORSE_ON_NEW_DATA : SAME_ON_NEW_DATA)], 'forge-tooltip');
  }

  private cards(): Card[] {
    return [...(this.currentCard === undefined ? [] : [this.currentCard]), ...this.pinnedCards];
  }

  private refreshPin(): void {
    const card = this.currentCard;
    if (card?.pin === undefined)
      return;
    // Hidden, not removed: the header keeps the icon's height.
    card.pin.style.visibility = this.pinnedCards.some((c) => isSameInputs(c.run, card.run)) ? 'hidden' : '';
  }

  /** Shows the cards of the form's data, Compare by with their metrics, each card's values and status, the order and
   * **Best**. */
  private render(): void {
    const visible = this.cards().filter((c) => c.key === this.key);
    for (const card of this.cards())
      ui.setDisplay(card.root, card.key === this.key);
    ui.setDisplay(this.group.root, visible.length > 0);
    if (visible.length === 0)
      return;
    const ids = METRIC_IDS.filter((id) => visible.some((c) => c.metricIds.includes(id)));
    if (!this.isMetricChosen || this.metric === undefined || !ids.includes(this.metric)) {
      this.isMetricChosen = false;
      this.metric = ModelCards.defaultMetricOf(visible[0].run, ids);
    }
    this.showCompareBy(ids);
    const metric = this.metric;
    if (metric === undefined)
      return;
    for (const card of visible)
      ModelCards.fill(card, metric);
    const direction = directionOf(metric);
    const rounded = new Map(this.cards().map((c): [Card, number | undefined] => {
      const v = c.run.metrics?.validation[metric];
      return [c, v === undefined ? undefined : roundTo3(v)];
    }));
    const rank = (c: Card) => {
      const v = rounded.get(c);
      return v === undefined ? -Number.MAX_VALUE : direction * v;
    };
    const ordered = [...(this.currentCard === undefined ? [] : [this.currentCard]),
      ...[...this.pinnedCards].sort((a, b) => rank(b) - rank(a))];
    if (ordered.some((c, i) => this.root.children[i] !== c.root))
      this.root.append(...ordered.map((c) => c.root));
    const shown = ordered.filter((c) => c.key === this.key);
    const marks = ModelCards.bestMarks(shown.map((c) => rounded.get(c)), direction);
    for (const [i, card] of shown.entries())
      ui.setDisplay(card.best, marks[i]);
  }

  /** Lists [ids] in **Compare by** with the chosen metric, without counting it as the user's pick. */
  private showCompareBy(ids: MetricId[]): void {
    const labels = ids.map((id) => METRIC_LABELS[id]);
    this.compareBy.notify = false;
    try {
      if (labels.length !== this.compareBy.items.length || labels.some((l, i) => l !== this.compareBy.items[i])) {
        this.compareBy.items = labels;
        TuneSettings.describeMetricOptions(this.compareBy.input);
      }
      const label = this.metric === undefined ? null : METRIC_LABELS[this.metric];
      // A set while the input fires its own change (the user's pick) throws on the platform.
      if (this.compareBy.value !== label)
        this.compareBy.value = label;
    } finally {
      this.compareBy.notify = true;
    }
  }

  /** The card's metric label, validation value and status for [metric]. */
  private static fill(card: Card, metric: MetricId): void {
    if (card.filled === metric)
      return;
    card.filled = metric;
    card.metric.textContent = METRIC_LABELS[metric];
    card.value.textContent = fixed(card.run.metrics?.validation[metric]);
    const overfitting = ModelCards.overfittingOfRun(card.run, metric);
    ui.setDisplay(card.status, overfitting !== undefined);
    if (overfitting === undefined)
      return;
    const {isOverfitting} = overfitting;
    card.status.classList.toggle('forge-card-ok', !isOverfitting);
    card.status.classList.toggle('forge-card-warn', isOverfitting);
    ui.empty(card.status);
    card.status.append(isOverfitting ? ui.span(['!'], 'forge-card-mark') : ui.iconFA('check'),
      ui.span([isOverfitting ? OVERFITTING : NO_OVERFITTING]));
  }

  private static overfittingOfRun({metrics}: SessionRun, metric: MetricId): Overfitting | undefined {
    const [train, validation] = [metrics?.train[metric], metrics?.validation[metric]];
    return train === undefined || validation === undefined ? undefined : overfittingOf(metric, train, validation);
  }

  /** A card of [run]; a click on a pinned card, but on its icon, restores its run. */
  private build(run: SessionRun, isPinned: boolean): Card {
    const root: HTMLElement = ui.element('div', 'forge-card');
    if (!isPinned)
      root.classList.add('forge-card-current');
    root.dataset.forgeRun = `${run.number}`;
    const icon = isPinned ? ui.icons.delete(() => this.unpin(run), UNPIN_TOOLTIP) :
      ui.icons.add(() => this.pin(run), PIN_TOOLTIP);
    icon.classList.add('forge-card-icon', isPinned ? 'forge-card-unpin' : 'forge-card-pin');
    const best = ui.span(['Best'], 'forge-card-best');
    ui.setDisplay(best, false);
    const metric = ui.span([], 'forge-card-metric');
    const value = ui.span([], 'forge-card-value');
    const score = ui.divH([metric, value], 'forge-card-score');
    const status = ui.divH([], 'forge-card-status');
    ui.tooltip.bind(metric, () => this.metric === undefined ? null : METRIC_DESCRIPTIONS[this.metric]);
    ui.tooltip.bind(value, () => this.metric === undefined ? null : ModelCards.metricsTooltip(run, this.metric));
    ui.tooltip.bind(status, () => this.metric === undefined ? null : ModelCards.statusTooltip(run, this.metric));
    const title = ui.divH([ui.span([run.engine.name], 'forge-card-title'), best], 'forge-card-title-row');
    const settings = ModelCards.settingsBlock(ModelCards.kindLine(run), ModelCards.settingsOf(run));
    root.append(ui.divH([ui.divV([title, ...settings], 'forge-card-body'), icon],
      'forge-card-header'), ui.divH([score, status], 'forge-card-main'));
    if (isPinned) {
      root.addEventListener('click', (e) => {
        if (e.target instanceof Element && !icon.contains(e.target))
          this.onRestore(run);
      });
    }
    return {run, key: ModelCards.keyOf(run.targetName, run.validation, run.seed), root,
      pin: isPinned ? undefined : icon, best, metric, value, status,
      metricIds: run.metrics === undefined ? [] : metricIdsOf({metrics: run.metrics})};
  }

  /** The metric the cards compare by until the user picks one: the one [run]'s search optimized, else Balanced
   * accuracy, else R2, else the first of [ids]. */
  private static defaultMetricOf(run: SessionRun, ids: MetricId[]): MetricId | undefined {
    return [run.search?.metric, defaultMetric('classification'), defaultMetric('regression')]
      .find((id): id is MetricId => id !== undefined && ids.includes(id)) ?? ids[0];
  }

  private static captionOf(p: DG.Property): string {
    return DG.StringUtils.camelCaseToSentence(p.caption || p.name);
  }
}

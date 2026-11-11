import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {visibleHyperparametersOf} from '../engines/engine';
import {METRIC_DESCRIPTIONS, METRIC_LABELS, metricFindings} from '../metrics/metrics';
import {ValidationSettings} from '../training/k-fold';
import {defaultMetric, SearchMode} from '../training/search';
import {findingsOf, METRIC_HEADERS, metricIdsOf, metricTooltip, valueColumns} from './metrics-table';
import {isSameInputs, SessionRun} from './runs-panel';

const PERFORMANCE = 'Performance';
const SETTINGS = 'Settings';
const PIN_TOOLTIP = 'Pin to compare with other models';
const UNPIN_TOOLTIP = 'Remove from the comparison';
const NO_SETTINGS = 'No settings';
/** The card's caption after a search: what it chose. */
const CARD_CAPTIONS: Record<SearchMode, string> = {best: 'Auto-tuned model', tune: 'Auto-tuned settings'};
/** The metric rows a card shows until Show more. */
const SHOWN_METRICS = 3;

interface Card {
  run: SessionRun;
  key: string;
  root: HTMLElement;
  tabs: DG.TabControl;
  sub: rxjs.Subscription;
  /** The current card's "+", hidden while a pinned card has the run's inputs. */
  pin?: HTMLElement;
}

/** The cards under the form: the current model (the newest completed training of the table) first, then the pinned
 * ones, newest first; each with its Performance and Settings tabs, switched together, and Performance's own Show more /
 * Show less. A card is shown only while the form's target, validation and split seed are its own ({@link keyOf}).
 * Hidden without a card. */
export class ModelCards {
  readonly root = ui.divV([], 'forge-cards');
  private currentCard: Card | undefined;
  private pinnedCards: Card[] = [];
  /** The tab every card shows; kept across a Table change. */
  private tab = PERFORMANCE;
  private key: string | undefined;

  constructor(private readonly onRestore: (run: SessionRun) => void) {
    ui.setDisplay(this.root, false);
  }

  get current(): SessionRun | undefined {
    return this.currentCard?.run;
  }

  /** The runs of the pinned cards, newest first. */
  get pinned(): SessionRun[] {
    return this.pinnedCards.map((c) => c.run);
  }

  /** The tabs of the card [root]; undefined for another element. */
  tabsOf(root: HTMLElement): DG.TabControl | undefined {
    return this.cards().find((c) => c.root === root)?.tabs;
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
      ModelCards.dispose(this.currentCard);
    this.root.prepend(card.root);
    this.currentCard = card;
    this.refreshPin();
    this.show();
  }

  /** Pins [run] right under the current card; nothing while a pinned card has its inputs. */
  pin(run: SessionRun): void {
    if (this.pinnedCards.some((c) => isSameInputs(c.run, run)))
      return;
    const card = this.build(run, true);
    this.pinnedCards.unshift(card);
    this.currentCard!.root.after(card.root);
    this.refreshPin();
    this.show();
  }

  /** Removes the pinned card of [run]. */
  unpin(run: SessionRun): void {
    const index = this.pinnedCards.findIndex((c) => c.run === run);
    if (index < 0)
      return;
    ModelCards.dispose(this.pinnedCards[index]);
    this.pinnedCards.splice(index, 1);
    this.refreshPin();
    this.show();
  }

  /** Shows the cards trained on [key], the form's data ({@link keyOf}), and hides the others; undefined hides all. */
  update(key: string | undefined): void {
    this.key = key;
    this.show();
  }

  /** Removes every card: a Table change starts a new list. */
  clear(): void {
    for (const card of this.cards())
      ModelCards.dispose(card);
    this.currentCard = undefined;
    this.pinnedCards = [];
    this.show();
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

  private show(): void {
    const cards = this.cards();
    for (const card of cards)
      ui.setDisplay(card.root, card.key === this.key);
    ui.setDisplay(this.root, cards.length > 0);
  }

  /** A card of [run]: the method and the icon ("+" on the current card, the trash on a pinned one), after a search what
   * it chose ({@link CARD_CAPTIONS}), then the tabs. A click on a pinned card, but on its tab headers or its icon,
   * restores its run. */
  private build(run: SessionRun, isPinned: boolean): Card {
    // Not a ui.div: ui.css gives a box in a ui-div 400 x 300 px, and the tabs are as tall as the shown pane.
    const root: HTMLElement = ui.element('div', 'forge-card');
    if (!isPinned)
      root.classList.add('forge-card-current');
    root.dataset.forgeRun = `${run.number}`;
    const icon = isPinned ? ui.icons.delete(() => this.unpin(run), UNPIN_TOOLTIP) :
      ui.icons.add(() => this.pin(run), PIN_TOOLTIP);
    icon.classList.add('forge-card-icon', isPinned ? 'forge-card-unpin' : 'forge-card-pin');
    const header = ui.divH([ui.divText(run.engine.name, 'forge-card-title'), icon], 'forge-card-header');
    const caption = run.search === undefined ? [] :
      [ui.divText(CARD_CAPTIONS[run.search.mode], 'forge-card-caption')];
    const performance = ModelCards.performance(run);
    const tabs = ui.tabControl({[PERFORMANCE]: performance, [SETTINGS]: () => ModelCards.settingsTable(run)});
    tabs.currentPane = tabs.getPane(this.tab);
    const sub = tabs.onTabChanged.subscribe(() => this.syncTabs(tabs));
    root.append(header, ...caption, tabs.root);
    if (isPinned) {
      root.addEventListener('click', (e) => {
        if (e.target instanceof Element && !tabs.header.contains(e.target) && !icon.contains(e.target) &&
          e.target.closest('.forge-card-more') === null)
          this.onRestore(run);
      });
    }
    return {run, key: ModelCards.keyOf(run.targetName, run.validation, run.seed), root, tabs, sub,
      pin: isPinned ? undefined : icon};
  }

  /** The Performance pane of [run]: its metrics table, the rows after the first {@link SHOWN_METRICS} hidden, and with
   * more metrics than that Show more, which shows them and turns into Show less, and back. */
  private static performance(run: SessionRun): HTMLElement {
    const table = ModelCards.performanceTable(run);
    if (!(table instanceof HTMLTableElement) || table.rows.length <= SHOWN_METRICS + 1)
      return ui.divV([table]);
    const extra = Array.from(table.rows).slice(SHOWN_METRICS + 1);
    let showsAll = false;
    const more = ui.link('', () => {
      ui.tooltip.hide();
      showsAll = !showsAll;
      show();
    }, undefined, 'forge-card-more');
    ui.tooltip.bind(more, () => showsAll ? 'Show fewer metrics' : 'Show all metrics');
    const show = () => {
      for (const row of extra)
        ui.setDisplay(row, showsAll);
      more.textContent = showsAll ? 'Show less' : 'Show more';
    };
    show();
    return ui.divV([table, more]);
  }

  /** A click on a tab of [source] switches every other card to it. */
  private syncTabs(source: DG.TabControl): void {
    this.tab = source.currentPane.name;
    for (const card of this.cards()) {
      if (card.tabs.currentPane.name !== this.tab)
        card.tabs.currentPane = card.tabs.getPane(this.tab);
    }
  }

  private static dispose(card: Card): void {
    card.sub.unsubscribe();
    ui.remove(card.root);
  }

  /** The run's metrics as the Performance grid shows them: the same value columns, three decimals, uncolored; hover a
   * header, a metric or a value as there (a value's tooltip names its findings). The first row is the metric the run's
   * search optimized, else the task's default one; then the others in the grid's order. */
  private static performanceTable({metrics, search}: SessionRun): HTMLElement {
    if (metrics === undefined)
      return ui.divText('No metrics', 'forge-note');
    const columns = valueColumns({metrics});
    const findings = metricFindings(metrics);
    const ids = metricIdsOf({metrics});
    const first = [search?.metric, defaultMetric('classification'), defaultMetric('regression')]
      .find((id) => id !== undefined && ids.includes(id)) ?? ids[0];
    const ordered = [first, ...ids.filter((id) => id !== first)];
    const table = ui.table(ordered, (id) => [
      ModelCards.described(ui.span([METRIC_LABELS[id]]), METRIC_DESCRIPTIONS[id]),
      ...columns.map(([column, values]) => {
        const value = values[id];
        return ModelCards.described(ui.span([value === undefined ? '' : value.toFixed(3)]),
          metricTooltip(value, column, findingsOf(findings, id, column)));
      }),
    ], ['Metric', ...columns.map(([column]) => column)]);
    table.classList.add('forge-card-table');
    for (const cell of Array.from(table.rows[0].cells))
      ModelCards.described(cell, METRIC_HEADERS[cell.textContent ?? ''] ?? '');
    return table;
  }

  /** Every setting of the method as trained, captioned as its Manual input (hover a row for its meaning), in a table
   * styled as Performance's. */
  private static settingsTable(run: SessionRun): HTMLElement {
    const properties = visibleHyperparametersOf(run.engine);
    if (properties.length === 0)
      return ui.divV([ui.divText(NO_SETTINGS, 'forge-note')], 'forge-card-settings');
    const table = ui.table(properties, (p) => [ui.span([DG.StringUtils.camelCaseToSentence(p.caption || p.name)]),
      ui.span([String(run.hyperparameters[p.name] ?? '')])], ['Setting', 'Value']);
    table.classList.add('forge-card-table');
    for (const [i, p] of properties.entries())
      ModelCards.described(table.rows[i + 1], p.description ?? '');
    return ui.divV([table], 'forge-card-settings');
  }

  private static described<T extends HTMLElement>(element: T, text: string): T {
    if (text !== '')
      ui.tooltip.bind(element, text);
    return element;
  }
}

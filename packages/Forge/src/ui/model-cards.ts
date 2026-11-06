import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {hyperparametersOf} from '../engines/engine';
import {METRIC_DESCRIPTIONS, METRIC_LABELS, metricFindings} from '../metrics/metrics';
import {ValidationSettings} from '../training/k-fold';
import {findingsOf, holdOutLine, METRIC_HEADERS, metricIdsOf, metricTooltip, valueColumns} from './metrics-table';
import {ForgeModelHandler} from './model-handler';
import {isSameInputs, RunsPanel, SessionRun} from './runs-panel';
import {STATUS_TOOLTIPS} from './tune-settings';

const PERFORMANCE = 'Performance';
const SETTINGS = 'Settings';
const PIN_TOOLTIP = 'Pin to compare with other models';
const UNPIN_TOOLTIP = 'Remove from the comparison';
const NO_SETTINGS = 'No settings';

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
 * ones, newest first; each with its Performance and Settings tabs, switched together. A card is shown only
 * while the form's target, validation and split seed are its own ({@link keyOf}). Hidden without a card. */
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

  /** A card of [run]: the method and the icon ("+" on the current card, the trash on a pinned one), the line of the
   * features and the preparation, then the tabs. A click on a pinned card, but on its tab headers or its icon,
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
    const line = ui.divText(ModelCards.lineText(run), 'forge-card-line');
    ui.tooltip.bind(line, () => ModelCards.lineTooltip(run));
    const tabs = ui.tabControl({[PERFORMANCE]: () => ModelCards.performanceTable(run),
      [SETTINGS]: () => ModelCards.settingsTable(run)});
    tabs.currentPane = tabs.getPane(this.tab);
    const sub = tabs.onTabChanged.subscribe(() => this.syncTabs(tabs));
    root.append(header, line, tabs.root);
    if (isPinned) {
      root.addEventListener('click', (e) => {
        if (e.target instanceof Node && !tabs.header.contains(e.target) && !icon.contains(e.target))
          this.onRestore(run);
      });
    }
    return {run, key: ModelCards.keyOf(run.targetName, run.validation, run.seed), root, tabs, sub,
      pin: isPinned ? undefined : icon};
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

  private static lineText(run: SessionRun): string {
    return [ForgeModelHandler.featureCount(run.featureNames.length), ...RunsPanel.preparationItems(run)].join(' · ');
  }

  /** The feature names, the preparation and, after a search with a hold out, the rows it held out. */
  private static lineTooltip(run: SessionRun): HTMLElement {
    const holdOut = run.search?.holdOut;
    const lines = [`Features: ${run.featureNames.join(', ')}`,
      `Preparation: ${RunsPanel.preparationItems(run).join(', ') || 'none'}`,
      ...(holdOut === undefined ? [] : [holdOutLine(holdOut.rows)])];
    return ui.divV(lines.map((text) => ui.divText(text)), 'forge-tooltip');
  }

  /** The run's metrics as the Performance grid shows them: the same rows and value columns, three decimals, uncolored;
   * hover a header, a metric or a value as there (a value's tooltip names its findings). */
  private static performanceTable({metrics}: SessionRun): HTMLElement {
    if (metrics === undefined)
      return ui.divText('No metrics', 'forge-note');
    const columns = valueColumns({metrics});
    const findings = metricFindings(metrics);
    const table = ui.table(metricIdsOf({metrics}), (id) => [
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

  /** The settings status recorded with the run (hover it for its meaning and the search that found the settings), then
   * every setting of the method as trained, captioned as its Manual input (hover a row for its meaning), in a table
   * styled as Performance's. */
  private static settingsTable(run: SessionRun): HTMLElement {
    const properties = hyperparametersOf(run.engine);
    const status = ui.span([run.settingsStatus], 'forge-tune-status-text');
    ui.tooltip.bind(status, () =>
      [STATUS_TOOLTIPS[run.settingsStatus], ...RunsPanel.searchLines(run.search)].join('. '));
    if (properties.length === 0)
      return ui.divV([status, ui.divText(NO_SETTINGS, 'forge-note')], 'forge-card-settings');
    const table = ui.table(properties, (p) => [ui.span([DG.StringUtils.camelCaseToSentence(p.caption || p.name)]),
      ui.span([String(run.hyperparameters[p.name] ?? '')])], ['Setting', 'Value']);
    table.classList.add('forge-card-table');
    for (const [i, p] of properties.entries())
      ModelCards.described(table.rows[i + 1], p.description ?? '');
    return ui.divV([status, table], 'forge-card-settings');
  }

  private static described<T extends HTMLElement>(element: T, text: string): T {
    if (text !== '')
      ui.tooltip.bind(element, text);
    return element;
  }
}

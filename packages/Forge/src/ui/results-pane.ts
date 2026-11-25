import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {onFrame} from '../preparation/shared-frame';
import {evaluationInsights, Insight, summaryInsight} from '../training/insights';
import {MetricsTable} from './metrics-table';
import {reportError} from './report-error';
import {PanelKey, resultKindOf, ResultsRow, resultsLayout} from './results-layout';
import {buildShown, markCurves, Mismatches, PANEL_TEXTS, refillShown, releaseLinkedFrame, releaseShown, ResultsData,
  resultsPanel, resultsTitle, Shown, updateWrongPredictions} from './results-viewers';
import {bulletList} from './rich-tooltip';
import '../../css/forge.css';

/** A slot's preset width per unit of its weight: the splitter derives the slots' growth from the widths it finds once
 * the row is in the DOM, so it keeps the ratio of the weights. */
const SLOT_WIDTH = 100;

/** A mounted row: [root] is the `ui.splitH` element, [slots] its panel boxes with their [weights], [height] its inline
 * height the last time it was shown. */
interface RowRecord {
  signature: string;
  root: HTMLElement;
  slots: HTMLElement[];
  weights: number[];
  height: string;
}

/** A block under the rows: [block] in its box [root], which keeps the height the user dragged. */
interface BottomRecord {
  block: HTMLElement;
  root: HTMLElement;
  height: string;
}

/** The results of the Train view and the evaluation. Before a result: the **Performance** panel (its header, a hint)
 * and **Insights & Tips** under it. With a result: rows of panels in the platform's splitters (`results-layout.ts`),
 * the two panels among them (**Insights & Tips** never hidden); the metrics grid kept in place while the metric rows
 * stay; the viewers of `results-viewers.ts` over one frame of the user's columns and Forge's own (no copies), linked
 * to the table by the hovered row. A result of the same data keeps the viewers and refills the model's columns; a row
 * whose panels stay the same keeps the sizes the user gave it. The blocks given for the bottom (the Train view's Search
 * and Runs) follow the rows in the same vertical splitter while they are shown, with or without a result. */
export class ResultsPane {
  readonly root: HTMLElement;
  private readonly insightsList = ui.element('ul');
  private readonly insightsPanel: HTMLElement;
  private readonly metricsHost = ui.div([], 'forge-results-body');
  private readonly metricsPanel: HTMLElement;
  /** The two panels while there is no result. */
  private readonly hintHost: HTMLElement;
  private readonly rowsHost = ui.divV([], 'forge-results-rows');
  private rows: RowRecord[] = [];
  private readonly bottom: BottomRecord[];
  /** The boxes in the vertical splitter, in order. */
  private sequence: HTMLElement[] = [];
  /** The panels of the mounted rows. */
  private mounted = new Set<PanelKey>();
  /** The splitters measure their boxes whenever their size changes, and a hidden view measures 0 px: the rows get back
   * the heights they had while shown, the slots their growth. A workaround of the platform: `splitV`'s `handleResize`
   * (`js-api/ui.ts:1947-1959`) lacks the `document.contains` guard of `splitH`'s, so a hidden view collapses every
   * row to 0 px; remove it once js-api has the guard. */
  private readonly sizeKeeper = new MutationObserver(() => this.keepSizes());
  private metricsTable: MetricsTable | undefined;
  private dataInsights: Insight[] = [];
  /** The insights of the shown result and its metric findings. */
  private resultLines: Insight[] = [];
  /** The line shown alone when there is no other: the result's {@link summaryInsight}. */
  private fallbackLine: Insight | undefined;
  private shown: Shown | undefined;
  /** Counts the shows and hints: a show that a newer call outdated while `visualize` ran mounts nothing. */
  private version = 0;

  /** [status] is shown next to the **Performance** title (the Train view's loader); [bottom] are the blocks under the
   * rows with their initial heights in px, mounted while shown ({@link layout}). */
  constructor(status?: HTMLElement, bottom: {block: HTMLElement; height: number}[] = []) {
    this.bottom = bottom.map(({block, height}) => {
      const root = ui.box(block);
      root.style.height = `${height}px`;
      return {block, root, height: root.style.height};
    });
    const header = resultsTitle('metrics');
    header.classList.add('forge-results-header');
    if (status !== undefined)
      header.append(status);
    this.metricsPanel = ui.divV([header, this.metricsHost], 'forge-results-panel forge-results-metrics');
    this.insightsPanel = ui.divV([resultsTitle('insights'), this.insightsList], 'forge-results-panel forge-insights');
    this.hintHost = ui.divV([this.metricsPanel, this.insightsPanel]);
    this.root = ui.divV([this.hintHost, this.rowsHost], 'forge-results');
    this.fillInsights();
  }

  /** Mounts the rows and the shown bottom blocks in one vertical splitter, rebuilt only when the sequence of boxes
   * changes; a hidden block leaves it. */
  layout(): void {
    const shown = this.bottom.filter((b) => b.block.style.display !== 'none');
    const roots = [...this.rows.map((r) => r.root), ...shown.map((b) => b.root)];
    if (roots.length === this.sequence.length && roots.every((root, i) => root === this.sequence[i]))
      return;
    this.observe();
    this.sequence = roots;
    // Every box out of the old splitter before it is disposed: a hidden block is kept for later.
    for (const b of this.bottom)
      b.root.remove();
    // A lone box in a vertical splitter would get half the height: a single box is mounted as is.
    const host = roots.length > 1 ? ui.splitV(roots, null, true) : roots[0];
    // A lone row may sit in the old splitter, which emptying disposes.
    host?.remove();
    ui.empty(this.rowsHost);
    if (host !== undefined)
      this.rowsHost.append(host);
  }

  /** The size keeper watches the rows, their slots and the bottom boxes. */
  private observe(): void {
    this.sizeKeeper.disconnect();
    for (const element of [...this.rows.flatMap((row) => [row.root, ...row.slots]), ...this.bottom.map((b) => b.root)])
      this.sizeKeeper.observe(element, {attributes: true, attributeFilter: ['style']});
  }

  /** The shared frames of the shown result: the viewers' and the wrong predictions'. */
  get frames(): readonly DG.DataFrame[] {
    const shown = this.shown;
    return shown === undefined ? [] : [shown.frame, ...(shown.wrong === undefined ? [] : [shown.wrong.frame])];
  }

  /** [insights] of the data, shown above the result ones until the next call. */
  setDataInsights(insights: Insight[]): void {
    this.dataInsights = insights;
    this.fillInsights();
  }

  /** [text] in the **Performance** panel instead of the result ('' for nothing); the data insights stay under it. The
   * metrics grid is kept for the next result. */
  showHint(text: string): void {
    this.version++;
    this.unmount();
    this.resultLines = [];
    this.fallbackLine = undefined;
    this.fillInsights();
    this.metricsTable?.root.remove();
    ui.empty(this.metricsHost);
    if (text !== '')
      this.metricsHost.append(ui.divText(text));
  }

  /** Shows [data]: the metrics grid at once (the same one while the metric rows stay) and the result insights; the
   * rows of panels once the method's own views are ready, mounted together. */
  async show(data: ResultsData): Promise<void> {
    const version = ++this.version;
    const {summary, target, prediction, positiveClass} = data;
    if (this.metricsTable === undefined || !this.metricsTable.update(summary)) {
      this.dropMetricsTable();
      this.metricsTable = new MetricsTable(summary);
    }
    if (this.metricsTable.root.parentElement !== this.metricsHost) {
      ui.empty(this.metricsHost);
      this.metricsHost.append(this.metricsTable.root);
    }
    this.resultLines = evaluationInsights({values: summary.evaluation?.values ?? summary.metrics.train, target,
      prediction, positiveClass}, summary.metrics);
    this.fallbackLine = summaryInsight(summary.metrics, summary.evaluation?.values);
    this.fillInsights();
    const elements = await this.visualizeElements(data);
    if (version !== this.version) {
      ui.empty(ui.div(elements));
      return;
    }
    const kind = resultKindOf(target);
    const old = this.shown;
    const isRefilled = old !== undefined && refillShown(old, data, kind);
    const shown = isRefilled ? old : buildShown(data, kind);
    let dropped: Mismatches | undefined;
    try {
      markCurves(shown, data);
      dropped = updateWrongPredictions(shown);
    } catch (e) {
      if (!isRefilled)
        releaseShown(shown);
      throw e;
    }
    const panels = new Map<PanelKey, HTMLElement>([...shown.panels, ['metrics', this.metricsPanel],
      ['insights', this.insightsPanel]]);
    for (const [i, element] of elements.entries())
      panels.set(`visualize:${i}`, resultsPanel('visualize', element, data.engine.name));
    this.mountRows(resultsLayout(kind, new Set(panels.keys()), elements.length), panels);
    if (dropped !== undefined)
      releaseLinkedFrame(dropped);
    if (!isRefilled && old !== undefined)
      releaseShown(old);
    this.shown = shown;
    this.fillInsights();
  }

  /** Everything out: the rows, their frames and subscriptions, the grid, the insights. */
  clear(): void {
    this.dataInsights = [];
    this.showHint('');
    this.dropMetricsTable();
  }

  /** The insights, each with a tooltip of its numbers and, last and in grey, a pointer to the mounted panel it is
   * about; a viewer panel's insights also become the viewer's description. Without any, the one line reading a metric
   * ({@link fallbackLine}); the panel is hidden only while empty under the hint, never with a result. */
  private fillInsights(): void {
    ui.empty(this.insightsList);
    const texts = new Set<string>();
    const lines = [...this.dataInsights, ...this.resultLines].filter(({text}) => !texts.has(text) && texts.add(text));
    const shown = lines.length > 0 || this.fallbackLine === undefined ? lines : [this.fallbackLine];
    for (const {text, detail, items, panel} of shown) {
      const li = ui.element('li');
      li.append(ui.divText(text));
      const pointer = panel !== undefined && this.mounted.has(panel) ? [`See ${PANEL_TEXTS[panel][0]}.`] : [];
      if (detail !== undefined || items !== undefined || pointer.length > 0) {
        ui.tooltip.bind(li, () => ui.divV([...(detail ?? '').split('\n').filter((line) => line !== '')
          .map((line) => ui.divText(line)), ...(items === undefined ? [] : [bulletList(items.map((item) => [item]))]),
        ...pointer.map((line) => ui.divText(line, 'forge-tooltip-minor'))]));
      }
      this.insightsList.append(li);
    }
    ui.setDisplay(this.insightsPanel, shown.length > 0 || this.mounted.has('insights'));
    for (const [key, viewer] of this.shown?.viewers ?? []) {
      const description = shown.filter(({panel}) => panel === key && panel !== 'metrics')
        .map(({text}) => ResultsPane.markdownText(text)).join(' ');
      if ((viewer.props.description ?? '') !== description)
        viewer.setOptions({description, descriptionVisibilityMode: 'Always'});
    }
  }

  /** [text] as literal Markdown (a viewer's description is Markdown): no markup from the names in it. */
  private static markdownText(text: string): string {
    return text.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/[\\`*_[\]]/g, '\\$&');
  }

  /** The metrics grid detached and disposed. */
  private dropMetricsTable(): void {
    if (this.metricsTable !== undefined)
      ui.empty(ui.div([this.metricsTable.root]));
    this.metricsTable = undefined;
  }

  /** The two kept panels back under the hint, the rows removed (`ui.empty` detaches their viewers; the bottom blocks
   * stay) and no longer observed, the frames given back. */
  private unmount(): void {
    this.hintHost.append(this.metricsPanel, this.insightsPanel);
    for (const row of this.rows)
      ui.empty(row.root);
    this.rows = [];
    this.layout();
    this.mounted = new Set();
    if (this.shown !== undefined)
      releaseShown(this.shown);
    this.shown = undefined;
  }

  /** The content of [slot] removed: the two kept panels go back under the hint, the [kept] ones are detached for
   * their new slots, the rest are disposed. */
  private emptySlot(slot: HTMLElement, kept: ReadonlySet<HTMLElement>): void {
    for (const child of Array.from(slot.children)) {
      if (child === this.metricsPanel || child === this.insightsPanel)
        this.hintHost.append(child);
      else if (child instanceof HTMLElement && kept.has(child))
        child.remove();
    }
    ui.empty(slot);
  }

  /** Mounts [layout]: a row with the same panels and height as a mounted one keeps its elements (and the sizes the
   * user dragged) and gets the [panels] it does not hold yet; another row is built at its initial sizes, the slots
   * sharing the width in the ratio of the row's weights. The vertical splitter is rebuilt only when the sequence of
   * rows changes. */
  private mountRows(layout: ResultsRow[], panels: Map<PanelKey, HTMLElement>): void {
    const wanted = new Set(panels.values());
    const kept = new Set<RowRecord>();
    const rows = layout.map(({height, panels: keys, weights}) => {
      const signature = `${keys.join('|')}@${height}`;
      const contents = keys.map((key) => panels.get(key) ?? ui.div([]));
      const reused = this.rows.find((r) => r.signature === signature && !kept.has(r));
      if (reused !== undefined) {
        kept.add(reused);
        for (const [i, slot] of reused.slots.entries()) {
          if (slot.firstElementChild === contents[i])
            continue;
          this.emptySlot(slot, wanted);
          slot.append(contents[i]);
        }
        return reused;
      }
      // The splitter's dividers bind the slots themselves, so a reused slot takes new content.
      const {root, slots} = ResultsPane.weightedRow(contents, weights);
      root.style.height = `${height}px`;
      return {signature, root, slots, weights, height: root.style.height};
    });
    for (const row of this.rows.filter((r) => !kept.has(r))) {
      for (const slot of row.slots)
        this.emptySlot(slot, wanted);
      ui.empty(row.root);
    }
    this.rows = rows;
    this.mounted = new Set(layout.flatMap((row) => row.panels));
    this.layout();
  }

  /** While the rows and the bottom blocks are shown their heights are noted (a drag, a resize); while hidden, a changed
   * height goes back. A row whose slots all lost their growth (measured while hidden) gets its weights back. */
  private keepSizes(): void {
    const isShown = this.rowsHost.getClientRects().length > 0;
    for (const box of [...this.rows, ...this.bottom]) {
      if (isShown)
        box.height = box.root.style.height;
      else if (box.root.style.height !== box.height)
        box.root.style.height = box.height;
    }
    for (const row of this.rows) {
      if (row.slots.every((slot) => !(Number(slot.style.flexGrow) > 0)))
        ResultsPane.grow(row);
    }
  }

  /** The slots of [row] share it in the ratio of its weights. */
  private static grow({slots, weights}: Pick<RowRecord, 'slots' | 'weights'>): void {
    for (const [i, slot] of slots.entries())
      slot.style.flexGrow = `${weights[i]}`;
  }

  /** A row of results (`forge-results-row`): a horizontal splitter of [contents], each in its box (of [boxClass]), the
   * boxes sharing the width in the ratio of [weights]: preset widths, from which the splitter derives their growth,
   * then the growth set after it. */
  static weightedRow(contents: HTMLElement[], weights: number[], boxClass?: string):
    {root: HTMLElement; slots: HTMLElement[]} {
    const slots = contents.map((content, i) => {
      const slot = ui.box(null, boxClass);
      slot.append(content);
      slot.style.width = `${weights[i] * SLOT_WIDTH}px`;
      return slot;
    });
    const root = ui.splitH(slots, null, true);
    root.classList.add('forge-results-row');
    ResultsPane.grow({slots, weights});
    return {root, slots};
  }

  /** The method's own views; a failure is reported and gives none. */
  private async visualizeElements({engine, columns, target, prediction, blob}: ResultsData): Promise<HTMLElement[]> {
    const visualize = engine.functions.visualize;
    if (visualize === undefined)
      return [];
    let result: unknown;
    try {
      result = await onFrame(columns, (df) => visualize.apply({df, targetColumn: target, predictColumn: prediction,
        model: blob}));
    } catch (e) {
      reportError(e, true);
      return [];
    }
    return (Array.isArray(result) ? result : [result]).filter((e): e is HTMLElement => e instanceof HTMLElement);
  }
}

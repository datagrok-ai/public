/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import * as echarts from 'echarts';
import 'echarts-wordcloud';

import $ from 'cash-dom';
import {Observable, Subject} from 'rxjs';

import {ERROR_CLASS, MessageHandler, unsubscribeAll} from '../../utils/utils';
import {wordCloudStatus} from './word-cloud-status';


const MAX_UNIQUE_CATEGORIES_NUMBER = 500;

@grok.decorators.viewer({
  name: 'Word cloud',
  description: 'Creates a word cloud viewer',
  icon: 'icons/wordcloud-viewer.svg',
  toolbox: true,
})
export class WordCloudViewer extends DG.JsViewer {
  wordColumnName: string;
  shape: string;
  minTextSize: number;
  maxTextSize: any;
  minRotationDegree: number;
  maxRotationDegree: number;
  rotationStep: number;
  gridSize: number;
  drawOutOfBound: boolean;
  fontFamily: string;
  bold: boolean;
  strColumns: DG.Column[];
  initialized: boolean;
  chart: any; //echarts.EChartsType
  private _counts = new Map<string, number>();
  private _error: string | null = null;
  private _renderPending = 0;
  private _renderTimer: any = null;
  private _onRendered = new Subject<void>();

  /** Fires after every render pass — what automation settles on together with `isRenderPending`. */
  get onRendered(): Observable<void> {return this._onRendered;}

  get isRenderPending(): boolean {return this._renderPending > 0;}

  /** Word to rows over the viewer's own filter, as the last render counted them; empty while
   * `renderError` is set. */
  get wordCounts(): Map<string, number> {return this._counts;}

  /** The message the viewer put in place of the cloud, `null` while the cloud is drawn. */
  get renderError(): string | null {return this._error;}

  getWidgetStatus(): DG.IWidgetStatus {return wordCloudStatus(this);}

  constructor() {
    super();

    this.wordColumnName = this.string('wordColumnName', '', {columnTypeFilter: DG.COLUMN_TYPE.STRING});

    this.shape = this.string('shape', 'circle', {
      choices: ['circle', 'diamond', 'triangle-forward', 'triangle', 'pentagon', 'star'],
    });

    this.minTextSize = this.int('minTextSize', 14);
    this.maxTextSize = this.int('maxTextSize', 100);

    this.minRotationDegree = this.int('minRotationDegree', -30);
    this.maxRotationDegree = this.int('maxRotationDegree', 30);
    this.rotationStep = this.int('rotationStep', 5, {min: 1});

    this.gridSize = this.int('gridSize', 8);

    this.drawOutOfBound = this.bool('drawOutOfBound', true);

    this.fontFamily = this.string('fontFamily', 'sans-serif', {choices: ['sans-serif', 'serif', 'monospace']});

    this.bold = this.bool('bold', true);

    this.strColumns = [];
    this.initialized = false;
  }

  init() {
    this.initialized = true;
  }

  _testColumns() {
    const columns = this.dataFrame.columns.toList();
    const strCols = columns.filter((col) => col.type === DG.TYPE.STRING);
    return strCols.length >= 1;
  }

  onTableAttached() {
    unsubscribeAll(this.subs);
    this.addSubs();

    this.init();

    const columns = this.dataFrame.columns.toList();
    this.strColumns = columns.filter((col) => col.type === DG.TYPE.STRING);

    if (this._testColumns())
      this.wordColumnName = this.strColumns.filter((col) => col.categories.length <= MAX_UNIQUE_CATEGORIES_NUMBER && col.categories.length > 1)[0]?.name ?? '';

    this.render();
  }

  addSubs() {
    // the flag has to go up when the change arrives, not 50 ms later when the render runs
    for (const stream of [this.dataFrame.filter.onChanged, ui.onSizeChanged(this.root)]) {
      this.subs.push(stream.subscribe((_: any) => this._renderPending = 1));
      this.subs.push(DG.debounce(stream, 50).subscribe((_: any) => this.render()));
    }
  }

  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);
    if (this.initialized && this._testColumns())
      this.render();
  }

  onSourceRowsChanged() {
    this.render();
  }

  detach() {
    if (this._renderTimer !== null)
      clearTimeout(this._renderTimer);
    this._renderTimer = null;
    this._renderPending = 0;
    super.detach();
  }

  /** The message replaces the cloud, so the words the previous frame drew are no longer on screen —
   * `render` returns before re-creating the chart, and `chart` would still hold their geometry. */
  private showError(message: string) {
    this._error = message;
    this._counts = new Map();
    MessageHandler._showMessage(this.root, message, ERROR_CLASS);
    this.renderFinished();
  }

  private renderFinished(attempt = 0) {
    this._renderTimer = null;
    if (this._renderPending === 0)
      return;
    // The counts are in place before `setOption`, so a status read between the two would report
    // four words over no picture. A cloud that owes a canvas is not rendered yet; under load
    // zrender can take more than the one frame the happy path needs.
    if (this._error === null && this.wordColumnName && this.root.querySelector('canvas') === null && attempt < 20) {
      this._renderTimer = setTimeout(() => requestAnimationFrame(() => this.renderFinished(attempt + 1)));
      return;
    }
    this._renderPending = 0;
    this._onRendered.next();
  }

  render() {
    this._renderPending = 1;
    if (!this._testColumns()) {
      this.showError('Not enough data to produce the result.');
      return;
    }
    if (!this.wordColumnName || this.dataFrame.getCol(this.wordColumnName).categories.length > MAX_UNIQUE_CATEGORIES_NUMBER) {
      this.showError('The Word cloud viewer requires categorical column with 500 or fewer unique categories');
      return;
    }

    this._error = null;
    $(this.root).empty();

    if (this.wordColumnName === null || this.wordColumnName === '') {
      this._counts = new Map();
      // A viewer that has just been added has no column until the property default lands, and it
      // renders again when it does. Declaring this frame finished would let a settle read the empty
      // one as the answer — the status would then look exactly like the message state. The timer is
      // the release bound, for a column the user cleared on purpose and that will never arrive.
      this._renderTimer = setTimeout(() => requestAnimationFrame(() => this.renderFinished()), 300);
      return;
    }

    const margin = {top: 10, right: 10, bottom: 10, left: 10};
    const width = this.root.parentElement!.clientWidth - margin.left - margin.right;
    const height = this.root.parentElement!.clientHeight - margin.top - margin.bottom;
    const strColumn = this.dataFrame.getCol(this.wordColumnName);
    const table = this.dataFrame;

    const counts = new Map<string, number>();
    for (const i of this.filter.getSelectedIndexes()) {
      const word = strColumn.get(i);
      counts.set(word, (counts.get(word) ?? 0) + 1);
    }
    this._counts = counts;
    const data = Array.from(counts, ([name, value]) => ({
      name: name,
      value: value,
      textStyle: {
        color: DG.Color.toHtml(DG.Color.getCategoryColor(strColumn, name)),
      },
    }));

    if (this.chart !== undefined)
      this.chart.dispose();

    this.chart = echarts.init(<HTMLDivElement | HTMLCanvasElement> this.root);

    this.chart.setOption({
      width: width + margin.left + margin.right,
      height: height + margin.top + margin.bottom,
      series: [{
        type: 'wordCloud',
        shape: this.shape,
        left: 'center',
        top: 'center',
        width: `${width}`,
        height: `${height}`,
        right: null,
        bottom: null,
        sizeRange: [this.minTextSize, this.maxTextSize],
        gridSize: this.gridSize,
        rotationRange: [this.minRotationDegree, this.maxRotationDegree],
        rotationStep: this.rotationStep,
        drawOutOfBound: this.drawOutOfBound,
        textStyle: {
          fontFamily: this.fontFamily,
          fontWeight: this.bold ? 'bold' : 'normal',
        },
        emphasis: {
          focus: 'self',
          textStyle: {
            shadowBlur: 10,
            shadowColor: '#333',
          },
        },
        data: data,
      }],
    });

    this.chart
      .on('mouseover', (d: any) => ui.tooltip.showRowGroup(table, (i) => {
        return d.name === strColumn.get(i);
      }, d.event.event.x + 10, d.event.event.y + 10))
      .on('mouseout', () => ui.tooltip.hide())
      .on('mousedown', (d: any) => {
        table.selection.handleClick((i) => {
          return d.name === strColumn.get(i);
        }, d.event.event);
      });

    // `layoutAnimation` is absent from the wordcloud series defaults, so the layout helper takes its
    // synchronous branch and lays every word out in the macrotask `setOption` queued; zrender paints
    // them on the frame after that. echarts' own `finished` fires before either.
    this._renderTimer = setTimeout(() => requestAnimationFrame(() => this.renderFinished()));
  }
}

/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {mergeChartOptions, mergeSeries, getOrCreateParsedChartData, mergeProperties, chartDataId} from './fit-chart-data';
import {FitChartData, fitChartDataProperties, IFitChartData, IFitChartOptions, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {debounce, tap} from 'rxjs/operators';
import {interval, merge, Observable, Subject} from 'rxjs';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {FitChartCellRenderer} from './fit-renderer';
import {chartTooltip} from './fit-interaction';
import {drawnLegend, LEGEND_CORNERS} from './fit-legend';

const ERROR_CLASS = 'd4-viewer-error';
/** The properties the viewer tells "set" from "unset" by its field. */
const FIELD_SET_PROPERTIES = ['logX', 'logY', 'allowXZeroes', 'mergeCellSeries', 'showColumnLabel'];

@grok.decorators.viewer({
  name: 'MultiCurveViewer',
  description: 'A viewer that superimposes multiple in-cell curves on one chart',
  icon: 'icons/multi-curve-viewer.png',
  trellisable: true
})
export class MultiCurveViewer extends DG.JsViewer {
  [index: string]: any;
  canvas: HTMLCanvasElement;
  gridCellWidget: DG.GridCellWidget;
  curvesColumnNames?: string[] = [];
  legendColumnName?: string;
  showSelectedRowsCurves: boolean = true;
  showCurrentRowCurve: boolean = true;
  showMouseOverRowCurve: boolean = true;
  mergeColumnSeries: boolean = false;
  showOutliers?: boolean = true;
  rows: number[] = [];
  data: IFitChartData = new FitChartData();
  fitIdentities: (string | undefined)[] = [];
  logX?: boolean;
  logY?: boolean;
  allowXZeroes?: boolean;
  mergeCellSeries?: boolean;
  showColumnLabel?: boolean;
  private renderPending = false;
  private rendered = new Subject<void>();

  /** Fires after every render - what automation settles on with {@link isRenderPending}. */
  get onRendered(): Observable<void> { return this.rendered; }

  /** Whether a change of the current, mouse-over or selected rows is waiting for its debounced render. */
  get isRenderPending(): boolean { return this.renderPending; }

  /** What a trellis plot puts on its control panel. */
  get trellisProperties(): string[] { return ['mergeColumnSeries', 'logX', 'logY']; }

  private isInTrellis(): boolean {
    let curRoot = this.root;
    let i = 0;
    const maxDepth = 7;
    while (curRoot && i < maxDepth) {
      if (curRoot.classList.contains('d4-trellis-plot-cell'))
        return true;
      curRoot = curRoot.parentElement!;
      i++;
    }
    return false;
  }

  constructor() {
    super();
    this.gridCellWidget = DG.GridCellWidget.fromGridCell(this.createGridCell(''));
    this.canvas = this.gridCellWidget.canvas;
    this.root.append(this.gridCellWidget.canvas);
    ui.tools.handleResize(this.root, (w: number, h: number) => {
      this.canvas.width = w;
      this.canvas.height = h;
      this.render();
    });
    // the grid routes these to the cell renderer; here the viewer owns the canvas and does it itself
    this.canvas.addEventListener('mousemove', (e: MouseEvent) => {
      if ((this.data.series?.length ?? 0) > 0 &&
        !chartTooltip(this.data, this.chartBounds(), e, {x: e.offsetX, y: e.offsetY}))
        ui.tooltip.hide();
    });
    this.canvas.addEventListener('mouseleave', () => ui.tooltip.hide());

    this.curvesColumnNames = this.addProperty('curvesColumnNames', DG.TYPE.COLUMN_LIST, [], {semType: FitConstants.FIT_SEM_TYPE});
    this.legendColumnName = this.addProperty('legendColumnName', DG.TYPE.STRING, '', {description: 'Column to be used for curves names'});

    this.showSelectedRowsCurves = this.bool('showSelectedRowsCurves', true, {description: 'Adds curves from the selected rows'});
    this.showCurrentRowCurve = this.bool('showCurrentRowCurve', true);
    this.showMouseOverRowCurve = this.bool('showMouseOverRowCurve', true);
    this.mergeColumnSeries = this.bool('mergeColumnSeries', false);
    this.showOutliers = this.bool('showOutliers', true);

    for (const p of fitChartDataProperties) {
      const name = p.name === 'mergeSeries' ? 'mergeCellSeries' : p.name;
      const value = this.addProperty(name, p.propertyType, p.defaultValue, {...p.options, showSlider: false});
      // the panel reads the field, so a default the render applies on its own has to be in it: Show Legend
      // read unchecked over a legend that was drawn
      if (value != null && !FIELD_SET_PROPERTIES.includes(name))
        this[name] = value;
    }
  }

  static fromChartData(chartData: IFitChartData): MultiCurveViewer {
    const viewer = new MultiCurveViewer();
    viewer.data = chartData;
    viewer.render();
    return viewer;
  }

  createGridCell(value: string): DG.GridCell {
    const gridCell = DG.GridCell.fromValue(value);
    gridCell.cellType = FitConstants.FIT_CELL_TYPE;
    return gridCell;
  }

  applyViewerProperties(): void {
    mergeProperties(fitChartDataProperties, this, this.data.chartOptions);
  }

  mergeViewerChartOptions(): void {
    for (const p of fitChartDataProperties) {
      const isDefaultValueChanged = FIELD_SET_PROPERTIES.includes(p.name === 'mergeSeries' ?
        'mergeCellSeries' : p.name) ? this[p.name === 'mergeSeries' ? 'mergeCellSeries' : p.name] !== undefined :
        this.props.get(p.name) !== undefined && this.props.get(p.name) !== p.defaultValue;
      if (isDefaultValueChanged) {
        if (FIELD_SET_PROPERTIES.includes(p.name === 'mergeSeries' ? 'mergeCellSeries' : p.name))
          this.data.chartOptions![p.name as keyof IFitChartOptions] = this[p.name];
        else
          this.data.chartOptions![p.name as keyof IFitChartOptions] = this.props.get(p.name) as unknown as any;
      }
    }
  }

  createChartData(): void {
    if (this.curvesColumnNames?.length === 0)
      return;
    this.rows.length = 0;
    let selectionStart = -1;
    if (this.isInTrellis()) { this.rows.push(...this.dataFrame.filter.getSelectedIndexes()); } else {
      if (this.showCurrentRowCurve && this.dataFrame.currentRowIdx !== -1)
        this.rows.push(this.dataFrame.currentRowIdx);
      if (this.showMouseOverRowCurve && this.dataFrame.mouseOverRowIdx !== -1)
        this.rows.push(this.dataFrame.mouseOverRowIdx);
      if (this.showSelectedRowsCurves) {
        selectionStart = this.rows.length;
        this.rows.push(...this.dataFrame.selection.clone().and(this.dataFrame.filter).getSelectedIndexes());
      }
    }

    this.data = new FitChartData();
    //const _grid = this.isInTrellis() ? grok.shell.tableView(this.dataFrame.name).grid : this.tableView?.grid!;
    const mergeCellSeries = this.props.get('mergeCellSeries') as unknown as boolean;
    const chartOptions: IFitChartOptions[] = [];
    this.fitIdentities = [];
    for (const colName of this.curvesColumnNames!) {
      const series = [];
      const identities: (string | undefined)[] = [];
      for (const i of new Set(this.rows)) {
        const tableCell = this.dataFrame.cell(i, colName);
        if (!tableCell || !tableCell.value)
          continue;
        const cellCurves = getOrCreateParsedChartData(tableCell);
        if (!cellCurves.series || cellCurves.series.length === 0)
          continue;
        // the parsed data is cached and shared with the grid, so the viewer works on its own copies -
        // merging into it used to rewrite the cell itself
        let cellSeries: IFitSeries[] = cellCurves.series.map((series) => ({...series, columnName: tableCell.column.name}));
        const legendCol = !this.legendColumnName || !this.dataFrame.col(this.legendColumnName) ? null : this.dataFrame.col(this.legendColumnName);
        if (legendCol) {
          for (const series of cellSeries)
            series.auxLegendName = cellSeries.length > 1 ? `${legendCol.get(i)} ${series.name}` : `${legendCol.get(i)}`;
        }
        const currentChartOptions = cellCurves.chartOptions;
        if (currentChartOptions != null)
          chartOptions[chartOptions.length] = currentChartOptions;
        if (mergeCellSeries) {
          const mergedSeries = mergeSeries(cellSeries)!;
          if (currentChartOptions?.title !== undefined && currentChartOptions?.title !== '')
            mergedSeries.name = currentChartOptions?.title;
          cellSeries = [mergedSeries];
        }
        // the cell these came from, so the render can reuse the fit instead of redoing it every time
        const cellId = chartDataId(cellCurves);
        identities.push(...cellSeries.map((_, j) => `${cellId} || idx: ${mergeCellSeries ? 'merged' : j}`));
        series.push(...cellSeries);
      }
      if (series.length === 0)
        continue;
      // a deep copy: the points are still the cell's own objects, and rendering mutates them
      this.data.series?.push(...JSON.parse(JSON.stringify(this.mergeColumnSeries ? [mergeSeries(series)!] : series)));
      // a column merge draws from many cells at once, so no one cell identifies its fit
      this.fitIdentities.push(...(this.mergeColumnSeries ? [undefined] : identities));
    }
    this.data.chartOptions = mergeChartOptions(chartOptions, this.data.series);
    this.data.chartOptions.useAuxLegendNames = true;
    this.mergeViewerChartOptions();
    // assigned, not set: the panel shows the name the data carries, and '' stays a name of none
    for (const name of ['title', 'xAxisName', 'yAxisName'])
      if (this[name] === undefined && this.data.chartOptions[name as keyof IFitChartOptions] !== undefined)
        this[name] = this.data.chartOptions[name as keyof IFitChartOptions];
    this.data.series?.forEach((series, i) => {
      if (!series.fitLineColor && !series.pointColor) {
        series.pointColor = DG.Color.toHtml(DG.Color.getCategoricalColor(this.data.series?.length! > 20 ? 0 : i));
        series.fitLineColor = DG.Color.toHtml(DG.Color.getCategoricalColor(this.data.series?.length! > 20 ? 0 : i));
      }
      series.showOutliers = this.showOutliers;
      series.showCurveConfidenceInterval = false;
      series.droplines = [];
      if (this.data.series?.length! > 20) {
        series.showPoints = '';
        series.lineStyle = 'solid';
      }
      if (this.data.series?.length! > 10 && this.data.series?.length! < 100 && i >= selectionStart) {
        const color = DG.Color.fromHtml(series.fitLineColor!);
        series.fitLineColor = `rgba(${DG.Color.r(color)}, ${DG.Color.g(color)}, ${DG.Color.b(color)}, 0.2)`;
      }
    });
  }

  onPropertyChanged(property: DG.Property | null): void {
    if (property?.name === 'curvesColumnNames') {
      this.props.set('showColumnLabel', this.curvesColumnNames?.length! > 1 as unknown as object);
      this.showColumnLabel = this.curvesColumnNames?.length! > 1;
    }
    if (FIELD_SET_PROPERTIES.includes(property?.name === 'mergeSeries' ? 'mergeCellSeries' : property?.name ?? ''))
      this[property?.name === 'mergeSeries' ? 'mergeCellSeries' : property?.name ?? ''] = this[property?.name ?? ''];
    this.applyViewerProperties();
    this.createChartData();
    this.render();
  }

  onTableAttached(): void {
    const grid = this.tableView?.grid!;
    const fitCol = this.dataFrame.columns.bySemType(FitConstants.FIT_SEM_TYPE);
    if (fitCol !== null)
      this.curvesColumnNames = [fitCol.name];

    this.subs.push(merge(...[this.dataFrame.onCurrentCellChanged, ...(grid ? [grid.onCellMouseEnter] : []), this.dataFrame.onSelectionChanged])
      .pipe(tap(() => this.renderPending = true), debounce((_) => interval(50)))
      .subscribe((_) => {
        this.renderPending = false;
        if (this.dataFrame) {
          this.createChartData();
          this.render();
        }
      }));
    this.createChartData();
    if (!this.isInTrellis())
      this.render();
  }

  _showErrorMessage(msg: string) {
    this.root.appendChild(ui.divText(msg, 'd4-viewer-error'));
  }

  _removeErrorMessage() {
    const divTextElement = this.root.getElementsByClassName(ERROR_CLASS)[0];
    if (divTextElement)
      this.root.removeChild(divTextElement);
  }

  render(): void {
    const g = this.canvas.getContext('2d')!;
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);
    this._removeErrorMessage();

    if (this.curvesColumnNames?.length === 0)
      this._showErrorMessage('The MultiCurveViewer viewer requires a minimum of 1 curves column.');
    else if (this.data?.series?.length === 0)
      this._showErrorMessage('No data to show.');
    else {
      // drawn from this very object rather than through a serialized cell, so that hovering can find the
      // legend the render laid out
      new FitChartCellRenderer().renderCurves(g, this.chartBounds(), this.data, undefined, this.fitIdentities);
    }
    this.rendered.next();
  }

  private chartBounds(): DG.Rect {
    return new DG.Rect(0, 0, this.canvas.width, this.canvas.height);
  }

  /** What the viewer shows, read off the last render: the curves, and the legend as it was laid out -
   * a hit area per row it drew and one for `+N more`, in CSS pixels of the canvas. */
  getWidgetStatus(): DG.IWidgetStatus {
    const status = super.getWidgetStatus();
    const drawn = this.curvesColumnNames?.length === 0 ? 0 : this.data?.series?.length ?? 0;
    const legend = drawn === 0 ? null : drawnLegend(this.data, this.chartBounds());
    const values: {[name: string]: number | string | boolean} = {...status.values, 'curves': drawn,
      'legend shown': legend !== null};
    const hitAreas = {...status.hitAreas};
    if (drawn > 0) {
      // the axes as the render lays them out, taking in the curves as well as the points (GROK-20661)
      const bounds = new FitChartCellRenderer().layOutChart(this.chartBounds(), this.data, undefined, this.fitIdentities)
        .dataBounds;
      const num = (x: number) => Number(x.toPrecision(4));
      values['x range'] = `${num(bounds.minX)} .. ${num(bounds.maxX)}`;
      values['y range'] = `${num(bounds.minY)} .. ${num(bounds.maxY)}`;
      // each axis at its own scale: concentrations in molar span 1e-9 .. 1e-5
      const outside = (v: number, min: number, max: number) => {
        const slack = 1e-9 * Math.max(Math.abs(max - min), Math.abs(min), Math.abs(max));
        return v < min - slack || v > max + slack;
      };
      values['points outside the axes'] = (this.data.series ?? []).reduce((n, series) => n + series.points.filter((p) =>
        outside(p.x, bounds.minX, bounds.maxX) || outside(p.y, bounds.minY, bounds.maxY)).length, 0);
    }
    if (legend) {
      const scale = this.canvas.width > 0 ? this.canvas.getBoundingClientRect().width / this.canvas.width : 1;
      const css = (r: DG.Rect) => ({x: r.x * scale, y: r.y * scale, width: r.width * scale, height: r.height * scale});
      const shortened: number[] = [];
      for (const [i, row] of legend.rows.entries()) {
        const whole = row.shown === row.entry.text;
        hitAreas[`${whole ? '' : 'shortened '}legend row ${i + 1}`] = css(row.box);
        values[`legend row ${i + 1}`] = row.shown;
        if (!whole)
          shortened.push(i + 1);
      }
      if (legend.moreBox)
        hitAreas['legend more'] = css(legend.moreBox);
      values['legend corner'] = LEGEND_CORNERS[legend.corner];
      values['legend rows'] = legend.rows.length;
      values['legend shortened rows'] = shortened.join(', ');
      values['legend hidden'] = legend.entries.length - legend.rows.length;
    }
    return {...status, parts: {...status.parts, canvas: this.canvas}, hitAreas, values};
  }
}

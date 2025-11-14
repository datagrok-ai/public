/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {mergeChartOptions, mergeSeries} from './fit-renderer';
import {getOrCreateParsedChartData, mergeProperties} from './fit-renderer';
import {FitChartData, fitChartDataProperties, IFitChartData, IFitChartOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {debounce} from 'rxjs/operators';
import {interval, merge} from 'rxjs';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';

const ERROR_CLASS = 'd4-viewer-error';

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
  logX?: boolean;
  logY?: boolean;
  allowXZeroes?: boolean;
  mergeCellSeries?: boolean;
  showColumnLabel?: boolean;


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

    this.curvesColumnNames = this.addProperty('curvesColumnNames', DG.TYPE.COLUMN_LIST, [], {semType: FitConstants.FIT_SEM_TYPE});
    this.legendColumnName = this.addProperty('legendColumnName', DG.TYPE.STRING, '', {description: 'Column to be used for curves names'});

    this.showSelectedRowsCurves = this.bool('showSelectedRowsCurves', true, {description: 'Adds curves from the selected rows'});
    this.showCurrentRowCurve = this.bool('showCurrentRowCurve', true);
    this.showMouseOverRowCurve = this.bool('showMouseOverRowCurve', true);
    this.mergeColumnSeries = this.bool('mergeColumnSeries', false);
    this.showOutliers = this.bool('showOutliers', true);

    for (const p of fitChartDataProperties)
      this.addProperty(p.name === 'mergeSeries' ? 'mergeCellSeries' : p.name, p.propertyType, p.defaultValue, {...p.options, showSlider: false});
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
      const isDefaultValueChanged = ['logX', 'logY', 'allowXZeroes', 'mergeCellSeries', 'showColumnLabel'].includes(p.name === 'mergeSeries' ?
        'mergeCellSeries' : p.name) ? this[p.name === 'mergeSeries' ? 'mergeCellSeries' : p.name] !== undefined :
        this.props.get(p.name) !== undefined && this.props.get(p.name) !== p.defaultValue;
      if (isDefaultValueChanged) {
        if (['title', 'xAxisName', 'yAxisName'].includes(p.name) && this.props.get(p.name) as unknown as string === '')
          continue;
        else if (['logX', 'logY', 'allowXZeroes', 'mergeCellSeries', 'showColumnLabel'].includes(p.name === 'mergeSeries' ? 'mergeCellSeries' : p.name))
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
    for (const colName of this.curvesColumnNames!) {
      const series = [];
      for (const i of new Set(this.rows)) {
        const tableCell = this.dataFrame.cell(i, colName);
        if (!tableCell || !tableCell.value)
          continue;
        const cellCurves = getOrCreateParsedChartData(tableCell);
        if (!cellCurves.series || cellCurves.series.length === 0)
          continue;
        cellCurves.series.forEach((series) => series.columnName = tableCell.column.name);
        const legendCol = !this.legendColumnName || !this.dataFrame.col(this.legendColumnName) ? null : this.dataFrame.col(this.legendColumnName);
        for (let j = 0; j < cellCurves.series.length; j++) {
          const series = cellCurves.series[j];
          if (legendCol)
            series.auxLegendName = cellCurves.series.length > 1 ? `${legendCol.get(i)} ${series.name}` : `${legendCol.get(i)}`;
        }
        const currentChartOptions = cellCurves.chartOptions;
        if (currentChartOptions != null)
          chartOptions[chartOptions.length] = currentChartOptions;
        if (mergeCellSeries) {
          const mergedSeries = mergeSeries(cellCurves.series!)!;
          if (currentChartOptions?.title !== undefined && currentChartOptions?.title !== '')
            mergedSeries.name = currentChartOptions?.title;
          cellCurves.series = [mergedSeries];
        }
        series.push(...cellCurves.series!);
      }
      if (this.mergeColumnSeries)
        this.data.series?.push(mergeSeries(series)!);
      else
        this.data.series?.push(...JSON.parse(JSON.stringify(series)));
    }
    this.data.chartOptions = mergeChartOptions(chartOptions);
    this.data.chartOptions.useAuxLegendNames = true;
    this.mergeViewerChartOptions();
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
    if (['logX', 'logY', 'allowXZeroes', 'mergeCellSeries', 'showColumnLabel'].includes(property?.name === 'mergeSeries' ? 'mergeCellSeries' : property?.name ?? ''))
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

    merge(...[this.dataFrame.onCurrentCellChanged, ...(grid ? [grid.onCellMouseEnter] : []), this.dataFrame.onSelectionChanged])
      .pipe(debounce((_) => interval(50)))
      .subscribe((_) => {
        if (this.dataFrame) {
          this.createChartData();
          this.render();
        }
      });
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

    if (this.curvesColumnNames?.length === 0) {
      this._showErrorMessage('The MultiCurveViewer viewer requires a minimum of 1 curves column.');
      return;
    }
    if (this.data?.series?.length === 0) {
      this._showErrorMessage('No data to show.');
      return;
    }
    this.gridCellWidget.gridCell = this.createGridCell(JSON.stringify(this.data));
  }
}

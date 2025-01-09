import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {CellRenderViewer} from '@datagrok-libraries/utils/src/viewers/cell-render-viewer';
import {FitChartCellRenderer, mergeChartOptions, mergeSeries} from './fit-renderer';
import {getOrCreateParsedChartData, mergeProperties} from './fit-renderer';
import {FitChartData, fitChartDataProperties, IFitChartData, IFitChartOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';

import {debounce} from 'rxjs/operators';
import {interval, merge} from 'rxjs';
import {FitConstants} from './const';


const ERROR_CLASS = 'd4-viewer-error';

@grok.decorators.viewer({
  name: 'MultiCurveViewer',
  description: 'A viewer that superimposes multiple in-cell curves on one chart',
  icon: 'icons/multi-curve-viewer.png',
})
export class MultiCurveViewer extends CellRenderViewer<FitChartCellRenderer> {
  [index: string]: any;
  curvesColumnNames?: string[] = [];
  showSelectedRowsCurves: boolean = false;
  showCurrentRowCurve: boolean = true;
  showMouseOverRowCurve: boolean = true;
  mergeColumnSeries: boolean = false;
  rows: number[] = [];
  data: IFitChartData = new FitChartData();
  logX?: boolean;
  logY?: boolean;
  allowXZeroes?: boolean;
  mergeCellSeries?: boolean;
  showColumnLabel?: boolean;

  constructor() {
    super(new FitChartCellRenderer());

    this.curvesColumnNames = this.addProperty('curvesColumnNames', DG.TYPE.COLUMN_LIST, [], {semType: FitConstants.FIT_SEM_TYPE});

    this.showSelectedRowsCurves = this.bool('showSelectedRowsCurves', false, { description: 'Adds curves from the selected rows'});
    this.showCurrentRowCurve = this.bool('showCurrentRowCurve', true);
    this.showMouseOverRowCurve = this.bool('showMouseOverRowCurve', true);
    this.mergeColumnSeries = this.bool('mergeColumnSeries', false);

    for (const p of fitChartDataProperties)
      this.addProperty(p.name === 'mergeSeries' ? 'mergeCellSeries' : p.name, p.propertyType, p.defaultValue, p.options);
  }

  static fromChartData(chartData: IFitChartData): MultiCurveViewer {
    let viewer = new MultiCurveViewer();
    viewer.data = chartData;
    viewer.render();
    return viewer;
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
    if (this.showCurrentRowCurve && this.dataFrame.currentRowIdx !== -1)
      this.rows.push(this.dataFrame.currentRowIdx);
    if (this.showMouseOverRowCurve && this.dataFrame.mouseOverRowIdx !== -1)
      this.rows.push(this.dataFrame.mouseOverRowIdx);
    if (this.showSelectedRowsCurves && this.dataFrame.mouseOverRowIdx !== -1) {
      selectionStart = this.rows.length;
      this.rows.push(...this.dataFrame.selection.getSelectedIndexes());
    }

    this.data = new FitChartData();
    const grid = this.tableView?.grid!;
    const mergeCellSeries = this.props.get('mergeCellSeries') as unknown as boolean;
    const chartOptions: IFitChartOptions[] = [];
    for (const colName of this.curvesColumnNames!) {
      const series = [];
      for (const i of new Set(this.rows)) {
        const gridCell = grid.cell(colName, grid.tableRowToGrid(i));
        if (gridCell.cell.value === '')
          continue;
        const cellCurves = getOrCreateParsedChartData(gridCell);
        cellCurves.series?.forEach((series) => series.columnName = gridCell.cell.column.name);
        const currentChartOptions = cellCurves.chartOptions;
        if (currentChartOptions !== undefined && currentChartOptions !== null)
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
        this.data.series?.push(...series);
    }
    this.data.chartOptions = mergeChartOptions(chartOptions);
    this.mergeViewerChartOptions();
    this.data.series?.forEach((series, i) => {
      series.pointColor = DG.Color.toHtml(DG.Color.getCategoricalColor(this.data.series?.length! > 20 ? 0 : i));
      series.fitLineColor = DG.Color.toHtml(DG.Color.getCategoricalColor(this.data.series?.length! > 20 ? 0 : i));
      series.showCurveConfidenceInterval = false;
      series.droplines = [];
      if (this.data.series?.length! > 20) {
        series.showPoints = '';
        series.lineStyle = 'solid';
      }
      if (this.data.series?.length! > 10 && this.data.series?.length! < 100 && i >= selectionStart) {
        const color = DG.Color.fromHtml(series.fitLineColor);
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

    merge(this.dataFrame.onCurrentCellChanged, grid.onCellMouseEnter, this.dataFrame.onSelectionChanged)
      .pipe(debounce(_ => interval(50)))
      .subscribe(_ => {
        this.createChartData();
        this.render();
      });
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
    const g = this.canvas.getContext('2d')!
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);

    if (this.curvesColumnNames?.length === 0) {
      this._showErrorMessage('The MultiCurveViewer viewer requires a minimum of 1 curves column.');
      return;
    }
    this._removeErrorMessage();
    this.renderer.renderCurves(g, FitChartCellRenderer.inflateScreenBounds(new DG.Rect(0, 0, this.canvas.width, this.canvas.height)), this.data);
  }
}

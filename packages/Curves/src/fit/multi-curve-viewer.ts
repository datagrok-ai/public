import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {CellRenderViewer} from './cell-render-viewer';
import {FitChartCellRenderer} from './fit-renderer';
import {getChartData, mergeProperties} from './fit-renderer';
import {FIT_SEM_TYPE, FitChartData, fitChartDataProperties, IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';

import {debounce} from 'rxjs/operators';
import {interval, merge} from 'rxjs';


export class MultiCurveViewer extends CellRenderViewer<FitChartCellRenderer> {
  curvesColumnNames?: string[] = [];
  showSelectedRowsCurves: boolean = false;
  showCurrentRowCurve: boolean = true;
  showMouseOverRowCurve: boolean = true;
  rows: number[] = [];
  data: IFitChartData = new FitChartData();

  constructor() {
    super(new FitChartCellRenderer());

    this.curvesColumnNames = this.addProperty('curvesColumnNames', DG.TYPE.COLUMN_LIST);

    this.showSelectedRowsCurves = this.bool('showSelectedRowsCurves', false, { description: 'Adds curves from the selected rows'});
    this.showCurrentRowCurve = this.bool('showCurrentRowCurve', true);
    this.showMouseOverRowCurve = this.bool('showMouseOverRowCurve', true);

    for (const p of fitChartDataProperties)
      this.addProperty(p.name, p.propertyType, p.defaultValue, p.options);
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

  createChartData(): void {
    if (this.curvesColumnNames?.length === 0)
      return;
    this.rows.length = 0;
    if (this.showCurrentRowCurve && this.dataFrame.currentRowIdx !== -1)
      this.rows.push(this.dataFrame.currentRowIdx);
    if (this.showMouseOverRowCurve && this.dataFrame.mouseOverRowIdx !== -1)
      this.rows.push(this.dataFrame.mouseOverRowIdx);
    if (this.showSelectedRowsCurves && this.dataFrame.mouseOverRowIdx !== -1)
      this.rows.push(...this.dataFrame.selection.getSelectedIndexes());

    this.data = new FitChartData();
    const grid = this.tableView?.grid!;
    this.data.chartOptions!.showColumnLabel = this.props.get('showColumnLabel') as unknown as boolean;
    for (const colName of this.curvesColumnNames!)
      for (let i of this.rows) {
        const gridCell = grid.cell(colName, grid.tableRowToGrid(i));
        const cellCurves = getChartData(gridCell);
        cellCurves.series?.forEach((series) => series.columnName = gridCell.cell.column.name);
        this.data.series?.push(...cellCurves.series!);
      }
  }

  onPropertyChanged(property: DG.Property | null): void {
    this.applyViewerProperties();
    this.createChartData();
    this.render();
  }

  onTableAttached(): void {
    const grid = this.tableView?.grid!;
    const fitCol = this.dataFrame.columns.bySemType(FIT_SEM_TYPE);
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

  render(): void {
    if (this.curvesColumnNames?.length === 0) {
      this._showErrorMessage('The MultiCurveViewer viewer requires a minimum of 1 curves column.');
      return;
    }

    const g = this.canvas.getContext('2d')!
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);
    this.renderer.renderCurves(g, 0, 0, this.canvas.width, this.canvas.height, this.data);
  }
}

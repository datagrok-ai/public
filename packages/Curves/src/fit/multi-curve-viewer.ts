import * as DG from 'datagrok-api/dg';
import {CellRenderViewer} from "./cell-render-viewer";
import {FitChartCellRenderer} from "./fit-renderer";
import {
  FIT_SEM_TYPE,
  FitChartData,
  fitChartDataProperties,
  getChartData,
  IFitChartData,
  mergeProperties
} from "./fit-data";
import {debounce, filter} from 'rxjs/operators';
import {interval, merge} from 'rxjs';


export class MultiCurveViewer extends CellRenderViewer<FitChartCellRenderer> {
  curvesColumnName?: string;
  showSelectedRowsCurves: boolean = false;
  showCurrentRowCurve: boolean = true;
  showMouseOverRowCurve: boolean = true;
  rows: number[] = [];
  data: IFitChartData = new FitChartData();

  constructor() {
    super(new FitChartCellRenderer());

    this.curvesColumnName = this.column('curves');

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

  applyViewerProperties() {
    mergeProperties(fitChartDataProperties, this, this.data.chartOptions);
  }

  createChartData() {
    this.rows.length = 0;
    if (this.showCurrentRowCurve && this.dataFrame.currentRowIdx != -1)
      this.rows.push(this.dataFrame.currentRowIdx);
    if (this.showMouseOverRowCurve && this.dataFrame.mouseOverRowIdx != -1)
      this.rows.push(this.dataFrame.mouseOverRowIdx);
    if (this.showSelectedRowsCurves && this.dataFrame.mouseOverRowIdx != -1)
      this.rows.push(...this.dataFrame.selection.getSelectedIndexes());

    this.data = new FitChartData();
    const grid = this.tableView?.grid!;
    for (let i of this.rows) {
      const cellCurves = getChartData(grid.cell(this.curvesColumnName!, grid.tableRowToGrid(i)));
      this.data.series?.push(...cellCurves.series!);
    }
  }

  onPropertyChanged(property: DG.Property | null) {
    this.applyViewerProperties();
    this.render();
  }

  onTableAttached() {
    const grid = this.tableView?.grid!;
    this.curvesColumnName ??= this.dataFrame.columns.bySemType(FIT_SEM_TYPE)?.name;

    merge(this.dataFrame.onCurrentCellChanged, grid.onCellMouseEnter, this.dataFrame.onSelectionChanged)
      .pipe(debounce(_ => interval(50)))
      .subscribe(_ => {
        this.createChartData();
        this.render();
      });
  }

  render() {
    const g = this.canvas.getContext('2d')!
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);
    this.renderer.renderCurves(g, 0, 0, this.canvas.width, this.canvas.height, this.data);
  }
}

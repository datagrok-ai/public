import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { Grid } from 'datagrok-api/dg';

const CHART_TYPES = {
  histogram: DG.VIEWER.HISTOGRAM,
  barchart: DG.VIEWER.BAR_CHART,
  linechart: DG.VIEWER.LINE_CHART,
  piechart: DG.VIEWER.PIE_CHART,
};
const CHART_TYPE_TAG = 'chart_type';
export class SubstituentAnalysisViewer extends DG.JsViewer {
  initialized: boolean = false;
  name = 'substituent';
  grouppingCols: string[] = [];
  colsToAnalyze: string[] = [];
  totalCols: string[] = [];
  grouppedDf: DG.DataFrame | undefined = undefined;
  grouppedGridDiv = ui.div();
  mainView = ui.splitV([]);

  constructor() {
    super();
  }

  init(): void {
    this.initialized = true;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async onTableAttached(): Promise<void> {
    this.init();
    const grouppingColsDiv = ui.div();
    const colsToAnalyzeDiv = ui.div();
    this.totalCols = this.dataFrame.columns.names().concat(['']);
    this.updateColumnChoices(this.grouppingCols, 'Group by', grouppingColsDiv);
    this.updateColumnChoices(this.colsToAnalyze, 'Values', colsToAnalyzeDiv);
    this.mainView.append(ui.box(ui.panel([grouppingColsDiv], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    this.mainView.append(ui.box(ui.panel([colsToAnalyzeDiv], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    this.mainView.append(this.grouppedGridDiv);
    this.root.append(this.mainView);
  }

  updateColumnChoices(selectedCols: string[], choicesInputsName: string, choicesDiv: HTMLDivElement) {
    const colsChoicesDiv = ui.divH([]);
    colsChoicesDiv.append(ui.divText(
      choicesInputsName, {style: {fontWeight: 'bold', paddingTop: '10px', paddingRight: '10px'}}));
    for (let i = 0; i < selectedCols.length + 1; i++) {
      const selectedValue = selectedCols.length === i ? '' : selectedCols[i];
      const groupChoiceInput = ui.choiceInput('', selectedValue, this.totalCols);
      ui.tooltip.bind(groupChoiceInput.root, () => selectedValue);
      groupChoiceInput.onChanged(() => {
        if (groupChoiceInput.value === '')
          selectedCols.splice(i, 1);
        else {
          selectedCols.length === i ? selectedCols.push(groupChoiceInput.value!) :
            selectedCols[i] = groupChoiceInput.value!;
        }
        this.updateColumnChoices(selectedCols, choicesInputsName, choicesDiv);
        this.updateGrid();
      });
      groupChoiceInput.input.style.width = '100px';
      colsChoicesDiv.append(groupChoiceInput.root);
    }
    ui.empty(choicesDiv);
    choicesDiv.append(colsChoicesDiv);
  }

  updateGrid() {
    this.grouppedDf = this.dataFrame.groupBy(this.grouppingCols).aggregate();
    const grid = this.grouppedDf.plot.grid();
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    for (const col of this.grouppingCols)
      grid.col(col)!.width = 150;

    for (const col of this.colsToAnalyze) {
      this.grouppedDf.columns.addNewString(col).init((i) => '');
      grid.col(col)!.cellType = 'html';
      grid.col(col)!.width = 150;
    }

    grid.onCellPrepare((gc) => {
      if (gc.isTableCell && gc.cell.value != null && this.colsToAnalyze.includes(gc.gridColumn.name)) {
        const condition: {[key: string]: any} = {};
        for (const col of this.grouppingCols)
          condition[col] = this.grouppedDf!.get(col, gc.tableRowIndex!);
        const colList = [];
        this.grouppingCols.forEach((col) => colList.push(this.dataFrame.col(col)!));
        colList.push(this.dataFrame.col(gc.gridColumn.name)!);
        const df = DG.DataFrame.fromColumns(colList);
        df.columns.addNewString('group').init((j) => {
          for (const group of this.grouppingCols) {
            if (df.get(group, j) !== this.grouppedDf!.get(group, gc.tableRowIndex!))
              return 'other';
          }
          return 'current_group';
        });
        let chartType = this.dataFrame.col(gc.gridColumn.name)!.getTag(CHART_TYPE_TAG);
        if (!chartType) {
          chartType = Object.keys(CHART_TYPES)[0];
          this.dataFrame.col(gc.gridColumn.name)!.setTag(CHART_TYPE_TAG, Object.keys(CHART_TYPES)[0]);
        }
        const viewer = DG.Viewer.fromType((CHART_TYPES as any)[chartType], df);
        gc.element = viewer.root;
        df.rows.match(condition).filter();
      }
    });

    grid.onCellClick.subscribe((gc) => {
      if (!gc.isTableCell && this.colsToAnalyze.includes(gc.gridColumn.name))
        this.createPropertyPanel(gc.gridColumn.name, this.dataFrame.col(gc.gridColumn.name)!.getTag(CHART_TYPE_TAG), grid);
    });

    ui.empty(this.grouppedGridDiv);
    this.grouppedGridDiv.append(grid.root);
  }

  createPropertyPanel(colName: string, chartType: string, grid: DG.Grid) {
    const chartChoice = ui.choiceInput('Chart', chartType, Object.keys(CHART_TYPES));
    chartChoice.onChanged(() => {
      this.dataFrame.col(colName)!.setTag(CHART_TYPE_TAG, chartChoice.value!);
      grid.invalidate();
    });
    grok.shell.o = chartChoice.root;
  }
}

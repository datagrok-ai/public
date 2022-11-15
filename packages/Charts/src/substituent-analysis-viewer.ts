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

const AGGREGATE_TYPES = {
  'min': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {return query.min(colName, resColName);},
  'max': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {return query.max(colName, resColName);},
  'avg': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {return query.avg(colName, resColName);},
  'count': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {return query.count(resColName);},
};

interface IAnalyzedColumn {
  colName: string, // name of analyzed column
  type: string, // aggregation or chart
  typeName: string //name of exact aggr method or exact chart
}

const AGGR_TYPE = 'aggr';
const CHART_TYPE = 'chart';
export class SubstituentAnalysisViewer extends DG.JsViewer {
  initialized: boolean = false;
  name = 'substituent';
  grouppingCols: string[] = [];
  colsToAnalyze: IAnalyzedColumn[] = [];
  totalCols: string[] = [];
  grouppedDf: DG.DataFrame | undefined = undefined;
  grouppedGridDiv = ui.div();
  mainView = ui.splitV([]);
  parentViewers: {[key: string]: DG.Viewer} = {};
  viewersStorage: {[key: string]: {[key: number]: DG.Viewer}} = {};

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
    this.totalCols = this.dataFrame.columns.names().concat(['']);
    this.grouppingCols.push(this.totalCols[0]);
    this.updateColumnChoices(this.grouppingCols, 'Group by', grouppingColsDiv);
    this.mainView.append(ui.box(ui.panel([grouppingColsDiv], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    const addColToAnalyze = ui.icons.add(() => {this.createAddColumnDialog();}, 'Add column to analyze');
    this.mainView.append(ui.box(ui.panel([addColToAnalyze], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    this.mainView.append(this.grouppedGridDiv);
    this.root.append(this.mainView);
    this.updateGrid();
  }

  createAddColumnDialog() {
    const columnInput = ui.columnInput('Column', this.dataFrame, this.dataFrame.columns.byIndex(0));
    columnInput.input.style.width = '100px';

    const aggrCheckbox = ui.boolInput('Aggregate', true);
    const chartCheckbox = ui.boolInput('Chart', false);
    const checkboxesDiv = ui.divH([
      aggrCheckbox.root,
      chartCheckbox.root,
    ]);

    const aggrTypesChoice = ui.choiceInput('', Object.keys(AGGREGATE_TYPES)[0], Object.keys(AGGREGATE_TYPES));
    const chartTypesChoice = ui.choiceInput('', Object.keys(CHART_TYPES)[0], Object.keys(CHART_TYPES));
    const columnTypeDiv = ui.div();
    columnTypeDiv.append(aggrTypesChoice.root);

    this.registerCheckBoxChanged(aggrCheckbox, [chartCheckbox], aggrTypesChoice, columnTypeDiv);
    this.registerCheckBoxChanged(chartCheckbox, [aggrCheckbox], chartTypesChoice, columnTypeDiv);


    ui.dialog('Add column')
      .add(ui.narrowForm([
        columnInput,
        //@ts-ignore
        checkboxesDiv,
        //@ts-ignore
        columnTypeDiv,
      ]))
      .onOK(() => {
        if ([aggrCheckbox, chartCheckbox].every((it) => !it.value)) {
          grok.shell.error(`Chart/statistic hasn't been selected`);
          return;
        } else {
          const selectedType = aggrCheckbox.value ? AGGR_TYPE : CHART_TYPE;
          const selectedFunc = aggrCheckbox.value ? aggrTypesChoice.value : chartTypesChoice.value;
          if (!this.colExists(columnInput.value!.name, selectedType, selectedFunc!))
            this.colsToAnalyze.push({colName: columnInput.value!.name, type: selectedType, typeName: selectedFunc!});
          else
            grok.shell.warning('Column has already been selected');
        }
        this.updateGrid();
      })
      .show();
  }

  colExists(colName: string, type: string, typeName: string) {
    return this.colsToAnalyze
      .filter((it) => it.colName === colName && it.type === type && it.typeName === typeName).length > 0;
  }


  registerCheckBoxChanged(currentCheckbox: DG.InputBase,
    otherCheckboxes: DG.InputBase[], columnTypeChoice: DG.InputBase, colTypeDiv: HTMLDivElement) {
    currentCheckbox.root.onclick = () => {
      if (currentCheckbox.value) {
        otherCheckboxes.forEach((it) => it.value = false);
        columnTypeChoice.input.style.width = '100px';
        ui.empty(colTypeDiv);
        colTypeDiv.append(columnTypeChoice.root);
      } else
        ui.empty(colTypeDiv);
    };
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
    let grouppedDfQuery = this.dataFrame.groupBy(this.grouppingCols);

    for (const col of this.colsToAnalyze) {
      if (col.type === AGGR_TYPE)
        grouppedDfQuery = (AGGREGATE_TYPES as any)[col.typeName](grouppedDfQuery, col.colName);
      else
        grouppedDfQuery = (AGGREGATE_TYPES as any)['count'](grouppedDfQuery, col.typeName, `${col.colName}_${col.typeName}`);
    }

    this.grouppedDf = grouppedDfQuery.aggregate();
    const grid = this.grouppedDf.plot.grid();
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    for (const col of this.grouppingCols)
      grid.col(col)!.width = 150;

    for (const col of this.colsToAnalyze) {
      if (col.type === CHART_TYPE)
        grid.col(`${col.colName}_${col.typeName}`)!.cellType = 'html';
    }

    grid.onCellPrepare((gc) => {
      const colToAnalyze = this.colsToAnalyze.length ?
        this.colsToAnalyze.filter((it) => gc.gridColumn.name === `${it.colName}_${it.typeName}`) : [];
      if (gc.isTableCell && colToAnalyze.length) {
        const condition = this.createFilterCondition(gc.tableRowIndex!);
        let df: DG.DataFrame;
        if (!this.parentViewers[gc.gridColumn.name]) {
          df = this.createViewerDf(colToAnalyze[0].colName, gc.tableRowIndex!);
          this.parentViewers[gc.gridColumn.name] = DG.Viewer.fromType((CHART_TYPES as any)[colToAnalyze[0].typeName], df);
        }
        if (!this.viewersStorage[gc.gridColumn.name])
          this.viewersStorage[gc.gridColumn.name] = {};
        else {
          if (!this.viewersStorage[gc.gridColumn.name][gc.gridRow]) {
            df ??= this.createViewerDf(colToAnalyze[0].colName, gc.tableRowIndex!);
            const viewer = DG.Viewer.fromType((CHART_TYPES as any)[colToAnalyze[0].typeName], df);
            this.viewersStorage[gc.gridColumn.name][gc.gridRow] = viewer;
            //@ts-ignore
            viewer.copyViewersLook(this.parentViewers[gc.gridColumn.name]);
            //@ts-ignore
            this.parentViewers[gc.gridColumn.name].onDartPropertyChanged
              .subscribe(() => {
              //@ts-ignore
                viewer.copyViewersLook(this.parentViewers[gc.gridColumn.name]);
                df.rows.match(condition).filter();
              });
          }
        }
        gc.element = this.viewersStorage[gc.gridColumn.name][gc.gridRow].root;
        this.viewersStorage[gc.gridColumn.name][gc.gridRow].dataFrame.rows.match(condition).filter();
      }
    });

    grid.onCellClick.subscribe((gc) => {
      const colToAnalyze = this.colsToAnalyze.filter((it) => gc.gridColumn.name === `${it.colName}_${it.typeName}`);
      if (!gc.isTableCell && colToAnalyze.length)
        grok.shell.o = this.parentViewers[gc.gridColumn.name];
    });

    ui.empty(this.grouppedGridDiv);
    this.grouppedGridDiv.append(grid.root);
  }

  createFilterCondition(idx: number): { [key: string]: any } {
    const condition: { [key: string]: any } = {};
    for (const col of this.grouppingCols)
      condition[col] = this.grouppedDf!.get(col, idx);
    return condition;
  }

  createViewerDf(colToAnalyzeName: string, idx: number): DG.DataFrame {
    const colList = [];
    this.grouppingCols.forEach((col) => colList.push(this.dataFrame.col(col)!));
    colList.push(this.dataFrame.col(colToAnalyzeName)!);
    const df = DG.DataFrame.fromColumns(colList);
    df.columns.addNewString('group').init((j) => {
      for (const group of this.grouppingCols) {
        if (df.get(group, j) !== this.grouppedDf!.get(group, idx!))
          return 'other';
      }
      return 'current_group';
    });
    return df;
  }
}

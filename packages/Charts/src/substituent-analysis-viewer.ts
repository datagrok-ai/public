import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { Grid } from 'datagrok-api/dg';
import {tTest} from '@datagrok-libraries/statistics/src/tests';

const AGGR_TYPE = 'Aggregate';
const CHART_TYPE = 'Chart';
const STAT_TYPE = 'Statistics';

const COL_TYPES = {
  [AGGR_TYPE]: {
    'min' : (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder  => { return query.min(colName, resColName) },
    'max': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => { return query.max(colName, resColName) },
    'avg': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => { return query.avg(colName, resColName) },
    'count': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => { return query.count(resColName) },
  },
  [CHART_TYPE]: {
    histogram: {viewer: DG.VIEWER.HISTOGRAM, params: {split: 'group'}},
    barchart: {viewer: DG.VIEWER.BAR_CHART, params: {split: 'group', showCategorySelector: false, showValueAxis: false}},
    piechart: {viewer: DG.VIEWER.PIE_CHART, params: {category: 'group'}},
  },
  [STAT_TYPE]: {
    'T-test': {}
  }
}


interface IAnalyzedColumn {
  colName: string, // name of analyzed column
  type: string, // aggregation or chart
  typeName: string //name of exact aggr method or exact chart
}

export class SubstituentAnalysisViewer extends DG.JsViewer {
  initialized: boolean = false;
  name = 'substituent';
  groupByColumns: string[];
  analyzedColumns: IAnalyzedColumn[];
  totalColumns: string[] | undefined = undefined;
  grouppedDf: DG.DataFrame | undefined = undefined;
  parentViewers: {[key: string]: DG.Viewer};
  viewersStorage: {[key: string]: {[key: number]: DG.Viewer}} = {};
  grouppingColsDiv = ui.div();
  grouppedGridDiv = ui.div();
  mainView = ui.splitV([]);
  grid: DG.Grid | undefined = undefined;

  constructor() {
    super();
    this.groupByColumns = this.stringList('groupByColumns', undefined);
    this.analyzedColumns = this.addProperty('analyzedColumns', 'object', [], {'userEditable': false});
    this.parentViewers = this.addProperty('parentViewers', 'object', {}, {'userEditable': false});
  }

  init(): void {
    this.initialized = true;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async onTableAttached(): Promise<void> {
    this.init();
    this.totalColumns = this.dataFrame.columns.names().concat(['']);
    this.groupByColumns ??= [this.totalColumns[0]];
    this.updateColumnChoices(this.groupByColumns, 'Group by', this.grouppingColsDiv);
    this.mainView.append(ui.box(ui.panel([this.grouppingColsDiv], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    const addColToAnalyze = ui.icons.add(() => {this.createAddColumnDialog()}, 'Add column to analyze');
    this.mainView.append(ui.box(ui.panel([addColToAnalyze], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    this.mainView.append(this.grouppedGridDiv);
    this.root.append(this.mainView);
    this.updateGrid();
  }

  onPropertyChanged(p: DG.Property) {
    if (p?.name === 'groupByColumns') {
      this.updateColumnChoices(this.groupByColumns, 'Group by', this.grouppingColsDiv);
      this.updateGrid();
    }      
  }

  createAddColumnDialog(){
    const columnInput = ui.columnInput('Column', this.dataFrame, this.dataFrame.columns.byIndex(0));
    columnInput.input.style.width = '100px';

    const colTypeInput = ui.choiceInput('Column Type', Object.keys(COL_TYPES)[0], Object.keys(COL_TYPES));

    const aggrTypesChoice = ui.choiceInput('', Object.keys(COL_TYPES[AGGR_TYPE])[0], Object.keys(COL_TYPES[AGGR_TYPE]));
    const chartTypesChoice = ui.choiceInput('', Object.keys(COL_TYPES[CHART_TYPE])[0], Object.keys(COL_TYPES[CHART_TYPE]));
    const statTypesChoice = ui.choiceInput('', Object.keys(COL_TYPES[STAT_TYPE])[0], Object.keys(COL_TYPES[STAT_TYPE]));
    const columnTypeDiv = ui.div();
    columnTypeDiv.append(aggrTypesChoice.root);

    let currentColType = AGGR_TYPE;
    let currentFuncChoice = aggrTypesChoice;

    function updateColTypeDiv(choiceInput: DG.InputBase){
      currentFuncChoice = choiceInput;
      columnTypeDiv.append(choiceInput.root);
    }

    colTypeInput.onChanged(() => {
      currentColType = colTypeInput.value!;
      ui.empty(columnTypeDiv);
      switch (colTypeInput.value) {
        case AGGR_TYPE: {
          updateColTypeDiv(aggrTypesChoice);
          break;
        }
        case CHART_TYPE: {
          updateColTypeDiv(chartTypesChoice);
          break;
        }
        case STAT_TYPE: {
          updateColTypeDiv(statTypesChoice);
          break;
        }
      }
    })

    ui.dialog('Add column')
      .add(ui.narrowForm([
        columnInput,
        colTypeInput,
        //@ts-ignore
        columnTypeDiv,
      ]))
      .onOK(() => {
        this.checkColExistsAndadd(columnInput.value!.name, currentColType, currentFuncChoice.value!);
      })
      .show();
  }

  checkColExistsAndadd(colName: string, type: string, typeName: string){
    if(this.analyzedColumns.filter((it) => it.colName === colName && it.type === type && it.typeName === typeName).length)
      grok.shell.warning('Column already exists');
    else {
      this.addCalculatedCol(colName, type, typeName);
    }
  }

  addCalculatedCol(colName: string, type: string, typeName: string) {
    let grouppedDfQuery = this.dataFrame.groupBy(this.groupByColumns);
    if (type === AGGR_TYPE) {
      grouppedDfQuery = (COL_TYPES[AGGR_TYPE] as any)[typeName](grouppedDfQuery, colName);
    }
    else if (type === CHART_TYPE){
      grouppedDfQuery = (COL_TYPES[AGGR_TYPE] as any)['count'](grouppedDfQuery, colName, `${colName}_${typeName}`);
    } else {
      grouppedDfQuery = (COL_TYPES[AGGR_TYPE] as any)['min'](grouppedDfQuery, colName, `pValue(${colName})`);
    }
    try {
      const grouppedDf = grouppedDfQuery.aggregate();
      const colToReturnName = grouppedDf.columns.names().filter(it => !this.groupByColumns.includes(it))[0];
      this.grouppedDf?.columns.add(grouppedDf.col(colToReturnName)!);
      this.analyzedColumns.push({colName: colName, type: type, typeName: typeName});
      this.updateGridColsWidth();
    } catch {
      grok.shell.error(`Incorrect column type`);
    }
  }

  updateColumnChoices(selectedCols: string[], choicesInputsName: string, choicesDiv: HTMLDivElement) {
    const colsChoicesDiv = ui.divH([]);
    colsChoicesDiv.append(ui.divText(
      choicesInputsName, {style: {fontWeight: 'bold', paddingTop: '10px', paddingRight: '10px'}}));
    for (let i = 0; i < selectedCols.length + 1; i++) {
      const selectedValue = selectedCols.length === i ? '' : selectedCols[i];
      const groupChoiceInput = ui.choiceInput('', selectedValue, this.totalColumns!);
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
    let grouppedDfQuery = this.dataFrame.groupBy(this.groupByColumns);

    for (const col of this.analyzedColumns) {
      if (col.type === AGGR_TYPE) {
        grouppedDfQuery = (COL_TYPES[AGGR_TYPE] as any)[col.typeName](grouppedDfQuery, col.colName);
      }
      else if (col.type === CHART_TYPE) {
        grouppedDfQuery = (COL_TYPES[AGGR_TYPE] as any)['count'](grouppedDfQuery, col.typeName, `${col.colName}_${col.typeName}`);
      } else {
        grouppedDfQuery = (COL_TYPES[AGGR_TYPE] as any)['count'](grouppedDfQuery, col.typeName, `pValue(${col.colName})`);
      }
    }

    this.grouppedDf = grouppedDfQuery.aggregate();
    this.grid = this.grouppedDf.plot.grid();
    this.grid.root.style.width = '100%';
    this.grid.root.style.height = '100%';

    for (const col of this.analyzedColumns) {
      if (col.type === CHART_TYPE) {
        this.grid.col(`${col.colName}_${col.typeName}`)!.cellType = 'html';
      }
    }

    this.grid.onCellPrepare((gc) => {
      const chartCol = this.analyzedColumns.length ?
        this.analyzedColumns.filter(it => gc.gridColumn.name === `${it.colName}_${it.typeName}`) : [];
      const tTestCol = this.analyzedColumns.length ?
        this.analyzedColumns.filter(it => gc.gridColumn.name === `pValue(${it.colName})`) : [];
      if (gc.isTableCell && (chartCol.length || tTestCol.length)) {
        if (chartCol.length) {
          let df: DG.DataFrame;
          if (!this.parentViewers[gc.gridColumn.name]) {
            df = this.createViewerDf(chartCol[0].colName, gc.tableRowIndex!);
            const parentViewer = DG.Viewer.fromType((COL_TYPES[CHART_TYPE] as any)[chartCol[0].typeName].viewer, df);
            //@ts-ignore
            parentViewer.toCompactLook();
            parentViewer.setOptions((COL_TYPES[CHART_TYPE] as any)[chartCol[0].typeName].params);
            this.parentViewers[gc.gridColumn.name] = parentViewer;         
          }
          if (!this.viewersStorage[gc.gridColumn.name])
            this.viewersStorage[gc.gridColumn.name] = {};
          else {
            if (!this.viewersStorage[gc.gridColumn.name][gc.gridRow]) {
              df ??= this.createViewerDf(chartCol[0].colName, gc.tableRowIndex!);
              const viewer = DG.Viewer.fromType((COL_TYPES[CHART_TYPE] as any)[chartCol[0].typeName].viewer, df);
              this.viewersStorage[gc.gridColumn.name][gc.gridRow] = viewer;
              //@ts-ignore
              viewer.copyViewersLook(this.parentViewers[gc.gridColumn.name]);
              //@ts-ignore
              this.parentViewers[gc.gridColumn.name].onDartPropertyChanged
              .subscribe(() => {
                //@ts-ignore
                viewer.copyViewersLook(this.parentViewers[gc.gridColumn.name]);
              });
            }
          }
          gc.element = this.viewersStorage[gc.gridColumn.name][gc.gridRow].root;
        } else {
          const pValue = this.performTTest(tTestCol[0].colName, gc.tableRowIndex!);
          gc.cell.value = pValue;
        }
      }
    });

    this.grid.onCellClick.subscribe((gc) => {
      const colToAnalyze = this.analyzedColumns.filter(it => gc.gridColumn.name === `${it.colName}_${it.typeName}`);
      if (!gc.isTableCell && colToAnalyze.length)
        grok.shell.o = this.parentViewers[gc.gridColumn.name];
    });

    ui.empty(this.grouppedGridDiv);
    this.grouppedGridDiv.append(this.grid.root);
    this.updateGridColsWidth();  
  }

  updateGridColsWidth() {
    setTimeout(() => {
      this.groupByColumns.forEach(col => {
        this.grid!.col(col)!.width = 100;
      });
      this.analyzedColumns.forEach(col => {
        if (col.type === AGGR_TYPE) {
          col.typeName === 'count' ? this.grid!.col('count')!.width = 50 : this.grid!.col(`${col.typeName}(${col.colName})`)!.width = 50;
        }
        else if (col.type === CHART_TYPE){
          this.grid!.col(`${col.colName}_${col.typeName}`)!.width = 100;
        }
      })
      this.grid!.setOptions({rowHeight: 70});
    }, 200);
  }

  createFilterCondition(idx: number): { [key: string]: any } {
    const condition: { [key: string]: any } = {};
    for (const col of this.groupByColumns)
      condition[col] = this.grouppedDf!.get(col, idx);
    return condition;
  }

  performTTest(colToAnalyzeName: string, idx: number): number {
    const df = this.extractGroupAndAnalyzedColFromInitialDf(colToAnalyzeName);
    const currentGroup = [];
    const otherGroup = [];
    const rawData = df.col(colToAnalyzeName)!.getRawData();
    const rowCount = df.rowCount;
    for (let i = 0; i < rowCount; i++) {
      for (const group of this.groupByColumns) {
        if (df.get(group, i) !== this.grouppedDf!.get(group, idx!))
          currentGroup.push(rawData[i]);
        else
          otherGroup.push(rawData[i])
      }
    };
    const res = tTest(currentGroup, otherGroup);
    return res['p-value'];
  }

  extractGroupAndAnalyzedColFromInitialDf(colToAnalyzeName: string): DG.DataFrame {
    const colList = [];
    this.groupByColumns.forEach((col) => colList.push(this.dataFrame.col(col)!));
    colList.push(this.dataFrame.col(colToAnalyzeName)!);
    const df = DG.DataFrame.fromColumns(colList);
    return df;
  }

  createViewerDf(colToAnalyzeName: string, idx: number): DG.DataFrame {
    const df = this.extractGroupAndAnalyzedColFromInitialDf(colToAnalyzeName);
    df.columns.addNewString('group').init((j) => {
      for (const group of this.groupByColumns) {
        if (df.get(group, j) !== this.grouppedDf!.get(group, idx!))
          return 'other';
      }
      return 'current_group';
    });
    return df;
  }

}

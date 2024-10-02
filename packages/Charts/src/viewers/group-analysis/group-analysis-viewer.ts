import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {tTest} from '@datagrok-libraries/statistics/src/tests';

const AGGR_TYPE = 'Aggregate';
const CHART_TYPE = 'Chart';
const STAT_TYPE = 'Statistics';
// const ROW_HEIGHT = 70;

const COL_TYPES = {
  [AGGR_TYPE]: {
    'min': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.min(colName, resColName);
    },
    'max': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.max(colName, resColName);
    },
    'avg': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.avg(colName, resColName);
    },
    'count': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.count(resColName);
    },
    'med': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.med(colName, resColName);
    },
    'first': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.first(colName, resColName);
    },
    'key': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.key(colName, resColName);
    },
    'missingValueCount': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.missingValueCount(colName, resColName);
    },
    'pivot': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.pivot(colName, resColName);
    },
    'q1': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.q1(colName, resColName);
    },
    'q2': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.q2(colName, resColName);
    },
    'q3': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.q3(colName, resColName);
    },
    'stdev': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.stdev(colName, resColName);
    },
    'sum': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.sum(colName, resColName);
    },
    'uniqueCount': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.uniqueCount(colName, resColName);
    },
    'valueCount': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.valueCount(colName, resColName);
    },
    'variance': (query: DG.GroupByBuilder, colName: string, resColName?: string): DG.GroupByBuilder => {
      return query.variance(colName, resColName);
    },
  },
  [CHART_TYPE]: {
    histogram: {viewer: DG.VIEWER.HISTOGRAM, params: {split: 'group', marginTop: 5, marginBottom: 5}},
    barchart: {viewer: DG.VIEWER.BAR_CHART, params: {split: 'group', showCategorySelector: false,
      showValueAxis: false}},
    piechart: {viewer: DG.VIEWER.PIE_CHART, params: {category: 'group'}},
  },
  [STAT_TYPE]: {
    'T-test': {},
  },
};

const COL_WIDTH: {[key: string]: number}= {
  [AGGR_TYPE]: 50,
  [CHART_TYPE]: 100,
  [STAT_TYPE]: 70,
};


interface IAnalyzedColumn {
  colName: string, // name of analyzed column
  type: string, // aggregation or chart
  typeName: string //name of exact aggr method or exact chart
  gridColName?: string
}


@grok.decorators.viewer({
  name: 'Group Analysis',
  description: 'Creates a group analysis viewer',
  icon: 'icons/groupanalysis-viewer.svg',
})
export class GroupAnalysisViewer extends DG.JsViewer {
  initialized: boolean = false;
  name = 'group';
  groupByColumnNames: string[];
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
    this.groupByColumnNames = this.addProperty('groupByColumnNames', DG.TYPE.COLUMN_LIST);
    this.analyzedColumns = [];
    this.parentViewers = {};
  }

  init(): void {
    this.initialized = true;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async onTableAttached(): Promise<void> {
    this.init();
    this.initChartEventListeners();
    this.totalColumns = this.dataFrame.columns.names().concat(['']);
    this.groupByColumnNames ??= [this.totalColumns[0]];
    this.updateColumnChoices(this.groupByColumnNames, 'Group by', this.grouppingColsDiv);
    this.mainView.append(
      ui.box(ui.panel([this.grouppingColsDiv], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    const addColToAnalyze = ui.icons.add(() => {this.createAddColumnDialog();}, 'Add column to analyze');
    this.mainView.append(ui.box(ui.panel([addColToAnalyze], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    this.mainView.append(this.grouppedGridDiv);
    this.root.append(this.mainView);
    this.updateGrid();
  }

  initChartEventListeners() {
    this.dataFrame.onRowsFiltered.subscribe((_) => {
      this.updateGrid();
    });
  }

  onPropertyChanged(p: DG.Property) {
    if (p?.name === 'groupByColumnNames') {
      this.updateColumnChoices(this.groupByColumnNames, 'Group by', this.grouppingColsDiv);
      this.updateGrid();
    }
  }

  _testColumns() {
    return this.dataFrame.columns.length >= 1;
  }

  createAddColumnDialog() {
    const columnInput = ui.input.column('Column', {table: this.dataFrame, value: this.dataFrame.columns.byIndex(0)});
    columnInput.input.style.width = '100px';

    const colTypeInput = ui.input.choice('Column type', {value: Object.keys(COL_TYPES)[0], items: Object.keys(COL_TYPES)});

    const aggrTypesChoice = ui.input.choice('Function', {value: Object.keys(COL_TYPES[AGGR_TYPE])[0],
      items: Object.keys(COL_TYPES[AGGR_TYPE])});
    const chartTypesChoice = ui.input.choice('Chart', {value: Object.keys(COL_TYPES[CHART_TYPE])[0],
      items: Object.keys(COL_TYPES[CHART_TYPE])});
    const statTypesChoice = ui.input.choice('Statistic', {value: Object.keys(COL_TYPES[STAT_TYPE])[0],
      items: Object.keys(COL_TYPES[STAT_TYPE])});
    const columnTypeDiv = ui.div();
    columnTypeDiv.append(aggrTypesChoice.root);

    let currentColType = AGGR_TYPE;
    let currentFuncChoice = aggrTypesChoice;

    function updateColTypeDiv(choiceInput: DG.InputBase) {
      // @ts-ignore
      currentFuncChoice = choiceInput;
      columnTypeDiv.append(choiceInput.root);
    }

    colTypeInput.onChanged.subscribe((value) => {
      currentColType = value!;
      ui.empty(columnTypeDiv);
      switch (value) {
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
    });

    ui.dialog('Add column')
      .add(ui.form([
        columnInput,
        colTypeInput,
        //@ts-ignore
        columnTypeDiv,
      ]))
      .onOK(() => {
        this.checkColExistsAndAdd(columnInput.value!.name, currentColType, currentFuncChoice.value!);
      })
      .show();
  }

  checkColExistsAndAdd(colName: string, type: string, typeName: string) {
    if (this.analyzedColumns.filter((it) => it.colName === colName &&
      it.type === type && it.typeName === typeName).length)
      grok.shell.warning('Column already exists');
    else {
      const columnToAdd: IAnalyzedColumn = {colName: colName, type: type, typeName: typeName};
      const col = this.getCalculatedCol(columnToAdd);
      if (col) {
        columnToAdd.gridColName = col.name;
        this.analyzedColumns.push(columnToAdd);
        this.addColumnAndSetWidth(col, type);
      }
    }
  }

  addColumnAndSetWidth(col: DG.Column, type: string) {
    this.grouppedDf?.columns.add(col);
    setTimeout(() => {this.grid!.col(col!.name)!.width = COL_WIDTH[type];}, 50);
  }

  getCalculatedCol(column: IAnalyzedColumn): DG.Column | null {
    let col: DG.Column | null = null;
    if (column.type === AGGR_TYPE) {
      const res = this.getAggregateCols([column]);
      if (res) col = res[0];
    } else if (column.type === CHART_TYPE)
      col = this.getChartCol(column, this.grouppedDf!.rowCount);
    else
      col = this.getStatisticsCol(column, this.grouppedDf!.rowCount);
    return col;
  }

  getAggregateCols(columnList: IAnalyzedColumn[]): DG.Column[] {
    let grouppedDfQuery = this.dataFrame.groupBy(this.groupByColumnNames).whereRowMask(this.dataFrame.filter);
    for (const col of columnList)
      grouppedDfQuery = (COL_TYPES[AGGR_TYPE] as any)[col.typeName](grouppedDfQuery, col.colName);
    try {
      const grouppedDf = grouppedDfQuery.aggregate();
      const colsToReturnNames = grouppedDf.columns.names().filter((it) => !this.groupByColumnNames.includes(it));
      const colList: DG.Column[] = colsToReturnNames.map((colName) => grouppedDf.col(colName)!);
      return colList;
    } catch {
      grok.shell.error(`Incorrect column type`);
      return [];
    }
  }

  getChartCol(column: IAnalyzedColumn, length: number) {
    return DG.Column.string(`${column.colName}_${column.typeName}`, length);
  }

  getStatisticsCol(column: IAnalyzedColumn, length: number) {
    const col = DG.Column.float(`pValue(${column.colName}`, length).init((i) => this.performTTest(column.colName, i));
    col.meta.format = '#.0000000';
    return col;
  }

  updateColumnChoices(selectedCols: string[], choicesInputsName: string, choicesDiv: HTMLDivElement) {
    const colsChoicesDiv = ui.divH([]);
    colsChoicesDiv.append(ui.divText(
      choicesInputsName, {style: {fontWeight: 'bold', paddingTop: '10px', paddingRight: '10px'}}));
    for (let i = 0; i < selectedCols.length + 1; i++) {
      const selectedValue = selectedCols.length === i ? '' : selectedCols[i];
      const groupChoiceInput = ui.input.choice('', {value: selectedValue, items: this.totalColumns!});
      ui.tooltip.bind(groupChoiceInput.root, () => selectedValue);
      groupChoiceInput.onChanged.subscribe((value) => {
        if (value === '')
          selectedCols.splice(i, 1);
        else {
          selectedCols.length === i ? selectedCols.push(value!) :
            selectedCols[i] = value!;
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
    this.grouppedDf = this.dataFrame.groupBy(this.groupByColumnNames).whereRowMask(this.dataFrame.filter).aggregate();
    this.grid = this.grouppedDf.plot.grid();
    const aggregateCols = this.analyzedColumns.filter((it) => it.type === AGGR_TYPE);
    const chartAndStatCols = this.analyzedColumns.filter((it) => it.type !== AGGR_TYPE);

    const calculatedAggrCols = this.getAggregateCols(aggregateCols);
    for (const col of calculatedAggrCols)
      this.addColumnAndSetWidth(col, AGGR_TYPE);

    for (const col of chartAndStatCols) {
      const calculatedCol = this.getCalculatedCol(col);
      if (calculatedCol)
        this.addColumnAndSetWidth(calculatedCol, col.type);
    }

    this.grid.columns.setOrder(this.groupByColumnNames.concat(this.analyzedColumns.map((it) => it.gridColName!)));
    this.grid.root.style.width = '100%';
    this.grid.root.style.height = '100%';

    this.grid.onCellPrepare((gc) => {
      const chartCol = this.analyzedColumns.length ?
        this.analyzedColumns.filter((it) => gc.gridColumn.name === `${it.colName}_${it.typeName}`) : [];
      if (gc.isTableCell && (chartCol.length)) {
        if (chartCol.length) {
          let df: DG.DataFrame;
          if (!this.parentViewers[gc.gridColumn.name]) {
            df = this.createViewerDf(chartCol[0].colName, gc.tableRowIndex!);
            const parentViewer = DG.Viewer.fromType((COL_TYPES[CHART_TYPE] as any)[chartCol[0].typeName].viewer, df);
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
              viewer.copyViewersLook(this.parentViewers[gc.gridColumn.name]);
              this.parentViewers[gc.gridColumn.name].onDartPropertyChanged
                .subscribe(() => {
                  viewer.copyViewersLook(this.parentViewers[gc.gridColumn.name]);
                });
            }
          }
          if (!(Object.keys(this.viewersStorage[gc.gridColumn.name]).length === 0))
            gc.element = this.viewersStorage[gc.gridColumn.name][gc.gridRow].root;
        }
      }
    });

    this.grid.onCellClick.subscribe((gc) => {
      const colToAnalyze = this.analyzedColumns.filter((it) => gc.gridColumn.name === `${it.colName}_${it.typeName}`);
      if (!gc.isTableCell && colToAnalyze.length)
        grok.shell.o = this.parentViewers[gc.gridColumn.name];
    });

    ui.empty(this.grouppedGridDiv);
    this.grouppedGridDiv.append(this.grid.root);
  }


  createFilterCondition(idx: number): { [key: string]: any } {
    const condition: { [key: string]: any } = {};
    for (const col of this.groupByColumnNames)
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
      let otherGroupFound = false;
      for (const group of this.groupByColumnNames) {
        if (df.get(group, i) !== this.grouppedDf!.get(group, idx!))
          otherGroupFound = true;
      }
      otherGroupFound ? otherGroup.push(rawData[i]) : currentGroup.push(rawData[i]);
    };
    const res = tTest(currentGroup, otherGroup);
    return res['p-value'];
  }

  extractGroupAndAnalyzedColFromInitialDf(colToAnalyzeName: string): DG.DataFrame {
    const colList = [];
    const filteredDf = this.dataFrame.groupBy(this.dataFrame.columns.names()).
      whereRowMask(this.dataFrame.filter).aggregate();
    this.groupByColumnNames.forEach((col) => colList.push(filteredDf.col(col)!));
    colList.push(filteredDf.col(colToAnalyzeName)!);
    const df = DG.DataFrame.fromColumns(colList);
    return df;
  }

  createViewerDf(colToAnalyzeName: string, idx: number): DG.DataFrame {
    const df = this.extractGroupAndAnalyzedColFromInitialDf(colToAnalyzeName);
    df.columns.addNewString('group').init((j) => {
      for (const group of this.groupByColumnNames) {
        if (df.get(group, j) !== this.grouppedDf!.get(group, idx!))
          return 'other';
      }
      return 'current_group';
    });
    return df;
  }
}
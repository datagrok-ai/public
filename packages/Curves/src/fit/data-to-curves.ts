/* eslint-disable max-len */

import {FIT_FUNCTION_4PL_REGRESSION, FitChartData, FitMarkerType, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import wu from 'wu';

export function dataToCurvesUI() {
  const tv = grok.shell.tv;
  if (!tv || !tv.dataFrame) {
    grok.shell.warning('No open tableview');
    return;
  }

  const df = tv.dataFrame;

  const assayColInput = ui.input.column('Assay Column', {table: df, nullable: true, filter: (c) => c.isCategorical, tooltipText: 'Column with assay names'});
  const concentrationColInput = ui.input.column('Concentration Column', {table: df, nullable: false, filter: (c) => c.isNumerical, tooltipText: 'Column with concentration values'});
  const readoutColInput = ui.input.column('Readout Column', {table: df, nullable: false, filter: (c) => c.isNumerical, tooltipText: 'Column with readout values'});
  const excludeOutliersInput = ui.input.column('Outliers Column', {table: df, nullable: true, filter: (c) => (c.isCategorical || c.type === DG.TYPE.BOOL) && c.categories.length === 2, tooltipText: 'Column with outliers to exclude from the analysis'});
  const compoundIDColInput = ui.input.column('Compound ID Column', {table: df, nullable: true, tooltipText: 'Column with compound IDs'});
  const batchIDColInput = ui.input.column('Batch ID Column', {table: df, nullable: false, tooltipText: 'Column with batch IDs. this column along with the Assay column will be used to group the data into curves'});
  const runIDColInput = ui.input.column('Run ID Column', {table: df, nullable: true, tooltipText: 'Column with run IDs. this column along with the Assay column will be used to group the data into curves'});

  assayColInput.value = wu(df.columns.categorical).find((c) => c.name.toLowerCase().includes('assay')) ?? null;
  concentrationColInput.value = df.columns.firstWhere((c) => c.isNumerical && c.type !== DG.COLUMN_TYPE.DATE_TIME && c.name.toLowerCase().includes('concentration')) ?? null;
  readoutColInput.value = df.columns.firstWhere((c) => c.isNumerical && c.type !== DG.COLUMN_TYPE.DATE_TIME && (c.name.toLowerCase().includes('inhibition'))) ?? null;
  excludeOutliersInput.value = df.columns.firstWhere((c) => (c.isCategorical || c.type === DG.TYPE.BOOL) && c.categories.length === 2 && (c.name.toLowerCase().includes('outlier') || (c.name.toLowerCase().includes('include') && c.name.toLowerCase().includes('exclude')))) ?? null;
  compoundIDColInput.value = df.columns.firstWhere((c) => c.isCategorical && c.name.toLowerCase().includes('compound') && c.name.toLowerCase().includes('id')) ?? null;
  batchIDColInput.value = df.columns.firstWhere((c) => c.name.toLowerCase().includes('batch') && c.name.toLowerCase().includes('id')) ?? null;
  runIDColInput.value = df.columns.firstWhere((c) => (c.name.toLowerCase().includes('run') || c.name.toLowerCase().includes('set')) && c.name.toLowerCase().includes('id')) ?? null;

  const onOK = async () => {
    if (!concentrationColInput.value || !readoutColInput.value || !batchIDColInput.value) {
      grok.shell.warning('Please fill all required fields');
      return;
    }

    const resultObj: {[assay: string]: {[batchID: string]: {[runID: string]: {x: number[], y: number[], outliers: boolean[], compoundId: any}}}} = {};

    const assayCol = assayColInput.value;
    const concentrationCol = concentrationColInput.value!;
    const readoutCol = readoutColInput.value!;
    const excludeOutliersCol = excludeOutliersInput.value;
    const compoundIDCol = compoundIDColInput.value;
    const batchIDCol = batchIDColInput.value!;
    const runIDCol = runIDColInput.value;

    const getAssay = assayCol ? (row: number) => assayCol.get(row) : () => 'All';
    const otlierIndexes = excludeOutliersCol ? excludeOutliersCol.getRawData() : null;
    const isOutlier = excludeOutliersCol ? (excludeOutliersCol.type == DG.COLUMN_TYPE.BOOL ? (row: number) => !!otlierIndexes![row] : (row: number) => !otlierIndexes![row] ) : (_row: number) => false;
    const getConcentration = (row: number) => concentrationCol.get(row);
    const getReadout = (row: number) => readoutCol.get(row);
    const getBatchID = (row: number) => batchIDCol.get(row);
    const getCompoundID = compoundIDCol ? (row: number) => compoundIDCol.get(row) : () => null;
    const getRunID = runIDCol ? (row: number) => runIDCol.get(row) : () => null;
    const markerTypes: FitMarkerType[] = ['circle', 'triangle bottom', 'star', 'cross border', 'diamond', 'square', 'triangle left', 'triangle right', 'triangle top', 'asterisk'];
    let markerTypeIndex = 0; // will keep track of the markers
    const runIdToMarkerType: {[runID: string]: FitMarkerType} = {};
    const getMarkerType = (runID: string) => {
      if (runID in runIdToMarkerType) { return runIdToMarkerType[runID]; } else {
        const markerType = markerTypes[markerTypeIndex];
        runIdToMarkerType[runID] = markerType;
        markerTypeIndex = (markerTypeIndex + 1) % markerTypes.length;
        return markerType;
      }
    };

    for (let i = 0; i < df.rowCount; i++) {
      const assay = getAssay(i);
      const batchID = getBatchID(i);

      const concentration = getConcentration(i);
      const readout = getReadout(i);
      const outlier = isOutlier(i);
      if (!resultObj[assay])
        resultObj[assay] = {};

      const assayObj = resultObj[assay];
      if (!assayObj[batchID])
        assayObj[batchID] = {};

      const assayBatchObj = assayObj[batchID];
      // {x: [], y: [], outliers: [], compoundId: compoundID};
      // use batchID as a key if runID is not provided
      const runID = getRunID(i) ?? batchID;
      if (!assayBatchObj[runID]) {
        const compoundID = getCompoundID(i);
        assayBatchObj[runID] = {x: [], y: [], outliers: [], compoundId: compoundID};
      }
      const assayBatchRunObj = assayBatchObj[runID];
      assayBatchRunObj.x.push(concentration);
      assayBatchRunObj.y.push(readout);
      assayBatchRunObj.outliers.push(outlier);
    }

    // create fit series

    const assays: string[] = [];
    const batchIDs: string[] = [];
    const compoundIDs: string[] = [];
    const seriesData: Omit<FitChartData, 'seriesOptions'>[] = [];

    for (const assay in resultObj) {
      for (const batchID in resultObj[assay]) {
        const assayBatchObj = resultObj[assay][batchID];
        for (const runID in assayBatchObj) {
          assays.push(assay);
          batchIDs.push(batchID);
          const compoundID = Object.values(assayBatchObj)[0]?.compoundId;
          compoundIDs.push(compoundID);
          const x = assayBatchObj[runID].x;
          const y = assayBatchObj[runID].y;
          const outliers = assayBatchObj[runID].outliers;
          const markerType = getMarkerType(runID);
          const s: IFitSeries = {
            fit: undefined, fitFunction: FIT_FUNCTION_4PL_REGRESSION, clickToToggle: true, droplines: ['IC50'], name: `${runID}`,
            points: x.map((xv, i) => ({x: xv, y: y[i], outlier: outliers[i], marker: markerType, size: 8})).sort((a, b) => a.x - b.x),
          };
          const fitData: Omit<FitChartData, 'seriesOptions'> = {
            chartOptions: {
              xAxisName: concentrationCol.name,
              yAxisName: readoutCol.name,
              logX: true,
              title: `${assay} - ${batchID} - ${runID}`,
            },
            series: [s],
          };
          seriesData.push(fitData);
        }
      }
    }

    const assaysColumn = DG.Column.fromStrings('Assay', assays);
    const batchIDsColumn = DG.Column.fromStrings('Batch ID', batchIDs);
    const compoundIDsColumn = DG.Column.fromStrings('Compound ID', compoundIDs);
    const seriesColumn = DG.Column.fromStrings('Series', seriesData.map((s) => JSON.stringify(s)));
    const allResultsColumn = DG.Column.fromStrings('All Results', seriesData.map((s) => s.series[0]?.points.map((p) => `(${p.x}, ${p.y})`).join(', ')));
    const runIdColumn = DG.Column.fromStrings('Run ID', seriesData.map((s) => s.series[0]!.name ?? ''));
    seriesColumn.semType = 'fit';
    seriesColumn.setTag('cell.renderer', 'fit');

    const resDF = DG.DataFrame.fromColumns([seriesColumn, assaysColumn, batchIDsColumn, compoundIDsColumn, runIdColumn, allResultsColumn]);
    resDF.name = 'Fitted Curves';

    const actualStatNames = {
      'interceptX': 'IC50',
      'top': 'Max',
      'bottom': 'Min',
      'slope': 'Hill',
      'auc': 'AUC',
    };

    Object.entries(actualStatNames).forEach(([statName, alias]) => {
      const params = {df: resDF, colName: seriesColumn.name, propName: statName, seriesName: 'series 0', seriesNumber: 0, newColName: alias};
      DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(params).callSync({processed: false});
    });
    if (actualStatNames['interceptX'])
      resDF.col(actualStatNames['interceptX']) && (resDF.col(actualStatNames['interceptX'])!.meta.format = 'scientific');

    return resDF;
  };

  ui.dialog('Data to Curves')
    .add(assayColInput)
    .add(batchIDColInput)
    .add(runIDColInput)
    .add(concentrationColInput)
    .add(readoutColInput)
    .add(compoundIDColInput)
    .add(excludeOutliersInput)
    .onOK(async () => {
      const df = await onOK();
      if (df) {
        const tv = grok.shell.addTableView(df);
        const trellis = tv.trellisPlot({xColumnNames: [], yColumnNames: ['Batch ID'], viewerType: 'MultiCurveViewer',
          showControlPanel: false, showXLabels: false, showYLabels: false, showXSelectors: false, showYSelectors: true, packCategories: false,
          onClick: 'Select'
        });
        tv.dockManager.dock(trellis, DG.DOCK_TYPE.TOP, null, 'Fitted Curves', 0.5);
        const selectionGrid = tv.addViewer(DG.VIEWER.GRID, {rowSource: 'Selected', selectedRowsColor: DG.Color.white});
        tv.dockManager.dock(selectionGrid, DG.DOCK_TYPE.FILL, tv.dockManager.findNode(tv.grid.root), 'Selected Curves');
        selectionGrid.props.title = 'Selected Curves';
        tv.dockManager.dock(tv.grid, DG.DOCK_TYPE.FILL, tv.dockManager.findNode(selectionGrid.root));
        tv.grid.props.title = 'Fitted Curves';

        // add the pivot table

        const pivot = tv.addViewer(DG.VIEWER.PIVOT_TABLE, {
          rowSource: 'Selected', pivotColumnNames: [],
          groupByColumnNames: ['Assay'], aggregateColumnNames: ['IC50', 'AUC', 'Hill', 'Min', 'Max'],
          ggregateAggTypes: ['geomean', 'avg', 'geomean', 'avg', 'avg'], showHeader: false
        });
        pivot.props.title = 'Seleted Statistics';
        tv.dockManager.dock(pivot, DG.DOCK_TYPE.TOP, tv.dockManager.findNode(trellis.root), 'Selected Statistics', 0.1);
      }
    })
    .show()
    .history(() => {
      const obj = {
        assayCol: assayColInput.value?.name,
        batchIDCol: batchIDColInput.value?.name,
        concentrationCol: concentrationColInput.value?.name,
        readoutCol: readoutColInput.value?.name,
        compoundIDCol: compoundIDColInput.value?.name,
        excludeOutliersCol: excludeOutliersInput.value?.name,
        runIDCol: runIDColInput.value?.name,
      };
      return Object.fromEntries(Object.entries(obj).filter(([_, v]) => !!v));
    }, (v) => {
      assayColInput.value = v.assayCol ? df.columns.byName(v.assayCol) : null;
      batchIDColInput.value = v.batchIDCol ? df.columns.byName(v.batchIDCol) : null;
      concentrationColInput.value = df.columns.byName(v.concentrationCol);
      readoutColInput.value = df.columns.byName(v.readoutCol);
      compoundIDColInput.value = v.compoundIDCol ? df.columns.byName(v.compoundIDCol) : null;
      excludeOutliersInput.value = v.excludeOutliersCol ? df.columns.byName(v.excludeOutliersCol) : null;
      runIDColInput.value = v.runIDCol ? df.columns.byName(v.runIDCol) : null;
    });
}

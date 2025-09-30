/* eslint-disable no-irregular-whitespace */
/* eslint-disable max-len */
import {AnalysisOptions, IAnalysisWidgetCoordinator, PlateWidget} from './plate-widget';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {mapFromRow} from './utils';
import {IPlateWellFilter, Plate, randomizeTableId} from './plate';
import {FIT_FUNCTION_4PL_REGRESSION, FitCurve, FitMarkerType, IFitChartData, IFitPoint} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '../fit/const';
import {FitCellOutlierToggleArgs, FitChartCellRenderer, setOutlier} from '../fit/fit-renderer';
import {savePlate, AnalysisProperty, createAnalysisRun, PlateProperty, getOrCreateProperty, saveAnalysisResult} from '../plates/plates-crud';
import {BaseAnalysisView} from './base-analysis-view';


export const drcAnalysisProperties: AnalysisProperty[] = [
  {name: 'Curve', type: DG.TYPE.STRING},
  {name: 'IC50', type: DG.TYPE.FLOAT},
  {name: 'Hill Slope', type: DG.TYPE.FLOAT},
  {name: 'R Squared', type: DG.TYPE.FLOAT},
  {name: 'Min', type: DG.TYPE.FLOAT},
  {name: 'Max', type: DG.TYPE.FLOAT},
  {name: 'AUC', type: DG.TYPE.FLOAT},
];

function _getOptionsForExcel(plate: Plate, defaultOptions: AnalysisOptions): AnalysisOptions {
  const finalOptions = {...defaultOptions};
  const aliases = {
    concentrationName: ['conc', 'dose', 'conc.', 'dilution', 'concentrations'],
    valueName: ['response', 'value', 'signal', 'raw data', 'raw signal', 'raw', 'response', 'values'],
    roleName: ['role', 'plate layout', 'plate positions', 'group', 'compound', 'well'],
  };

  Object.entries(aliases).forEach(([key, value]) => {
    const colName = finalOptions[key as keyof typeof finalOptions] as string;
    if (plate.data.columns.contains(colName))
      return;
    const alias = value.find((v) => plate.data.columns.contains(v));
    if (alias)
      (finalOptions[key as keyof typeof finalOptions] as string) = alias;
  });

  if (plate.data.columns.contains(finalOptions.roleName) && !finalOptions.controlColumns.every((cc) => plate.data.col(finalOptions.roleName)!.categories.includes(cc))) {
    const controlCols = plate.data.col(finalOptions.roleName)!.categories.filter((cat) => cat.toLowerCase().includes('control'));
    finalOptions.controlColumns = controlCols.slice(0, 2);
  }

  return finalOptions;
}

function _getOptionsForCsv(plate: Plate, defaultOptions: AnalysisOptions): AnalysisOptions {
  const finalOptions = {...defaultOptions};

  const csvColumnMapping = {
    roleName: 'SampleID', // Column identifying the compound/sample for grouping.
    concentrationName: 'Concentration', // Column with dose/concentration values.
    valueName: 'Activity', // Column with the measured response/activity.
  };

  Object.assign(finalOptions, csvColumnMapping);

  if (plate.data.columns.contains(finalOptions.roleName)) {
    const roleCol = plate.data.col(finalOptions.roleName)!;
    const controlCats = roleCol.categories.filter((c) => c && c.toLowerCase().includes('control'));
    if (controlCats.length >= 2)
      finalOptions.controlColumns = controlCats.slice(0, 2).sort();
  }

  return finalOptions;
}


export class PlateDrcAnalysis {
  static analysisView(plate: Plate, options?: Partial<AnalysisOptions>, layout: 'csv' | 'excel' = 'csv'): PlateWidget | null {
    const defaultOptions: AnalysisOptions = {
      roleName: 'layout', concentrationName: 'concentration', valueName: 'readout',
      normalize: true, controlColumns: ['High Control', 'Low Control'], autoFilterOutliers: true,
      categorizeFormula: '${rSquared} > 0.8 && ${Hill} > 0.25 && ${Max} > 80 && ${Max} < 120',
      statisticsColumns: ['Min', 'Max', 'rSquared', 'Hill', 'IC50', 'AUC'],
    };

    let actOptions: AnalysisOptions & {normalizedColName?: string};
    if (layout === 'csv')
      actOptions = {..._getOptionsForCsv(plate, defaultOptions), ...options};
    else
      actOptions = {..._getOptionsForExcel(plate, defaultOptions), ...options};

    if (!plate.data.columns.contains(actOptions.roleName) ||
        !plate.data.columns.contains(actOptions.concentrationName) ||
        !plate.data.columns.contains(actOptions.valueName))
      return null;


    const pw = PlateWidget.detailedView(plate);
    const gridRoot = pw.grid.root;
    const detailsRoot = pw.detailsDiv;
    pw.plateActionsDiv = ui.div();
    pw.plateDetailsDiv = ui.div();
    detailsRoot!.prepend(pw.plateDetailsDiv);
    detailsRoot!.append(pw.plateActionsDiv);

    if (!gridRoot || !detailsRoot)
      return pw;

    const container = ui.divH([gridRoot, detailsRoot]);
    pw.tabsContainer.appendChild(container);
    detailsRoot.style.flex = '0 0 min(250px, 25%)';
    gridRoot.style.flexGrow = '1';
    gridRoot.style.removeProperty('width');

    const drFilterOptions: IPlateWellFilter = {exclude: {[actOptions.roleName]: actOptions.controlColumns}};
    let normed = false;
    if (actOptions.normalize && actOptions.controlColumns.length === 2 && plate.data.col(actOptions.roleName)?.categories.includes(actOptions.controlColumns[0]) && plate.data.col(actOptions.roleName)?.categories.includes(actOptions.controlColumns[1])) {
      const [lStats, hStats] = actOptions.controlColumns.map((colName) => {
        return plate.getStatistics(actOptions.valueName, ['mean', 'std'], {match: {[actOptions.roleName]: colName}});
      }).sort((a, b) => a.mean - b.mean);

      if (hStats.mean !== lStats.mean) {
        defaultOptions.plateStatistics = {
          'Z Prime': (_plate) => 1 - (3 * (hStats.std + lStats.std) / Math.abs(hStats.mean - lStats.mean)),
          'Signal to background': (_plate) => hStats.mean / lStats.mean,
        };
        if (!actOptions.plateStatistics)
          actOptions.plateStatistics = defaultOptions.plateStatistics;
        actOptions.normalizedColName = plate.normalize(actOptions.valueName, (v) => (hStats.mean - v) / (hStats.mean - lStats.mean) * 100).name;
        normed = true;
        const plateStatMap: Record<string, string> = {};
        Object.entries(actOptions.plateStatistics).forEach(([statName, statFunc]) => {
          const value = statFunc(plate);
          plateStatMap[statName] = typeof value === 'number' ? DG.format(value, '#0.0000') : value;
        });
        const statTable = ui.tableFromMap(plateStatMap);
        pw.plateDetailsDiv!.appendChild(statTable);

        //         if (actOptions.autoFilterOutliers)
        //           plate.markOutliersWhere(actOptions.normalizedColName!, (v) => v > 106 || v < -6, drFilterOptions);
        if (actOptions.autoFilterOutliers)
          plate.markOutliersWhere(actOptions.normalizedColName!, (v) => v > 106 || v < -6, drFilterOptions, 'auto-qc');
      }
    }

    const series = plate.doseResponseSeries({...drFilterOptions, value: normed ? actOptions.normalizedColName! : actOptions.valueName, concentration: actOptions.concentrationName, groupBy: actOptions.roleName});
    const seriesVals = Object.entries(series);
    if (seriesVals.length === 0) return pw;

    const minMax = normed ? {minY: -10, maxY: 110} : {};
    const roleCol = DG.Column.string(actOptions.roleName, seriesVals.length);
    const curveCol = DG.Column.string('Curve', seriesVals.length);

    curveCol.init((i) => {
      const currentSeries = seriesVals[i][1];
      const seriesData = {
        points: currentSeries.points,
        name: currentSeries.name,
        fitFunction: FIT_FUNCTION_4PL_REGRESSION,
        showPoints: 'points',
        clickToToggle: true,
        droplines: ['IC50'],
      };

      return JSON.stringify({
        chartOptions: {
          xAxisName: actOptions.concentrationName,
          yAxisName: normed ? `norm(${actOptions.valueName})` : actOptions.valueName,
          logX: true,
          title: `${seriesVals[i][0]}`,
          ...minMax,
        },
        series: [seriesData],
      });
    });

    roleCol.init((i) => seriesVals[i][0]);

    const df = DG.DataFrame.fromColumns([roleCol, curveCol]);
    df.name = plate.data.name;
    df.id = randomizeTableId();
    curveCol.semType ='fit';
    const curvesGrid = df.plot.grid();

    const statisticsAliases = {
      'rSquared': ['rsquared', 'r2', 'r squared', 'r^2'],
      'slope': ['slope', 'hill', 'steepness', 'hill slope'],
      'bottom': ['bottom', 'min', 'minimum', 'miny', 'min y'],
      'top': ['top', 'max', 'maximum', 'maxy', 'max y'],
      'interceptX': ['interceptx', 'intercept x', 'ic50', 'ic 50', 'ic-50', 'ic_50'],
      'auc': ['auc', 'area under the curve', 'area under curve', 'area'],
    };

    const actualStatNames: Record<string, string> = {};
    actOptions.statisticsColumns.forEach((stat) => {
      const alias = Object.entries(statisticsAliases).find(([_, aliases]) => aliases.includes(stat.toLowerCase()));
      if (alias)
        actualStatNames[alias[0]] = stat;
    });

    Object.entries(actualStatNames).forEach(([statName, alias]) => {
      const params = {table: df, colName: curveCol.name, propName: statName, seriesNumber: 0};
      const col: DG.Column = (DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(params).callSync({processed: false})).getOutputParamValue();
      col.name = alias;
    });

    if (actualStatNames['interceptX'])
      df.col(actualStatNames['interceptX']) && (df.col(actualStatNames['interceptX'])!.meta.format = 'scientific');

    const criteriaCol = df.columns.addNewString(df.columns.getUnusedName('Criteria'));
    criteriaCol.applyFormula(`if(${actOptions.categorizeFormula}, "Qualified", "Fails Criteria")`, DG.COLUMN_TYPE.STRING);
    criteriaCol.meta.colors.setCategorical({'Fails Criteria': 4294922560, 'Qualified': 4283477800});

    curvesGrid.root.style.width = '100%';
    pw.root.style.display = 'flex';
    pw.root.style.flexDirection = 'column';
    pw.root.style.height = '100%';

    pw.root.appendChild(curvesGrid.root);

    let prevSelection: {seriesIndex: number, pointIndex: number, markerType: FitMarkerType, markerSize: number, markerColor: string, curvesGridCell?: DG.GridCell} | null = null;

    pw.mapFromRowFunc = (row) => mapFromRow(row, (rowIdx, checkBoxState) => {
      // NEW: Use source tracking for checkbox-originated outlier markings
      plate._markOutlierWithSource(rowIdx, checkBoxState, 'user-checkbox');
      if (prevSelection && prevSelection.curvesGridCell)
        setOutlier(prevSelection.curvesGridCell, {x: 0, y: 0, outlier: !checkBoxState}, 0, prevSelection.pointIndex);
      pw.grid.invalidate();
    });


    function clearPreviousSelection() {
      try {
        if (!prevSelection || (prevSelection.seriesIndex ?? -1) < 0)
          return;

        const series = curveCol.get(prevSelection.seriesIndex);
        if (!series)
          return;
        const parsed = JSON.parse(series);
        const points: IFitPoint[] = parsed.series[0]?.points;
        if (points && points.length && points.length > prevSelection.pointIndex) {
          points[prevSelection.pointIndex].marker = prevSelection.markerType;
          points[prevSelection.pointIndex].size = prevSelection.markerSize;
          points[prevSelection.pointIndex].color = prevSelection.markerColor;
          curveCol.set(prevSelection.seriesIndex, JSON.stringify(parsed));
        }
      } catch (e) {
        console.error(e);
      } finally {
        prevSelection = null;
      }
    }

    pw.subs.push(pw.grid.onCurrentCellChanged.subscribe((gc) => {
      clearPreviousSelection();
      if (gc?.gridRow == null || gc?.gridRow == -1 || !gc?.gridColumn || gc?.gridColumn.idx == 0)
        return;
      const row = pw.plate._idx(gc.gridRow, gc.gridColumn.idx - 1);
      if (row == undefined || row < 0)
        return;
      const catValue = pw.plate.data.get(actOptions.roleName, row)?.toLowerCase();
      if (!catValue)
        return;
      const seriesIndex = seriesVals.findIndex(([serName, _]) => serName?.toLowerCase() === catValue);
      if (seriesIndex < 0)
        return;
      const conscentration: number = pw.plate.data.get(actOptions.concentrationName, row);
      const value: number = pw.plate.data.get(normed ? actOptions.normalizedColName! : actOptions.valueName, row);

      const pointInSeriesIndex: number = seriesVals[seriesIndex][1].points.findIndex((p) => p.x === conscentration && p.y === value);
      if (pointInSeriesIndex < 0)
        return;

      prevSelection = {seriesIndex, pointIndex: pointInSeriesIndex, markerType: seriesVals[seriesIndex][1].points[pointInSeriesIndex].marker ?? DG.MARKER_TYPE.CIRCLE,
        markerSize: seriesVals[seriesIndex][1].points[pointInSeriesIndex].size ?? FitConstants.POINT_PX_SIZE, markerColor: seriesVals[seriesIndex][1].points[pointInSeriesIndex].color ?? DG.Color.toHtml(DG.Color.getCategoricalColor(0))};

      const curveJSON = JSON.parse(curveCol.get(seriesIndex)!);
      const points: IFitPoint[] = curveJSON.series[0]?.points;
      if (points && points.length && points.length > pointInSeriesIndex) {
        points[pointInSeriesIndex].marker = DG.MARKER_TYPE.SQUARE;
        points[pointInSeriesIndex].size = FitConstants.POINT_PX_SIZE * 2;
        points[pointInSeriesIndex].color = DG.Color.toHtml(DG.Color.green);
        curveCol.set(seriesIndex, JSON.stringify(curveJSON));
      }

      const curveGridRow = curvesGrid.tableRowToGrid(seriesIndex);
      prevSelection.curvesGridCell = curvesGrid.cell(curveCol.name, curveGridRow);
      prevSelection.curvesGridCell && curvesGrid.scrollToCell(curveCol.name, curveGridRow);
    }));

    pw.subs.push(grok.events.onCustomEvent('fit-cell-outlier-toggle').subscribe((args: FitCellOutlierToggleArgs) => {
      if (!args || !args.gridCell || !args.series || args.pointIdx == null || args.gridCell.cell.column !== curveCol)
        return;


      const point: IFitPoint = args.series.points[args.pointIdx];

      if (point.meta !== null) {
        plate._markOutlier(point.meta, !!point.outlier);
        pw.grid.invalidate();
      }
    }));


    const s = DG.debounce(curvesGrid.onAfterDrawContent, 300).subscribe(() => {
      s.unsubscribe();
      curvesGrid.col(curveCol.name) && (curvesGrid.col(curveCol.name)!.width = 400);
      curvesGrid.props.rowHeight = 200;
    });

    if (actOptions.submitAction) {
      const btn = ui.button('Save to ELN', () => {
        savePlate(plate).then((_) => grok.shell.info('Plate saved'));
      });
      pw.plateActionsDiv!.appendChild(btn);
    }
    return pw;
  }

  static createAnalysisViewWithMapping(
    plate: Plate,
    currentMappings: Map<string, string>,
    onMap: (target: string, source: string) => void,
    onUndo: (target: string) => void,
    plateWidget: PlateWidget
  ): HTMLElement {
    const analysisView = new BaseAnalysisView(
      plate,
      {
        analysisName: 'Dose Response Curve',
        requiredFields: [
          {name: 'Activity', required: true, description: 'Response/activity values'},
          {name: 'Concentration', required: true, description: 'Concentration/dose values'},
          {name: 'SampleID', required: true, description: 'Sample/compound identifiers'}
        ],
        createResultsView: (plate, mappings) => {
          return this.createCurvesGrid(plate, plateWidget, mappings, 'csv');
        }
      },
      currentMappings,
      onMap,
      onUndo
    );

    return analysisView.getRoot();
  }


  static createCurvesGrid(
    plate: Plate,
    plateWidget: PlateWidget,
    mappingsOrOptions?: Partial<AnalysisOptions> | Map<string, string>,
    layout: 'csv' | 'excel' = 'csv'
  ): HTMLElement | null {
    const defaultOptions: AnalysisOptions = {
      roleName: 'layout', concentrationName: 'concentration', valueName: 'readout',
      normalize: true, controlColumns: ['High Control', 'Low Control'], autoFilterOutliers: true,
      categorizeFormula: '${rSquared} > 0.8 && ${Hill} > 0.25 && ${Max} > 80 && ${Max} < 120',
      statisticsColumns: ['Min', 'Max', 'rSquared', 'Hill', 'IC50', 'AUC'],
    };

    let actOptions: AnalysisOptions & {normalizedColName?: string};

    // Handle both mappings and options
    if (mappingsOrOptions instanceof Map) {
    // It's a mappings Map - use the mapped column names
      const mappings = mappingsOrOptions;
      if (layout === 'csv') {
        actOptions = {
          ..._getOptionsForCsv(plate, defaultOptions),
          roleName: mappings.get('SampleID') || defaultOptions.roleName,
          concentrationName: mappings.get('Concentration') || defaultOptions.concentrationName,
          valueName: mappings.get('Activity') || defaultOptions.valueName,
        };
      } else {
        actOptions = {
          ..._getOptionsForExcel(plate, defaultOptions),
          roleName: mappings.get('SampleID') || defaultOptions.roleName,
          concentrationName: mappings.get('Concentration') || defaultOptions.concentrationName,
          valueName: mappings.get('Activity') || defaultOptions.valueName,
        };
      }
    } else {
    // It's options - use existing logic
      const options = mappingsOrOptions || {};
      if (layout === 'csv')
        actOptions = {..._getOptionsForCsv(plate, defaultOptions), ...options};
      else
        actOptions = {..._getOptionsForExcel(plate, defaultOptions), ...options};
    }

    // Rest of the method remains the same...
    if (!plate.data.columns.contains(actOptions.roleName) ||
    !plate.data.columns.contains(actOptions.concentrationName) ||
    !plate.data.columns.contains(actOptions.valueName))
      return null;


    const drFilterOptions: IPlateWellFilter = {exclude: {[actOptions.roleName]: actOptions.controlColumns}};
    let normed = false;
    if (actOptions.normalize && actOptions.controlColumns.length === 2 && plate.data.col(actOptions.roleName)?.categories.includes(actOptions.controlColumns[0]) && plate.data.col(actOptions.roleName)?.categories.includes(actOptions.controlColumns[1])) {
      const [lStats, hStats] = actOptions.controlColumns.map((colName) => plate.getStatistics(actOptions.valueName, ['mean', 'std'], {match: {[actOptions.roleName]: colName}})).sort((a, b) => a.mean - b.mean);
      if (hStats.mean !== lStats.mean) {
        actOptions.normalizedColName = plate.normalize(actOptions.valueName, (v) => (hStats.mean - v) / (hStats.mean - lStats.mean) * 100).name;
        normed = true;
      }
    }

    const series = plate.doseResponseSeries({...drFilterOptions, value: normed ? actOptions.normalizedColName! : actOptions.valueName, concentration: actOptions.concentrationName, groupBy: actOptions.roleName});
    const seriesVals = Object.entries(series);
    if (seriesVals.length === 0) return ui.divText('No dose-response series found to display.');

    const hasAtLeastOneCurve = seriesVals.some(([_, s]) => s.points.length > 1);
    if (!hasAtLeastOneCurve)
      return ui.divText('Data is available, but no compound has multiple data points to fit a curve. Please check your data layout.', 'info-message');


    const minMax = normed ? {minY: -10, maxY: 110} : {};
    const roleCol = DG.Column.string(actOptions.roleName, seriesVals.length);
    const curveCol = DG.Column.string('Curve', seriesVals.length);

    curveCol.init((i) => {
      const currentSeries = seriesVals[i][1];
      const seriesData = {
        points: currentSeries.points, name: currentSeries.name, fitFunction: FIT_FUNCTION_4PL_REGRESSION,
        clickToToggle: true, droplines: ['IC50'], showPoints: 'points',
      };
      return JSON.stringify({
        chartOptions: {
          xAxisName: actOptions.concentrationName, yAxisName: normed ? `norm(${actOptions.valueName})` : actOptions.valueName,
          logX: true, title: `${seriesVals[i][0]}`, ...minMax,
        },
        series: [seriesData],
      });
    });
    roleCol.init((i) => seriesVals[i][0]);

    const df = DG.DataFrame.fromColumns([roleCol, curveCol]);
    df.name = plate.data.name;
    df.id = randomizeTableId();
    curveCol.semType ='fit';
    const curvesGrid = df.plot.grid();

    const statisticsAliases = {
      'rSquared': ['rsquared', 'r2', 'r squared', 'r^2'],
      'slope': ['slope', 'hill', 'steepness', 'hill slope'],
      'bottom': ['bottom', 'min', 'minimum', 'miny', 'min y'],
      'top': ['top', 'max', 'maximum', 'maxy', 'max y'],
      'interceptX': ['interceptx', 'intercept x', 'ic50', 'ic 50', 'ic-50', 'ic_50'],
      'auc': ['auc', 'area under the curve', 'area under curve', 'area'],
    };

    const actualStatNames: Record<string, string> = {};
    actOptions.statisticsColumns.forEach((stat) => {
      const alias = Object.entries(statisticsAliases).find(([_, aliases]) => aliases.includes(stat.toLowerCase()));
      if (alias)
        actualStatNames[alias[0]] = stat;
    });

    Object.entries(actualStatNames).forEach(([statName, alias]) => {
      const params = {table: df, colName: curveCol.name, propName: statName, seriesNumber: 0};
      const col: DG.Column = (DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(params).callSync({processed: false})).getOutputParamValue();
      col.name = alias;
    });

    if (actualStatNames['interceptX'])
      df.col(actualStatNames['interceptX'])!.meta.format = 'scientific';


    let prevSelection: {seriesIndex: number, pointIndex: number, markerType: FitMarkerType, markerSize: number, markerColor: string, curvesGridCell?: DG.GridCell} | null = null;

    function clearPreviousSelection() {
      try {
        if (!prevSelection || (prevSelection.seriesIndex ?? -1) < 0) return;
        const s = curveCol.get(prevSelection.seriesIndex);
        if (!s) return;
        const parsed = JSON.parse(s);
        const points: IFitPoint[] = parsed.series[0]?.points;
        if (points && points.length > prevSelection.pointIndex) {
          points[prevSelection.pointIndex].marker = prevSelection.markerType;
          points[prevSelection.pointIndex].size = prevSelection.markerSize;
          points[prevSelection.pointIndex].color = prevSelection.markerColor;
          curveCol.set(prevSelection.seriesIndex, JSON.stringify(parsed));
        }
      } finally {
        prevSelection = null;
      }
    }

    plateWidget.subs.push(plateWidget.grid.onCurrentCellChanged.subscribe((gc) => {
      clearPreviousSelection();
      if (!gc.isTableCell) return;

      const row = plate._idx(gc.gridRow, gc.gridColumn.idx - 1);
      const catValue = plate.data.get(actOptions.roleName, row)?.toLowerCase();
      if (!catValue) return;

      const seriesIndex = seriesVals.findIndex(([serName, _]) => serName?.toLowerCase() === catValue);
      if (seriesIndex < 0) return;

      const concentration: number = plate.data.get(actOptions.concentrationName, row);
      const value: number = plate.data.get(normed ? actOptions.normalizedColName! : actOptions.valueName, row);
      const pointInSeriesIndex: number = seriesVals[seriesIndex][1].points.findIndex((p) => p.x === concentration && p.y === value);
      if (pointInSeriesIndex < 0) return;

      const currentPoints = seriesVals[seriesIndex][1].points;
      prevSelection = {seriesIndex, pointIndex: pointInSeriesIndex,
        markerType: currentPoints[pointInSeriesIndex].marker ?? DG.MARKER_TYPE.CIRCLE,
        markerSize: currentPoints[pointInSeriesIndex].size ?? FitConstants.POINT_PX_SIZE,
        markerColor: currentPoints[pointInSeriesIndex].color ?? DG.Color.toHtml(DG.Color.getCategoricalColor(0))
      };

      const curveJSON = JSON.parse(curveCol.get(seriesIndex)!);
      const points: IFitPoint[] = curveJSON.series[0]?.points;
      if (points?.length > pointInSeriesIndex) {
        points[pointInSeriesIndex].marker = DG.MARKER_TYPE.SQUARE;
        points[pointInSeriesIndex].size = FitConstants.POINT_PX_SIZE * 2;
        points[pointInSeriesIndex].color = DG.Color.toHtml(DG.Color.green);
        curveCol.set(seriesIndex, JSON.stringify(curveJSON));
      }

      const curveGridRow = curvesGrid.tableRowToGrid(seriesIndex);
      prevSelection.curvesGridCell = curvesGrid.cell(curveCol.name, curveGridRow);
      prevSelection.curvesGridCell?.grid.scrollToCell(curveCol.name, curveGridRow);
    }));

    plateWidget.subs.push(grok.events.onCustomEvent('fit-cell-outlier-toggle').subscribe((args: FitCellOutlierToggleArgs) => {
      if (!args || !args.gridCell || !args.series || args.pointIdx == null || args.gridCell.cell.column !== curveCol) return;
      const point: IFitPoint = args.series.points[args.pointIdx];
      if (point.meta !== null) {
        plate._markOutlier(point.meta, !!point.outlier);
        plateWidget.grid.invalidate();
      }
    }));
    const s = DG.debounce(curvesGrid.onAfterDrawContent, 300).subscribe(() => {
      s.unsubscribe();
      const col = curvesGrid.col(curveCol.name);
      if (col)
        col.width = 400;

      curvesGrid.props.rowHeight = 200;
      curvesGrid.root.style.width = '100%';
    });
    const coordinator = new DrcAnalysisCoordinator(
      plate,
      plateWidget,
      curvesGrid,
      curveCol,
      seriesVals,
      actOptions,
      drFilterOptions,
      normed
    );

    const saveButton = ui.button('SAVE RESULTS', async () => {
      // Ensure the plate has an ID before proceeding.
      if (!plate.id) {
        grok.shell.warning('Please use the main "CREATE" button to save the plate before saving analysis results.');
        return;
      }

      // This is the results dataframe we've already built.
      const resultsDf = df;

      // These are the parameters for the analysis run.
      const analysisParams = {
        roleName: actOptions.roleName,
        concentrationName: actOptions.concentrationName,
        valueName: actOptions.valueName,
        normalize: actOptions.normalize,
      };

      const roleColumnName = actOptions.roleName;

      // This now calls the NEWLY REFACTORED saveDrcAnalysisResults function.
      // The function signature is the same, but its internal logic now correctly
      // creates a run and saves to the analysis_results table.
      await saveDrcAnalysisResults(plate, resultsDf, analysisParams, roleColumnName);
    });

    saveButton.style.marginTop = '8px';

    const container = ui.divV([
      curvesGrid.root,
      ui.div([saveButton], {style: {display: 'flex', justifyContent: 'flex-end', paddingRight: '4px'}})
    ], 'drc-grid-container');
    container.style.width = '100%';
    container.style.height = '100%';
    container.style.display = 'flex';
    container.style.flexDirection = 'column';

    curvesGrid.root.style.width = '100%';
    curvesGrid.root.style.flexGrow = '1';

    // 3. Return the new container
    return container;
  }
}

class DrcAnalysisCoordinator implements IAnalysisWidgetCoordinator {
  constructor(
    private plate: Plate,
    private plateWidget: PlateWidget,
    private curvesGrid: DG.Grid,
    private curveCol: DG.Column,
    private seriesVals: Array<[string, any]>,
    private actOptions: AnalysisOptions & {normalizedColName?: string},
    private drFilterOptions: IPlateWellFilter,
    private normed: boolean
  ) {}

  onPlateDataChanged(changeType: 'outlier' | 'data' | 'layer', details: any): void {
    if (changeType === 'outlier') {
      // START of CHANGE
      this.regenerateCurveData();
      // END of CHANGE
    }
  }


  refreshAnalysisView(): void {
    // Full refresh - not needed for outlier changes
    console.log('Full analysis view refresh requested');
  }

  private regenerateCurveData(): void {
    // Regenerate the series data with current outlier states
    const series = this.plate.doseResponseSeries({
      ...this.drFilterOptions,
      value: this.normed ? this.actOptions.normalizedColName! : this.actOptions.valueName,
      concentration: this.actOptions.concentrationName,
      groupBy: this.actOptions.roleName
    });

    const newSeriesVals = Object.entries(series);

    // Update the curve column with fresh data that includes current outlier states
    for (let i = 0; i < Math.min(newSeriesVals.length, this.curveCol.length); i++) {
      const currentSeries = newSeriesVals[i][1];
      const outlierCount = currentSeries.points.filter((p) => p.outlier).length;

      const seriesData = {
        points: currentSeries.points,
        name: currentSeries.name,
        fitFunction: FIT_FUNCTION_4PL_REGRESSION,
        clickToToggle: true,
        droplines: ['IC50'],
        showPoints: 'points',
      };

      const minMax = this.normed ? {minY: -10, maxY: 110} : {};

      const curveJson = JSON.stringify({
        chartOptions: {
          xAxisName: this.actOptions.concentrationName,
          yAxisName: this.normed ? `norm(${this.actOptions.valueName})` : this.actOptions.valueName,
          logX: true,
          title: `${newSeriesVals[i][0]}`,
          ...minMax,
        },
        series: [seriesData],
      });

      // Force column to recognize change by using notify=true
      this.curveCol.set(i, curveJson, true);
    }

    // Update our stored series values
    this.seriesVals.length = 0;
    this.seriesVals.push(...newSeriesVals);
  }
}

// Add this new function to the file.

async function saveDrcAnalysisResults(
  plate: Plate,
  resultsDf: DG.DataFrame,
  analysisParams: object, // Mappings, etc.
  roleColumnName: string
): Promise<void> {
  // --- DEBUGGING STEP 0: Add console log to confirm function is called ---
  console.log('%c--- saveDrcAnalysisResults initiated ---', 'color: blue; font-weight: bold;');
  if (!plate.id) {
    grok.shell.error('Plate must be saved before saving analysis results (plate.id is missing).');
    return;
  }

  const pi = DG.TaskBarProgressIndicator.create('Saving DRC results...');
  try {
    const groups = resultsDf.col(roleColumnName)?.categories;
    // --- DEBUGGING STEP 1: Validate the groups ---
    if (!groups || groups.length === 0) {
      const errorMsg = `Could not find groups in column "${roleColumnName}". Cannot save.`;
      console.error(errorMsg);
      grok.shell.error(errorMsg);
      pi.close();
      return;
    }
    console.log('DEBUG: Found groups:', groups);

    const runId = await createAnalysisRun(plate.id, 'DRC', groups);
    // --- DEBUGGING STEP 2: Validate the runId ---
    if (typeof runId !== 'number' || runId <= 0) {
      const errorMsg = `Failed to create a valid analysis run. Received runId: ${runId}`;
      console.error(errorMsg);
      grok.shell.error(errorMsg);
      pi.close();
      return;
    }
    console.log(`DEBUG: Created AnalysisRun with ID: ${runId}`);

    const columnToDbPropName: Record<string, string> = {
      'Curve': 'Curve', 'IC50': 'IC50', 'Hill': 'Hill Slope',
      'rSquared': 'R Squared', 'Min': 'Min', 'Max': 'Max', 'AUC': 'AUC',
    };

    const propertyCache = new Map<string, PlateProperty>();
    for (const dbPropName of Object.values(columnToDbPropName)) {
      const propDef = drcAnalysisProperties.find((p) => p.name === dbPropName);
      if (propDef)
        propertyCache.set(dbPropName, await getOrCreateProperty(propDef.name, propDef.type));
    }
    // --- DEBUGGING STEP 3: Check the created property cache ---
    console.log('DEBUG: Property Cache created:', Object.fromEntries(propertyCache.entries()));
    if (propertyCache.size === 0) {
      grok.shell.error('Property cache is empty. Check property name definitions.');
      pi.close();
      return;
    }


    // --- DEBUGGING STEP 4: Log the first row to be processed ---
    console.log(`DEBUG: Starting to loop through ${resultsDf.rowCount} result rows.`);
    let isFirstRow = true;

    for (const row of resultsDf.rows) {
      const groupKey = row.get(roleColumnName);
      if (!groupKey) continue;
      const groupCombination = [groupKey];

      if (isFirstRow)
        console.log(`DEBUG: Processing first row for group: "${groupKey}"`);


      for (const [columnName, dbPropName] of Object.entries(columnToDbPropName)) {
        if (resultsDf.columns.contains(columnName)) {
          const propObject = propertyCache.get(dbPropName);
          if (!propObject) continue;

          const value = row.get(columnName);
          if (value === null || value === undefined || (typeof value === 'number' && !isFinite(value)))
            continue;

          if (isFirstRow)
            console.log(`  - Preparing to save property "${dbPropName}" (from column "${columnName}") with value:`, value);


          // This is the actual DB call
          await saveAnalysisResult({
            runId: runId,
            propertyId: propObject.id,
            propertyType: propObject.type,
            value: value,
            groupCombination: groupCombination,
          });
        }
      }
      isFirstRow = false; // Only log details for the first row to avoid spam
    }
    grok.shell.info(`Saved ${resultsDf.rowCount} dose-response curves for plate ${plate.barcode}.`);
  } catch (e) {
    // --- DEBUGGING STEP 5: ENHANCED ERROR HANDLING ---
    // This is the most important part. We will force the error to be visible.
    console.error('--- ERROR CAUGHT IN saveDrcAnalysisResults ---', e);
    grok.shell.error('Failed to save DRC results. Check console for details.');
    grok.shell.error(e instanceof Error ? e.message : String(e));
  } finally {
    pi.close();
  }
}

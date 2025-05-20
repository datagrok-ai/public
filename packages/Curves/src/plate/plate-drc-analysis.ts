/* eslint-disable max-len */
import { AnalysisOptions } from "./plate-widget";

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {mapFromRow} from './utils';
import {IPlateWellFilter, Plate, randomizeTableId} from './plate';
import {FIT_FUNCTION_4PL_REGRESSION, FIT_FUNCTION_SIGMOID, FitMarkerType, IFitPoint} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import { FitConstants } from '../fit/const';
import {FitCellOutlierToggleArgs, setOutlier} from '../fit/fit-renderer';
import { _package } from '../package';
import {savePlate} from "../plates/plates-crud";
import { PlateWidget } from './plate-widget';


export class PlateDrcAnalysis {
  static analysisView(plate: Plate, options?: Partial<AnalysisOptions>): PlateWidget {
    const pw = PlateWidget.detailedView(plate.data);
    const gridRoot = pw.grid.root;
    const detailsRoot = pw.detailsDiv;
    pw.plateActionsDiv = ui.div();
    pw.plateDetailsDiv = ui.div();
    detailsRoot!.prepend(pw.plateDetailsDiv);
    detailsRoot!.append(pw.plateActionsDiv);
    const defaultOptions: AnalysisOptions = {
      roleName: 'layout', concentrationName: 'concentration', valueName: 'readout',
      normalize: true, controlColumns: ['High Control', 'Low Control'], autoFilterOutliers: true,
      categorizeFormula: '${rSquared} > 0.8 && ${Hill} > 0.25 && ${Max} > 80 && ${Max} < 120', statisticsColumns: ['Min', 'Max', 'rSquared', 'Hill', 'IC50', 'AUC'],
    };

    // system of aliases for the column names
    const aliases = {
      concentrationName: ['conc', 'dose', 'conc.', 'dilution', 'concentrations'],
      valueName: ['response', 'value', 'signal', 'raw data', 'raw signal', 'raw', 'response', 'values'],
      roleName: ['role', 'plate layout', 'plate positions', 'group', 'compound', 'well']
    };
    Object.entries(aliases).forEach(([key , value]) => {
      const colName = defaultOptions[key as keyof typeof defaultOptions] as string;
      if (plate.data.columns.contains(colName))
        return;
      const alias = value.find((v) => plate.data.columns.contains(v));
      if (alias)
        (defaultOptions[key as keyof typeof defaultOptions] as string) = alias;
    })

    // find control columns
    if (plate.data.columns.contains(defaultOptions.roleName) && !defaultOptions.controlColumns.every((cc) => plate.data.col(defaultOptions.roleName)!.categories.includes(cc))) {
      //  if control columns are not provided, try to find them in the data
      const controlCols = plate.data.col(defaultOptions.roleName)!.categories.filter((cat) => cat.toLowerCase().includes('control'));
      // take maximum of 2 control columns
      defaultOptions.controlColumns = controlCols.slice(0, 2);
    }

    const actOptions: AnalysisOptions & {normalizedColName?: string} = {...defaultOptions, ...options};

    const statisticsAliases = {
      'rSquared': ['rsquared', 'r2', 'r squared', 'r^2'],
      'slope': ['slope', 'hill', 'steepness', 'hill slope'],
      'bottom': ['bottom', 'min', 'minimum', 'miny', 'min y'],
      'top': ['top', 'max', 'maximum', 'maxy', 'max y'],
      'interceptX': ['interceptx', 'intercept x','ic50', 'ic 50', 'ic-50', 'ic_50', 'ic 50 value', 'ic50 value', 'ic50value', 'ec50', 'ec 50', 'ec-50', 'ec_50', 'ec 50 value', 'ec50 value', 'ec50value'],
      'auc': ['auc', 'area under the curve', 'area under curve', 'area'],
    };

    const actualStatNames: Record<string, string> = {}; // will hold actual statistic name to its alias column name
    actOptions.statisticsColumns.forEach((stat) => {
      const alias = Object.entries(statisticsAliases).find(([_, aliases]) => aliases.includes(stat.toLowerCase()));
      if (alias)
        actualStatNames[alias[0]] = stat;
    });


    if (!gridRoot || !detailsRoot)
      return pw;
    // place them horizontally
    const container = ui.divH([gridRoot, detailsRoot]);
    pw.root.appendChild(container);
    detailsRoot.style.flex = '0 0 min(250px, 25%)';
    gridRoot.style.flexGrow = '1';
    gridRoot.style.removeProperty('width');

    const drFilterOptions: IPlateWellFilter = {exclude: {[actOptions.roleName]: actOptions.controlColumns}};
    let normed = false;
    if (actOptions.normalize && actOptions.controlColumns.length === 2) {

      const [lStats, hStats] = actOptions.controlColumns.map((colName) => {
        return plate.getStatistics(actOptions.valueName, ['mean', 'std'], {match: {[actOptions.roleName]: colName}});
      }).sort((a,b) => a.mean - b.mean);

      defaultOptions.plateStatistics = {
        'Z Prime': (_plate) => 1 - (3 * (hStats.std + lStats.std) / Math.abs(hStats.mean - lStats.mean)),
        'Signal to background': (_plate) => hStats.mean / lStats.mean,
      }
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
      // mark outliers that go outside the bounds

      if (actOptions.autoFilterOutliers) {
        plate.markOutliersWhere(actOptions.normalizedColName!, (v) => v > 106 || v < -6, drFilterOptions);
      }
    }

    const series = plate.doseResponseSeries({...drFilterOptions, value: normed ? actOptions.normalizedColName! : actOptions.valueName , concentration: actOptions.concentrationName, groupBy: actOptions.roleName});
    const seriesVals = Object.entries(series);
    const minMax = normed ? {minY: -10, maxY: 110} : 0;
    const roleCol = DG.Column.string(actOptions.roleName, seriesVals.length);
    const curveCol = DG.Column.string('Curve', seriesVals.length);
    curveCol.init((i) => JSON.stringify(
      {
                  "chartOptions": {
                    "xAxisName": actOptions.concentrationName,
                    "yAxisName": normed ? `norm(${actOptions.valueName})` : actOptions.valueName,
                    "logX": true,
                    "title": `${seriesVals[i][0]}`,
                    ...minMax
                  },
                  // TODO: change to 4PL regression once fixed for normed data
                series: [{...seriesVals[i][1], fit: undefined, fitFunction: FIT_FUNCTION_4PL_REGRESSION, clickToToggle: true, droplines: ['IC50'], name: seriesVals[i][0]}]
      }

      ));
    roleCol.init((i) => seriesVals[i][0]);

    const df = DG.DataFrame.fromColumns([roleCol, curveCol]);
    df.name = pw.plateData.name;
    df.id = randomizeTableId();

    curveCol.semType ='fit';
    const curvesGrid = df.plot.grid();

    // add statistics columns
    Object.entries(actualStatNames).forEach(([statName, alias]) => {
      const params = {table: df, colName: curveCol.name, propName: statName, seriesNumber: 0};
      const col: DG.Column = ( DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(params).callSync({processed: false})).getOutputParamValue();
      col.name = alias;
    });
    if (actualStatNames['interceptX'])
      df.col(actualStatNames['interceptX']) && (df.col(actualStatNames['interceptX'])!.meta.format = 'scientific');

    //categorize the curves

    const criteriaCol = df.columns.addNewString(df.columns.getUnusedName('Criteria'));
    criteriaCol.applyFormula(`if(${actOptions.categorizeFormula}, "Qualified", "Fails Criteria")`, 'string');
    criteriaCol.meta.colors.setCategorical({ "Fails Criteria":4294922560, "Qualified":4283477800 });

    curvesGrid.root.style.width = '100%';
    pw.root.style.display = 'flex';
    pw.root.style.flexDirection = 'column';
    pw.root.style.height = '100%';

    pw.root.appendChild(curvesGrid.root);

    // when selecting a cell on plate, go to appropriate curve and mark the corresponding point with different marker
    // remember the previous change, so that we can revert it
    let prevSelection: {seriesIndex: number, pointIndex: number, markerType: FitMarkerType, markerSize: number, markerColor: string, curvesGridCell?: DG.GridCell} | null = null;

    pw.mapFromRowFunc = (row) => mapFromRow(row, (rowIdx, checkBoxState) => {
      plate._markOutlier(rowIdx, checkBoxState);
      if (prevSelection && prevSelection.curvesGridCell)
      setOutlier(prevSelection.curvesGridCell, {x: 0, y: 0, outlier: !checkBoxState}, 0, prevSelection.pointIndex);
      pw.grid.invalidate();
    });

    function clearPreviousSelection() {
      try {
        if (!prevSelection || (prevSelection.seriesIndex ?? -1) < 0) {
          return;
        }
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

    // highlight selected point in the curve
    pw.subs.push(pw.grid.onCurrentCellChanged.subscribe((gc) => {
      clearPreviousSelection();
      if (gc?.gridRow == null || gc?.gridRow == -1 || !gc?.gridColumn || gc?.gridColumn.idx == 0)
        return;
      const row = pw.dataRow(gc);
      if (row == undefined || row < 0)
        return;
      const catValue = pw.plateData.get(actOptions.roleName, row)?.toLowerCase();
      if (!catValue)
        return;
      const seriesIndex = seriesVals.findIndex(([serName, _]) => serName?.toLowerCase() === catValue);
      if (seriesIndex < 0)
        return;
      const conscentration: number = pw.plateData.get(actOptions.concentrationName, row);
      const value: number = pw.plateData.get(normed ? actOptions.normalizedColName! : actOptions.valueName, row);

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

      // scroll in grid
      const curveGridRow = curvesGrid.tableRowToGrid(seriesIndex);
      prevSelection.curvesGridCell = curvesGrid.cell(curveCol.name, curveGridRow);
      prevSelection.curvesGridCell && curvesGrid.scrollToCell(curveCol.name, curveGridRow);
    }));

    // mark outliers in the original plate if it is switched from curve manually

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
        //actOptions.submitAction!(plate, df);
      });
      pw.plateActionsDiv!.appendChild(btn);
    }
    return pw;
  }

}

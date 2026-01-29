/* eslint-disable max-len */
//@ts-ignore
export * from './package.g';

/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {FitGridCellHandler, calculateSeriesStats, getChartDataAggrStats} from './fit/fit-grid-cell-handler';
import {getOrCreateParsedChartData, substituteZeroes} from './fit/fit-renderer';
import {curveDemo} from './fit/fit-demo';
import {convertXMLToIFitChartData} from './fit/fit-parser';
import {LogOptions} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {FitStatistics} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';

// import {PlateWidget} from './plate/plate-widget';

import * as api from './package-api';
import {convertDataToCurves, dataToCurvesUI} from './fit/data-to-curves';

export const _package = new DG.Package();


const SOURCE_COLUMN_TAG = '.sourceColumn';
const SERIES_NUMBER_TAG = '.seriesNumber';
const SERIES_AGGREGATION_TAG = '.seriesAggregation';
const STATISTICS_TAG = '.statistics';

export class Sync {
  private static _currentPromise: Promise<any> = Promise.resolve();
  public static async runWhenDone<T>(func: () => Promise<T>): Promise<T> {
    Sync._currentPromise = Sync._currentPromise.then(async () => { try { return await func(); } catch (e) { _package.logger.error(e); } });
    return Sync._currentPromise;
  }
  // the number at the end is the column version
}


export class PackageFunctions {
  @grok.decorators.demo({
    name: 'Curve fitting',
    description: 'Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points',
    meta: {demoPath: 'Curves | Curve Fitting'},
    test: {test: 'curveFitDemo()', wait: '2000'}
  })
  static async curveFitDemo(): Promise<void> {
    await curveDemo();
  }

  // @grok.decorators.func({
  //   name: 'Assay Plates',
  //   description: 'Assasy plates with concentration, layout and readout data',
  //   meta: {demoPath: 'Curves | Assay Plates'},
  // })
  // static async assayPlatesDemo(): Promise<void> {
  //   const plateFile = (await grok.dapi.files.list('System:DemoFiles/hts/xlsx_plates'))[0];
  //   grok.shell.addView(await PackageFunctions.previewPlateXlsx(plateFile) as DG.ViewBase);
  // }

  @grok.decorators.init()
  static _initCurves(): void {
    DG.ObjectHandler.register(new FitGridCellHandler());
  }

  @grok.decorators.func()
  static async dataToCurves(df: DG.DataFrame, concentrationCol: DG.Column, readoutCol: DG.Column, batchIDCol: DG.Column, assayCol: DG.Column,
    runIDCol: DG.Column, compoundIDCol: DG.Column, targetEntityCol: DG.Column, @grok.decorators.param({options: {nullable: true}})excludeOutliersCol?: DG.Column,
    // rest is parent level data
    @grok.decorators.param({options: {nullable: true}})parentTable?: DG.DataFrame, // these inputs need to be string and resolved here bellow, because this function is used in datasync, otherwise context is lost
    @grok.decorators.param({options: {nullable: true}})fitParamColumns?: string[],
    @grok.decorators.param({options: {nullable: true}})reportedIC50Column?: string,
    @grok.decorators.param({options: {nullable: true}})reportedQualifiedIC50Column?: string,
    @grok.decorators.param({options: {nullable: true}})experimentIDColumn?: string, @grok.decorators.param({options: {nullable: true}})qualifierColumn?: string,
    @grok.decorators.param({options: {nullable: true}})additionalColumns?: string[],
    @grok.decorators.param({options: {nullable: true}})wellLevelJoinCol?: string,
    @grok.decorators.param({options: {nullable: true}})parentLevelJoinCol?: string,
    @grok.decorators.param({options: {nullable: true, optional: true}})wellLevelAdditionalColumns?: string[]
  ): Promise<DG.DataFrame> {
    const pt = parentTable;
    const joinInfo = pt && wellLevelJoinCol && parentLevelJoinCol && df.col(wellLevelJoinCol) && pt.col(parentLevelJoinCol) ? {
      wellLevelCol: df.col(wellLevelJoinCol)!,
      parentLevelCol: pt.col(parentLevelJoinCol)!,
    } : undefined;
    const wellLevelAdditionalColumnsAct = (wellLevelAdditionalColumns ?? []).map((c) => df.col(c)).filter((c) => c != null) as DG.Column[];
    // this needs to work with datasync so we use wide format
    return convertDataToCurves(df, concentrationCol, readoutCol, batchIDCol, assayCol, runIDCol, compoundIDCol, targetEntityCol, excludeOutliersCol, {
      table: pt,
      fitParamColumns: (fitParamColumns ?? []).map((c) => pt?.col(c)).filter((c) => c != null) as DG.Column[],
      reportedIC50Column: reportedIC50Column ? pt?.col(reportedIC50Column) ?? undefined : undefined,
      reportedQualifiedIC50Column: reportedQualifiedIC50Column ? pt?.col(reportedQualifiedIC50Column) ?? undefined : undefined,
      experimentIDColumn: experimentIDColumn ? pt?.col(experimentIDColumn) ?? undefined : undefined,
      qualifierColumn: qualifierColumn ? pt?.col(qualifierColumn) ?? undefined : undefined,
      additionalColumns: (additionalColumns ?? []).map((c) => pt?.col(c)).filter((c) => c != null) as DG.Column[],
    },
    joinInfo, wellLevelAdditionalColumnsAct
    );
  }

  @grok.decorators.func({'top-menu': 'Data | Curves | Data to Curves', 'outputs': [{'name': 'result', 'type': 'dynamic'}]})
  static async dataToCurvesTopMenu() {
    dataToCurvesUI();
  }

  @grok.decorators.func({meta: {vectorFunc: 'true', role: 'transform'}})
  static addStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, @grok.decorators.param({type: 'int'}) seriesNumber: number): DG.Column {
    const df = table;
    const col = df.col(colName)!;
    const sourceColName = col.name;
    const initialName = df.columns.getUnusedName(`${colName} ${seriesNumber + 1} ${propName}`);

    //const grid = DG.Viewer.grid(df);
    const column = DG.Column.float(initialName, col.length);

    column.tags[SOURCE_COLUMN_TAG] = col.name;
    column.tags[SERIES_NUMBER_TAG] = seriesNumber?.toString();
    column.tags[STATISTICS_TAG] = propName;

    column
      .init((i) => {
        //const gridCell = DG.GridCell.fromColumnRow(grid, colName, grid.tableRowToGrid(i));
        const cell = df.cell(i, sourceColName);
        if (!cell || !cell.value)
          return null;
        const chartData = cell.column.getTag(FitConstants.TAG_FIT_CHART_FORMAT) === FitConstants.TAG_FIT_CHART_FORMAT_3DX ?
          convertXMLToIFitChartData(cell.value) : getOrCreateParsedChartData(cell, true); // false because there is no dataframe
        if (chartData.series![seriesNumber] === undefined || chartData.series![seriesNumber].points.every((p) => p.outlier))
          return null;
        if (chartData.chartOptions?.allowXZeroes && chartData.chartOptions?.logX &&
          chartData.series?.some((series) => series.points.some((p) => p.x === 0)))
          substituteZeroes(chartData);
        const chartLogOptions: LogOptions = {logX: chartData.chartOptions?.logX, logY: chartData.chartOptions?.logY};
        const fitResult = calculateSeriesStats(chartData.series![seriesNumber], seriesNumber, chartLogOptions, cell, true);
        return fitResult[propName as keyof FitStatistics];
      });

    df.columns.insert(column, df.columns.names().indexOf(colName) + 1);
    return column;
  }

  @grok.decorators.func({meta: {vectorFunc: 'true', role: 'transform'}})
  static addAggrStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, aggrType: string): DG.Column {
    const df = table;
    const col = df.col(colName)!;
    const nName = `${colName} ${aggrType} ${propName}`;
    const column = DG.Column.float(df.columns.getUnusedName(nName), col?.length);

    column.tags[SOURCE_COLUMN_TAG] = colName;
    column.tags[SERIES_AGGREGATION_TAG] = aggrType;
    column.tags[STATISTICS_TAG] = propName;
    column
      .init((i) => {
        const cell = df.cell(i, colName);
        if (!cell || !cell.value)
          return null;
        const chartData =
          cell.column.getTag(FitConstants.TAG_FIT_CHART_FORMAT) === FitConstants.TAG_FIT_CHART_FORMAT_3DX ?
            convertXMLToIFitChartData(cell.value) : getOrCreateParsedChartData(cell);
        if (chartData.series?.every((series) => series.points.every((p) => p.outlier)))
          return null;
        if (chartData.chartOptions?.allowXZeroes && chartData.chartOptions?.logX &&
          chartData.series?.some((series) => series.points.some((p) => p.x === 0)))
          substituteZeroes(chartData);
        const fitResult = getChartDataAggrStats(chartData, aggrType, cell);
        return fitResult[propName as keyof FitStatistics];
      });
    df.columns.insert(column, df.columns.names().indexOf(colName) + 1);
    return column;
  }
}

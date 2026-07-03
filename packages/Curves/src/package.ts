/* eslint-disable max-len */
//@ts-ignore
export * from './package.g';

/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {FitGridCellHandler, calculateSeriesStats, getChartDataAggrStats} from './fit/fit-grid-cell-handler';
import {getOrCreateParsedChartData, substituteZeroes} from './fit/fit-renderer';
import {assayCurvesDemo, curveDemo} from './fit/fit-demo';
import {convertXmlCurveToJson} from './fit/converters/xml-converter';
import {convertCompactDrToJson} from './fit/converters/compact-dr-converter';
import {convertPzfxToJson} from './fit/converters/pzfx-converter';
import {registerCurveConverter, initExternalConverters} from './fit/curve-converter';
import {LogOptions} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {FitStatistics} from '@datagrok-libraries/statistics/src/fit/fit-curve';

// import {PlateWidget} from './plate/plate-widget';

import * as api from './package-api';
import {convertDataToCurves, dataToCurvesUI} from './fit/data-to-curves';
import {parsePzfxXml, pzfxTableToFitChartData, pzfxToFitDataFrame, pzfxTableToDataFrame} from './formats/pzfx/pzfx-parser';

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

  @grok.decorators.demo({
    name: 'Assay Curves',
    description: 'Dashboard with curves for multiple compounds, assays and targets',
    meta: {demoPath: 'Curves | Assay Curves'},
  })
  static async assayCurveFitDemo(): Promise<void> {
    await assayCurvesDemo();
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
  static async _initCurves(): Promise<void> {
    DG.ObjectHandler.register(new FitGridCellHandler());
    // Register local converters by format name (avoids circular dependency with DG.Func.find)
    registerCurveConverter('3dx', convertXmlCurveToJson);
    registerCurveConverter('compact-dr', convertCompactDrToJson);
    registerCurveConverter('pzfx', convertPzfxToJson);
    // Discover converters from external packages (skips already-registered local ones)
    await initExternalConverters();
  }

  @grok.decorators.func({
    name: 'Fit Dose-Response Curves',
    description: 'Group well-level assay data by compound, assay, target, and run, then fit a dose-response curve per group.',
  })
  static async dataToCurves(df: DG.DataFrame,
    @grok.decorators.param({options: {description: 'Concentration (dose) column'}}) concentrationCol: DG.Column,
    @grok.decorators.param({options: {description: 'Readout (response) column'}}) readoutCol: DG.Column,
    @grok.decorators.param({options: {description: 'Batch identifier column'}}) batchIDCol: DG.Column,
    @grok.decorators.param({options: {description: 'Assay name column'}}) assayCol: DG.Column,
    @grok.decorators.param({options: {description: 'Run identifier column'}}) runIDCol: DG.Column,
    @grok.decorators.param({options: {description: 'Compound identifier column'}}) compoundIDCol: DG.Column,
    @grok.decorators.param({options: {description: 'Target entity column'}}) targetEntityCol: DG.Column,
    @grok.decorators.param({options: {nullable: true, description: 'Boolean column marking points to exclude as outliers'}})excludeOutliersCol?: DG.Column,
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

  @grok.decorators.func({
    name: 'Add Curve Statistic Column',
    description: 'Extract a fit statistic (e.g. IC50, AUC, R²) from a specific curve series into a new column.',
    meta: {vectorFunc: 'true', role: 'transform'},
  })
  static addStatisticsColumn(table: DG.DataFrame,
    @grok.decorators.param({options: {description: 'Name of the curve column to read'}}) colName: string,
    @grok.decorators.param({options: {description: 'Fit statistic to extract (e.g. IC50, AUC, R²)'}}) propName: string,
    @grok.decorators.param({type: 'int', options: {description: 'Zero-based index of the curve series'}}) seriesNumber: number): DG.Column {
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
        const chartData = getOrCreateParsedChartData(cell, true);
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

  @grok.decorators.func({
    name: 'Add Aggregated Curve Statistic Column',
    description: 'Aggregate a fit statistic across all series of a curve into a new column.',
    meta: {vectorFunc: 'true', role: 'transform'},
  })
  static addAggrStatisticsColumn(table: DG.DataFrame,
    @grok.decorators.param({options: {description: 'Name of the curve column to read'}}) colName: string,
    @grok.decorators.param({options: {description: 'Fit statistic to aggregate (e.g. IC50, AUC, R²)'}}) propName: string,
    @grok.decorators.param({options: {description: 'Aggregation type applied across series (e.g. avg, min, max)'}}) aggrType: string): DG.Column {
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
        const chartData = getOrCreateParsedChartData(cell);
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

  @grok.decorators.func({description: 'Returns XML 3DX curve converter function', meta: {role: 'curveConverter', curveFormat: '3dx'}})
  static convertXmlCurveToJsonFunc(): (value: string) => string {
    return convertXmlCurveToJson;
  }

  @grok.decorators.func({description: 'Returns compact dose-response JSON converter function', meta: {role: 'curveConverter', curveFormat: 'compact-dr'}})
  static convertCompactDrToJsonFunc(): (value: string) => string {
    return convertCompactDrToJson;
  }

  @grok.decorators.func({description: 'Returns PZFX curve converter function', meta: {role: 'curveConverter', curveFormat: 'pzfx'}})
  static convertPzfxToJsonFunc(): (value: string) => string {
    return convertPzfxToJson;
  }

  @grok.decorators.fileViewer({fileViewer: 'pzfx'})
  static async previewPzfx(file: DG.FileInfo): Promise<DG.View> {
    const view = DG.View.create();
    view.name = file.name;

    const text = await file.readAsString();
    const tables = parsePzfxXml(text);
    const xyTables = tables.filter((t) => t.tableType === 'XY');

    if (xyTables.length === 0) {
      const nonXyDfs = tables.map(pzfxTableToDataFrame);
      for (const df of nonXyDfs)
        view.append(DG.Viewer.grid(df).root);
      return view;
    }

    const df = pzfxToFitDataFrame(xyTables);
    const grid = DG.Viewer.grid(df);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    view.append(grid.root);

    return view;
  }

  @grok.decorators.fileHandler({
    ext: 'pzfx',
    description: 'Open a GraphPad Prism (.pzfx) file as data tables, fitting XY curve tables.',
  })
  static pzfxFileHandler(@grok.decorators.param({type: 'list', options: {description: 'Raw bytes of the .pzfx file'}}) bytes: Uint8Array): DG.DataFrame[] {
    const text = new TextDecoder().decode(new Uint8Array(bytes));
    const tables = parsePzfxXml(text);
    const results: DG.DataFrame[] = [];

    const xyTables = tables.filter((t) => t.tableType === 'XY');
    if (xyTables.length > 0)
      results.push(pzfxToFitDataFrame(xyTables));

    for (const table of tables.filter((t) => t.tableType !== 'XY'))
      results.push(pzfxTableToDataFrame(table));

    return results;
  }
}

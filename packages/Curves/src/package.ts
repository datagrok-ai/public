/* eslint-disable max-len */
//@ts-ignore
export * from './package.g';

/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {FitGridCellHandler} from './fit/fit-grid-cell-handler';
import {calculateSeriesFit, getChartDataAggrStats, curveStatisticAt, curveAggrStatisticAt} from './fit/fit-statistics';
import {getOrCreateParsedChartData, substituteZeroes} from './fit/fit-chart-data';
import {assayCurvesDemo, curveDemo} from './fit/fit-demo';
import {convertXmlCurveToJson} from './fit/converters/xml-converter';
import {convertCompactDrToJson} from './fit/converters/compact-dr-converter';
import {convertPzfxToJson} from './fit/converters/pzfx-converter';
import {registerCurveConverter, initExternalConverters, parseCellValue} from './fit/curve-converter';
import {FitStatistics, IFitChartData, LogOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {getStatistic} from '@datagrok-libraries/statistics/src/fit/fit-engine';

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
    /* The declared `type` must stay `column`: an `options.type` on its own
     * becomes the parameter's OWN type in the generated annotation, turning the
     * column slot into a plain numerical/categorical value. (And never spell an
     * `input:` annotation out inside a comment — the server scans package.ts for
     * those lines and would try to parse it as a real parameter.) */
    // All seven are dereferenced unconditionally, so they are `nullable: false`
    // — a column parameter defaults to nullable, which read as "optional" and
    // let a half-configured node run.
    @grok.decorators.param({type: 'column', options: {type: 'numerical', nullable: false, description: 'Concentration (dose) column'}}) concentrationCol: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'numerical', nullable: false, description: 'Readout (response) column'}}) readoutCol: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'categorical', nullable: false, description: 'Batch identifier column'}}) batchIDCol: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'categorical', nullable: false, description: 'Assay name column'}}) assayCol: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'categorical', nullable: false, description: 'Run identifier column'}}) runIDCol: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'categorical', nullable: false, description: 'Compound identifier column'}}) compoundIDCol: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'categorical', nullable: false, description: 'Target entity column'}}) targetEntityCol: DG.Column,
    @grok.decorators.param({type: 'column', options: {nullable: true, description: 'Boolean column marking points to exclude as outliers'}})excludeOutliersCol?: DG.Column,
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
    name: 'curveStatistic',
    description: 'Extract a fit statistic (e.g. IC50, AUC, R²) from a curve series into a calculated column.',
    meta: {vectorFunc: 'true'},
    outputs: [{name: 'result', type: 'column', options: {action: 'join(table)'}}],
  })
  static curveStatistic(table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'fit', nullable: false, description: 'Curve column to read'}}) curveColumn: DG.Column,
    @grok.decorators.param({options: {nullable: false, description: 'Fit statistic to extract (e.g. ic50, auc, rSquared)'}}) propName: string,
    @grok.decorators.param({type: 'int', options: {nullable: false, description: 'Zero-based index of the curve series'}}) seriesNumber: number): DG.Column {
    // stable name: recalculation matches the result back by name, and AddNewColumn makes it unique on insert
    const result = DG.Column.float(`${curveColumn.name} ${seriesNumber + 1} ${propName}`, curveColumn.length);
    result.init((i) => curveStatisticAt(curveColumn, i, propName, seriesNumber, table));
    return result;
  }

  @grok.decorators.func({
    name: 'curveAggrStatistic',
    description: 'Aggregate a fit statistic across all series of a curve into a calculated column.',
    meta: {vectorFunc: 'true'},
    outputs: [{name: 'result', type: 'column', options: {action: 'join(table)'}}],
  })
  static curveAggrStatistic(table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'fit', nullable: false, description: 'Curve column to read'}}) curveColumn: DG.Column,
    @grok.decorators.param({options: {nullable: false, description: 'Fit statistic to aggregate'}}) propName: string,
    @grok.decorators.param({options: {choices: ['avg', 'med', 'min', 'max', 'q1', 'q2', 'q3'], initialValue: 'med', nullable: false, description: 'Aggregation applied across the series of each curve'}}) aggrType: string): DG.Column {
    // stable name: recalculation matches the result back by name, and AddNewColumn makes it unique on insert
    const result = DG.Column.float(`${curveColumn.name} ${aggrType} ${propName}`, curveColumn.length);
    result.init((i) => curveAggrStatisticAt(curveColumn, i, propName, aggrType, table));
    return result;
  }

  @grok.decorators.func({
    name: 'Add Curve Statistic Column',
    description: 'Extract a fit statistic (e.g. IC50, AUC, R²) from a specific curve series into a new column.',
    meta: {vectorFunc: 'true', role: 'transform'},
  })
  static addStatisticsColumn(table: DG.DataFrame,
    @grok.decorators.param({options: {description: 'Name of the curve column to read'}}) colName: string,
    // Literal strings, not `statisticsProperties.map(...)` — the func generator
    // only reads literal arrays and emits nothing for a computed one.
    @grok.decorators.param({options: {choices: ['rSquared', 'auc', 'interceptX', 'interceptY', 'slope', 'top', 'bottom'], initialValue: 'interceptX', description: 'Fit statistic to extract. interceptX is IC50, top and bottom are max/min Y'}}) propName: string,
    @grok.decorators.param({type: 'int', options: {initialValue: '0', description: 'Zero-based index of the curve series'}}) seriesNumber: number): DG.Column {
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
        // resolve by name so both current statistic names and legacy ones from saved projects work
        const fitResult = calculateSeriesFit(chartData.series![seriesNumber], seriesNumber, chartLogOptions, cell, true);
        return getStatistic(fitResult, propName) ?? null;
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
    @grok.decorators.param({options: {choices: ['rSquared', 'auc', 'interceptX', 'interceptY', 'slope', 'top', 'bottom'], initialValue: 'interceptX', description: 'Fit statistic to aggregate. interceptX is IC50, top and bottom are max/min Y'}}) propName: string,
    @grok.decorators.param({options: {choices: ['min', 'max', 'avg', 'med', 'q1', 'q2', 'q3'], initialValue: 'med', description: 'Aggregation applied across the series of each curve'}}) aggrType: string): DG.Column {
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

  /* The two functions above address the curve column by NAME, because their
   * `colName` string is what the "+" buttons in the Fit pane and the Data to
   * Curves pipeline pass, and both are recorded as `role: transform` steps that
   * must keep replaying. On a pipeline canvas a name string means no column
   * picker and no `fit` filter, so these twins take a real column slot and
   * delegate to the same implementation. */

  @grok.decorators.func({
    name: 'Add Curve Statistic',
    description: 'Extracts a fit statistic from one series of a curve column into a new column.',
    outputs: [{name: 'result', type: 'column'}],
    meta: {vectorFunc: 'true', role: 'transform'},
  })
  static addCurveStatistic(
    @grok.decorators.param({options: {caption: 'Table', nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'fit', caption: 'Curves', nullable: false, description: 'Column of fitted curves'}}) curvesCol: DG.Column,
    @grok.decorators.param({type: 'string', options: {caption: 'Statistic', nullable: false, choices: ['rSquared', 'auc', 'interceptX', 'interceptY', 'slope', 'top', 'bottom'], initialValue: 'interceptX', description: 'interceptX is IC50, top and bottom are max/min Y'}}) statistic: string = 'interceptX',
    @grok.decorators.param({type: 'int', options: {caption: 'Series', nullable: false, initialValue: '0', description: 'Zero-based index of the curve series'}}) seriesNumber: number = 0): DG.Column {
    return PackageFunctions.addStatisticsColumn(table, curvesCol.name, statistic, seriesNumber);
  }

  @grok.decorators.func({
    name: 'Add Aggregated Curve Statistic',
    description: 'Aggregates a fit statistic across all series of a curve column into a new column.',
    outputs: [{name: 'result', type: 'column'}],
    meta: {vectorFunc: 'true', role: 'transform'},
  })
  static addAggrCurveStatistic(
    @grok.decorators.param({options: {caption: 'Table', nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'fit', caption: 'Curves', nullable: false, description: 'Column of fitted curves'}}) curvesCol: DG.Column,
    @grok.decorators.param({type: 'string', options: {caption: 'Statistic', nullable: false, choices: ['rSquared', 'auc', 'interceptX', 'interceptY', 'slope', 'top', 'bottom'], initialValue: 'interceptX', description: 'interceptX is IC50, top and bottom are max/min Y'}}) statistic: string = 'interceptX',
    @grok.decorators.param({type: 'string', options: {caption: 'Aggregation', nullable: false, choices: ['med', 'avg', 'min', 'max', 'sum', 'stdev', 'variance', 'q1', 'q2', 'q3'], initialValue: 'med', description: 'Applied across the series of each curve'}}) aggregation: string = 'med'): DG.Column {
    return PackageFunctions.addAggrStatisticsColumn(table, curvesCol.name, statistic, aggregation);
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

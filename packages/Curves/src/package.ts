import {_PlateGridCellRenderer} from './package.g';
import {_MultiCurveViewer} from './package.g';
import {_FitChartCellRenderer} from './package.g';

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
import {FitConstants} from './fit/const';
import {PlateCellHandler} from "./plate/plate-cell-renderer";
import {FitSeries} from '@datagrok-libraries/statistics/src/fit/new-fit-API';
import {Plate} from './plate/plate';
import {getPlatesFolderPreview, PlateWidget} from './plate/plate-widget';
//@ts-ignore
import * as jStat from 'jstat';
import {PlateReader} from "./plate/plate-reader";
import {initPlatesAppTree, platesAppView} from "./plates/plates_app";
import {__createDummyPlateData, initPlates, savePlate} from "./plates/plates_crud";

export const _package = new DG.Package();
const SOURCE_COLUMN_TAG = '.sourceColumn';
const SERIES_NUMBER_TAG = '.seriesNumber';
const SERIES_AGGREGATION_TAG = '.seriesAggregation';
const STATISTICS_TAG = '.statistics';



export class PackageFunctions {

  @grok.decorators.demo({
    name: 'Curve fitting',
    description: 'Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points',
    meta: { demoPath: 'Curves | Curve Fitting' },
    test: { test: 'curveFitDemo()', wait: '2000' }
  })
  static async curveFitDemo(): Promise<void> {
    await curveDemo();
  }

  @grok.decorators.func({
    name: 'Assay Plates',
    description: 'Assasy plates with concentration, layout and readout data',
    meta: { demoPath: 'Curves | Assay Plates' },
  })
  static async assayPlatesDemo(): Promise<void> {
    const plateFile = (await grok.dapi.files.list('System:DemoFiles/hts/xlsx_plates'))[0];
    grok.shell.addView(await PackageFunctions.previewPlateXlsx(plateFile) as DG.ViewBase);
  }

  @grok.decorators.init({})
  static _initCurves(): void {
    DG.ObjectHandler.register(new FitGridCellHandler());
    DG.ObjectHandler.register(new PlateCellHandler());
  }

  @grok.decorators.func({ tags: ['Transform'] })
  static addStatisticsColumn(df: DG.DataFrame, colName: string, propName: string, seriesName: string, seriesNumber: number, newColName?: string): void {
    const grid = DG.Viewer.grid(df);
    const chartColumn = grid.col(colName)!;
    const column = DG.Column.float(newColName ?? `${colName} ${seriesName} ${propName}`, chartColumn.column?.length);
    column.tags[SOURCE_COLUMN_TAG] = colName;
    column.tags[SERIES_NUMBER_TAG] = seriesNumber;
    column.tags[STATISTICS_TAG] = propName;

    column
      .init((i) => {
        const gridCell = DG.GridCell.fromColumnRow(grid, colName, grid.tableRowToGrid(i));
        if (gridCell.cell.value === '')
          return null;
        const chartData = gridCell.cell.column.getTag(FitConstants.TAG_FIT_CHART_FORMAT) === FitConstants.TAG_FIT_CHART_FORMAT_3DX ?
          convertXMLToIFitChartData(gridCell.cell.value) : getOrCreateParsedChartData(gridCell);
        if (chartData.series![seriesNumber] === undefined || chartData.series![seriesNumber].points.every((p) => p.outlier))
          return null;
        if (chartData.chartOptions?.allowXZeroes && chartData.chartOptions?.logX &&
          chartData.series?.some((series) => series.points.some((p) => p.x === 0)))
          substituteZeroes(chartData);
        const chartLogOptions: LogOptions = { logX: chartData.chartOptions?.logX, logY: chartData.chartOptions?.logY };
        const fitResult = calculateSeriesStats(chartData.series![seriesNumber], seriesNumber, chartLogOptions, gridCell);
        return fitResult[propName as keyof FitStatistics];
      });
    df.columns.insert(column, chartColumn.idx);
  }

  @grok.decorators.func({ tags: ['Transform'] })
  static addAggrStatisticsColumn(df: DG.DataFrame, colName: string, propName: string, aggrType: string): void {
    const grid = DG.Viewer.grid(df);
    const chartColumn = grid.col(colName)!;
    const column = DG.Column.float(`${colName} ${aggrType} ${propName}`, chartColumn.column?.length);
    column.tags[SOURCE_COLUMN_TAG] = colName;
    column.tags[SERIES_AGGREGATION_TAG] = aggrType;
    column.tags[STATISTICS_TAG] = propName;

    column
      .init((i) => {
        const gridCell = DG.GridCell.fromColumnRow(grid, colName, grid.tableRowToGrid(i));
        if (gridCell.cell.value === '')
          return null;
        const chartData = gridCell.cell.column.getTag(FitConstants.TAG_FIT_CHART_FORMAT) === FitConstants.TAG_FIT_CHART_FORMAT_3DX ?
          convertXMLToIFitChartData(gridCell.cell.value) : getOrCreateParsedChartData(gridCell);
        if (chartData.series?.every((series) => series.points.every((p) => p.outlier)))
          return null;
        if (chartData.chartOptions?.allowXZeroes && chartData.chartOptions?.logX &&
          chartData.series?.some((series) => series.points.some((p) => p.x === 0)))
          substituteZeroes(chartData);
        const fitResult = getChartDataAggrStats(chartData, aggrType, gridCell);
        return fitResult[propName as keyof FitStatistics];
      });
    df.columns.insert(column, chartColumn.idx);
  }
  
  @grok.decorators.folderViewer({})
  static async platesFolderPreview(folder: DG.FileInfo, files: DG.FileInfo[]): Promise<DG.Widget | DG.ViewBase | undefined> {
    const nameLowerCase = folder.name?.toLowerCase();
    if (!nameLowerCase?.includes('plate'))
      return undefined;
    return getPlatesFolderPreview(files);
  }

  @grok.decorators.fileViewer({  fileViewer: 'txt', fileViewerCheck: 'Curves:checkFileIsPlate' })
  static previewPlate(file: DG.FileInfo): DG.View {
    const view = DG.View.create();
    view.name = file.name;
    file.readAsString().then((content) => {
      const plate = PlateReader.read(content);
      if (plate !== null) {
        view.root.appendChild(PlateWidget.fromPlate(plate).root);
      }
    });
    return view;
  }
  
  @grok.decorators.fileHandler({ ext: 'txt', fileViewerCheck: 'Curves:checkFileIsPlate' })
  static async importPlate(fileContent: string): Promise<DG.DataFrame[]> {
    const plate = PlateReader.read(fileContent);
    const view = DG.View.create();
    if (plate !== null) {
      view.root.appendChild(PlateWidget.fromPlate(plate).root);
    }
    view.name = 'Plate';
    grok.shell.addView(view);
    return [];
  }
  
  @grok.decorators.fileHandler({ outputs:[], ext: 'xlsx', fileViewerCheck: 'Curves:checkExcelIsPlate' })
  static async importPlateXlsx(fileContent: Uint8Array) : Promise<any[]> {
    const view = DG.View.create();
    const plate = await PackageFunctions.parseExcelPlate(fileContent);
    view.root.appendChild(PlateWidget.analysisView(plate).root);
    view.name = 'Plate';
    grok.shell.addView(view);
    return [];
  }
  
  @grok.decorators.fileViewer({ name: 'viewPlateXlsx', fileViewer: 'xlsx', fileViewerCheck: 'Curves:checkExcelIsPlate' })
  static async previewPlateXlsx(file: DG.FileInfo): Promise<DG.View> {
    const view = DG.View.create();
    view.name = file.friendlyName;
    const plate = await PackageFunctions.parseExcelPlate(await file.readAsBytes());
    view.root.appendChild(PlateWidget.analysisView(plate).root);
    return view;
  }
  
  @grok.decorators.func({ })
  static async checkExcelIsPlate(content: Uint8Array): Promise<boolean> {
    try {
      if (content.length > 1_000_000) // haven't really seen a plate file larger than 1MB
        return false;
      const plate = await PackageFunctions.parseExcelPlate(content);
      return plate !== null;
    } catch (e) {
      return false;
    }
  }

  static async parseExcelPlate(content: string | Uint8Array, name?: string) {
    if (typeof content === 'string') {
      const blob = new Blob([content], { type: 'application/octet-binary' });
      const buf = await blob.arrayBuffer();
      const plate = await Plate.fromExcel(new Uint8Array(buf), name);
      return plate;
    } else
      return await Plate.fromExcel(content, name);
  }
  
  @grok.decorators.func({ })
  static checkFileIsPlate(content: string): boolean {
    if (content.length > 1_000_000)
      return false;
    return PlateReader.getReader(content) != null;
  }

  @grok.decorators.app({ name: 'Browse', browsePath: 'Plates'})
  static platesApp() : any {
    return platesAppView();
  }
  
  @grok.decorators.treeBrowser({ })
  static async platesAppTreeBrowser(treeNode: DG.TreeViewGroup) : Promise<void> {
    await initPlatesAppTree(treeNode);
  }
  
  @grok.decorators.func({ })
  static async getPlateByBarcode(barcode: string): Promise<Plate> {
    await initPlates();
    const df: DG.DataFrame = await grok.functions.call('Curves:getWellValuesByBarcode', { barcode: barcode });
    return Plate.fromDbDataFrame(df);
  }

  @grok.decorators.func({ })
  static async createDummyPlateData() : Promise<void> {
    await __createDummyPlateData();
  }
}


//name: Browse
//tags: app
//meta.browsePath: Plates
//output: view result
export function platesApp() {
  return PackageFunctions.platesApp();
}

//input: dynamic treeNode
//input: view browseView
export async function platesAppTreeBrowser(treeNode: DG.TreeViewGroup) {
  await PackageFunctions.platesAppTreeBrowser(treeNode);
}

//name: getPlateByBarcode
//input: string barcode
//output: dynamic plate
export async function getPlateByBarcode(barcode: string): Promise<Plate> {
  return PackageFunctions.getPlateByBarcode(barcode);
}

//name: createDummyPlateData
export async function createDummyPlateData() {
  await PackageFunctions.createDummyPlateData();
}


export * from './package.g';export {_FitChartCellRenderer};
export {_MultiCurveViewer};
export {_PlateGridCellRenderer};

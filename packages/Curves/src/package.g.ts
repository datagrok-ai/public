import {PackageFunctions} from './package';
import {MultiCurveViewer} from './fit/multi-curve-viewer';
import {FitChartCellRenderer} from './fit/fit-renderer';
import * as DG from 'datagrok-api/dg';
//name: Fit
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: fit
//meta.virtual: true
export function _FitChartCellRenderer() {
  return new FitChartCellRenderer();
}

//name: MultiCurveViewer
//description: A viewer that superimposes multiple in-cell curves on one chart
//output: viewer result
//meta.role: viewer
//meta.icon: icons/multi-curve-viewer.png
//meta.trellisable: true
export function _MultiCurveViewer() {
  return new MultiCurveViewer();
}


//name: Curve fitting
//description: Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points
//meta.demoPath: Curves | Curve Fitting
//test: curveFitDemo() //wait: 2000 
export async function curveFitDemo() : Promise<void> {
  await PackageFunctions.curveFitDemo();
}

//name: Assay Curves
//description: Dashboard with curves for multiple compounds, assays and targets
//meta.demoPath: Curves | Assay Curves
export async function assayCurveFitDemo() : Promise<void> {
  await PackageFunctions.assayCurveFitDemo();
}

//meta.role: init
export async function _initCurves() : Promise<void> {
  await PackageFunctions._initCurves();
}

//name: Fit Dose-Response Curves
//description: Group well-level assay data by compound, assay, target, and run, then fit a dose-response curve per group.
//input: dataframe df 
//input: column concentrationCol { description: Concentration (dose) column }
//input: column readoutCol { description: Readout (response) column }
//input: column batchIDCol { description: Batch identifier column }
//input: column assayCol { description: Assay name column }
//input: column runIDCol { description: Run identifier column }
//input: column compoundIDCol { description: Compound identifier column }
//input: column targetEntityCol { description: Target entity column }
//input: column excludeOutliersCol { nullable: true; description: Boolean column marking points to exclude as outliers }
//input: dataframe parentTable { nullable: true }
//input: list<string> fitParamColumns { nullable: true }
//input: string reportedIC50Column { nullable: true }
//input: string reportedQualifiedIC50Column { nullable: true }
//input: string experimentIDColumn { nullable: true }
//input: string qualifierColumn { nullable: true }
//input: list<string> additionalColumns { nullable: true }
//input: string wellLevelJoinCol { nullable: true }
//input: string parentLevelJoinCol { nullable: true }
//input: list<string> wellLevelAdditionalColumns { nullable: true; optional: true }
//output: dataframe result
export async function dataToCurves(df: DG.DataFrame, concentrationCol: DG.Column, readoutCol: DG.Column, batchIDCol: DG.Column, assayCol: DG.Column, runIDCol: DG.Column, compoundIDCol: DG.Column, targetEntityCol: DG.Column, excludeOutliersCol?: DG.Column, parentTable?: DG.DataFrame, fitParamColumns?: string[], reportedIC50Column?: string, reportedQualifiedIC50Column?: string, experimentIDColumn?: string, qualifierColumn?: string, additionalColumns?: string[], wellLevelJoinCol?: string, parentLevelJoinCol?: string, wellLevelAdditionalColumns?: string[]) : Promise<any> {
  return await PackageFunctions.dataToCurves(df, concentrationCol, readoutCol, batchIDCol, assayCol, runIDCol, compoundIDCol, targetEntityCol, excludeOutliersCol, parentTable, fitParamColumns, reportedIC50Column, reportedQualifiedIC50Column, experimentIDColumn, qualifierColumn, additionalColumns, wellLevelJoinCol, parentLevelJoinCol, wellLevelAdditionalColumns);
}

//output: dynamic result
//top-menu: Data | Curves | Data to Curves
export async function dataToCurvesTopMenu() {
  return await PackageFunctions.dataToCurvesTopMenu();
}

//name: Add Curve Statistic Column
//description: Extract a fit statistic (e.g. IC50, AUC, R²) from a specific curve series into a new column.
//input: dataframe table 
//input: string colName { description: Name of the curve column to read }
//input: string propName { description: Fit statistic to extract (e.g. IC50, AUC, R²) }
//input: int seriesNumber { description: Zero-based index of the curve series }
//output: column result
//meta.vectorFunc: true
//meta.role: transform
export function addStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, seriesNumber: number) : any {
  return PackageFunctions.addStatisticsColumn(table, colName, propName, seriesNumber);
}

//name: Add Aggregated Curve Statistic Column
//description: Aggregate a fit statistic across all series of a curve into a new column.
//input: dataframe table 
//input: string colName { description: Name of the curve column to read }
//input: string propName { description: Fit statistic to aggregate (e.g. IC50, AUC, R²) }
//input: string aggrType { description: Aggregation type applied across series (e.g. avg, min, max) }
//output: column result
//meta.vectorFunc: true
//meta.role: transform
export function addAggrStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, aggrType: string) : any {
  return PackageFunctions.addAggrStatisticsColumn(table, colName, propName, aggrType);
}

//description: Returns XML 3DX curve converter function
//output: dynamic result
//meta.role: curveConverter
//meta.curveFormat: 3dx
export function convertXmlCurveToJsonFunc() : any {
  return PackageFunctions.convertXmlCurveToJsonFunc();
}

//description: Returns compact dose-response JSON converter function
//output: dynamic result
//meta.role: curveConverter
//meta.curveFormat: compact-dr
export function convertCompactDrToJsonFunc() : any {
  return PackageFunctions.convertCompactDrToJsonFunc();
}

//description: Returns PZFX curve converter function
//output: dynamic result
//meta.role: curveConverter
//meta.curveFormat: pzfx
export function convertPzfxToJsonFunc() : any {
  return PackageFunctions.convertPzfxToJsonFunc();
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: pzfx
export async function previewPzfx(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewPzfx(file);
}

//description: Open a GraphPad Prism (.pzfx) file as data tables, fitting XY curve tables.
//input: list bytes { description: Raw bytes of the .pzfx file }
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: pzfx
export function pzfxFileHandler(bytes: Uint8Array) : any {
  return PackageFunctions.pzfxFileHandler(bytes);
}

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

//meta.role: init
export function _initCurves() : void {
  PackageFunctions._initCurves();
}

//input: dataframe df 
//input: column concentrationCol 
//input: column readoutCol 
//input: column batchIDCol 
//input: column assayCol 
//input: column runIDCol 
//input: column compoundIDCol 
//input: column targetEntityCol 
//input: column excludeOutliersCol { nullable: true }
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

//input: dataframe table 
//input: string colName 
//input: string propName 
//input: int seriesNumber 
//output: column result
//meta.vectorFunc: true
//meta.role: transform
export function addStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, seriesNumber: number) : any {
  return PackageFunctions.addStatisticsColumn(table, colName, propName, seriesNumber);
}

//input: dataframe table 
//input: string colName 
//input: string propName 
//input: string aggrType 
//output: column result
//meta.vectorFunc: true
//meta.role: transform
export function addAggrStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, aggrType: string) : any {
  return PackageFunctions.addAggrStatisticsColumn(table, colName, propName, aggrType);
}

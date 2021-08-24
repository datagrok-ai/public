import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { ALT, AP, AST, BILIRUBIN, TREATMENT_ARM } from '../constants';
import { addTreatmentArm } from './utils';

export function createMaxValuesData(dataframe, aggregatedColName, filerValue){ 
	let condition = `LBTEST = ${filerValue}`; 
	let grouped = dataframe.groupBy(['USUBJID','LBORRES', 'LBORNRHI'])
	  .where(condition)
	  .aggregate(); 
	grouped.columns.addNewFloat(aggregatedColName)
      .init((i) => parseFloat(grouped.get('LBORRES', i))/parseFloat(grouped.get('LBORNRHI', i)) );
    grouped = grouped.groupBy(['USUBJID'])
	  .max(aggregatedColName)
	  .aggregate();
    grouped.getCol(`max(${aggregatedColName})`).name = aggregatedColName;
    return grouped;
}


export function createHysLawDataframe(lb: DG.DataFrame, dm: DG.DataFrame) {
    let alt = createMaxValuesData(lb, ALT, 'Alanine Aminotransferase');
    let bln = createMaxValuesData(lb, BILIRUBIN, 'Bilirubin');
    let ast = createMaxValuesData(lb, AST, 'Aspartate Aminotransferase');
    let ap = createMaxValuesData(lb, AP, 'Alkaline Phosphatase');

    let joined = grok.data.joinTables(bln, alt, [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', BILIRUBIN ], [ ALT ], DG.JOIN_TYPE.LEFT, false);
    joined = grok.data.joinTables(joined, ast, [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', ALT, BILIRUBIN ], [ AST ], DG.JOIN_TYPE.LEFT, false);
    joined = grok.data.joinTables(joined, ap, [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', ALT, BILIRUBIN, AST ], [ AP ], DG.JOIN_TYPE.LEFT, false);
    let withTreatmentArm = addTreatmentArm(joined, dm, [ 'USUBJID', ALT, BILIRUBIN, AST, AP ]);
    return withTreatmentArm;
}


export function createFilteredFloatValuesDataframe(df: DG.DataFrame & any, condition: string, groupCols: string[], newFloatCol: string, colToTransform: string) {
  let filteredData = df.groupBy(groupCols)
  .where(condition)
  .aggregate();
  filteredData.columns.addNewFloat(newFloatCol).init((i) => parseFloat(filteredData.get(colToTransform, i)));
  return filteredData;
}


export function createBaselineEndpointDataframe(lb: DG.DataFrame,
  dm: DG.DataFrame,
  labValue: string, 
  blVisit: string, 
  epVisit: string, 
  visitCol: string,
  blNumColumn: string,
  epNumColumn: string) {
  let condition = `LBTEST = ${labValue} and ${visitCol} IN`;
  let filteredDataBaseline = createFilteredFloatValuesDataframe(lb, `${condition} (${blVisit})`, [ 'USUBJID', 'LBORRES', 'LBORNRLO', 'LBORNRHI', visitCol ], blNumColumn, 'LBORRES');
  let filteredDataEndpoint = createFilteredFloatValuesDataframe(lb, `${condition} (${epVisit})`, [ 'USUBJID', 'LBORRES', 'LBORNRLO', 'LBORNRHI', visitCol ], epNumColumn, 'LBORRES');
  let joined = grok.data.joinTables(filteredDataBaseline, filteredDataEndpoint,
    [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', blNumColumn, 'LBORNRLO', 'LBORNRHI' ], [ epNumColumn ],
    DG.JOIN_TYPE.LEFT, false);
  let withTreatmentArm = addTreatmentArm(joined, dm, [ 'USUBJID', blNumColumn, epNumColumn, 'LBORNRLO', 'LBORNRHI' ]);  
  return withTreatmentArm;
  
}


export function createLabValuesByVisitDataframe(lb: DG.DataFrame, dm: DG.DataFrame,
  labValue: string,
  treatmentArm: string,
  labValueNumCol: string,
  visitCol: string) {
  let condition = `LBTEST = ${labValue}`;
  let joined =  grok.data.joinTables(lb, dm, [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', 'LBORRES', 'LBTEST', visitCol ],
  [ TREATMENT_ARM ], DG.JOIN_TYPE.LEFT, false);
  let filtered =  createFilteredFloatValuesDataframe(joined, condition, [ 'USUBJID', 'LBORRES', TREATMENT_ARM, visitCol ], labValueNumCol, 'LBORRES');
  filtered.rows.filter((row) => row.visitdy != DG.INT_NULL && row.actarm == treatmentArm);
  return filtered;
}


export function createKaplanMeierDataframe(){
  return grok.shell.table('kaplan_meier_data');
}


export function createUniqueCountDataframe(df: DG.DataFrame, groupCol: string[], uniqueCountCol: string, renameCol: string){
  let grouped = df.groupBy(groupCol)
    .uniqueCount(uniqueCountCol)
    .aggregate();
    grouped.getCol(`unique(${uniqueCountCol})`).name = renameCol;
    return grouped;
}


export function createAERiskAssessmentDataframe(dm: DG.DataFrame, ae: DG.DataFrame, ex: DG.DataFrame) {

  let subjArm = createUniqueCountDataframe(ex, [ 'EXTRT' ], 'USUBJID', 'TOTALSUBJ');
  let joinedAeEX = grok.data.joinTables(ae, ex, [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', 'AETERM' ], [ 'EXTRT' ], DG.JOIN_TYPE.LEFT, false);
  let subjAE = createUniqueCountDataframe(joinedAeEX, [ 'AETERM', 'EXTRT' ], 'USUBJID', 'AECOUNT');

  let tj = grok.data.joinTables(subjAE, subjArm, [ 'EXTRT' ], [ 'EXTRT' ], [ 'AETERM', 'EXTRT', 'AECOUNT' ], [ 'TOTALSUBJ' ], DG.JOIN_TYPE.LEFT, false);

  tj.columns.addNewFloat('PERCENT')
    .init((i) => parseFloat(tj.get('AECOUNT', i)) / parseFloat(tj.get('TOTALSUBJ', i)));

  let aeRiskRatioPlacebo = tj.groupBy([ 'AETERM', 'EXTRT', 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ])
    .where('EXTRT = PLACEBO')
    .aggregate();

  let aeRiskRatioActive = tj.groupBy([ 'AETERM', 'EXTRT', 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ])
    .where('EXTRT = XANOMELINE')
    .aggregate();


  let tj2 = grok.data.joinTables(aeRiskRatioActive, aeRiskRatioPlacebo, [ 'AETERM' ], [ 'AETERM' ],
    [ 'AETERM', 'EXTRT', 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ],
    [ 'AETERM', 'EXTRT', 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ], DG.JOIN_TYPE.OUTER, false);



  let column1 = tj2.columns.byName('null.AETERM');
  let percent1 = tj2.columns.byName('null.PERCENT');
  let exposed1 = tj2.columns.byName('null.AECOUNT');
  let column2 = tj2.columns.byName('null.AETERM (2)');
  let percent2 = tj2.columns.byName('null.PERCENT (2)');
  let exposed2 = tj2.columns.byName('null.AECOUNT (2)');
  let totalExposed1;
  let totalExposed2;

  let rowCount = tj2.rowCount;

  for (let i = 0; i < rowCount; i++) {
    if (column1.isNone(i)) {
      const value = column2.get(i);
      column1.set(i, value)
      percent1.set(i, 0);
      exposed1.set(i, 0);
    } else {
      totalExposed1 = parseFloat(tj2.get('null.TOTALSUBJ', i))
    }
  }


  for (let i = 0; i < rowCount; i++) {
    if (column2.isNone(i)) {
      const value = column1.get(i);
      column2.set(i, value)
      percent2.set(i, 0);
      exposed2.set(i, 0);
    } else {
      totalExposed2 = parseFloat(tj2.get('null.TOTALSUBJ (2)', i))
    }
  }


  tj2.columns.addNewFloat('RELATIVE RISK')
    .init((i) => parseFloat(tj2.get('null.PERCENT', i)) / parseFloat(tj2.get('null.PERCENT (2)', i)));

  tj2.columns.addNewFloat('RISK DIFF')
    .init((i) => parseFloat(tj2.get('null.PERCENT', i)) - parseFloat(tj2.get('null.PERCENT (2)', i)));


  tj2.columns.addNewFloat('ODDS RATIO').init((i) => (parseFloat(tj2.get('null.AECOUNT', i)) * (totalExposed2 - parseFloat(tj2.get('null.AECOUNT (2)', i)))) /
    (parseFloat(tj2.get('null.AECOUNT (2)', i)) * (totalExposed1 - parseFloat(tj2.get('null.AECOUNT', i)))));

}



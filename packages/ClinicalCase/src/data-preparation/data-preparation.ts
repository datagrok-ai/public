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



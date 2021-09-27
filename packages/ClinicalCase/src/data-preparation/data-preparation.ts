import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { ALT, AP, AST, BILIRUBIN, TREATMENT_ARM } from '../constants';
import { addTreatmentArm, dateDifferenceInDays, filterBooleanColumn, filterNulls } from './utils';
import { study } from '../clinical-study';

export function createMaxValuesDataForHysLaw(dataframe, aggregatedColName, filerValue){ 
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
    let alt = createMaxValuesDataForHysLaw(lb, ALT, 'Alanine Aminotransferase');
    let bln = createMaxValuesDataForHysLaw(lb, BILIRUBIN, 'Bilirubin');
    let ast = createMaxValuesDataForHysLaw(lb, AST, 'Aspartate Aminotransferase');
    let ap = createMaxValuesDataForHysLaw(lb, AP, 'Alkaline Phosphatase');

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


export function addColumnWithDrugPlusDosage(df: DG.DataFrame, drugCol: string, dosageCol: string, unitsCol: string, newCol: string){
  df.columns.addNewString(newCol)
    .init((i) => `${df.get(drugCol, i).toString()} ${df.get(dosageCol, i).toString()}${df.get(unitsCol, i).toString()}` );
  return df;
}

export function createSurvivalData(endpoint: string, SDTMendpoint: string, covariates: string[]) {
  let dm = study.domains.dm.clone();
  filterNulls(dm, 'RFENDTC');
  if (SDTMendpoint === 'AESTDTC') {
    const ae = study.domains.ae.clone();
    if(endpoint == 'HOSPITALIZATION'){
      filterBooleanColumn(ae, 'AESHOSP', false);
    }
    const condition = endpoint == 'DRUG RELATED AE' ? 'AEREL not in (NONE, NOT RELATED)' : 'AESEV = SEVERE';
    const aeGrouped = ae.groupBy([ 'USUBJID' ]).
      min('AESEQ').
      where(condition).
      aggregate();
    const aeJoined = grok.data.joinTables(ae, aeGrouped, [ 'USUBJID', 'AESEQ' ],
      [ 'USUBJID', 'min(AESEQ)' ], [ 'USUBJID', 'AESTDTC' ], [ 'min(AESEQ)' ], DG.JOIN_TYPE.LEFT, false);
    filterNulls(aeJoined, 'min(AESEQ)');
    dm = grok.data.joinTables(dm, aeJoined, [ 'USUBJID' ], [ 'USUBJID' ], dm.columns.names(), [ 'AESTDTC' ], DG.JOIN_TYPE.LEFT, false);
  }
  dm.columns.addNewInt('time')
    .init((i) => getSurvivalTime(dm.columns.byName(SDTMendpoint), dm.columns.byName('RFSTDTC'), dm.columns.byName('RFENDTC'), i));
  dm.columns.addNewInt('status')
    .init((i) => getSurvivalStatus(dm.columns.byName(SDTMendpoint), i));
  return dm.groupBy([ 'USUBJID', 'time', 'status' ].concat(covariates)).aggregate();
}

export function getSurvivalTime(eventColumn: DG.Column, startColumn: DG.Column, endColumn: DG.Column, i: number) {
  if (eventColumn.isNone(i)) {
    return dateDifferenceInDays(startColumn.get(i).toString(), endColumn.get(i).toString());
  } else {
    return dateDifferenceInDays(startColumn.get(i).toString(), eventColumn.get(i).toString());
  }
}

export function getSurvivalStatus(eventColumn: DG.Column, i: number) {
  return eventColumn.isNone(i) ? 0 : 1;
}

export function createAERiskAssessmentDataframe(ae: DG.DataFrame, ex: DG.DataFrame) {

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

  tj2.getCol(`null.AETERM`).name = 'AETERM';

  return tj2.groupBy([ 'AETERM', 'RELATIVE RISK', 'RISK DIFF', 'ODDS RATIO'])
    .aggregate();

}


export function labDynamicComparedToBaseline(dataframe, baselineVisitNum: string, newColName: string) {
  const dfName = dataframe.name;
  let grouped = dataframe.groupBy([ 'USUBJID', 'LBTEST', 'LBORRES' ])
    .where(`VISITDY = ${baselineVisitNum}`)
    .aggregate();
  grouped.getCol(`LBORRES`).name = 'BL_LBORRES';
  dataframe.join(grouped,
    [ 'USUBJID', 'LBTEST' ], [ 'USUBJID', 'LBTEST' ],
    dataframe.columns.names(), [ 'BL_LBORRES' ], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => parseFloat(dataframe.get('LBORRES', i)) / parseFloat(dataframe.get('BL_LBORRES', i)));
  dataframe.name = dfName;
}

export function labDynamicComparedToMinMax(dataframe, newColName: string) {
  const dfName = dataframe.name;
  dataframe.columns.addNewFloat('LBORRES_FLOAT')
    .init((i) => parseFloat(dataframe.get('LBORRES', i)));
  let groupedMax = dataframe.groupBy([ 'LBTEST'])
    .max('LBORRES_FLOAT')
    .aggregate();
  let groupedMin = dataframe.groupBy([ 'LBTEST'])
    .min('LBORRES_FLOAT')
    .aggregate();
  dataframe.join(groupedMax, [ 'LBTEST' ], [ 'LBTEST' ], dataframe.columns.names(), [ 'max(LBORRES_FLOAT)' ], DG.JOIN_TYPE.LEFT, true)
  .join(groupedMin,[ 'LBTEST' ], [ 'LBTEST' ], dataframe.columns.names(), [ 'min(LBORRES_FLOAT)' ], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => 
    (parseFloat(dataframe.get('LBORRES', i)) - parseFloat(dataframe.get('min(LBORRES_FLOAT)', i))) / 
    (parseFloat(dataframe.get('max(LBORRES_FLOAT)', i)) - parseFloat(dataframe.get('min(LBORRES_FLOAT)', i))));
  dataframe.name = dfName;
}


export function cumulativeEnrollemntByDay(df: DG.DataFrame, dateCol: string, subjIDCol: string, cumCol: string) {
  const subjsPerDay = df.groupBy([ dateCol ]).uniqueCount(subjIDCol).aggregate();
  const refStartCol = subjsPerDay.col(dateCol);
  for (let i = 0, rowCount = subjsPerDay.rowCount; i < rowCount; i++) {
    if (refStartCol.isNone(i))
      subjsPerDay.rows.removeAt(i);
  }
  subjsPerDay.columns.addNewInt(cumCol)
    .init((i) => {
      if (i > 0) {
        return subjsPerDay.get(cumCol, i - 1) + subjsPerDay.get(`unique(${subjIDCol})`, i)
      }
      return subjsPerDay.get(`unique(${subjIDCol})`, i);
    });
  subjsPerDay.columns.remove(`unique(${subjIDCol})`)
  return subjsPerDay;
}



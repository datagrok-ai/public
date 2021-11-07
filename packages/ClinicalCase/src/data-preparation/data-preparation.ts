import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { ALT, AP, AST, BILIRUBIN } from '../constants';
import { addDataFromDmDomain, dateDifferenceInDays, filterBooleanColumn, filterNulls, getUniqueValues } from './utils';
import { study } from '../clinical-study';
import { AE_CAUSALITY, AE_REQ_HOSP, AE_SEQ, AE_SEVERITY, AE_START_DATE, AE_TERM, INV_DRUG_NAME, LAB_HI_LIM_N, LAB_LO_LIM_N, LAB_RES_N, LAB_TEST, SUBJECT_ID, SUBJ_REF_ENDT, SUBJ_REF_STDT, TREATMENT_ARM } from '../columns-constants';

export function createMaxValuesDataForHysLaw(dataframe, aggregatedColName, filerValue){ 
	let condition = `${LAB_TEST} = ${filerValue}`; 
	let grouped = dataframe.groupBy([SUBJECT_ID,LAB_RES_N, LAB_HI_LIM_N])
	  .where(condition)
	  .aggregate(); 
	grouped.columns.addNewFloat(aggregatedColName)
      .init((i) => grouped.get(LAB_RES_N, i)/grouped.get(LAB_HI_LIM_N, i));
    grouped = grouped.groupBy([SUBJECT_ID])
	  .max(aggregatedColName)
	  .aggregate();
    grouped.getCol(`max(${aggregatedColName})`).name = aggregatedColName;
    return grouped;
}


export function createHysLawDataframe(lb: DG.DataFrame, dm: DG.DataFrame) {
    let alt = createMaxValuesDataForHysLaw(lb, ALT, 'Alanine Aminotransferase');
    let bln = createMaxValuesDataForHysLaw(lb, BILIRUBIN, 'Bilirubin');
    let ast = createMaxValuesDataForHysLaw(lb, AST, 'Aspartate Aminotransferase');
    //let ap = createMaxValuesDataForHysLaw(lb, AP, 'Alkaline Phosphatase');

    let joined = grok.data.joinTables(bln, alt, [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, BILIRUBIN ], [ ALT ], DG.JOIN_TYPE.LEFT, false);
    joined = grok.data.joinTables(joined, ast, [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, ALT, BILIRUBIN ], [ AST ], DG.JOIN_TYPE.LEFT, false);
    //joined = grok.data.joinTables(joined, ap, [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, ALT, BILIRUBIN, AST ], [ AP ], DG.JOIN_TYPE.LEFT, false);
    let withTreatmentArm = addDataFromDmDomain(joined, dm, [ SUBJECT_ID, ALT, BILIRUBIN, AST ], [TREATMENT_ARM]);
    return withTreatmentArm;
}


export function createFilteredFloatValuesDataframe(df: DG.DataFrame & any, condition: string, groupCols: string[], newFloatCol: string, colToTransform: string) {
  let filteredData = df.groupBy(groupCols)
  .where(condition)
  .aggregate();
  filteredData.columns.addNewFloat(newFloatCol).init((i) => parseFloat(filteredData.get(colToTransform, i)));
  return filteredData;
}

export function createFilteredDataframe(df: DG.DataFrame & any, condition: string, groupCols: string[], newFloatCol: string, colToTransform: string) {
  let filteredData = df.groupBy(groupCols)
  .where(condition)
  .aggregate();
  filteredData.columns.addNewFloat(newFloatCol).init((i) => parseFloat(filteredData.get(colToTransform, i)));
  return filteredData;
}


export function createBaselineEndpointDataframe(lb: DG.DataFrame,
  dm: DG.DataFrame,
  columnToExtractFromDm: string[],
  labValue: string, 
  blVisit: string, 
  epVisit: string, 
  visitCol: string,
  blNumColumn: string,
  epNumColumn = null) {
  let condition = `${LAB_TEST} = ${labValue} and ${visitCol} IN`;
  let filteredDataBaseline = createFilteredDataframe(lb, `${condition} (${blVisit})`, [ SUBJECT_ID, LAB_RES_N, LAB_LO_LIM_N, LAB_HI_LIM_N, visitCol, LAB_TEST ], blNumColumn, LAB_RES_N);
  let finalDf;
  let columnsToEXtract = [];
  if(epVisit){
    let filteredDataEndpoint = createFilteredDataframe(lb, `${condition} (${epVisit})`, [ SUBJECT_ID, LAB_RES_N, LAB_LO_LIM_N, LAB_HI_LIM_N, visitCol, LAB_TEST ], epNumColumn, LAB_RES_N);
    finalDf = grok.data.joinTables(filteredDataBaseline, filteredDataEndpoint,
    [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, blNumColumn, LAB_LO_LIM_N, LAB_HI_LIM_N, LAB_TEST], [ epNumColumn ],
    DG.JOIN_TYPE.LEFT, false);
    columnsToEXtract = [ SUBJECT_ID, blNumColumn, epNumColumn, LAB_LO_LIM_N, LAB_HI_LIM_N, LAB_TEST ];
  } else {
    finalDf = filteredDataBaseline;
    columnsToEXtract = [ SUBJECT_ID, blNumColumn, LAB_LO_LIM_N, LAB_HI_LIM_N, LAB_TEST ];
  }
  let withDmData = addDataFromDmDomain(finalDf, dm, columnsToEXtract, columnToExtractFromDm);  
  return withDmData;
}


export function createLabValuesByVisitDataframe(lb: DG.DataFrame, dm: DG.DataFrame,
  labValue: string,
  treatmentArm: string,
  labValueNumCol: string,
  visitCol: string) {
  let condition = `${LAB_TEST} = ${labValue}`;
  let joined =  grok.data.joinTables(lb, dm, [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, LAB_RES_N, LAB_TEST, visitCol ],
  [ TREATMENT_ARM ], DG.JOIN_TYPE.LEFT, false);
  let filtered =  createFilteredDataframe(joined, condition, [ SUBJECT_ID, LAB_RES_N, TREATMENT_ARM, visitCol ], labValueNumCol, LAB_RES_N);
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
  filterNulls(dm, SUBJ_REF_ENDT); // according to SDTMIG v3.3 RFENDTC is null for screen failures or unassigned subjects
  if (SDTMendpoint === AE_START_DATE) {
    const ae = study.domains.ae.clone();
    if(endpoint == 'HOSPITALIZATION'){
      filterBooleanColumn(ae, AE_REQ_HOSP, false);
    }
    const condition = endpoint == 'DRUG RELATED AE' ? `${AE_CAUSALITY} not in (NONE, NOT RELATED)` : `${AE_SEVERITY} = SEVERE`;
    const aeGrouped = ae.groupBy([ SUBJECT_ID ]).
      min(AE_SEQ).
      where(condition).
      aggregate();
    const aeJoined = grok.data.joinTables(ae, aeGrouped, [ SUBJECT_ID, AE_SEQ ],
      [ SUBJECT_ID, `min(${AE_SEQ})` ], [ SUBJECT_ID, AE_START_DATE ], [ `min(${AE_SEQ})` ], DG.JOIN_TYPE.LEFT, false);
    filterNulls(aeJoined, `min(${AE_SEQ})`);
    dm = grok.data.joinTables(dm, aeJoined, [ SUBJECT_ID ], [ SUBJECT_ID ], dm.columns.names(), [ AE_START_DATE ], DG.JOIN_TYPE.LEFT, false);
  }
  dm.columns.addNewInt('time')
    .init((i) => getSurvivalTime(dm.columns.byName(SDTMendpoint), dm.columns.byName(SUBJ_REF_STDT), dm.columns.byName(SUBJ_REF_ENDT), i));
  dm.columns.addNewInt('status')
    .init((i) => getSurvivalStatus(dm.columns.byName(SDTMendpoint), i));
  return dm.groupBy([ SUBJECT_ID, 'time', 'status' ].concat(covariates)).aggregate();
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

export function createAERiskAssessmentDataframe(ae: DG.DataFrame, ex: DG.DataFrame, placeboArm: string, activeArm: string) {

  let subjArm = createUniqueCountDataframe(ex, [ INV_DRUG_NAME ], SUBJECT_ID, 'TOTALSUBJ');
  let joinedAeEX = grok.data.joinTables(ae, ex, [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, AE_TERM ], [ INV_DRUG_NAME ], DG.JOIN_TYPE.LEFT, false);
  let subjAE = createUniqueCountDataframe(joinedAeEX, [ AE_TERM, INV_DRUG_NAME ], SUBJECT_ID, 'AECOUNT');

  let tj = grok.data.joinTables(subjAE, subjArm, [ INV_DRUG_NAME ], [ INV_DRUG_NAME ], [ AE_TERM, INV_DRUG_NAME, 'AECOUNT' ], [ 'TOTALSUBJ' ], DG.JOIN_TYPE.LEFT, false);

  tj.columns.addNewFloat('PERCENT')
    .init((i) => parseFloat(tj.get('AECOUNT', i)) / parseFloat(tj.get('TOTALSUBJ', i)));

  let aeRiskRatioPlacebo = tj.groupBy([ AE_TERM, INV_DRUG_NAME, 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ])
    .where(`${INV_DRUG_NAME} = ${placeboArm}`)
    .aggregate();

  let aeRiskRatioActive = tj.groupBy([ AE_TERM, INV_DRUG_NAME, 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ])
    .where(`${INV_DRUG_NAME} = ${activeArm}`)
    .aggregate();


  let tj2 = grok.data.joinTables(aeRiskRatioActive, aeRiskRatioPlacebo, [ AE_TERM ], [ AE_TERM ],
    [ AE_TERM, INV_DRUG_NAME, 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ],
    [ AE_TERM, INV_DRUG_NAME, 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ], DG.JOIN_TYPE.OUTER, false);



  let column1 = tj2.columns.byName(`null.${AE_TERM}`);
  let percent1 = tj2.columns.byName('null.PERCENT');
  let exposed1 = tj2.columns.byName('null.AECOUNT');
  let column2 = tj2.columns.byName(`null.${AE_TERM} (2)`);
  let percent2 = tj2.columns.byName('null.PERCENT (2)');
  let exposed2 = tj2.columns.byName('null.AECOUNT (2)');
  let totalExposed1;
  let totalExposed2;

  let rowCount = tj2.rowCount;

  for (let i = 0; i < rowCount; i++) {
    if (column1.isNone(i)) {
      const value = column2.get(i);
      column1.set(i, value, false)
      percent1.set(i, 0, false);
      exposed1.set(i, 0, false);
    } else {
      totalExposed1 = parseFloat(tj2.get('null.TOTALSUBJ', i))
    }
  }


  for (let i = 0; i < rowCount; i++) {
    if (column2.isNone(i)) {
      const value = column1.get(i);
      column2.set(i, value, false)
      percent2.set(i, 0, false);
      exposed2.set(i, 0, false);
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

  tj2.getCol(`null.${AE_TERM}`).name = AE_TERM;

  return tj2.groupBy([ AE_TERM, 'RELATIVE RISK', 'RISK DIFF', 'ODDS RATIO'])
    .aggregate();

}


export function labDynamicComparedToBaseline(dataframe, baselineVisitNum: string, blVisitColumn: string, newColName: string, renameNewCol: boolean) {
  const dfName = dataframe.name;
  let grouped = dataframe.groupBy([ SUBJECT_ID, LAB_TEST, LAB_RES_N ])
    .where(`${blVisitColumn} = ${baselineVisitNum}`)
    .aggregate();
  grouped.getCol(LAB_RES_N).name = `BL_${LAB_RES_N}`;
  dataframe.join(grouped,
    [ SUBJECT_ID, LAB_TEST ], [ SUBJECT_ID, LAB_TEST ],
    dataframe.columns.names(), [ `BL_${LAB_RES_N}` ], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => {
      if(dataframe.getCol(LAB_RES_N).isNone(i) || dataframe.getCol(`BL_${LAB_RES_N}`).isNone(i)) {
        return null;
      } else {
        return (dataframe.get(LAB_RES_N, i) - dataframe.get(`BL_${LAB_RES_N}`, i)) / dataframe.get(`BL_${LAB_RES_N}`, i);
      }
    });
    if(renameNewCol){
      dataframe.columns.remove(LAB_RES_N);
      dataframe.getCol(newColName).name = LAB_RES_N;
    }
  dataframe.name = dfName;
}

export function labDynamicComparedToMinMax(dataframe, newColName: string) {
  const dfName = dataframe.name;
  let groupedMax = dataframe.groupBy([ LAB_TEST])
    .max(LAB_RES_N)
    .aggregate();
  let groupedMin = dataframe.groupBy([ LAB_TEST])
    .min(LAB_RES_N)
    .aggregate();
  dataframe.join(groupedMax, [ LAB_TEST ], [ LAB_TEST ], dataframe.columns.names(), [ `max(${LAB_RES_N})` ], DG.JOIN_TYPE.LEFT, true)
  .join(groupedMin,[ LAB_TEST ], [ LAB_TEST ], dataframe.columns.names(), [ `min(${LAB_RES_N})` ], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => 
    (dataframe.get(LAB_RES_N, i) - dataframe.get(`min(${LAB_RES_N})`, i)) / 
    (dataframe.get(`max(${LAB_RES_N})`, i) - dataframe.get(`min(${LAB_RES_N})`, i)));
  dataframe.name = dfName;
}


export function labDynamicRelatedToRef(df: DG.DataFrame, newColName: string){
  df.columns.addNewFloat(newColName)
    .init((i) => {
      const val = df.get(LAB_RES_N, i);
      const min = df.get(LAB_LO_LIM_N, i);
      const max = df.get(LAB_HI_LIM_N, i);
      return val >= max ? val/max : val <= min ? 0 - min/val : ((val-min)/(max-min) - 0.5 )*2;
    });
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

export function getSubjectDmData(subjId: string, columnsToReturn: string[]){
  const dmData = study.domains.dm
  .groupBy(study.domains.dm.columns.names())
  .where(`${SUBJECT_ID} = ${subjId}`)
  .aggregate();
  const dmDict = {};
  columnsToReturn.forEach(it => dmDict[it] = dmData.get(it, 0));
  return dmDict;
}


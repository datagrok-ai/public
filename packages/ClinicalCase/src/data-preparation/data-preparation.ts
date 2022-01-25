import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { ALT, AP, AST, BILIRUBIN } from '../constants';
import { addDataFromDmDomain, dataframeContentToRow, dateDifferenceInDays, filterBooleanColumn, filterNulls, getUniqueValues } from './utils';
import { study } from '../clinical-study';
import { AE_CAUSALITY, AE_REQ_HOSP, AE_SEQ, AE_SEVERITY, AE_START_DATE, AE_TERM, INV_DRUG_NAME, LAB_HI_LIM_N, LAB_LO_LIM_N, LAB_RES_N, LAB_TEST, SUBJECT_ID, SUBJ_REF_ENDT, SUBJ_REF_STDT, TREATMENT_ARM } from '../columns-constants';
var { jStat } = require('jstat')


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


export function createHysLawDataframe(lb: DG.DataFrame, dm: DG.DataFrame, altName: string, astName: string, blnName: string) {
    let alt = createMaxValuesDataForHysLaw(lb, ALT, altName);
    let bln = createMaxValuesDataForHysLaw(lb, BILIRUBIN, blnName);
    let ast = createMaxValuesDataForHysLaw(lb, AST, astName);
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


export function createBaselineEndpointDataframe(df: DG.DataFrame,
  dm: DG.DataFrame,
  columnToExtractFromDm: string[],
  testCol: string,
  resCol: string,
  additionalCols: string[],
  value: string, 
  blVisit: string, 
  epVisit: string, 
  visitCol: string,
  blNumColumn: string,
  epNumColumn = null,) {
  let condition = `${testCol} = ${value} and ${visitCol} IN`;
  let filteredDataBaseline = createFilteredDataframe(df, `${condition} (${blVisit})`, [ SUBJECT_ID, resCol, visitCol, testCol ].concat(additionalCols), blNumColumn, resCol);
  let finalDf;
  let columnsToEXtract = [];
  if(epVisit){
    let filteredDataEndpoint = createFilteredDataframe(df, `${condition} (${epVisit})`, [ SUBJECT_ID, resCol, visitCol, testCol ].concat(additionalCols), epNumColumn, resCol);
    finalDf = grok.data.joinTables(filteredDataBaseline, filteredDataEndpoint,
    [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, blNumColumn, testCol].concat(additionalCols), [ epNumColumn ],
    DG.JOIN_TYPE.LEFT, false);
    columnsToEXtract = [ SUBJECT_ID, blNumColumn, epNumColumn, testCol ].concat(additionalCols);
  } else {
    finalDf = filteredDataBaseline;
    columnsToEXtract = [ SUBJECT_ID, blNumColumn, testCol ].concat(additionalCols);
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
  filtered.rows.filter((row) => row.visitdy !== DG.INT_NULL && row.actarm === treatmentArm);
  return filtered;
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

export function createSurvivalData(dmDf: DG.DataFrame, aeDF: DG.DataFrame, endpoint: string, SDTMendpoint: string, covariates: string[]) {
  let dm = dmDf.clone();
  filterNulls(dm, SUBJ_REF_ENDT); // according to SDTMIG v3.3 RFENDTC is null for screen failures or unassigned subjects
  if (SDTMendpoint === AE_START_DATE) {
    const ae = aeDF.clone();
    if(endpoint == 'HOSPITALIZATION'){
      filterBooleanColumn(ae, AE_REQ_HOSP, false);
    }
    const condition = endpoint == 'DRUG RELATED AE' ? `${AE_CAUSALITY} in (PROBABLE, POSSIBLE, RELATED, UNLIKELY RELATED, POSSIBLY RELATED, RELATED)` : `${AE_SEVERITY} = SEVERE`;
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

  let numberSubjPlaceboArm =  subjArm.groupBy(subjArm.columns.names())
     .where(`${INV_DRUG_NAME} = ${placeboArm}`)
    .aggregate()
    .get('TOTALSUBJ', 0);
  let numberSubjActiveArm =  subjArm.groupBy(subjArm.columns.names())
     .where(`${INV_DRUG_NAME} = ${activeArm}`)
    .aggregate()
    .get('TOTALSUBJ', 0);
  let totalSubj = numberSubjPlaceboArm + numberSubjActiveArm;

  let joinedAeEX = grok.data.joinTables(ae, ex, [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, AE_TERM ], [ INV_DRUG_NAME ], DG.JOIN_TYPE.LEFT, false);
  let subjAE = createUniqueCountDataframe(joinedAeEX, [ AE_TERM, INV_DRUG_NAME ], SUBJECT_ID, 'TOTALSUBJ_WITH_AE');

  let tj = grok.data.joinTables(subjAE, subjArm, [ INV_DRUG_NAME ], [ INV_DRUG_NAME ], [ AE_TERM, INV_DRUG_NAME, 'TOTALSUBJ_WITH_AE' ], [ 'TOTALSUBJ' ], DG.JOIN_TYPE.LEFT, false);

  tj.columns.addNewFloat('PERCENT')
    .init((i) => parseFloat(tj.get('TOTALSUBJ_WITH_AE', i)) / parseFloat(tj.get('TOTALSUBJ', i)));

  let aeRiskRatioPlacebo = tj.groupBy([ AE_TERM, INV_DRUG_NAME, 'TOTALSUBJ_WITH_AE', 'TOTALSUBJ', 'PERCENT' ])
    .where(`${INV_DRUG_NAME} = ${placeboArm}`)
    .aggregate();

  let aeRiskRatioActive = tj.groupBy([ AE_TERM, INV_DRUG_NAME, 'TOTALSUBJ_WITH_AE', 'TOTALSUBJ', 'PERCENT' ])
    .where(`${INV_DRUG_NAME} = ${activeArm}`)
    .aggregate();


  let tj2 = grok.data.joinTables(aeRiskRatioActive, aeRiskRatioPlacebo, [ AE_TERM ], [ AE_TERM ],
    [ AE_TERM, INV_DRUG_NAME, 'TOTALSUBJ_WITH_AE', 'TOTALSUBJ', 'PERCENT' ],
    [ AE_TERM, INV_DRUG_NAME, 'TOTALSUBJ_WITH_AE', 'TOTALSUBJ', 'PERCENT' ], DG.JOIN_TYPE.OUTER, false);



  let column1 = tj2.columns.byName(`null.${AE_TERM}`);
  let percent1 = tj2.columns.byName('null.PERCENT');
  let exposed1 = tj2.columns.byName('null.TOTALSUBJ_WITH_AE');
  let column2 = tj2.columns.byName(`null.${AE_TERM} (2)`);
  let percent2 = tj2.columns.byName('null.PERCENT (2)');
  let exposed2 = tj2.columns.byName('null.TOTALSUBJ_WITH_AE (2)');
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

  tj2.columns.addNewFloat('RR')
    .init((i) => {
return (tj2.get('null.TOTALSUBJ_WITH_AE', i)*numberSubjPlaceboArm)/(tj2.get('null.TOTALSUBJ_WITH_AE (2)', i)*numberSubjActiveArm);
    });

  tj2.columns.addNewFloat('RISK DIFF')
    .init((i) => parseFloat(tj2.get('null.PERCENT', i)) - parseFloat(tj2.get('null.PERCENT (2)', i)));

    tj2.columns.addNewFloat('RD')
    .init((i) => tj2.get('null.TOTALSUBJ_WITH_AE', i)/numberSubjActiveArm - tj2.get('null.TOTALSUBJ_WITH_AE (2)', i)/numberSubjPlaceboArm);


  tj2.columns.addNewFloat('ODDS RATIO').init((i) => (parseFloat(tj2.get('null.TOTALSUBJ_WITH_AE', i)) * (totalExposed2 - parseFloat(tj2.get('null.TOTALSUBJ_WITH_AE (2)', i)))) /
    (parseFloat(tj2.get('null.TOTALSUBJ_WITH_AE (2)', i)) * (totalExposed1 - parseFloat(tj2.get('null.TOTALSUBJ_WITH_AE', i)))));
    
    tj2.columns.addNewFloat('OR').init((i) => {
      const TE = tj2.get('null.TOTALSUBJ_WITH_AE', i);
        const PE = tj2.get('null.TOTALSUBJ_WITH_AE (2)', i);
        return  (TE*(numberSubjPlaceboArm - PE))/(PE*(numberSubjActiveArm - TE))
      });

    tj2.columns.addNewFloat('Z').init((i) => {
      const TE = tj2.get('null.TOTALSUBJ_WITH_AE', i);
      const TN = numberSubjActiveArm - TE;
      const PE = tj2.get('null.TOTALSUBJ_WITH_AE (2)', i);
      const PN = numberSubjPlaceboArm - TE;
      const N = totalSubj;

      const Z = (Math.pow(TE - (numberSubjActiveArm*(TE+PE))/N, 2)*Math.pow(N, 2)*(N-1))/
      (numberSubjActiveArm*numberSubjPlaceboArm*(TE+PE)*(TN+PN));
      return Z;
    });

    tj2.columns.addNewFloat('p-value').init((i) => {
      const Z = tj2.get('Z', i);
      const f = jStat.chisquare.cdf( Z, 1 );
      return 1-f < f ? 1-f : f;
    });

    tj2.columns.addNewFloat('-log 10(p-value)').init((i) => (-Math.log10(tj2.get('p-value', i))));


  tj2.getCol(`null.${AE_TERM}`).name = AE_TERM;

  return tj2;

  /* return tj2.groupBy([ AE_TERM, 'RELATIVE RISK', 'RISK DIFF', 'ODDS RATIO'])
    .aggregate(); */

}


export function dynamicComparedToBaseline(
  dataframe: DG.DataFrame, 
  testCol: string,
  resCol: string,
  baselineVisitNum: string, 
  blVisitColumn: string, 
  newColName: string, 
  renameNewCol: boolean) {
  const dfName = dataframe.name;
  let grouped = dataframe.groupBy([ SUBJECT_ID, testCol, resCol ])
    .where(`${blVisitColumn} = ${baselineVisitNum}`)
    .aggregate();
  grouped.getCol(resCol).name = `BL_${resCol}`;
  dataframe.join(grouped,
    [ SUBJECT_ID, testCol ], [ SUBJECT_ID, testCol ],
    dataframe.columns.names(), [ `BL_${resCol}` ], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => {
      if(dataframe.getCol(resCol).isNone(i) || dataframe.getCol(`BL_${resCol}`).isNone(i)) {
        return null;
      } else {
        return (dataframe.get(resCol, i) - dataframe.get(`BL_${resCol}`, i)) / dataframe.get(`BL_${resCol}`, i);
      }
    });
    if(renameNewCol){
      dataframe.columns.remove(resCol);
      dataframe.getCol(newColName).name = resCol;
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


export function labDynamicRelatedToRef(df: DG.DataFrame, newColName: string) {
  df.columns.addNewFloat(newColName)
    .init((i) => {
      const val = df.get(LAB_RES_N, i);
      const min = df.get(LAB_LO_LIM_N, i);
      const max = df.get(LAB_HI_LIM_N, i);
      if (!val || !min || !max) {
        return null;
      } else {
        return val >= max ? val / max : val <= min ? 0 - min / val : ((val - min) / (max - min) - 0.5) * 2;
      }
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


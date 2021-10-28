import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { ALT, AP, AST, BILIRUBIN, SUBJECT_ID, TREATMENT_ARM } from '../constants';
import { addDataFromDmDomain, dateDifferenceInDays, filterBooleanColumn, filterNulls, getUniqueValues } from './utils';
import { study } from '../clinical-study';

export function createMaxValuesDataForHysLaw(dataframe, aggregatedColName, filerValue){ 
	let condition = `LBTEST = ${filerValue}`; 
	let grouped = dataframe.groupBy(['USUBJID','LBSTRESN', 'LBSTNRHI'])
	  .where(condition)
	  .aggregate(); 
	grouped.columns.addNewFloat(aggregatedColName)
      .init((i) => grouped.get('LBSTRESN', i)/grouped.get('LBSTNRHI', i));
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
    let withTreatmentArm = addDataFromDmDomain(joined, dm, [ 'USUBJID', ALT, BILIRUBIN, AST, AP ], [TREATMENT_ARM]);
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
  let condition = `LBTEST = ${labValue} and ${visitCol} IN`;
  let filteredDataBaseline = createFilteredDataframe(lb, `${condition} (${blVisit})`, [ 'USUBJID', 'LBSTRESN', 'LBSTNRLO', 'LBSTNRHI', visitCol, 'LBTEST' ], blNumColumn, 'LBSTRESN');
  let finalDf;
  let columnsToEXtract = [];
  if(epVisit){
    let filteredDataEndpoint = createFilteredDataframe(lb, `${condition} (${epVisit})`, [ 'USUBJID', 'LBSTRESN', 'LBSTNRLO', 'LBSTNRHI', visitCol, 'LBTEST' ], epNumColumn, 'LBSTRESN');
    finalDf = grok.data.joinTables(filteredDataBaseline, filteredDataEndpoint,
    [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', blNumColumn, 'LBSTNRLO', 'LBSTNRHI', 'LBTEST'], [ epNumColumn ],
    DG.JOIN_TYPE.LEFT, false);
    columnsToEXtract = [ 'USUBJID', blNumColumn, epNumColumn, 'LBSTNRLO', 'LBSTNRHI', 'LBTEST' ];
  } else {
    finalDf = filteredDataBaseline;
    columnsToEXtract = [ 'USUBJID', blNumColumn, 'LBSTNRLO', 'LBSTNRHI', 'LBTEST' ];
  }
  let withDmData = addDataFromDmDomain(finalDf, dm, columnsToEXtract, columnToExtractFromDm);  
  return withDmData;
}


export function createLabValuesByVisitDataframe(lb: DG.DataFrame, dm: DG.DataFrame,
  labValue: string,
  treatmentArm: string,
  labValueNumCol: string,
  visitCol: string) {
  let condition = `LBTEST = ${labValue}`;
  let joined =  grok.data.joinTables(lb, dm, [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', 'LBSTRESN', 'LBTEST', visitCol ],
  [ TREATMENT_ARM ], DG.JOIN_TYPE.LEFT, false);
  let filtered =  createFilteredDataframe(joined, condition, [ 'USUBJID', 'LBSTRESN', TREATMENT_ARM, visitCol ], labValueNumCol, 'LBSTRESN');
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

export function createAERiskAssessmentDataframe(ae: DG.DataFrame, ex: DG.DataFrame, placeboArm: string, activeArm: string) {

  let subjArm = createUniqueCountDataframe(ex, [ 'EXTRT' ], 'USUBJID', 'TOTALSUBJ');
  let joinedAeEX = grok.data.joinTables(ae, ex, [ 'USUBJID' ], [ 'USUBJID' ], [ 'USUBJID', 'AETERM' ], [ 'EXTRT' ], DG.JOIN_TYPE.LEFT, false);
  let subjAE = createUniqueCountDataframe(joinedAeEX, [ 'AETERM', 'EXTRT' ], 'USUBJID', 'AECOUNT');

  let tj = grok.data.joinTables(subjAE, subjArm, [ 'EXTRT' ], [ 'EXTRT' ], [ 'AETERM', 'EXTRT', 'AECOUNT' ], [ 'TOTALSUBJ' ], DG.JOIN_TYPE.LEFT, false);

  tj.columns.addNewFloat('PERCENT')
    .init((i) => parseFloat(tj.get('AECOUNT', i)) / parseFloat(tj.get('TOTALSUBJ', i)));

  let aeRiskRatioPlacebo = tj.groupBy([ 'AETERM', 'EXTRT', 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ])
    .where(`EXTRT = ${placeboArm}`)
    .aggregate();

  let aeRiskRatioActive = tj.groupBy([ 'AETERM', 'EXTRT', 'AECOUNT', 'TOTALSUBJ', 'PERCENT' ])
    .where(`EXTRT = ${activeArm}`)
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

  tj2.getCol(`null.AETERM`).name = 'AETERM';

  return tj2.groupBy([ 'AETERM', 'RELATIVE RISK', 'RISK DIFF', 'ODDS RATIO'])
    .aggregate();

}


export function labDynamicComparedToBaseline(dataframe, baselineVisitNum: string, blVisitColumn: string, newColName: string, renameNewCol: boolean) {
  const dfName = dataframe.name;
  let grouped = dataframe.groupBy([ 'USUBJID', 'LBTEST', 'LBSTRESN' ])
    .where(`${blVisitColumn} = ${baselineVisitNum}`)
    .aggregate();
  grouped.getCol(`LBSTRESN`).name = 'BL_LBSTRESN';
  dataframe.join(grouped,
    [ 'USUBJID', 'LBTEST' ], [ 'USUBJID', 'LBTEST' ],
    dataframe.columns.names(), [ 'BL_LBSTRESN' ], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => {
      if(dataframe.getCol('LBSTRESN').isNone(i) || dataframe.getCol('BL_LBSTRESN').isNone(i)) {
        return null;
      } else {
        return (dataframe.get('LBSTRESN', i) - dataframe.get('BL_LBSTRESN', i)) / dataframe.get('BL_LBSTRESN', i);
      }
    });
    if(renameNewCol){
      dataframe.columns.remove('LBSTRESN')
      dataframe.getCol(newColName).name = 'LBSTRESN';
    }
  dataframe.name = dfName;
}

export function labDynamicComparedToMinMax(dataframe, newColName: string) {
  const dfName = dataframe.name;
  let groupedMax = dataframe.groupBy([ 'LBTEST'])
    .max('LBSTRESN')
    .aggregate();
  let groupedMin = dataframe.groupBy([ 'LBTEST'])
    .min('LBSTRESN')
    .aggregate();
  dataframe.join(groupedMax, [ 'LBTEST' ], [ 'LBTEST' ], dataframe.columns.names(), [ 'max(LBSTRESN)' ], DG.JOIN_TYPE.LEFT, true)
  .join(groupedMin,[ 'LBTEST' ], [ 'LBTEST' ], dataframe.columns.names(), [ 'min(LBSTRESN)' ], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => 
    (dataframe.get('LBSTRESN', i) - dataframe.get('min(LBSTRESN)', i)) / 
    (dataframe.get('max(LBSTRESN)', i) - dataframe.get('min(LBSTRESN)', i)));
  dataframe.name = dfName;
}


export function labDynamicRelatedToRef(df: DG.DataFrame, newColName: string){
  df.columns.addNewFloat(newColName)
    .init((i) => {
      const val = df.get('LBSTRESN', i);
      const min = df.get('LBSTNRLO', i);
      const max = df.get('LBSTNRHI', i);
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


export function labDataForCorrelationMatrix(): DG.DataFrame {
  let subjIds = Array.from(getUniqueValues(study.domains.lb, 'USUBJID'));
  let t = DG.DataFrame.create();
  let newCols = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
  newCols.forEach(col => t.columns.addNewFloat(col));
  t.columns.addNewString('USUBJID');
  t.columns.addNewString('VISIT');
  let newRowIndex = -1;
  subjIds.forEach(it => {
    let lbs = study.domains.lb.groupBy(study.domains.lb.columns.names())
      .where(`USUBJID = ${it}`)
      .aggregate()
    let order = lbs.getSortedOrder([ 'VISITDY' ] as any);
    let lbdy = lbs.get('VISITDY', order[ 0 ])
    t.rows.addNew();
    newRowIndex +=1;
    order.forEach(index => {
      let newLbDy = lbs.get('VISITDY', index)
      if (lbdy !== newLbDy) {
        lbdy = newLbDy;
        newRowIndex += 1;
        t.rows.addNew();
      }
      t.set(lbs.get('LBTEST', index), newRowIndex, lbs.get('LBSTRESN', index))
      t.set('USUBJID', newRowIndex, it)
      t.set('VISIT', newRowIndex, lbs.get('VISIT', index))
    })
  })
  return t;
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


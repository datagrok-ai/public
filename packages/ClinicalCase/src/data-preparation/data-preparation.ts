import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { ACTIVE_ARM_POSTTFIX, AE_PERCENT, ALT, AP, AST, BILIRUBIN, NEG_LOG10_P_VALUE, ODDS_RATIO, PLACEBO_ARM_POSTTFIX, P_VALUE, RELATIVE_RISK, RISK_DIFFERENCE, SE_RD_WITH_SIGN_LEVEL, STANDARD_ERROR_RD } from '../constants';
import { addDataFromDmDomain, dataframeContentToRow, dateDifferenceInDays, filterBooleanColumn, filterNulls } from './utils';
import { study } from '../clinical-study';
import { AE_CAUSALITY, AE_REQ_HOSP, AE_SEQ, AE_SEVERITY, AE_START_DATE, AE_TERM, LAB_HI_LIM_N, LAB_LO_LIM_N, LAB_RES_N, LAB_TEST, SUBJECT_ID, SUBJ_REF_ENDT, SUBJ_REF_STDT, TREATMENT_ARM } from '../columns-constants';
var { jStat } = require('jstat')


export function createMaxValuesDataForHysLaw(dataframe, aggregatedColName, filerValue) {
  let condition = `${LAB_TEST} = ${filerValue}`;
  let grouped = dataframe.groupBy([SUBJECT_ID, LAB_RES_N, LAB_HI_LIM_N])
    .where(condition)
    .aggregate();
  grouped.columns.addNewFloat(aggregatedColName)
    .init((i) => grouped.get(LAB_RES_N, i) / grouped.get(LAB_HI_LIM_N, i));
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

  let joined = grok.data.joinTables(bln, alt, [SUBJECT_ID], [SUBJECT_ID], [SUBJECT_ID, BILIRUBIN], [ALT], DG.JOIN_TYPE.LEFT, false);
  joined = grok.data.joinTables(joined, ast, [SUBJECT_ID], [SUBJECT_ID], [SUBJECT_ID, ALT, BILIRUBIN], [AST], DG.JOIN_TYPE.LEFT, false);
  //joined = grok.data.joinTables(joined, ap, [ SUBJECT_ID ], [ SUBJECT_ID ], [ SUBJECT_ID, ALT, BILIRUBIN, AST ], [ AP ], DG.JOIN_TYPE.LEFT, false);
  let withTreatmentArm = addDataFromDmDomain(joined, dm, [SUBJECT_ID, ALT, BILIRUBIN, AST], [TREATMENT_ARM]);
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
  let filteredDataBaseline = createFilteredDataframe(df, `${condition} (${blVisit})`, [SUBJECT_ID, resCol, visitCol, testCol].concat(additionalCols), blNumColumn, resCol);
  let finalDf;
  let columnsToEXtract = [];
  if (epVisit) {
    let filteredDataEndpoint = createFilteredDataframe(df, `${condition} (${epVisit})`, [SUBJECT_ID, resCol, visitCol, testCol].concat(additionalCols), epNumColumn, resCol);
    finalDf = grok.data.joinTables(filteredDataBaseline, filteredDataEndpoint,
      [SUBJECT_ID], [SUBJECT_ID], [SUBJECT_ID, blNumColumn, testCol].concat(additionalCols), [epNumColumn],
      DG.JOIN_TYPE.LEFT, false);
    columnsToEXtract = [SUBJECT_ID, blNumColumn, epNumColumn, testCol].concat(additionalCols);
  } else {
    finalDf = filteredDataBaseline;
    columnsToEXtract = [SUBJECT_ID, blNumColumn, testCol].concat(additionalCols);
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
  let joined = grok.data.joinTables(lb, dm, [SUBJECT_ID], [SUBJECT_ID], [SUBJECT_ID, LAB_RES_N, LAB_TEST, visitCol],
    [TREATMENT_ARM], DG.JOIN_TYPE.LEFT, false);
  let filtered = createFilteredDataframe(joined, condition, [SUBJECT_ID, LAB_RES_N, TREATMENT_ARM, visitCol], labValueNumCol, LAB_RES_N);
  filtered.rows.filter((row) => row.visitdy !== DG.INT_NULL && row.actarm === treatmentArm);
  return filtered;
}


export function createUniqueCountDataframe(df: DG.DataFrame, groupCol: string[], uniqueCountCol: string, renameCol: string) {
  let grouped = df.groupBy(groupCol)
    .uniqueCount(uniqueCountCol)
    .aggregate();
  grouped.getCol(`unique(${uniqueCountCol})`).name = renameCol;
  return grouped;
}


export function addColumnWithDrugPlusDosage(df: DG.DataFrame, drugCol: string, dosageCol: string, unitsCol: string, newCol: string) {
  df.columns.addNewString(newCol)
    .init((i) => `${df.get(drugCol, i).toString()} ${df.get(dosageCol, i).toString()}${df.get(unitsCol, i).toString()}`);
  return df;
}

export function createSurvivalData(dmDf: DG.DataFrame, aeDF: DG.DataFrame, endpoint: string, SDTMendpoint: string, covariates: string[]) {
  let dm = dmDf.clone();
  filterNulls(dm, SUBJ_REF_ENDT); // according to SDTMIG v3.3 RFENDTC is null for screen failures or unassigned subjects
  if (SDTMendpoint === AE_START_DATE) {
    const ae = aeDF.clone();
    if (endpoint == 'HOSPITALIZATION') {
      filterBooleanColumn(ae, AE_REQ_HOSP, false);
    }
    const condition = endpoint == 'DRUG RELATED AE' ? `${AE_CAUSALITY} in (PROBABLE, POSSIBLE, RELATED, UNLIKELY RELATED, POSSIBLY RELATED, RELATED)` : `${AE_SEVERITY} = SEVERE`;
    const aeGrouped = ae.groupBy([SUBJECT_ID]).
      min(AE_SEQ).
      where(condition).
      aggregate();
    const aeJoined = grok.data.joinTables(ae, aeGrouped, [SUBJECT_ID, AE_SEQ],
      [SUBJECT_ID, `min(${AE_SEQ})`], [SUBJECT_ID, AE_START_DATE], [`min(${AE_SEQ})`], DG.JOIN_TYPE.LEFT, false);
    filterNulls(aeJoined, `min(${AE_SEQ})`);
    dm = grok.data.joinTables(dm, aeJoined, [SUBJECT_ID], [SUBJECT_ID], dm.columns.names(), [AE_START_DATE], DG.JOIN_TYPE.LEFT, false);
  }
  dm.columns.addNewInt('time')
    .init((i) => getSurvivalTime(dm.columns.byName(SDTMendpoint), dm.columns.byName(SUBJ_REF_STDT), dm.columns.byName(SUBJ_REF_ENDT), i));
  dm.columns.addNewInt('status')
    .init((i) => getSurvivalStatus(dm.columns.byName(SDTMendpoint), i));
  return dm.groupBy([SUBJECT_ID, 'time', 'status'].concat(covariates)).aggregate();
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

export function createAERiskAssessmentDataframe(ae: DG.DataFrame, dm: DG.DataFrame, placeboArm: string[], activeArm: string[], signLevel: number) {

  let TRTARM = 'TREATMENT_ARM';
  let TOTALSUBJ = 'TOTALSUBJ';
  let TOTALSUBJ_WITH_AE = 'TOTALSUBJ_WITH_AE';
  let PERCENT = 'PERCENT';

  dm.columns.addNewString(TRTARM).init((i) => {
    if (placeboArm.includes(dm.get(TREATMENT_ARM, i))) {
      return placeboArm.join(',');
    }
    if (activeArm.includes(dm.get(TREATMENT_ARM, i))) {
      return activeArm.join(',');
    }
    return dm.get(TREATMENT_ARM, i);
  })

  let subjArm = createUniqueCountDataframe(dm, [TRTARM], SUBJECT_ID, TOTALSUBJ);

  let numberSubjPlaceboArm = subjArm.groupBy(subjArm.columns.names())
    .where({ [TRTARM]: `${placeboArm.join(',')}` })
    .aggregate()
    .get(TOTALSUBJ, 0);
  let numberSubjActiveArm = subjArm.groupBy(subjArm.columns.names())
    .where({ [TRTARM]: `${activeArm.join(',')}` })
    .aggregate()
    .get(TOTALSUBJ, 0);
  let totalSubj = numberSubjPlaceboArm + numberSubjActiveArm;

  let joinedAeEX = grok.data.joinTables(ae, dm, [SUBJECT_ID], [SUBJECT_ID], [SUBJECT_ID, AE_TERM], [TRTARM], DG.JOIN_TYPE.LEFT, false);
  let subjAE = createUniqueCountDataframe(joinedAeEX, [AE_TERM, TRTARM], SUBJECT_ID, TOTALSUBJ_WITH_AE);

  let subjTotalAndWithAE = grok.data.joinTables(subjAE, subjArm, [TRTARM], [TRTARM], [AE_TERM, TRTARM, TOTALSUBJ_WITH_AE], [TOTALSUBJ], DG.JOIN_TYPE.LEFT, false);

  subjTotalAndWithAE.columns.addNewFloat(PERCENT)
    .init((i) => parseFloat(subjTotalAndWithAE.get(TOTALSUBJ_WITH_AE, i)) / parseFloat(subjTotalAndWithAE.get(TOTALSUBJ, i)));

  let aeRiskRatioPlacebo = subjTotalAndWithAE.groupBy([AE_TERM, TRTARM, TOTALSUBJ_WITH_AE, TOTALSUBJ, PERCENT])
    .where({ [TRTARM]: `${placeboArm.join(',')}` })
    .aggregate();
  aeRiskRatioPlacebo.name = `${placeboArm.join(',')}`;

  let aeRiskRatioActive = subjTotalAndWithAE.groupBy([AE_TERM, TRTARM, TOTALSUBJ_WITH_AE, TOTALSUBJ, PERCENT])
    .where({ [TRTARM]: `${activeArm.join(',')}` })
    .aggregate();
  aeRiskRatioActive.name = `${activeArm.join(',')}`;

  let renameCols = (df: DG.DataFrame, columnNames: string[], arm: string) => {
    columnNames.forEach(it => {
      df.getCol(it).name = `${it}_${arm}`;
    })
  }

  renameCols(aeRiskRatioPlacebo, [TRTARM, TOTALSUBJ_WITH_AE, TOTALSUBJ, PERCENT], PLACEBO_ARM_POSTTFIX)
  renameCols(aeRiskRatioActive, [TRTARM, TOTALSUBJ_WITH_AE, TOTALSUBJ, PERCENT], ACTIVE_ARM_POSTTFIX)

  let joinedFinal = grok.data.joinTables(aeRiskRatioActive, aeRiskRatioPlacebo, [AE_TERM], [AE_TERM],
    [AE_TERM, `${TRTARM}_${ACTIVE_ARM_POSTTFIX}`, `${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, `${TOTALSUBJ}_${ACTIVE_ARM_POSTTFIX}`, `${PERCENT}_${ACTIVE_ARM_POSTTFIX}`],
    [AE_TERM, `${TRTARM}_${PLACEBO_ARM_POSTTFIX}`, `${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, `${TOTALSUBJ}_${PLACEBO_ARM_POSTTFIX}`, `${PERCENT}_${PLACEBO_ARM_POSTTFIX}`],
    DG.JOIN_TYPE.OUTER, false);


  joinedFinal.getCol(`${activeArm.join(',')}.${AE_TERM}`).name = `${AE_TERM}_${ACTIVE_ARM_POSTTFIX}`;
  joinedFinal.getCol(`${placeboArm.join(',')}.${AE_TERM}`).name = `${AE_TERM}_${PLACEBO_ARM_POSTTFIX}`;


  let removeNulls = (rowCount, colToCheck, colWithWalue, percentCol, totalSubjWithAECol) => {
    for (let i = 0; i < rowCount; i++) {
      if (colToCheck.isNone(i)) {
        const value = colWithWalue.get(i);
        colToCheck.set(i, value, false)
        percentCol.set(i, 0, false);
        totalSubjWithAECol.set(i, 0, false);
      }
    }
  }

  removeNulls(joinedFinal.rowCount,
    joinedFinal.col(`${AE_TERM}_${ACTIVE_ARM_POSTTFIX}`),
    joinedFinal.col(`${AE_TERM}_${PLACEBO_ARM_POSTTFIX}`),
    joinedFinal.col(`${PERCENT}_${ACTIVE_ARM_POSTTFIX}`),
    joinedFinal.col(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`));

  removeNulls(joinedFinal.rowCount,
    joinedFinal.col(`${AE_TERM}_${PLACEBO_ARM_POSTTFIX}`),
    joinedFinal.col(`${AE_TERM}_${ACTIVE_ARM_POSTTFIX}`),
    joinedFinal.col(`${PERCENT}_${PLACEBO_ARM_POSTTFIX}`),
    joinedFinal.col(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`));

  joinedFinal.columns.addNewFloat(RELATIVE_RISK)
    .init((i) => (joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, i) * numberSubjPlaceboArm) / (joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, i) * numberSubjActiveArm));

  joinedFinal.columns.addNewFloat(RISK_DIFFERENCE)
    .init((i) => joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, i) / numberSubjActiveArm - joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, i) / numberSubjPlaceboArm);

  joinedFinal.columns.addNewFloat(ODDS_RATIO).init((i) => {
    const TE = joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, i);
    const PE = joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, i);
    return (TE * (numberSubjPlaceboArm - PE)) / (PE * (numberSubjActiveArm - TE))
  });

  joinedFinal.columns.addNewFloat('Z').init((i) => {
    const TE = joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, i);
    const TN = numberSubjActiveArm - TE;
    const PE = joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, i);
    const PN = numberSubjPlaceboArm - PE;
    const N = totalSubj;
    const Z = (Math.pow(TE - (numberSubjActiveArm * (TE + PE)) / N, 2) * Math.pow(N, 2) * (N - 1)) /
      (numberSubjActiveArm * numberSubjPlaceboArm * (TE + PE) * (TN + PN));
    return Z;
  });

  joinedFinal.columns.addNewFloat(P_VALUE).init((i) => {
    const Z = joinedFinal.get('Z', i);
    const f = jStat.chisquare.cdf(Z, 1);
    return 1 - f < f ? 1 - f : f;
  });

  joinedFinal.columns.addNewFloat(NEG_LOG10_P_VALUE).init((i) => (-Math.log10(joinedFinal.get(P_VALUE, i))));

  joinedFinal.columns.addNewFloat(STANDARD_ERROR_RD).init((i) => {
    const TE = joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, i);
    const TN = numberSubjActiveArm - TE;
    const PE = joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, i);
    const PN = numberSubjPlaceboArm - PE;
    const SE = Math.sqrt((TE * TN / Math.pow(numberSubjActiveArm, 3) + PE * PN / Math.pow(numberSubjPlaceboArm, 3)));
    return SE;
  });

  joinedFinal.columns.addNewFloat(SE_RD_WITH_SIGN_LEVEL).init((i) => {
    const SE = joinedFinal.get(STANDARD_ERROR_RD, i);
    const quantile = jStat.normal.inv(1 - signLevel / 2, 0, 1);
    return 2 * quantile * SE;
  });

  const aePercent = ae.groupBy([AE_TERM]).count().aggregate();
  const totalAeEvents = ae.rowCount;
  aePercent.columns.addNewFloat(AE_PERCENT).init((i) => {
    return aePercent.get(`count`, i) / totalAeEvents;
  });
  aePercent.columns.remove(`count`);

  joinedFinal = grok.data.joinTables(joinedFinal, aePercent, [`${AE_TERM}_${ACTIVE_ARM_POSTTFIX}`], [AE_TERM], joinedFinal.columns.names(), [AE_PERCENT], DG.JOIN_TYPE.LEFT, false);

  return joinedFinal;
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
  let grouped = dataframe.groupBy([SUBJECT_ID, testCol, resCol])
    .where(`${blVisitColumn} = ${baselineVisitNum}`)
    .aggregate();
  grouped.getCol(resCol).name = `BL_${resCol}`;
  dataframe.join(grouped,
    [SUBJECT_ID, testCol], [SUBJECT_ID, testCol],
    dataframe.columns.names(), [`BL_${resCol}`], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => {
      if (dataframe.getCol(resCol).isNone(i) || dataframe.getCol(`BL_${resCol}`).isNone(i)) {
        return null;
      } else {
        return (dataframe.get(resCol, i) - dataframe.get(`BL_${resCol}`, i)) / dataframe.get(`BL_${resCol}`, i);
      }
    });
  if (renameNewCol) {
    dataframe.columns.remove(resCol);
    dataframe.getCol(newColName).name = resCol;
  }
  dataframe.name = dfName;
}

export function labDynamicComparedToMinMax(dataframe, newColName: string) {
  const dfName = dataframe.name;
  let groupedMax = dataframe.groupBy([LAB_TEST])
    .max(LAB_RES_N)
    .aggregate();
  let groupedMin = dataframe.groupBy([LAB_TEST])
    .min(LAB_RES_N)
    .aggregate();
  dataframe.join(groupedMax, [LAB_TEST], [LAB_TEST], dataframe.columns.names(), [`max(${LAB_RES_N})`], DG.JOIN_TYPE.LEFT, true)
    .join(groupedMin, [LAB_TEST], [LAB_TEST], dataframe.columns.names(), [`min(${LAB_RES_N})`], DG.JOIN_TYPE.LEFT, true);
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
  const subjsPerDay = df.groupBy([dateCol]).uniqueCount(subjIDCol).aggregate();
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

export function getSubjectDmData(subjId: string, columnsToReturn: string[]) {
  if (study.domains.dm) {
    const dmData = study.domains.dm
      .groupBy(study.domains.dm.columns.names())
      .where(`${SUBJECT_ID} = ${subjId}`)
      .aggregate();
    const dmDict = {};
    columnsToReturn.forEach(it => {
      if (dmData.columns.names().includes(it)) {
        dmDict[it] = dmData.get(it, 0);
      }
    });
    return dmDict;
  }
  return {};
}


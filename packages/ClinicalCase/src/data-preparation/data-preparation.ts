import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ACTIVE_ARM_POSTTFIX, AE_PERCENT, ALT, AST, BILIRUBIN, DOMAINS_WITH_EVENT_START_END_DAYS,
  NEG_LOG10_P_VALUE, ODDS_RATIO, PLACEBO_ARM_POSTTFIX, P_VALUE, RELATIVE_RISK, RISK_DIFFERENCE,
  SE_RD_WITH_SIGN_LEVEL, STANDARD_ERROR_RD} from '../constants/constants';
import {addDataFromDmDomain, dateDifferenceInDays, filterBooleanColumn, filterNulls} from './utils';
import {ACT_TRT_ARM, AE_CAUSALITY, AE_REQ_HOSP, AE_SEQ, AE_SEVERITY, AE_START_DATE, LAB_DAY, LAB_HI_LIM_N,
  LAB_LO_LIM_N, LAB_RES_N, LAB_TEST, PLANNED_TRT_ARM, SUBJECT_ID, SUBJ_REF_ENDT, SUBJ_REF_STDT, VISIT_DAY,
  VISIT, VISIT_START_DATE, VISIT_DAY_STR, LAB_TEST_CAT} from '../constants/columns-constants';
import {studies} from '../utils/app-utils';
const {jStat} = require('jstat');


export function createMaxValuesDataForHysLaw(dataframe, aggregatedColName, filerValue) {
  let missingDataFlag = false;
  const condition = `${LAB_TEST} = ${filerValue}`;
  let grouped = dataframe.groupBy([SUBJECT_ID, LAB_RES_N, LAB_HI_LIM_N])
    .where(condition)
    .aggregate();
  grouped.columns.addNewFloat(aggregatedColName)
    .init((i) => {
      if (!grouped.col(LAB_RES_N).isNone(i) && !grouped.col(LAB_HI_LIM_N).isNone(i))
        return grouped.get(LAB_RES_N, i) / grouped.get(LAB_HI_LIM_N, i);
      else {
        missingDataFlag = true;
        return 0;
      }
    });
  grouped = grouped.groupBy([SUBJECT_ID])
    .max(aggregatedColName)
    .aggregate();
  grouped.getCol(`max(${aggregatedColName})`).name = aggregatedColName;
  return {df: grouped, missingData: missingDataFlag};
}


export function createHysLawDataframe(lb: DG.DataFrame, dm: DG.DataFrame, altName: string,
  astName: string, blnName: string, trtArm: string) {
  const {df: alt, missingData: missingDataAlt} = createMaxValuesDataForHysLaw(lb, ALT, altName);
  const {df: bln, missingData: missingDataBln} = createMaxValuesDataForHysLaw(lb, BILIRUBIN, blnName);
  const {df: ast, missingData: missingDataAst} = createMaxValuesDataForHysLaw(lb, AST, astName);
  //let ap = createMaxValuesDataForHysLaw(lb, AP, 'Alkaline Phosphatase');

  let joined = grok.data.joinTables(bln, alt, [SUBJECT_ID], [SUBJECT_ID],
    [SUBJECT_ID, BILIRUBIN], [ALT], DG.JOIN_TYPE.LEFT, false);
  joined = grok.data.joinTables(joined, ast, [SUBJECT_ID], [SUBJECT_ID],
    [SUBJECT_ID, ALT, BILIRUBIN], [AST], DG.JOIN_TYPE.LEFT, false);
  //joined = grok.data.joinTables(joined, ap, [ SUBJECT_ID ], [ SUBJECT_ID ],
  // [ SUBJECT_ID, ALT, BILIRUBIN, AST ], [ AP ], DG.JOIN_TYPE.LEFT, false);
  if (dm && dm.col(trtArm))
    joined = addDataFromDmDomain(joined, dm, [SUBJECT_ID, ALT, BILIRUBIN, AST], [trtArm]);

  if (missingDataAlt || missingDataAst || missingDataBln) {
    grok.shell.error(`Some values were missing in ${LAB_RES_N} and/or ${LAB_HI_LIM_N}.
      The results of Hy's law should be reviewed manually.`);
  }

  return joined;
}


export function createFilteredFloatValuesDataframe(df: DG.DataFrame & any, condition: string, groupCols: string[],
  newFloatCol: string, colToTransform: string) {
  const filteredData = df.groupBy(groupCols)
    .where(condition)
    .aggregate();
  filteredData.columns.addNewFloat(newFloatCol).init((i) => parseFloat(filteredData.get(colToTransform, i)));
  return filteredData;
}

export function createFilteredDataframe(df: DG.DataFrame & any, condition: string, groupCols: string[],
  newFloatCol: string, colToTransform: string) {
  const filteredData = df.groupBy(groupCols)
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
  epNumColumn = null) {
  const condition = `${testCol} = ${value} and ${visitCol} IN`;
  const filteredDataBaseline = createFilteredDataframe(df, `${condition} (${blVisit})`,
    [SUBJECT_ID, resCol, visitCol, testCol].concat(additionalCols), blNumColumn, resCol);
  let finalDf;
  let columnsToEXtract = [];
  if (epVisit) {
    const filteredDataEndpoint = createFilteredDataframe(df, `${condition} (${epVisit})`,
      [SUBJECT_ID, resCol, visitCol, testCol].concat(additionalCols), epNumColumn, resCol);
    finalDf = grok.data.joinTables(filteredDataBaseline, filteredDataEndpoint,
      [SUBJECT_ID], [SUBJECT_ID], [SUBJECT_ID, blNumColumn, testCol].concat(additionalCols), [epNumColumn],
      DG.JOIN_TYPE.LEFT, false);
    columnsToEXtract = [SUBJECT_ID, blNumColumn, epNumColumn, testCol].concat(additionalCols);
  } else {
    finalDf = filteredDataBaseline;
    columnsToEXtract = [SUBJECT_ID, blNumColumn, testCol].concat(additionalCols);
  }
  if (dm && columnToExtractFromDm.every((it) => dm.columns.names().includes(it)))
    finalDf = addDataFromDmDomain(finalDf, dm, columnsToEXtract, columnToExtractFromDm);

  return finalDf;
}


export function createLabValuesByVisitDataframe(lb: DG.DataFrame, dm: DG.DataFrame,
  labValue: string,
  treatmentArmCol: string,
  treatmentArm: string,
  labValueNumCol: string,
  visitCol: string) {
  const condition = `${LAB_TEST} = ${labValue}`;
  let filtered = createFilteredDataframe(lb, condition, [SUBJECT_ID, LAB_RES_N, visitCol], labValueNumCol, LAB_RES_N);
  if (dm && dm.col(treatmentArmCol)) {
    filtered = addDataFromDmDomain(filtered, dm, filtered.columns.names(), [treatmentArmCol]);
    if (treatmentArm)
      filtered.rows.filter((row) => row.visitdy !== DG.INT_NULL && row[treatmentArmCol] === treatmentArm);
  } else
    filtered.rows.filter((row) => row.visitdy !== DG.INT_NULL);

  return filtered;
}


export function createUniqueCountDataframe(df: DG.DataFrame, groupCol: string[],
  uniqueCountCol: string, renameCol: string) {
  const grouped = df.groupBy(groupCol)
    .uniqueCount(uniqueCountCol)
    .aggregate();
  grouped.getCol(`unique(${uniqueCountCol})`).name = renameCol;
  return grouped;
}


export function addColumnWithDrugPlusDosage(df: DG.DataFrame, drugCol: string, dosageCol: string,
  unitsCol: string, newCol: string) {
  df.columns.addNewString(newCol)
    .init((i) => `${df.col(drugCol) ? df.get(drugCol, i).toString() : ''} ${df.col(dosageCol) ?
      df.get(dosageCol, i).toString() : ''}${df.col(unitsCol) ? df.get(unitsCol, i).toString() : ''}`);
  return df;
}

export function createSurvivalData(dmDf: DG.DataFrame, aeDf: DG.DataFrame,
  endpoint: string, SDTMendpoint: string, covariates: string[]) {
  let dm = dmDf.clone();
  filterNulls(dm, SUBJ_REF_ENDT); // according to SDTMIG v3.3 RFENDTC is null for screen failures or unassigned subjects
  if (SDTMendpoint === AE_START_DATE) {
    const ae = aeDf.clone();
    if (endpoint == 'HOSPITALIZATION')
      filterBooleanColumn(ae, AE_REQ_HOSP, false);

    // eslint-disable-next-line max-len
    const condition = endpoint == 'DRUG RELATED AE' ? `${AE_CAUSALITY} in (PROBABLE, POSSIBLE, RELATED, UNLIKELY RELATED, POSSIBLY RELATED, RELATED)` : `${AE_SEVERITY} = SEVERE`;
    const aeGrouped = ae.groupBy([SUBJECT_ID]).
      min(AE_SEQ).
      where(condition).
      aggregate();
    const aeJoined = grok.data.joinTables(ae, aeGrouped, [SUBJECT_ID, AE_SEQ],
      [SUBJECT_ID, `min(${AE_SEQ})`], [SUBJECT_ID, AE_START_DATE], [`min(${AE_SEQ})`], DG.JOIN_TYPE.LEFT, false);
    filterNulls(aeJoined, `min(${AE_SEQ})`);
    dm = grok.data.joinTables(dm, aeJoined, [SUBJECT_ID], [SUBJECT_ID], dm.columns.names(),
      [AE_START_DATE], DG.JOIN_TYPE.LEFT, false);
  }
  dm.columns.addNewInt('time')
    .init((i) => getSurvivalTime(dm.columns.byName(SDTMendpoint),
      dm.columns.byName(SUBJ_REF_STDT), dm.columns.byName(SUBJ_REF_ENDT), i));
  dm.columns.addNewInt('status')
    .init((i) => getSurvivalStatus(dm.columns.byName(SDTMendpoint), i));
  return dm.groupBy([SUBJECT_ID, 'time', 'status'].concat(covariates)).aggregate();
}

export function getSurvivalTime(eventColumn: DG.Column, startColumn: DG.Column, endColumn: DG.Column, i: number) {
  if (eventColumn.isNone(i))
    return dateDifferenceInDays(startColumn.get(i).toString(), endColumn.get(i).toString());
  else
    return dateDifferenceInDays(startColumn.get(i).toString(), eventColumn.get(i).toString());
}

export function getSurvivalStatus(eventColumn: DG.Column, i: number) {
  return eventColumn.isNone(i) ? 0 : 1;
}

export function createAERiskAssessmentDataframe(ae: DG.DataFrame, dm: DG.DataFrame, trtArmCol: string, aeTermCol:
  string, placeboArm: string[], activeArm: string[], signLevel: number) {
  const TRTARM = 'TREATMENT_ARM';
  const TOTALSUBJ = 'TOTALSUBJ';
  const TOTALSUBJ_WITH_AE = 'TOTALSUBJ_WITH_AE';
  const PERCENT = 'PERCENT';

  dm.columns.addNewString(TRTARM).init((i) => {
    if (placeboArm.includes(dm.get(trtArmCol, i)))
      return placeboArm.join(',');

    if (activeArm.includes(dm.get(trtArmCol, i)))
      return activeArm.join(',');

    return dm.get(trtArmCol, i);
  });

  const subjArm = createUniqueCountDataframe(dm, [TRTARM], SUBJECT_ID, TOTALSUBJ);

  const numberSubjPlaceboArm = subjArm.groupBy(subjArm.columns.names())
    .where({[TRTARM]: `${placeboArm.join(',')}`})
    .aggregate()
    .get(TOTALSUBJ, 0);
  const numberSubjActiveArm = subjArm.groupBy(subjArm.columns.names())
    .where({[TRTARM]: `${activeArm.join(',')}`})
    .aggregate()
    .get(TOTALSUBJ, 0);
  const totalSubj = numberSubjPlaceboArm + numberSubjActiveArm;

  const joinedAeEX = grok.data.joinTables(ae, dm, [SUBJECT_ID], [SUBJECT_ID],
    [SUBJECT_ID, aeTermCol], [TRTARM], DG.JOIN_TYPE.LEFT, false);
  const subjAE = createUniqueCountDataframe(joinedAeEX, [aeTermCol, TRTARM], SUBJECT_ID, TOTALSUBJ_WITH_AE);

  const subjTotalAndWithAE = grok.data.joinTables(subjAE, subjArm, [TRTARM], [TRTARM],
    [aeTermCol, TRTARM, TOTALSUBJ_WITH_AE], [TOTALSUBJ], DG.JOIN_TYPE.LEFT, false);

  subjTotalAndWithAE.columns.addNewFloat(PERCENT)
    .init((i) => parseFloat(subjTotalAndWithAE.get(TOTALSUBJ_WITH_AE, i)) /
    parseFloat(subjTotalAndWithAE.get(TOTALSUBJ, i)));

  const aeRiskRatioPlacebo = subjTotalAndWithAE.groupBy([aeTermCol, TRTARM, TOTALSUBJ_WITH_AE, TOTALSUBJ, PERCENT])
    .where({[TRTARM]: `${placeboArm.join(',')}`})
    .aggregate();
  aeRiskRatioPlacebo.name = `${placeboArm.join(',')}`;

  const aeRiskRatioActive = subjTotalAndWithAE.groupBy([aeTermCol, TRTARM, TOTALSUBJ_WITH_AE, TOTALSUBJ, PERCENT])
    .where({[TRTARM]: `${activeArm.join(',')}`})
    .aggregate();
  aeRiskRatioActive.name = `${activeArm.join(',')}`;

  const renameCols = (df: DG.DataFrame, columnNames: string[], arm: string) => {
    columnNames.forEach((it) => {
      df.getCol(it).name = `${it}_${arm}`;
    });
  };

  renameCols(aeRiskRatioPlacebo, [TRTARM, TOTALSUBJ_WITH_AE, TOTALSUBJ, PERCENT], PLACEBO_ARM_POSTTFIX);
  renameCols(aeRiskRatioActive, [TRTARM, TOTALSUBJ_WITH_AE, TOTALSUBJ, PERCENT], ACTIVE_ARM_POSTTFIX);

  let joinedFinal = grok.data.joinTables(aeRiskRatioActive, aeRiskRatioPlacebo, [aeTermCol], [aeTermCol],
    [aeTermCol, `${TRTARM}_${ACTIVE_ARM_POSTTFIX}`, `${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`,
      `${TOTALSUBJ}_${ACTIVE_ARM_POSTTFIX}`, `${PERCENT}_${ACTIVE_ARM_POSTTFIX}`],
    [aeTermCol, `${TRTARM}_${PLACEBO_ARM_POSTTFIX}`, `${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`,
      `${TOTALSUBJ}_${PLACEBO_ARM_POSTTFIX}`, `${PERCENT}_${PLACEBO_ARM_POSTTFIX}`],
    DG.JOIN_TYPE.OUTER, false);


  const aeTermActCol = joinedFinal.col(`${activeArm.join(',')}.${aeTermCol}`);
  joinedFinal.col(aeTermActCol ? `${activeArm.join(',')}.${aeTermCol}` : aeTermCol).name =
    `${aeTermCol}_${ACTIVE_ARM_POSTTFIX}`;
  const aeTermPlaceboCol = joinedFinal.col(`${placeboArm.join(',')}.${aeTermCol}`);
  joinedFinal.col(aeTermPlaceboCol ? `${placeboArm.join(',')}.${aeTermCol}`: aeTermCol).name =
    `${aeTermCol}_${PLACEBO_ARM_POSTTFIX}`;


  const removeNulls = (rowCount, colToCheck, colWithWalue, percentCol, totalSubjWithAECol) => {
    for (let i = 0; i < rowCount; i++) {
      if (colToCheck.isNone(i)) {
        const value = colWithWalue.get(i);
        colToCheck.set(i, value, false);
        percentCol.set(i, 0, false);
        totalSubjWithAECol.set(i, 0, false);
      }
    }
  };

  removeNulls(joinedFinal.rowCount,
    joinedFinal.col(`${aeTermCol}_${ACTIVE_ARM_POSTTFIX}`),
    joinedFinal.col(`${aeTermCol}_${PLACEBO_ARM_POSTTFIX}`),
    joinedFinal.col(`${PERCENT}_${ACTIVE_ARM_POSTTFIX}`),
    joinedFinal.col(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`));

  removeNulls(joinedFinal.rowCount,
    joinedFinal.col(`${aeTermCol}_${PLACEBO_ARM_POSTTFIX}`),
    joinedFinal.col(`${aeTermCol}_${ACTIVE_ARM_POSTTFIX}`),
    joinedFinal.col(`${PERCENT}_${PLACEBO_ARM_POSTTFIX}`),
    joinedFinal.col(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`));

  joinedFinal.columns.addNewFloat(RELATIVE_RISK)
    .init((i) => (joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, i) * numberSubjPlaceboArm) /
    (joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, i) * numberSubjActiveArm));

  joinedFinal.columns.addNewFloat(RISK_DIFFERENCE)
    .init((i) => joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, i) /
    numberSubjActiveArm - joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, i) / numberSubjPlaceboArm);

  joinedFinal.columns.addNewFloat(ODDS_RATIO).init((i) => {
    const TE = joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${ACTIVE_ARM_POSTTFIX}`, i);
    const PE = joinedFinal.get(`${TOTALSUBJ_WITH_AE}_${PLACEBO_ARM_POSTTFIX}`, i);
    return (TE * (numberSubjPlaceboArm - PE)) / (PE * (numberSubjActiveArm - TE));
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

  const aePercent = ae.groupBy([aeTermCol]).count().aggregate();
  const totalAeEvents = ae.rowCount;
  aePercent.columns.addNewFloat(AE_PERCENT).init((i) => {
    return aePercent.get(`count`, i) / totalAeEvents;
  });
  aePercent.columns.remove(`count`);

  joinedFinal = grok.data.joinTables(joinedFinal, aePercent, [`${aeTermCol}_${ACTIVE_ARM_POSTTFIX}`],
    [aeTermCol], joinedFinal.columns.names(), [AE_PERCENT], DG.JOIN_TYPE.LEFT, false);

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
  const grouped = dataframe.groupBy([SUBJECT_ID, testCol, resCol])
    .where(`${blVisitColumn} = ${baselineVisitNum}`)
    .aggregate();
  grouped.getCol(resCol).name = `BL_${resCol}`;
  dataframe.join(grouped,
    [SUBJECT_ID, testCol], [SUBJECT_ID, testCol],
    dataframe.columns.names(), [`BL_${resCol}`], DG.JOIN_TYPE.LEFT, true);
  dataframe.columns.addNewFloat(newColName)
    .init((i) => {
      if (dataframe.getCol(resCol).isNone(i) || dataframe.getCol(`BL_${resCol}`).isNone(i))
        return null;
      else
        return (dataframe.get(resCol, i) - dataframe.get(`BL_${resCol}`, i)) / dataframe.get(`BL_${resCol}`, i);
    });
  if (renameNewCol) {
    dataframe.columns.remove(resCol);
    dataframe.getCol(newColName).name = resCol;
  }
  dataframe.name = dfName;
}

export function labDynamicComparedToMinMax(dataframe, newColName: string) {
  const dfName = dataframe.name;
  const groupedMax = dataframe.groupBy([LAB_TEST])
    .max(LAB_RES_N)
    .aggregate();
  const groupedMin = dataframe.groupBy([LAB_TEST])
    .min(LAB_RES_N)
    .aggregate();
  dataframe.join(groupedMax, [LAB_TEST], [LAB_TEST], dataframe.columns.names(),
    [`max(${LAB_RES_N})`], DG.JOIN_TYPE.LEFT, true)
    .join(groupedMin, [LAB_TEST], [LAB_TEST], dataframe.columns.names(),
      [`min(${LAB_RES_N})`], DG.JOIN_TYPE.LEFT, true);
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
      if (!val || !min || !max)
        return null;
      else
        return val >= max ? val / max : val <= min ? 0 - min / val : ((val - min) / (max - min) - 0.5) * 2;
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
      if (i > 0)
        return subjsPerDay.get(cumCol, i - 1) + subjsPerDay.get(`unique(${subjIDCol})`, i);

      return subjsPerDay.get(`unique(${subjIDCol})`, i);
    });
  subjsPerDay.columns.remove(`unique(${subjIDCol})`);
  return subjsPerDay;
}

export function getSubjectDmData(subjId: string, columnsToReturn: string[], studyId: string) {
  if (studies[studyId].domains.dm) {
    const dmData = studies[studyId].domains.dm
      .groupBy(studies[studyId].domains.dm.columns.names())
      .where(`${SUBJECT_ID} = ${subjId}`)
      .aggregate();
    const dmDict = {};
    columnsToReturn.forEach((it) => {
      if (dmData.columns.names().includes(it))
        dmDict[it] = dmData.get(it, 0);
    });
    return dmDict;
  }
  return {};
}


export function createEventStartEndDaysCol(studyId: string) {
  if (studies[studyId].domains.sv != null) {
    if (studies[studyId].domains.sv.col(VISIT_START_DATE) &&
      studies[studyId].domains.sv.col(VISIT_DAY) && studies[studyId].domains.sv.col(SUBJECT_ID)) {
      const baselineDates = getSubjectBaselineDates(studyId);
      DOMAINS_WITH_EVENT_START_END_DAYS.forEach((domain) => {
        if (studies[studyId].domains[domain]) {
          addCalculatedStudyDayColumn(domain, 'ST', baselineDates, studyId);
          addCalculatedStudyDayColumn(domain, 'EN', baselineDates, studyId);
        };
      });
    };
  };
}

export function addVisitDayFromTvDomain(studyId: string) {
  if (studies[studyId].domains.tv != null && studies[studyId].domains.tv.col(VISIT_DAY)) {
    const visitNamesAndDays = studies[studyId].domains.tv.groupBy([VISIT, VISIT_DAY]).aggregate();
    studies[studyId].domains.all().forEach((domain) => {
      if (domain.name !== 'tv' && domain.col(VISIT) && !domain.col(VISIT_DAY)) {
        const tableName = domain.name;
        grok.data.joinTables(domain, visitNamesAndDays, [VISIT],
          [VISIT], domain.columns.names(), [VISIT_DAY], DG.JOIN_TYPE.LEFT, true);
        domain.name = tableName;
      }
    });
  };
}

export function addCalculatedStudyDayColumn(domain: string, prefix: string,
  baselineDates: DG.DataFrame, studyId: string) {
  if (studies[studyId].domains[domain].col(`${domain.toUpperCase()}${prefix}DTC`)) {
    grok.data.joinTables(studies[studyId].domains[domain],
      baselineDates, [SUBJECT_ID], [SUBJECT_ID], studies[studyId].domains[domain].columns.names(),
      [VISIT_START_DATE], DG.JOIN_TYPE.LEFT, true);
    studies[studyId].domains[domain].columns.addNewInt(`${domain.toUpperCase()}${prefix}DY_CALCULATED`).init((i) => {
      const baselineDate = studies[studyId].domains[domain].get(VISIT_START_DATE, i); ;
      const startDate = studies[studyId].domains[domain].get(`${domain.toUpperCase()}${prefix}DTC`, i);
      let dateDiff = dateDifferenceInDays(baselineDate, startDate);
      if (dateDiff >= 0)
        dateDiff += 1;

      return dateDiff;
    });
    studies[studyId].domains[domain].columns.remove(VISIT_START_DATE);
    studies[studyId].domains[domain].name = domain;
  };
}

export function getSubjectBaselineDates(studyId: string) {
  const subjBaselineDates = studies[studyId].domains.sv.groupBy([SUBJECT_ID, VISIT_START_DATE, VISIT_DAY])
    .where(`${VISIT_DAY} = 1`)
    .aggregate();
  return subjBaselineDates;
}

export function createVisitDayStrCol(df: DG.DataFrame, visitColNamesDict?: {[key: string]: string}) {
  if (!df)
    return;

  const getVisitDayCol = () => {
    const visitDayCol = df.col(VISIT_DAY);
    //TODO!!! For some interval related tests, like body weight gain, there are two columns:
    // BGDY(start of interval) and BGENDY (end of interval) - need to decide which to use and in which cases
    const domainSpecificVisitDayCol = df.col(`${df.name.toUpperCase()}DY`);
    const colToUse = visitDayCol ?? domainSpecificVisitDayCol;
    return colToUse;
  };
  const visitDayCol = getVisitDayCol();
  if (!df.col(VISIT_DAY_STR) && visitDayCol) {
    //create categorical visit day column
    df.columns.addNewString(VISIT_DAY_STR)
      .init((i) => visitDayCol.isNone(i) ? undefined : visitDayCol.get(i).toString());
  }
  //set visit day column in case it is domain specific
  if (visitColNamesDict) {
    if (visitDayCol)
      visitColNamesDict[df.name] = visitDayCol.name;
    else
      delete visitColNamesDict[df.name];
  }
}

/**
 * Calculates baseline-related and post-baseline summary columns for LB domain:
 * Baseline columns:
 * - LB_BASELINE: The LBSTRESN value at the minimum LBDY for each subject and test
 * - LB_CHG: Change from baseline (LBSTRESN - LB_BASELINE)
 * - LB_PCT_CHG: Percent change from baseline (100 * (LBSTRESN - LB_BASELINE) / LB_BASELINE)
 * Post-baseline columns:
 * - MAX_POST_VALUE: max LBSTRESN post-baseline per subject/test
 * - MIN_PCT_CHG: min LB_PCT_CHG post-baseline per subject/test
 * - MAX_PCT_CHG: max LB_PCT_CHG post-baseline per subject/test
 * @param {DG.DataFrame} lbDomain - The LB domain DataFrame
 */
export function calculateLBBaselineColumns(lbDomain: DG.DataFrame): void {
  if (!lbDomain || lbDomain.name.toLowerCase() !== 'lb')
    return;

  const subjectIdCol = lbDomain.col(SUBJECT_ID);
  const labDayCol = lbDomain.col(LAB_DAY);
  const labCatCol = lbDomain.col(LAB_TEST_CAT);
  const labTestCol = lbDomain.col(LAB_TEST);
  const labResNCol = lbDomain.col(LAB_RES_N);

  // Check if required columns exist
  if (!subjectIdCol || !labDayCol || !labTestCol || !labResNCol)
    return;

  // Column names
  const baselineColName = 'LB_BASELINE';
  const chgColName = 'LB_CHG';
  const pctChgColName = 'LB_PCT_CHG';
  const maxPostValueColName = 'MAX_POST_VALUE';
  const minPctChgColName = 'MIN_PCT_CHG';
  const maxPctChgColName = 'MAX_PCT_CHG';

  const removeColIfExists = (colName: string) => {
    if (lbDomain.col(colName))
      lbDomain.columns.remove(colName);
  };
  // Remove existing columns if they exist (to allow recalculation)
  removeColIfExists(baselineColName);
  removeColIfExists(chgColName);
  removeColIfExists(pctChgColName);
  removeColIfExists(maxPostValueColName);
  removeColIfExists(minPctChgColName);
  removeColIfExists(maxPctChgColName);

  // Find baseline day (minimum LBDY) for each subject/test - calculate once
  const groupByCols = labCatCol ? [SUBJECT_ID, LAB_TEST_CAT, LAB_TEST] : [SUBJECT_ID, LAB_TEST];
  const minDayGrouped = lbDomain.groupBy(groupByCols)
    .min(LAB_DAY)
    .aggregate();

  // Join back to get the LBSTRESN value at the minimum day
  const baselineData = grok.data.joinTables(
    minDayGrouped,
    lbDomain,
    groupByCols.concat([`min(${LAB_DAY})`]),
    groupByCols.concat([LAB_DAY]),
    groupByCols.concat([`min(${LAB_DAY})`]),
    [LAB_RES_N],
    DG.JOIN_TYPE.LEFT,
    false,
  );

  // Build maps: baseline value and baseline day (calculated once, used for both)
  const baselineMap: {[key: string]: number} = {};
  const baselineDayMap: {[key: string]: number} = {};
  const baselineSubjCol = baselineData.col(SUBJECT_ID);
  const baselineCatCol = baselineData.col(LAB_TEST_CAT);
  const baselineTestCol = baselineData.col(LAB_TEST);
  const baselineResCol = baselineData.col(LAB_RES_N);
  const minDayColName = `min(${LAB_DAY})`;
  const baselineDayCol = baselineData.col(minDayColName);

  for (let i = 0; i < baselineData.rowCount; i++) {
    if (baselineResCol.isNone(i) || baselineDayCol.isNone(i))
      continue;

    const subjectId = baselineSubjCol.get(i);
    const cat = baselineCatCol ? baselineCatCol.get(i) : '';
    const test = baselineTestCol.get(i);
    const key = labCatCol ? `${subjectId}|${cat}|${test}` : `${subjectId}|${test}`;
    const baselineValue = baselineResCol.get(i);
    const baselineDay = baselineDayCol.get(i);

    // Store baseline value and day (if multiple rows match, first one wins)
    if (!baselineMap.hasOwnProperty(key)) {
      baselineMap[key] = baselineValue;
      baselineDayMap[key] = baselineDay;
    }
  }

  // Helper function to build key
  const buildKey = (subjectId: string, cat: string | null, test: string): string => {
    return labCatCol && cat !== null ? `${subjectId}|${cat}|${test}` : `${subjectId}|${test}`;
  };

  // Helper function to get subject/test info and build key for a row
  const getRowKey = (i: number): {subjectId: string, cat: string | null, test: string, key: string} | null => {
    if (labDayCol.isNone(i) || labResNCol.isNone(i))
      return null;

    const subjectId = subjectIdCol.get(i);
    const cat = labCatCol ? labCatCol.get(i) : null;
    const test = labTestCol.get(i);
    const key = buildKey(subjectId, cat, test);
    return {subjectId, cat, test, key};
  };

  // Collect post-baseline data while iterating (for post-baseline columns)
  interface PostBaselineData {
    maxValue: number | null;
    minPctChg: number | null;
    maxPctChg: number | null;
  }
  const postBaselineData: {[key: string]: PostBaselineData} = {};

  // Create baseline column
  const baselineCol = lbDomain.columns.addNewFloat(baselineColName);
  baselineCol.init((i) => {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      return null;
    return baselineMap[rowInfo.key] ?? null;
  });

  // Create change from baseline column
  const chgCol = lbDomain.columns.addNewFloat(chgColName);
  chgCol.init((i) => {
    if (labResNCol.isNone(i) || baselineCol.isNone(i))
      return null;

    const currentValue = labResNCol.get(i);
    const baselineValue = baselineCol.get(i);

    if (baselineValue === null || isNaN(baselineValue))
      return null;

    return currentValue - baselineValue;
  });

  // Create percent change from baseline column
  const pctChgCol = lbDomain.columns.addNewFloat(pctChgColName);
  pctChgCol.init((i) => {
    if (labResNCol.isNone(i) || baselineCol.isNone(i))
      return null;

    const currentValue = labResNCol.get(i);
    const baselineValue = baselineCol.get(i);

    if (baselineValue === null || isNaN(baselineValue) || baselineValue === 0)
      return null;

    return 100 * (currentValue - baselineValue) / baselineValue;
  });

  // Now collect post-baseline data in a single pass
  for (let i = 0; i < lbDomain.rowCount; i++) {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      continue;

    const baselineDay = baselineDayMap[rowInfo.key];

    if (baselineDay === undefined)
      continue;

    const currentDay = labDayCol.get(i);
    if (currentDay <= baselineDay)
      continue; // Skip baseline and pre-baseline rows

    const currentValue = labResNCol.get(i);
    const currentPctChg = pctChgCol && !pctChgCol.isNone(i) ? pctChgCol.get(i) : null;

    if (!postBaselineData[rowInfo.key]) {
      postBaselineData[rowInfo.key] = {
        maxValue: currentValue,
        minPctChg: currentPctChg,
        maxPctChg: currentPctChg,
      };
    } else {
      const data = postBaselineData[rowInfo.key];
      if (currentValue > (data.maxValue ?? -Infinity))
        data.maxValue = currentValue;
      if (currentPctChg !== null) {
        if (data.minPctChg === null || currentPctChg < data.minPctChg)
          data.minPctChg = currentPctChg;
        if (data.maxPctChg === null || currentPctChg > data.maxPctChg)
          data.maxPctChg = currentPctChg;
      }
    }
  }

  // Create MAX_POST_VALUE column
  const maxPostValueCol = lbDomain.columns.addNewFloat(maxPostValueColName);
  maxPostValueCol.init((i) => {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      return null;
    return postBaselineData[rowInfo.key]?.maxValue ?? null;
  });

  // Create MIN_PCT_CHG column
  const minPctChgCol = lbDomain.columns.addNewFloat(minPctChgColName);
  minPctChgCol.init((i) => {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      return null;
    return postBaselineData[rowInfo.key]?.minPctChg ?? null;
  });

  // Create MAX_PCT_CHG column
  const maxPctChgCol = lbDomain.columns.addNewFloat(maxPctChgColName);
  maxPctChgCol.init((i) => {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      return null;
    return postBaselineData[rowInfo.key]?.maxPctChg ?? null;
  });
}

/**
 * Calculates control-related columns for LB domain:
 * - CONTROL_MEAN: mean of LBSTRESN for concurrent control per test & visit
 * - DELTA_VS_CONTROL: VALUE – CONTROL_MEAN
 * - PCT_VS_CONTROL: 100 * (VALUE – CONTROL_MEAN) / CONTROL_MEAN
 * @param {DG.DataFrame} lbDomain - The LB domain DataFrame
 * @param {Array<{treatment: string, control: string}>} treatmentControlConfig - Treatment/control pairs
 */
export function calculateLBControlColumns(
  lbDomain: DG.DataFrame,
  treatmentControlConfig: {treatment: string, control: string}[],
): void {
  if (!lbDomain || lbDomain.name.toLowerCase() !== 'lb')
    return;

  if (!treatmentControlConfig || treatmentControlConfig.length === 0)
    return;

  const labTestCol = lbDomain.col(LAB_TEST);
  const labCatCol = lbDomain.col(LAB_TEST_CAT);
  const labResNCol = lbDomain.col(LAB_RES_N);
  // Get treatment arm column (ARM or ACTARM, should be joined from DM)
  const armCol = lbDomain.col(PLANNED_TRT_ARM) || lbDomain.col(ACT_TRT_ARM);
  // Get visit column (prefer VISITDY, then VISIT, then VISITNUM)
  const visitCol = lbDomain.col(LAB_DAY);

  // Check if required columns exist
  if (!labTestCol || !labResNCol || !armCol || !visitCol)
    return;

  // Column names
  const controlMeanColName = 'CONTROL_MEAN';
  const deltaVsControlColName = 'DELTA_VS_CONTROL';
  const pctVsControlColName = 'PCT_VS_CONTROL';

  // Remove existing columns if they exist (to allow recalculation)
  if (lbDomain.col(controlMeanColName))
    lbDomain.columns.remove(controlMeanColName);
  if (lbDomain.col(deltaVsControlColName))
    lbDomain.columns.remove(deltaVsControlColName);
  if (lbDomain.col(pctVsControlColName))
    lbDomain.columns.remove(pctVsControlColName);

  // Build a map of treatment -> control arm
  const treatmentToControlMap: {[treatment: string]: string} = {};
  const controlArms = new Set<string>();
  for (const config of treatmentControlConfig) {
    treatmentToControlMap[config.treatment] = config.control;
    controlArms.add(config.control);
  }

  // Build group columns for aggregation (test + optional category)
  const testGroupCols = labCatCol ? [LAB_TEST_CAT, LAB_TEST] : [LAB_TEST];
  const visitColName = visitCol.name;

  // Calculate mean LBSTRESN for control subjects grouped by test & visit
  // First, filter to only control subjects
  const controlMeanMap: {[key: string]: number} = {};

  // Group by test (+ optional category) and visit, calculate mean for control subjects
  for (const controlArm of controlArms) {
    // Use groupBy with where clause to filter control subjects
    const controlFiltered = lbDomain.groupBy(testGroupCols.concat([visitColName, LAB_RES_N]))
      .where(`${armCol.name} = ${controlArm}`)
      .aggregate();

    if (controlFiltered.rowCount === 0)
      continue;

    // Now group by test + visit and calculate mean
    const meanGrouped = controlFiltered.groupBy(testGroupCols.concat([visitColName]))
      .avg(LAB_RES_N)
      .aggregate();

    const avgColName = `avg(${LAB_RES_N})`;
    const meanResCol = meanGrouped.col(avgColName);
    const meanVisitCol = meanGrouped.col(visitColName);
    const meanTestCol = meanGrouped.col(LAB_TEST);
    const meanCatCol = labCatCol ? meanGrouped.col(LAB_TEST_CAT) : null;

    if (!meanResCol || !meanVisitCol || !meanTestCol)
      continue;

    for (let i = 0; i < meanGrouped.rowCount; i++) {
      if (meanResCol.isNone(i))
        continue;

      const test = meanTestCol.get(i);
      const cat = meanCatCol ? meanCatCol.get(i) : '';
      const visit = meanVisitCol.get(i);
      const meanValue = meanResCol.get(i);

      // Key format: controlArm|category|test|visit (category may be empty)
      const key = `${controlArm}|${cat}|${test}|${visit}`;
      controlMeanMap[key] = meanValue;
    }
  }

  // Create CONTROL_MEAN column
  const controlMeanCol = lbDomain.columns.addNewFloat(controlMeanColName);
  controlMeanCol.init((i) => {
    if (labResNCol.isNone(i))
      return null;

    const arm = armCol.get(i);
    // Find control arm for this subject's arm
    let controlArm: string | null = null;
    if (controlArms.has(arm)) {
      // Subject is in control group, use their own arm
      controlArm = arm;
    } else if (treatmentToControlMap[arm]) {
      // Subject is in treatment group, use mapped control
      controlArm = treatmentToControlMap[arm];
    }

    if (!controlArm)
      return null;

    const test = labTestCol.get(i);
    const cat = labCatCol ? labCatCol.get(i) : '';
    const visit = visitCol.get(i);
    const key = `${controlArm}|${cat}|${test}|${visit}`;

    return controlMeanMap[key] ?? null;
  });

  // Create DELTA_VS_CONTROL column
  const deltaVsControlCol = lbDomain.columns.addNewFloat(deltaVsControlColName);
  deltaVsControlCol.init((i) => {
    if (labResNCol.isNone(i) || controlMeanCol.isNone(i))
      return null;

    const currentValue = labResNCol.get(i);
    const controlMean = controlMeanCol.get(i);

    if (controlMean === null || isNaN(controlMean))
      return null;

    return currentValue - controlMean;
  });

  // Create PCT_VS_CONTROL column
  const pctVsControlCol = lbDomain.columns.addNewFloat(pctVsControlColName);
  pctVsControlCol.init((i) => {
    if (labResNCol.isNone(i) || controlMeanCol.isNone(i))
      return null;

    const currentValue = labResNCol.get(i);
    const controlMean = controlMeanCol.get(i);

    if (controlMean === null || isNaN(controlMean) || controlMean === 0)
      return null;

    return 100 * (currentValue - controlMean) / controlMean;
  });
}

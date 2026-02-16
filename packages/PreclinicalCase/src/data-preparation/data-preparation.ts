import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ACT_TRT_ARM, LAB_DAY, LAB_RES_N, LAB_TEST, LAB_TEST_CAT, LAB_LO_LIM_N,
  LAB_HI_LIM_N, PLANNED_TRT_ARM, SUBJECT_ID, VISIT_DAY, VISIT_DAY_STR,
  DOMAIN} from '../constants/columns-constants';
import { studies } from '../utils/app-utils';

export function createVisitDayStrCol(df: DG.DataFrame, visitColNamesDict?: {[key: string]: string}) {
  if (!df)
    return;

  const getVisitDayCol = () => {
    const visitDayCol = df.col(VISIT_DAY);
    const domainSpecificVisitDayCol = df.col(`${df.name.toUpperCase()}DY`);
    const colToUse = visitDayCol ?? domainSpecificVisitDayCol;
    return colToUse;
  };
  const visitDayCol = getVisitDayCol();
  if (!df.col(VISIT_DAY_STR) && visitDayCol) {
    df.columns.addNewString(VISIT_DAY_STR)
      .init((i) => visitDayCol.isNone(i) ? undefined : visitDayCol.get(i).toString());
  }
  if (visitColNamesDict) {
    if (visitDayCol)
      visitColNamesDict[df.name] = visitDayCol.name;
    else
      delete visitColNamesDict[df.name];
  }
}

export function calculateLBBaselineColumns(lbDomain: DG.DataFrame): void {
  if (!lbDomain || lbDomain.name.toLowerCase() !== 'lb')
    return;

  const subjectIdCol = lbDomain.col(SUBJECT_ID);
  const labDayCol = lbDomain.col(LAB_DAY);
  const labCatCol = lbDomain.col(LAB_TEST_CAT);
  const labTestCol = lbDomain.col(LAB_TEST);
  const labResNCol = lbDomain.col(LAB_RES_N);

  if (!subjectIdCol || !labDayCol || !labTestCol || !labResNCol)
    return;

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
  removeColIfExists(baselineColName);
  removeColIfExists(chgColName);
  removeColIfExists(pctChgColName);
  removeColIfExists(maxPostValueColName);
  removeColIfExists(minPctChgColName);
  removeColIfExists(maxPctChgColName);

  const groupByCols = labCatCol ? [SUBJECT_ID, LAB_TEST_CAT, LAB_TEST] : [SUBJECT_ID, LAB_TEST];
  const minDayGrouped = lbDomain.groupBy(groupByCols)
    .min(LAB_DAY)
    .aggregate();

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

  const baselineMap: {[key: string]: number} = {};
  const baselineDayMap: {[key: string]: number} = {};
  const baselineSubjCol = baselineData.col(SUBJECT_ID);
  const baselineCatCol = baselineData.col(LAB_TEST_CAT);
  const baselineTestCol = baselineData.col(LAB_TEST);
  const baselineResCol = baselineData.col(LAB_RES_N);
  const minDayColName = `min(${LAB_DAY})`;
  const baselineDayCol = baselineData.col(minDayColName);

  if (!baselineSubjCol || !baselineTestCol || !baselineResCol || !baselineDayCol)
    return;

  for (let i = 0; i < baselineData.rowCount; i++) {
    if (baselineResCol!.isNone(i) || baselineDayCol!.isNone(i))
      continue;

    const subjectId = baselineSubjCol.get(i);
    const cat = baselineCatCol ? baselineCatCol.get(i) : '';
    const test = baselineTestCol.get(i);
    const key = labCatCol ? `${subjectId}|${cat}|${test}` : `${subjectId}|${test}`;
    const baselineValue = baselineResCol.get(i);
    const baselineDay = baselineDayCol.get(i);

    if (!baselineMap.hasOwnProperty(key)) {
      baselineMap[key] = baselineValue;
      baselineDayMap[key] = baselineDay;
    }
  }

  const buildKey = (subjectId: string, cat: string | null, test: string): string => {
    return labCatCol && cat !== null ? `${subjectId}|${cat}|${test}` : `${subjectId}|${test}`;
  };

  const getRowKey = (i: number): {subjectId: string, cat: string | null, test: string, key: string} | null => {
    if (labDayCol.isNone(i) || labResNCol.isNone(i))
      return null;

    const subjectId = subjectIdCol.get(i);
    const cat = labCatCol ? labCatCol.get(i) : null;
    const test = labTestCol.get(i);
    const key = buildKey(subjectId, cat, test);
    return {subjectId, cat, test, key};
  };

  interface PostBaselineData {
    maxValue: number | null;
    minPctChg: number | null;
    maxPctChg: number | null;
  }
  const postBaselineData: {[key: string]: PostBaselineData} = {};

  const baselineCol = lbDomain.columns.addNewFloat(baselineColName);
  baselineCol.init((i) => {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      return null;
    return baselineMap[rowInfo.key] ?? null;
  });

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

  for (let i = 0; i < lbDomain.rowCount; i++) {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      continue;

    const baselineDay = baselineDayMap[rowInfo.key];
    if (baselineDay === undefined)
      continue;

    const currentDay = labDayCol.get(i);
    if (currentDay <= baselineDay)
      continue;

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

  const maxPostValueCol = lbDomain.columns.addNewFloat(maxPostValueColName);
  maxPostValueCol.init((i) => {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      return null;
    return postBaselineData[rowInfo.key]?.maxValue ?? null;
  });

  const minPctChgCol = lbDomain.columns.addNewFloat(minPctChgColName);
  minPctChgCol.init((i) => {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      return null;
    return postBaselineData[rowInfo.key]?.minPctChg ?? null;
  });

  const maxPctChgCol = lbDomain.columns.addNewFloat(maxPctChgColName);
  maxPctChgCol.init((i) => {
    const rowInfo = getRowKey(i);
    if (!rowInfo)
      return null;
    return postBaselineData[rowInfo.key]?.maxPctChg ?? null;
  });
}

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
  const armCol = lbDomain.col(PLANNED_TRT_ARM) || lbDomain.col(ACT_TRT_ARM);
  const visitCol = lbDomain.col(LAB_DAY);

  if (!labTestCol || !labResNCol || !armCol || !visitCol)
    return;

  const controlMeanColName = 'CONTROL_MEAN';
  const deltaVsControlColName = 'DELTA_VS_CONTROL';
  const pctVsControlColName = 'PCT_VS_CONTROL';

  if (lbDomain.col(controlMeanColName))
    lbDomain.columns.remove(controlMeanColName);
  if (lbDomain.col(deltaVsControlColName))
    lbDomain.columns.remove(deltaVsControlColName);
  if (lbDomain.col(pctVsControlColName))
    lbDomain.columns.remove(pctVsControlColName);

  const treatmentToControlMap: {[treatment: string]: string} = {};
  const controlArms = new Set<string>();
  for (const config of treatmentControlConfig) {
    treatmentToControlMap[config.treatment] = config.control;
    controlArms.add(config.control);
  }

  const testGroupCols = labCatCol ? [LAB_TEST_CAT, LAB_TEST] : [LAB_TEST];
  const visitColName = visitCol.name;

  const controlMeanMap: {[key: string]: number} = {};

  for (const controlArm of controlArms) {
    const controlFiltered = lbDomain.groupBy(testGroupCols.concat([visitColName, LAB_RES_N]))
      .where(`${armCol.name} = ${controlArm}`)
      .aggregate();

    if (controlFiltered.rowCount === 0)
      continue;

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

      const key = `${controlArm}|${cat}|${test}|${visit}`;
      controlMeanMap[key] = meanValue;
    }
  }

  const controlMeanCol = lbDomain.columns.addNewFloat(controlMeanColName);
  controlMeanCol.init((i) => {
    if (labResNCol.isNone(i))
      return null;

    const arm = armCol.get(i);
    let controlArm: string | null = null;
    if (controlArms.has(arm))
      controlArm = arm;
    else if (treatmentToControlMap[arm])
      controlArm = treatmentToControlMap[arm];

    if (!controlArm)
      return null;

    const test = labTestCol.get(i);
    const cat = labCatCol ? labCatCol.get(i) : '';
    const visit = visitCol.get(i);
    const key = `${controlArm}|${cat}|${test}|${visit}`;

    return controlMeanMap[key] ?? null;
  });

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

export function createAllMeasurementsDf(studyId: string): DG.DataFrame | null {
  //body weight gain measurements should be analyaed separately, since BGENDY should be considered intead of BGDY and
  //each body gain weight record corresponds to different time periods, for instance:
  //USUBJID,BGSTRESN,BGDY,BGENDY
//PC201708-4115,24,1,29
//PC201708-4115,57,1,92
//PC201708-4115,108,1,106
//PC201708-4115,14,29,57
//PC201708-4115,27,57,85
//PC201708-4115,-7,85,92
//PC201708-4115,51,92,106
  const excludeFromMeasurements = 'bg';
  let resDf: DG.DataFrame | null = null;
  for (const it of studies[studyId].domains.all()) {
    if (it.name === excludeFromMeasurements)
      continue;
    let requiredColumns = [SUBJECT_ID, DOMAIN, `${it.name.toUpperCase()}TEST`,
      `${it.name.toUpperCase()}STRESN`, `${it.name.toUpperCase()}STRESC`,
      `${it.name.toUpperCase()}STRESU`, `${it.name.toUpperCase()}DY`];
    if (requiredColumns.every((colName) => it.columns.names().includes(colName))) {
      createVisitDayStrCol(it);
      requiredColumns = requiredColumns.concat([VISIT_DAY_STR]);

      const catColName = `${it.name.toUpperCase()}CAT`;
      const specColName = `${it.name.toUpperCase()}SPEC`;
      const catCol = it.col(catColName);
      const specCol = it.col(specColName);
      const containsCatOrSpecCol = catCol || specCol;
      const fullTestNameCol = 'full_test_name';
      if (containsCatOrSpecCol) {
        it.columns.addNewString(fullTestNameCol)
          .init((i) => {
            const arr = [];
            if (catCol && !catCol.isNone(i))
              arr.push(catCol.get(i));
            if (specCol && !specCol.isNone(i))
              arr.push(specCol.get(i));
            return `${it.get(`${it.name.toUpperCase()}TEST`, i)} ${arr.length ? `(${arr.join(',')}) ` : ``}`;
          });
        requiredColumns.push(fullTestNameCol);
        requiredColumns = requiredColumns.filter((i) => i !== `${it.name.toUpperCase()}TEST`);
      }
      if (!catCol)
        it.columns.addNewString(catColName);
      if (!specCol)
        it.columns.addNewString(specColName);
      requiredColumns.push(`${it.name.toUpperCase()}CAT`);
      requiredColumns.push(`${it.name.toUpperCase()}SPEC`);

      // Per SENDIG v3.1.1: carry --BLFL (Baseline Flag) when available.
      // BLFL is Expected in LB, EG, CV, RE and Permissible in BW, PC, VS.
      const blflColName = `${it.name.toUpperCase()}BLFL`;
      const hasBlfl = it.col(blflColName) !== null;
      if (hasBlfl)
        requiredColumns.push(blflColName);

      const df = it.clone(null, requiredColumns);
      df.getCol(containsCatOrSpecCol ? fullTestNameCol : `${it.name.toUpperCase()}TEST`).name = 'test';
      df.getCol(`${it.name.toUpperCase()}CAT`).name = 'category';
      df.getCol(`${it.name.toUpperCase()}SPEC`).name = 'specimen';
      df.getCol(`${it.name.toUpperCase()}STRESN`).name = 'result';
      df.getCol(`${it.name.toUpperCase()}DY`).name = 'visit_day';
      df.getCol(`${it.name.toUpperCase()}STRESC`).name = 'string_result';
      df.getCol(`${it.name.toUpperCase()}STRESU`).name = 'units';
      if (hasBlfl)
        df.getCol(blflColName).name = 'blfl';
      else if (!df.col('blfl'))
        df.columns.addNewString('blfl');

      if (containsCatOrSpecCol)
        it.columns.remove(fullTestNameCol);

      if (!resDf)
        resDf = df;
      else
        resDf.append(df, true);
    }
  }

  if (resDf)
    calculateBaselineColumns(resDf);

  return resDf;
}

/** Per SENDIG v3.1.1, baseline is identified by --BLFL = "Y" (sponsor-defined).
 *  When BLFL is absent or empty, falls back to the earliest study day per subject/test.
 *  Baseline can be a derived average (--DRVFL = "Y" + --BLFL = "Y").
 *  Adds baseline, change, and pct_change columns to the unified measurements table. */
function calculateBaselineColumns(df: DG.DataFrame): void {
  const subjectCol = df.col(SUBJECT_ID);
  const testCol = df.col('test');
  const resultCol = df.col('result');
  const dayCol = df.col('visit_day');
  const blflCol = df.col('blfl');
  if (!subjectCol || !testCol || !resultCol || !dayCol)
    return;

  // Build baseline map: subject|test -> baseline value
  // Priority 1: rows where BLFL = "Y" (per SENDIG, the sponsor-defined baseline flag)
  // Priority 2: row with earliest study day (fallback when BLFL is not populated)
  const baselineMap: {[key: string]: number} = {};
  const fallbackMap: {[key: string]: {value: number, day: number}} = {};

  for (let i = 0; i < df.rowCount; i++) {
    if (resultCol.isNone(i) || subjectCol.isNone(i) || testCol.isNone(i))
      continue;

    const key = `${subjectCol.get(i)}|${testCol.get(i)}`;
    const value = resultCol.get(i);

    if (blflCol && !blflCol.isNone(i) && blflCol.get(i) === 'Y') {
      baselineMap[key] = value;
      continue;
    }

    if (!dayCol.isNone(i)) {
      const day = dayCol.get(i);
      if (!fallbackMap[key] || day < fallbackMap[key].day)
        fallbackMap[key] = {value, day};
    }
  }

  // Apply fallback for keys without explicit BLFL
  for (const key of Object.keys(fallbackMap)) {
    if (!baselineMap.hasOwnProperty(key))
      baselineMap[key] = fallbackMap[key].value;
  }

  const getKey = (i: number): string | null => {
    if (subjectCol.isNone(i) || testCol.isNone(i))
      return null;
    return `${subjectCol.get(i)}|${testCol.get(i)}`;
  };

  df.columns.addNewFloat('baseline').init((i) => {
    const key = getKey(i);
    if (!key || !baselineMap.hasOwnProperty(key))
      return null;
    return baselineMap[key];
  });

  const baselineCol = df.col('baseline')!;
  df.columns.addNewFloat('change').init((i) => {
    if (resultCol.isNone(i) || baselineCol.isNone(i))
      return null;
    return resultCol.get(i) - baselineCol.get(i);
  });

  df.columns.addNewFloat('pct_change').init((i) => {
    if (resultCol.isNone(i) || baselineCol.isNone(i))
      return null;
    const bl = baselineCol.get(i);
    if (bl === 0)
      return null;
    return 100 * (resultCol.get(i) - bl) / bl;
  });
}

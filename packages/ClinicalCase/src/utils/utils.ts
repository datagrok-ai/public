import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {scripts} from '../package-api';
import {ClinStudyConfig} from './types';
import {ACT_TRT_ARM, AGE, AGETXT, COL_HAS_ERRORS_POSTFIX, DEATH_DATE, ERRORS_POSTFIX,
  ETHNIC, HAS_VALIDATION_ERRORS_COL, PLANNED_TRT_ARM, RACE, SEX, SPECIES, SUBJ_REF_ENDT,
  SUBJ_REF_STDT, SUBJECT_ID, VIOLATED_RULES_COL} from '../constants/columns-constants';
import {studies} from './app-utils';

export function updateDivInnerHTML(div: HTMLElement, content: any) {
  div.innerHTML = '';
  div.append(content);
}

export function createFilters(df: DG.DataFrame) {
  return DG.Viewer.fromType('Filters', df, {
    'showContextMenu': false,
  });
}

export function removeExtension(filename: string): string {
  const lastDotIndex = filename.lastIndexOf('.');
  return lastDotIndex === -1 ? filename : filename.substring(0, lastDotIndex);
};

export async function readClinicalFile(file: DG.FileInfo): Promise<DG.DataFrame> {
  let df: DG.DataFrame | null = null;
  try {
    if (file.extension === 'xpt')
      df = await scripts.readSas(file);
    else
      df = await DG.DataFrame.fromCsv(await file.readAsString());
  } catch (e: any) {
    grok.shell.error(`Error loading ${file.name}: ${e?.message ?? e}`);
  }
  return df;
}

export function studyConfigToMap(studyConfig: ClinStudyConfig): {[key: string]: any} {
  const map: {[key: string]: any} = {};
  for (const key of Object.keys(studyConfig)) {
    if (key !== 'other' && key !== 'fieldsDefinitions') {
      let formattedKey = '';
      if (key === 'totalSubjects')
        formattedKey = 'Total subjects';
      else if (key === 'startDate')
        formattedKey = 'Start date';
      else if (key === 'endDate')
        formattedKey = 'End date';
      else
        formattedKey = key[0].toUpperCase() + key.slice(1);
      map[formattedKey] = studyConfig[key];
    }
  }
  const otherKeys = Object.keys(studyConfig.other);
  const keysToSkip: string[] = [];
  for (const key of Object.keys(studyConfig.other)) {
    if (keysToSkip.includes(key))
      continue;
    const keySplitted = key.split(' ');
    const unitsKey = otherKeys
      .filter((it) => it.toLowerCase() === `${keySplitted[0]} units`.toLocaleLowerCase() ||
          it.toLowerCase() === `${keySplitted[0]} unit`.toLowerCase());
    if (unitsKey.length) {
      map[keySplitted[0]] = `${studyConfig.other[key]} ${studyConfig.other[unitsKey[0]]}`;
      keysToSkip.push(unitsKey[0]);
    } else
      map[key] = studyConfig.other[key];
  }
  return map;
}

export function hideValidationColumns(tv: DG.TableView) {
  const colNames = tv.dataFrame.columns.names();
  for (const colName of colNames) {
    if (colName === VIOLATED_RULES_COL || colName.endsWith(COL_HAS_ERRORS_POSTFIX) || colName.endsWith(ERRORS_POSTFIX))
      tv.grid.col(colName).visible = false;
  }
}

export function addDomainFilters(tv: DG.TableView, studyId: string): DG.FilterGroup | null {
  const dmDomainColNames = studies[studyId].domains.dm.columns.names();
  const firstColumnsNames = [SUBJ_REF_ENDT, SUBJ_REF_STDT, DEATH_DATE, SPECIES, ETHNIC, RACE, SEX,
    AGETXT, AGE, PLANNED_TRT_ARM, ACT_TRT_ARM, SUBJECT_ID];
  const otherCols = dmDomainColNames
    .filter((it) => !firstColumnsNames.includes(it) && it !== HAS_VALIDATION_ERRORS_COL);
  const filterCols: {colName: string, filterType: DG.FILTER_TYPE}[] = [];
  const addFilterColNames = (columnNamesList: string[]) => {
    for (const colName of columnNamesList) {
      const col = studies[studyId].domains.dm.col(colName);
      if (studies[studyId].domains.dm.col(colName)) {
        //keep only filters containing more than one category
        if (col.type === DG.TYPE.STRING && col.categories.length === 1)
          continue;
        const filterType = [DG.TYPE.INT, DG.TYPE.FLOAT].includes(col.type as DG.TYPE) ? DG.FILTER_TYPE.HISTOGRAM :
          col.type === DG.TYPE.BOOL ? DG.FILTER_TYPE.BOOL_COLUMNS : DG.FILTER_TYPE.CATEGORICAL;
        filterCols.push({colName, filterType});
      }
    }
  };
  addFilterColNames(firstColumnsNames);
  addFilterColNames(otherCols);
  const fg = tv.getFiltersGroup({createDefaultFilters: false});
  if (filterCols.length) {
    for (const filter of filterCols) {
      fg.updateOrAdd({
        type: filter.filterType,
        column: filter.colName,
        columnName: filter.colName,
      });
    }
    return fg;
  }
  return null;
}

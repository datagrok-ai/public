import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {scripts} from '../package-api';
import {ClinStudyConfig} from './types';

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
    if (key !== 'other') {
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

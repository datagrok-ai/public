/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';

import {
  APP_PATH, PATTERN_APP_DATA_FILENAME, CODES_TO_HELM_DICT_FILENAME,
  CODES_TO_SYMBOLS_FILENAME, MONOMERS_WITH_PHOSPHATE_FILENAME
} from './const';
import {PatternAppData, CodeToSymbol, FormatToHELMDict} from './types';

const fileSource = new DG.FileSource(APP_PATH);

export let PATTERN_APP_DATA: PatternAppData;
export let CODES_TO_HELM_DICT: FormatToHELMDict;
export let CODES_TO_SYMBOLS_DICT: CodeToSymbol;
export let MONOMERS_WITH_PHOSPHATE: {[key: string]: string[]};

export async function loadJsonData(): Promise<void> {
  if (isAllJsonDataLoaded())
    return;

  const jsonFileNames = [
    PATTERN_APP_DATA_FILENAME, CODES_TO_HELM_DICT_FILENAME, CODES_TO_SYMBOLS_FILENAME, MONOMERS_WITH_PHOSPHATE_FILENAME
  ];

  [
    PATTERN_APP_DATA,
    CODES_TO_HELM_DICT,
    CODES_TO_SYMBOLS_DICT,
    MONOMERS_WITH_PHOSPHATE
  ] = await Promise.all(
    jsonFileNames.map((fileName) => loadAndParseJson(fileName))
  );
}

async function loadAndParseJson(filePath: string): Promise<any> {
  try {
    const content = await fileSource.readAsText(filePath);
    return JSON.parse(content);
  } catch (err) {
    console.error(`Error loading json from ${filePath}:`, err);
  }
}

function isAllJsonDataLoaded(): boolean {
  const data = [PATTERN_APP_DATA, CODES_TO_HELM_DICT, CODES_TO_SYMBOLS_DICT, MONOMERS_WITH_PHOSPHATE];

  return data.every((item) => item !== undefined);
}

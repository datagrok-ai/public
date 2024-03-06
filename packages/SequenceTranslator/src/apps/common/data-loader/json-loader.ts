/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {APP_PATH, AXOLABS_STYLE_FILENAME, CODES_TO_HELM_DICT_FILENAME, CODES_TO_SYMBOLS_FILENAME, MONOMERS_WITH_PHOSPHATE_FILENAME} from './const';
import {AxolabsStyle, FormatToHELMDict, CodeToSymbol} from './types';

const fileSource = new DG.FileSource(APP_PATH);

export let AXOLABS_STYLE_MAP: AxolabsStyle;
export let CODES_TO_HELM_DICT: FormatToHELMDict;
export let CODES_TO_SYMBOLS_DICT: CodeToSymbol;
export let MONOMERS_WITH_PHOSPHATE: {[key: string]: string[]};

export async function loadJsonData(): Promise<void> {
  if (isAllJsonDataLoaded())
    return;

  const jsonFileNames = [AXOLABS_STYLE_FILENAME, CODES_TO_HELM_DICT_FILENAME, CODES_TO_SYMBOLS_FILENAME, MONOMERS_WITH_PHOSPHATE_FILENAME];

  [
    AXOLABS_STYLE_MAP,
    CODES_TO_HELM_DICT,
    CODES_TO_SYMBOLS_DICT,
    MONOMERS_WITH_PHOSPHATE
  ] = await Promise.all(
    jsonFileNames.map(fileName => loadAndParseJson(fileName))
  );
}

async function loadAndParseJson(filePath: string): Promise<any> {
  try {
    const content = await fileSource.readAsText(filePath);
    return JSON.parse(content);
  } catch (err) {
    console.error(`Error loading json from ${filePath}:`, err);
    // throw err;
  }
}

function isAllJsonDataLoaded(): boolean {
  const data = [AXOLABS_STYLE_MAP, CODES_TO_HELM_DICT, CODES_TO_SYMBOLS_DICT, MONOMERS_WITH_PHOSPHATE];

  return data.every((item) => item !== undefined);
}

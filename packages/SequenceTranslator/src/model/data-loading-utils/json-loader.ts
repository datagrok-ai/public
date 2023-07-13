/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {APP_PATH, AXOLABS_STYLE_FILENAME, CODES_TO_HELM_DICT_FILENAME, CODES_TO_SYMBOLS_FILENAME, MONOMERS_WITH_PHOSPHATE_LINKERS} from './const';
import {AxolabsStyle, FormatToHELMDict, CodeToSymbol} from './types';

const fileSource = new DG.FileSource(APP_PATH);

export let axolabsStyleMap: AxolabsStyle;
export let codesToHelmDictionary: FormatToHELMDict;
export let codesToSymbolsDictionary: CodeToSymbol;
export let monomersWithPhosphateLinkers: {[key: string]: string[]};

export async function getJsonData(): Promise<void> {
  const data = [axolabsStyleMap, codesToHelmDictionary, codesToSymbolsDictionary, monomersWithPhosphateLinkers];

  if (data.every((item) => item !== undefined))
    return;

  axolabsStyleMap = await parse(AXOLABS_STYLE_FILENAME);
  codesToHelmDictionary = await parse(CODES_TO_HELM_DICT_FILENAME);
  codesToSymbolsDictionary = await parse(CODES_TO_SYMBOLS_FILENAME);
  monomersWithPhosphateLinkers = await parse(MONOMERS_WITH_PHOSPHATE_LINKERS);
}

async function parse(path: string): Promise<any> {
  let parsedJson: string;
  try {
    parsedJson = JSON.parse(await fileSource.readAsText(path))
  } catch (err: any) {
    const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
    throw new Error(`Error loading json from ${path}:` + errMsg);
  }
  return parsedJson;
}


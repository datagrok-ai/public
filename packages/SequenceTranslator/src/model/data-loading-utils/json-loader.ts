/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {APP_PATH, FORMAT_DICT_FILENAME, AXOLABS_STYLE_FILENAME} from './const';
import {FormatDict, AxolabsStyle} from './types';

const fileSource = new DG.FileSource(APP_PATH);
async function parse(path: string): Promise<any> {
  return JSON.parse(await fileSource.readAsText(path));
}

export let formatDictionary: FormatDict;
export let axolabsStyleMap: AxolabsStyle;

export async function getJsonData(): Promise<void> {
  try {
    formatDictionary = await parse(FORMAT_DICT_FILENAME);
    axolabsStyleMap = await parse(AXOLABS_STYLE_FILENAME);
  } catch (err: any) {
    const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
    throw new Error('ST: Loading format dictionary error: ' + errMsg);
  }
}

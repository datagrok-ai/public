/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type KeyToValue = {[key: string]: string};

export type Edges = {
  [key: string]: KeyToValue
}

export type Dict = {
  [sourceFormat: string]: {
    [targetFormat: string]: KeyToValue | Edges
  }
}

const DICT_PATH = 'System:AppData/SequenceTranslator';
const DEFAULT_DICT_FILENAME = 'format-converter-dict.json';

export class FormatDictionaryLoader {
  private constructor() { }

  private static instance?: FormatDictionaryLoader;
  private formatDictionary: Dict;

  static getInstance() {
    if (FormatDictionaryLoader.instance === undefined)
      FormatDictionaryLoader.instance = new FormatDictionaryLoader();
    return FormatDictionaryLoader.instance;
  }

  async init() {
    try {
      const fileSource = new DG.FileSource(DICT_PATH);
      const file = await fileSource.readAsText(DEFAULT_DICT_FILENAME);
      this.formatDictionary = JSON.parse(file);
    } catch (err: any) {
      const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
      throw new Error('ST: Loading format dictionary error: ' + errMsg);
    }
  }

  public getDictionary(): Dict {
    if (this.formatDictionary === undefined)
      throw new Error('ST: format dictionary not initialized');
    return this.formatDictionary as Dict;
  }
}

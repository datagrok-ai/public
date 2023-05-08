/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type KeyToValue = {[key: string]: string};

export type Edges = {
  [key: string]: KeyToValue
}

export type FormatDict = {
  [sourceFormat: string]: {
    [targetFormat: string]: KeyToValue | Edges
  }
}

export type AxolabsStyle = {
  [index: string]: {
    fullName: string,
    symbols: string[],
    color: string,
  }
};

const APP_PATH = 'System:AppData/SequenceTranslator';
const FORMAT_DICT_FILENAME = 'format-converter-dict.json';
const AXOLABS_STYLE_FILENAME = 'axolabs-style.json';

export class JsonLoader {
  private constructor() { }

  private static instance?: JsonLoader;
  private formatDictionary: FormatDict;
  private axolabsStyle: AxolabsStyle;

  static getInstance() {
    if (JsonLoader.instance === undefined)
      JsonLoader.instance = new JsonLoader();
    return JsonLoader.instance;
  }

  async init() {
    try {
      const fileSource = new DG.FileSource(APP_PATH);

      const formatDictFile = await fileSource.readAsText(FORMAT_DICT_FILENAME);
      this.formatDictionary = JSON.parse(formatDictFile);

      const axolabsStyleFile = await fileSource.readAsText(AXOLABS_STYLE_FILENAME);
      this.axolabsStyle = JSON.parse(axolabsStyleFile);
    } catch (err: any) {
      const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
      throw new Error('ST: Loading format dictionary error: ' + errMsg);
    }
  }

  public getFormatDictionary(): FormatDict {
    if (this.formatDictionary === undefined)
      throw new Error('ST: format dictionary not initialized');
    return this.formatDictionary as FormatDict;
  }

  public getAxolabsStyleDictionary(): AxolabsStyle {
    if (this.axolabsStyle === undefined)
      throw new Error('ST: axolabs style dictionary not initialized');
    return this.axolabsStyle as AxolabsStyle;
  }
}

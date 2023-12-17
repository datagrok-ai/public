/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {axolabsStyleMap} from '../../../model/data-loading-utils/json-loader';

export class DataManager {
  readonly baseChoices: string[] = Object.keys(axolabsStyleMap);
  readonly defaultBase: string = this.baseChoices[0];
}

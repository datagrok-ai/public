/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getSdfTab} from '../sdf-tab/sdf-tab';

export class SdfTabUI {
  constructor() {
    this._htmlDivElement = getSdfTab();
  }

  private _htmlDivElement;

  get htmlDivElement() { return this._htmlDivElement; }
}

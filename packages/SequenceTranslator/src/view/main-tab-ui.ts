/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DEFAULT_INPUT_MAIN_TAB} from './view-const';

// todo: port UI-related functionality here
import {getMainTab} from '../main-tab/main-tab';

export class MainTabUI {
  constructor(onInputChanged: (input: string) => void) {
    this._inputSequence = DEFAULT_INPUT_MAIN_TAB;
    this._htmlDivElement = ui.wait(() => getMainTab(onInputChanged));
  }

  private _htmlDivElement: HTMLDivElement;

  /** Sequence inserted on the Main tab to be translated  */
  private _inputSequence: string;

  get htmlDivElement() { return this._htmlDivElement; }

  get inputSequence() { return this._inputSequence; }

  // todo: add validation of the inserted sequence here!
  set inputSequence(newSequence: string) {
    this._inputSequence = newSequence;
  }
}

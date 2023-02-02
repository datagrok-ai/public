/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MAIN_TAB_DEFAULT_INPUT} from './view-const';

// todo: port UI-related functionality here
import {getMainTab} from '../main-tab/main-tab';

export class MainTabUI {
  constructor(onInputChanged: (input: string) => void) {
    this._inputSequence = MAIN_TAB_DEFAULT_INPUT;
    this._onInputChanged = onInputChanged;
  }

  /** Sequence inserted on the Main tab to be translated  */
  private _inputSequence: string;

  private _onInputChanged: (input: string) => void;

  get inputSequence() { return this._inputSequence; }

  /** Get the HTMLElement of the tab, async because of the internal
   * grok.functions.call() used to draw the sequence */
  async getHtmlElement(): Promise<HTMLDivElement> {
    return await getMainTab(this._onInputChanged);
  }

  // todo: add validation of the inserted sequence here!
  set inputSequence(newSequence: string) {
    this._inputSequence = newSequence;
  }
}

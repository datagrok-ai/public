/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class PatternManager {
  getLoadPatternInput(): DG.InputBase {
    // const loadPattern = ui.choiceInput('Load pattern', '', lstMy, (v: string) => parsePatternAndUpdateUi(v));
    const loadPattern = ui.choiceInput('Load pattern', '', [], () => {});
    return loadPattern;
  };
}

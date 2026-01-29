/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {expectTable as _expectTable} from '@datagrok-libraries/test/src/test';

export const _package = new DG.Package();
export * from './package.g';

export class PackageFunctions {
  @grok.decorators.func()
  static info() {
    grok.shell.info(_package.webRoot);
  }


  @grok.decorators.func()
  static expectTable(
    actual: DG.DataFrame,
    expected: DG.DataFrame): boolean {
    _expectTable(actual, expected);
    return true;
  }
}

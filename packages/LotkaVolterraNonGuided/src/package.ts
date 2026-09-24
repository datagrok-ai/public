/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';
import {runLotkaVolterra} from './app';

export const _package = new DG.Package();

//name: Lotka-Volterra (non-Guided)
//tags: app
export function LotkaVolterra(): void {
  runLotkaVolterra();
}

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

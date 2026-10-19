/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

export const _package = new DG.Package();

import {lotkaVolterraApp} from './lotka-volterra/app';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Lotka-Volterra Simulation (Guided)
//tags: app
export function lotkaVolterraSimulation(): void {
  lotkaVolterraApp(_package);
}

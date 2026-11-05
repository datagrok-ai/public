/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

import {levinsMetapopulationApp, setPackage} from './levins/app';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Levins Metapopulation Model
//tags: app
//description: Interactive ODE solver for the Levins model: simulation of occupied patch fraction p(t) dynamics
export function levinsMetapopulationModelApp(): void {
  setPackage(_package);
  levinsMetapopulationApp();
}

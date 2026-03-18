/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

export const _package = new DG.Package();

import {sirSeirEpidemicApp} from './sir-seir/app';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: SIR / SEIR Epidemic Simulation
//tags: app
//description: Interactive ODE simulation of SIR and SEIR epidemiological models
export function sirSeirEpidemicSimulation(): void {
  sirSeirEpidemicApp(_package);
}

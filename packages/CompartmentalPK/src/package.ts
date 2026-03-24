/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

export const _package = new DG.Package();

import {compartmentalPkApp} from './compartmental-pk/app';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Multi-Compartment PK Model
//description: Interactive two- and three-compartment pharmacokinetic model with IV bolus, IV infusion, and oral absorption input modes. Includes brute-force grid search parameter optimization.
//meta.role: app
//meta.browsePath: Compute | Simulations
export function compartmentalPkModelApp(): void {
  compartmentalPkApp(_package);
}

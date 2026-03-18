/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

export const _package = new DG.Package();

import {threeBodyProblemApp} from './cr3bp/app';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Circular Restricted Three-Body Problem
//tags: app
//description: Interactive simulation of the CR3BP in the co-rotating frame with orbits, Lagrange points, Jacobi constant, and zero-velocity curves
export function threeBodyProblemAppEntry(): void {
  threeBodyProblemApp(_package);
}

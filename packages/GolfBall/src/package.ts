/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

export const _package = new DG.Package();

import {golfBallApp} from './golf-ball/app';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Golf Ball Flight Simulator
//tags: app
//description: Interactive 2D golf ball flight simulation with drag, trajectory optimization, and sensitivity analysis
export function golfBallFlightSimulator(): void {
  golfBallApp(_package);
}

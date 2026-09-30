// Fitting methods registry.
//
// Worker-safe on purpose: the Diff Studio fitting worker picks optimizers from here,
// so nothing reachable from this module may touch grok/ui/DG at runtime.

import {singleObjective} from '@datagrok-libraries/sci-comp';
import {METHOD} from './constants';
import type {IOptimizer, Setting} from './optimizer-misc';
import {nelderMeadSettingsOpts, optimizeNM} from './optimizer-nelder-mead';
import {adaptOptimizer, adamSettingsOpts, buildAdamSettings, buildLbfgsbSettings, buildLbfgsSettings,
  buildPsoSettings, lbfgsbSettingsOpts, lbfgsSettingsOpts, psoSettingsOpts} from './optimizer-sci-comp';

export type OptimizerDescriptor = {
  optimizer: IOptimizer,
  settingsOpts: Map<string, Setting>,
  /** The standard path's worker pool runs only the sync codegen twin of Nelder-Mead */
  supportsWorker: boolean,
  /** The optimizer takes the box bounds of the varied inputs */
  wantsBounds: boolean,
};

export const OPTIMIZERS = new Map<METHOD, OptimizerDescriptor>([
  [METHOD.NELDER_MEAD, {
    optimizer: optimizeNM,
    settingsOpts: nelderMeadSettingsOpts,
    supportsWorker: true,
    wantsBounds: false,
  }],
  [METHOD.LBFGSB, {
    optimizer: adaptOptimizer(() => new singleObjective.LBFGSB(), buildLbfgsbSettings),
    settingsOpts: lbfgsbSettingsOpts,
    supportsWorker: false,
    wantsBounds: true,
  }],
  [METHOD.PSO, {
    optimizer: adaptOptimizer(() => new singleObjective.PSO(), buildPsoSettings),
    settingsOpts: psoSettingsOpts,
    supportsWorker: false,
    wantsBounds: true,
  }],
  [METHOD.LBFGS, {
    optimizer: adaptOptimizer(() => new singleObjective.LBFGS(), buildLbfgsSettings),
    settingsOpts: lbfgsSettingsOpts,
    supportsWorker: false,
    wantsBounds: true,
  }],
  [METHOD.ADAM, {
    optimizer: adaptOptimizer(() => new singleObjective.Adam(), buildAdamSettings),
    settingsOpts: adamSettingsOpts,
    supportsWorker: false,
    wantsBounds: true,
  }],
]);

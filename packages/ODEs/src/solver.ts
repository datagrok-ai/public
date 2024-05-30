// Solver of initial value problem

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ODEs, SolverOptions} from './solver-tools/solver-defs';
import {mrt} from './solver-tools/mrt-method';

/** Default solver of initial value problem. */
export const solveDefault = (odes: ODEs): DG.DataFrame => mrt(odes);

/** Customizable solver of initial value problem. */
export const solveIVP = (odes: ODEs, options?: Partial<SolverOptions>): DG.DataFrame => mrt(odes, options);

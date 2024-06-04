// Manager of callbacks

import {SolverOptions} from '../solver-defs';
import {Callback} from './callback-base';

/** Get callback corresponding to the options */
export function getCallback(options?: Partial<SolverOptions>): Callback | undefined {
  if (options === undefined)
    return undefined;

  return new Callback();
}

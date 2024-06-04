// Manager of callbacks

import {SolverOptions} from '../solver-defs';
import {Callback} from './callback-base';
import {TerminatorCallback} from './terminator-callback';

/** Get callback corresponding to the options */
export function getCallback(options?: Partial<SolverOptions>): Callback | undefined {
  if (options === undefined)
    return undefined;

  if (options.maxTimeMs !== undefined)
    return new TerminatorCallback(options.maxTimeMs);

  return undefined;
}

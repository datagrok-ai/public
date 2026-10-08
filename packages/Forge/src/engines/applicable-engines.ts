import * as DG from 'datagrok-api/dg';
import {releaseFrame, sharedFrame} from '../preparation/shared-frame';
import {Engine, isComplete} from './engine';
import {isApplicable} from './engine-calls';

export interface EngineFailure { engine: Engine; error: unknown }
export interface ApplicableEngines { applicable: Engine[]; failed: EngineFailure[] }

/** The complete [engines] whose `isApplicable` accepts the data, in the given order, asked in parallel; an engine
 * whose check throws is left out and returned in `failed`, for the caller to log. */
export async function applicableEngines(engines: Engine[], features: DG.Column[], target: DG.Column):
  Promise<ApplicableEngines> {
  const complete = engines.filter(isComplete);
  const frame = sharedFrame(features);
  try {
    const outcomes = await Promise.allSettled(complete.map((engine) => isApplicable(engine, frame, target)));
    const result: ApplicableEngines = {applicable: [], failed: []};
    for (const [i, outcome] of outcomes.entries()) {
      if (outcome.status === 'rejected')
        result.failed.push({engine: complete[i], error: outcome.reason});
      else if (outcome.value)
        result.applicable.push(complete[i]);
    }
    return result;
  } finally {
    releaseFrame(frame);
  }
}

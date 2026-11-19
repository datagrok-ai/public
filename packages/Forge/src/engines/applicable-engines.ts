import * as DG from 'datagrok-api/dg';
import {Engine, isComplete} from './engine';
import {CHECK_TIMEOUT_MS, isApplicable} from './engine-calls';

export interface EngineFailure { engine: Engine; error: unknown }
export interface ApplicableEngines { applicable: Engine[]; failed: EngineFailure[] }

/** The complete [engines] whose `isApplicable` accepts the data, in the given order, asked in parallel; an engine
 * whose check throws is left out and returned in `failed`, for the caller to log; a check that does not answer within
 * [timeoutMs] fails with `NoAnswerError`. */
export async function applicableEngines(engines: Engine[], features: DG.Column[], target: DG.Column,
  timeoutMs: number = CHECK_TIMEOUT_MS): Promise<ApplicableEngines> {
  const complete = engines.filter(isComplete);
  const outcomes = await Promise.allSettled(complete.map((engine) =>
    isApplicable(engine, features, target, timeoutMs)));
  const result: ApplicableEngines = {applicable: [], failed: []};
  for (const [i, outcome] of outcomes.entries()) {
    if (outcome.status === 'rejected')
      result.failed.push({engine: complete[i], error: outcome.reason});
    else if (outcome.value)
      result.applicable.push(complete[i]);
  }
  return result;
}

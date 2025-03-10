import * as DSL from '@datagrok/diff-grok';
import {DiffGrok} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';

export function getDiffGrok(ivp: DSL.IVP): DiffGrok | undefined {
  return {
    ivp: ivp,
    ivpWW: DSL.getIvp2WebWorker(ivp),
    pipelineCreator: DSL.getPipelineCreator(ivp),
  };
}

import * as DSL from '@datagrok/diff-grok';
import {DiffGrok} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';

export function getDiffGrok(ivp: DSL.IVP): DiffGrok | undefined {
  if ((ivp.loop !== null) || (ivp.updates !== null))
    return undefined;

  return {
    ivp: ivp,
    ivpWW: DSL.getIvp2WebWorker(ivp),
    pipelineCreator: new BasicModelPipelineCreator(ivp),
  };
}

class BasicModelPipelineCreator extends DSL.PipelineCreator {
  constructor(ivp: DSL.IVP) {
    super(ivp);
  }

  public getPipeline(inputs: Float64Array): DSL.Pipeline {
    return {
      wrappers: [
        {
          preproc: null,
          out: null,
          postproc: null,
        },
      ],
      out: DSL.getOutputCode(this.ivp),
    };
  }
} // BasicModelPipelineCreator

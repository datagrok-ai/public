import {ItemPathArray} from '../data/common-types';
import {pathJoin} from '../utils';
import {FuncallStateItem, PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed} from './config-processing-utils';
import {PipelineSelfRef, PipelineStepConfiguration} from './PipelineConfiguration';
import {PipelineInstanceState} from './PipelineInstance';

type TraverseItem = PipelineConfigurationProcessed | PipelineStepConfiguration<FuncallStateItem[]> | PipelineSelfRef;

function isPipelineStaticConfig(c: TraverseItem): c is PipelineConfigurationStaticProcessed {
  return !!((c as PipelineConfigurationStaticProcessed).type === 'static');
}

function isPipelineParallelConfig(c: TraverseItem): c is PipelineConfigurationParallelProcessed {
  return !!((c as PipelineConfigurationParallelProcessed).type === 'parallel');
}

function isPipelineSequentialConfig(c: TraverseItem): c is PipelineConfigurationSequentialProcessed {
  return !!((c as PipelineConfigurationSequentialProcessed).type === 'sequential');
}

function isPipelineStepConfig(c: TraverseItem): c is PipelineStepConfiguration<FuncallStateItem[]> {
  return !isPipelineStaticConfig(c) && !isPipelineParallelConfig(c) && isPipelineSequentialConfig(c) && !isPipelineSelfRef(c);
}

function isPipelineSelfRef(c: TraverseItem): c is PipelineSelfRef {
  return !!((c as PipelineSelfRef).type === 'selfRef');
}

export function getConfigByInstancePath(instancePath: ItemPathArray, config: PipelineConfigurationProcessed): TraverseItem {
  let node: TraverseItem = config;
  if (node.id !== config.id) {
    throw new Error(`Instance config id ${node.id} is not matching ${config.id}`);
  }
  const [, ...rest] = instancePath;
  // TODO
  for (const segment of rest) {
    if (isPipelineStaticConfig(node)) {

    } else if (isPipelineParallelConfig(node) || isPipelineSequentialConfig(node)) {

    } else if (isPipelineStepConfig(node)) {

    } else if (isPipelineSelfRef(node)) {

    }
  }
}

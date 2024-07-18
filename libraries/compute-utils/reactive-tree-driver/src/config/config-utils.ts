import {ItemPathArray} from '../data/common-types';
import {pathJoin} from '../utils';
import {PipelineConfigurationProcessed, TraverseItem, isPipelineParallelConfig, isPipelineSequentialConfig, isPipelineStaticConfig} from './config-processing-utils';

export function traverseConfig<T>(
  config: PipelineConfigurationProcessed,
  handler: (
    acc: T,
    conf: TraverseItem<PipelineConfigurationProcessed>,
    path: ItemPathArray,
    stop: () => void,
  ) => T,
  acc: T,
  path: ItemPathArray = [],
): T {
  const q = [{config, path}];
  let stop = false;
  const signal = () => stop = true;
  while (q.length) {
    const {config, path} = q.shift()!;
    acc = handler(acc, config, path, signal);
    if (stop)
      return acc;
    const items = getNextItems(config, path);
    q.push(...items);
  }
  return acc;
}

function getNextItems<C extends TraverseItem<PipelineConfigurationProcessed>>(config: C, path: ItemPathArray) {
  const nextItems = getNextConfigs(config);
  const items = nextItems.map((item) => ({...item, path: pathJoin(path, [config.id])}));
  return items;
}

function getNextConfigs<C extends TraverseItem<PipelineConfigurationProcessed>>(config: C) {
  if (isPipelineParallelConfig(config) || isPipelineSequentialConfig(config))
    return config.items.map((item) => ({config: item.config as C}));
  else if (isPipelineStaticConfig(config))
    return config.steps.map((config) => ({config: config as C}));

  return [];
}

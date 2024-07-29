import {ItemPathArray} from '../data/common-types';
import {buildTraverseD} from '../data/traversable';
import {FuncallStateItem, PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed} from './config-processing-utils';
import {PipelineSelfRef, PipelineStepConfiguration} from './PipelineConfiguration';

export type ConfigTraverseItem = PipelineConfigurationProcessed | PipelineStepConfiguration<FuncallStateItem[]> | PipelineSelfRef;

export type ConfigItem = PipelineConfigurationProcessed | PipelineStepConfiguration<FuncallStateItem[]>;

export function isPipelineStaticConfig(c: ConfigTraverseItem): c is PipelineConfigurationStaticProcessed {
  return !!((c as PipelineConfigurationStaticProcessed).type === 'static');
}

export function isPipelineParallelConfig(c: ConfigTraverseItem): c is PipelineConfigurationParallelProcessed {
  return !!((c as PipelineConfigurationParallelProcessed).type === 'parallel');
}

export function isPipelineSequentialConfig(c: ConfigTraverseItem): c is PipelineConfigurationSequentialProcessed {
  return !!((c as PipelineConfigurationSequentialProcessed).type === 'sequential');
}

export function isPipelineConfig(c: ConfigTraverseItem):  c is PipelineConfigurationProcessed {
  return isPipelineStaticConfig(c) ||isPipelineParallelConfig(c) || isPipelineSequentialConfig(c)
}

export function isPipelineStepConfig(c: ConfigTraverseItem): c is PipelineStepConfiguration<FuncallStateItem[]> {
  return !isPipelineConfig(c) && !isPipelineSelfRef(c);
}

export function isPipelineSelfRef(c: ConfigTraverseItem): c is PipelineSelfRef {
  return !!((c as PipelineSelfRef).type === 'selfRef');
}

export function buildRefMap(config: PipelineConfigurationProcessed): Map<string, PipelineConfigurationProcessed> {
  const res = new Map<string, PipelineConfigurationProcessed>();

  const getNextItem = (item: ConfigTraverseItem, path: ItemPathArray) => {
    if (isPipelineStepConfig(item)) {
      return [] as any;
    } else if (isPipelineSelfRef(item)) {
      const next = res.get(item.selfRef);
      if (next == null)
        throw new Error(`SelfRef ${item.selfRef} on path ${path.join('/')} not found`);
      const npath = [...path, next.id];
      return [next, npath] as const;
    }
    const npath = [...path, item.id];
    return [item, npath] as const;
  }

  const traverse = buildTraverseD([] as ItemPathArray, (item: ConfigTraverseItem, path: ItemPathArray) => {
    if (isPipelineSelfRef(item) || isPipelineStepConfig(item)) {
      return [] as [ConfigTraverseItem, ItemPathArray][];
    } else if(isPipelineStaticConfig(item)) {
      return item.steps.map((item) => getNextItem(item, path));
    }
    return item.stepTypes.map((item) => getNextItem(item, path));
  });

  traverse(config, (res, item, path) => {
    if (isPipelineConfig(item) && item.globalId) {
      if (res.has(item.globalId))
        throw new Error(`Duplicate globalId non-ref config ${item.globalId}, path ${path.join('/')}`);
      res.set(item.globalId, item);
    }
    return res;
  }, res);

  return res;
}

export function getConfigByInstancePath(instancePath: ItemPathArray, config: PipelineConfigurationProcessed, refMap: Map<string, PipelineConfigurationProcessed>) {
  let [startId, ...startPath] = instancePath;
  if (startId !== config.id)
    throw new Error(`Root segment ${startId} different id ${config.id}`);

  const findNextNode = (items: ConfigTraverseItem[], targetSegment: string, currentPath: ItemPathArray) => {
    for (const item of items) {
      if (isPipelineSelfRef(item)) {
        const nitem = refMap.get(item.selfRef)!;
        if (nitem.id === targetSegment)
          return nitem;
      } else {
        if (item.id === targetSegment)
          return item;
      }
    }
    throw new Error(`Segment ${targetSegment} not found on path ${currentPath}`);
  }

  const traverse = buildTraverseD([] as ItemPathArray, (item: ConfigTraverseItem, path: ItemPathArray, remainingPath?: ItemPathArray) => {
    const [segment, ...newRemainingPath] = remainingPath!;
    if (isPipelineStaticConfig(item)) {
      const node = findNextNode(item.steps, segment, path);
      return [[node, [...path, segment] as ItemPathArray, newRemainingPath] as const];
    } else if (isPipelineParallelConfig(item) || isPipelineSequentialConfig(item)) {
      const node = findNextNode(item.stepTypes, segment, path);
      return [[node, [...path, segment] as ItemPathArray, newRemainingPath] as const];
    } else if (isPipelineStepConfig(item)) {
      if (item.id !== segment)
        throw new Error(`Step segment ${segment} has different id ${item.id}`);
    } else if (isPipelineSelfRef(item)) {
      throw new Error(`Ref segment ${segment} was not mapped`);
    }
    if (newRemainingPath.length > 0)
      throw new Error(`Unmatched segments ${newRemainingPath}`);
    return [];
  }, startPath);

  const node = traverse(config, (_acc, node) => node, undefined as ConfigTraverseItem | undefined);
  return node as ConfigItem;
}

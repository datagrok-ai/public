import {ItemPathArray} from '../data/common-types';
import {buildTraverseD} from '../data/graph-traverse-utils';
import {addPipelineRef, containsPipelineRef, FuncCallIODescription, getPipelineRef, PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed, PipelineRefStore} from './config-processing-utils';
import {LinkIOParsed} from './LinkSpec';
import {PipelineSelfRef, PipelineStepConfiguration} from './PipelineConfiguration';


export type PipelineStepConfigurationProcessed = PipelineStepConfiguration<LinkIOParsed[], FuncCallIODescription[]>;
export type ConfigTraverseItem = PipelineConfigurationProcessed | PipelineStepConfigurationProcessed | PipelineSelfRef;

export type ConfigItem = PipelineConfigurationProcessed | PipelineStepConfigurationProcessed;

export function isPipelineStaticConfig(c: ConfigTraverseItem): c is PipelineConfigurationStaticProcessed {
  return !!((c as PipelineConfigurationStaticProcessed).type === 'static');
}

export function isPipelineParallelConfig(c: ConfigTraverseItem): c is PipelineConfigurationParallelProcessed {
  return !!((c as PipelineConfigurationParallelProcessed).type === 'parallel');
}

export function isPipelineSequentialConfig(c: ConfigTraverseItem): c is PipelineConfigurationSequentialProcessed {
  return !!((c as PipelineConfigurationSequentialProcessed).type === 'sequential');
}

export function isPipelineConfig(c: ConfigTraverseItem): c is PipelineConfigurationProcessed {
  return isPipelineStaticConfig(c) ||isPipelineParallelConfig(c) || isPipelineSequentialConfig(c);
}

export function isPipelineSelfRef(c: ConfigTraverseItem): c is PipelineSelfRef {
  return !!((c as PipelineSelfRef).type === 'selfRef');
}

export function isPipelineStepConfig(c: ConfigTraverseItem): c is PipelineStepConfigurationProcessed {
  return !isPipelineConfig(c) && !isPipelineSelfRef(c);
}

export type PipelineRefMap = PipelineRefStore<PipelineConfigurationProcessed>

export function buildRefMap(config: PipelineConfigurationProcessed): PipelineRefMap {
  const refMap: PipelineRefMap = new Map();

  function getNextItem(item: ConfigTraverseItem, path: ItemPathArray) {
    if (isPipelineStepConfig(item))
      return null;
    else if (isPipelineSelfRef(item)) {
      const next = getPipelineRef(refMap, item.nqName, item.version);
      if (next == null)
        throw new Error(`SelfRef nqName ${item.nqName} version ${item.version} on path ${path.join('/')} not found`);
      return null;
    }
    const npath = [...path, item.id];
    return [item, npath] as const;
  };

  const traverse = buildTraverseD([] as ItemPathArray, (item: ConfigTraverseItem, path: ItemPathArray) => {
    if (isPipelineSelfRef(item) || isPipelineStepConfig(item))
      return [] as [ConfigTraverseItem, ItemPathArray][];
    else if (isPipelineStaticConfig(item))
      return item.steps.map((item) => getNextItem(item, path)!).filter((x) => x);

    return item.stepTypes.map((item) => getNextItem(item, path)!).filter((x) => x);
  });

  traverse(config, (res, item) => {
    if (isPipelineConfig(item) && item.nqName) {
      if (!containsPipelineRef(res, item.nqName, item.version))
        addPipelineRef(res, item.nqName, item.version, item);
    }
    return res;
  }, refMap);

  return refMap;
}

export function getConfigByInstancePath(
  instancePath: ItemPathArray,
  config: PipelineConfigurationProcessed,
  refMap: PipelineRefMap,
) {
  if (instancePath.length == 0)
    return config;

  function findNextNode(
    items: ConfigTraverseItem[],
    targetSegment: string,
    currentPath: ItemPathArray,
    refMap: PipelineRefMap,
  ) {
    for (const item of items) {
      if (isPipelineSelfRef(item)) {
        const nitem = getPipelineRef(refMap, item.nqName, item.version)!;
        if (nitem.id === targetSegment)
          return nitem;
      } else {
        if (item.id === targetSegment)
          return item;
      }
    }
    throw new Error(`Segment ${targetSegment} not found on path ${currentPath.join('/')}`);
  }

  const traverse = buildTraverseD([] as ItemPathArray, (item: ConfigTraverseItem, path: ItemPathArray, remainingPath?: ItemPathArray) => {
    if (remainingPath?.length == 0)
      return [];
    const [segment, ...newRemainingPath] = remainingPath!;
    if (isPipelineStaticConfig(item)) {
      const node = findNextNode(item.steps, segment, path, refMap);
      return [[node, [...path, segment] as ItemPathArray, newRemainingPath] as const];
    } else if (isPipelineParallelConfig(item) || isPipelineSequentialConfig(item)) {
      const node = findNextNode(item.stepTypes, segment, path, refMap);
      return [[node, [...path, segment] as ItemPathArray, newRemainingPath] as const];
    } else if (isPipelineSelfRef(item))
      throw new Error(`Ref segment ${segment} was not mapped`);

    if (newRemainingPath.length > 0)
      throw new Error(`Unmatched segments ${newRemainingPath}`);
    return [];
  }, instancePath);

  const node = traverse(config, (_acc, node) => node, undefined as ConfigTraverseItem | undefined);
  return node as ConfigItem;
}

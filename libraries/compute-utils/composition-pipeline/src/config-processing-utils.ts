import cloneDeepWith from 'lodash.clonedeepwith';
import {PipelineConfiguration, PipelineCompositionConfiguration, ItemName, ItemPath} from './PipelineConfiguration';

//
// Internal config processing
//

export type NestedPipelineConfiguration = PipelineCompositionConfiguration & GraphNestedPipelines

export type CompositionGraphConfig = PipelineConfiguration | NestedPipelineConfiguration;

export type PipelineConfigVariants = PipelineConfiguration | PipelineCompositionConfiguration;

type GraphNestedPipelines = {
  nestedPipelines?: {
    [key: ItemName]: CompositionGraphConfig;
  };
};

export function isNestedPipelineConfig(config: CompositionGraphConfig): config is NestedPipelineConfiguration {
  return (config as any)?.nestedPipelines;
}

function isPipelineConfig(config: CompositionGraphConfig): config is PipelineConfiguration {
  return !(config as any)?.nestedPipelines;
}

export function cloneConfig<T>(config: T): T {
  return cloneDeepWith(config, (val) => {
    if (typeof val === 'function')
      return val;
  });
}

export type PathKey = string;
export function pathJoin(path: ItemPath, ...restPaths: ItemPath[]): ItemPath {
  return path.concat(...restPaths);
}

export function pathToKey(path: ItemPath): PathKey {
  return path.filter((x) => x).join('/');
}

export function keyToPath(key: PathKey) {
  return key.split('/').filter((x) => x);
}

export function getParentKey(key: PathKey): PathKey {
  const path = keyToPath(key);
  const ppath = path.slice(0, path.length - 1);
  return pathToKey(ppath);
}

export function getSuffix(key: PathKey, prefix: PathKey) {
  if (key.startsWith(prefix))
    return key.substring(prefix.length);
}

export function normalizePaths(paths: ItemPath | ItemPath[]): ItemPath[] {
  if (Array.isArray(paths[0]))
    return paths as ItemPath[];

  else
    return [paths as ItemPath];
}

export function traverseConfigPipelines<T>(
  graph: CompositionGraphConfig,
  nodeHandler: (
    acc: T,
    node: CompositionGraphConfig,
    path: ItemPath) => T,
  acc: T,
) {
  const stk: [{
    node: CompositionGraphConfig;
    path: ItemPath;
  }] = [{
    node: graph,
    path: [] as ItemPath,
  }];

  const queuedNested = new Set<object>();

  while (stk.length) {
    const {node, path} = stk.pop()!;
    if (isNestedPipelineConfig(node)) {
      if (queuedNested.has(node))
        acc = nodeHandler(acc, node, path);
      else {
        stk.push({node: node, path});
        for (const nextNode of Object.values(node.nestedPipelines ?? {}))
          stk.push({node: nextNode, path: [...path, node.id]});
        queuedNested.add(node);
      }
    } else if (isPipelineConfig(node))
      acc = nodeHandler(acc, node, path);
  }

  return acc;
}

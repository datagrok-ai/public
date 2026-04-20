import * as DG from 'datagrok-api/dg';
import {Observable, defer, of, merge, from} from 'rxjs';
import {map, mapTo, toArray, concatMap} from 'rxjs/operators';
import {NodePath, BaseTree, TreeNode} from '../data/BaseTree';
import {getPipelineRef, PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {isFuncCallSerializedState, PipelineInstanceConfig, PipelineSerializedState} from '../config/PipelineInstance';
import {buildTraverseD} from '../data/graph-traverse-utils';
import {buildRefMap, ConfigTraverseItem, getConfigByInstancePath, isPipelineDynamicConfig, isPipelineSelfRef, isPipelineStaticConfig, isPipelineStepConfig, PipelineRefMap, PipelineStepConfigurationProcessed} from '../config/config-utils';
import {FuncCallAdapter, FuncCallMockAdapter} from './FuncCallAdapters';
import {loadFuncCall, loadInstanceState, makeFuncCall} from './funccall-utils';
import {DynamicPipelineNode, FuncCallNode, isFuncCallNode, StateTreeNode, StaticPipelineNode} from './StateTreeNodes';
import {indexFromEnd} from '../utils';
import {DriverLogger} from '../data/Logger';
import {StateTree} from './StateTree';

const MAX_CONCURRENT_LOADS = 5;

export function loadStateTree({
  dbId,
  config,
  isReadonly = false,
  defaultValidators = false,
  mockMode = false,
}: {
  dbId: string,
  config: PipelineConfigurationProcessed,
  isReadonly?: boolean,
  defaultValidators?: boolean
  mockMode?: boolean
}): Observable<StateTree> {
  return defer(async () => {
    const [state, _] = await loadInstanceState(dbId);
    const tree = fromPipelineInstanceState({state, config, isReadonly, defaultValidators, mockMode});
    return tree;
  });
}

export function fromPipelineConfig({
  config,
  startNode = config,
  startPath = [],
  startState,
  isReadonly = false,
  defaultValidators = false,
  mockMode = false,
  logger,
} : {
  config: PipelineConfigurationProcessed;
  startNode?: PipelineConfigurationProcessed;
  startPath?: Readonly<NodePath>;
  startState?: StateTree;
  isReadonly?: boolean;
  defaultValidators?: boolean
  mockMode?: boolean;
  logger?: DriverLogger
}): StateTree {
  const refMap = buildRefMap(config);

  // TODO: initial infinite cycles detection
  const traverse = buildTraverseD(startPath, (data: ConfigTraverseItem, path) => {
    if (isPipelineDynamicConfig(data)) {
      const items = (data?.initialSteps ?? []).map((step, idx) => {
        const item = data.stepTypes.find((t) => {
          if (isPipelineSelfRef(t)) {
            const derefedStep = getPipelineRef(refMap, t.nqName, t.version);
            return derefedStep?.id === step.id;
          }
          return t.id === step.id;
        });
        if (!item)
          throw new Error(`Node ${step.id} not found on path ${JSON.stringify(path)}`);
        const nextPath = [...path, {id: step.id, idx}];
        const stepItem = {...step, ...item};
        return [stepItem, nextPath] as const;
      });
      return items;
    } else if (isPipelineStaticConfig(data))
      return data.steps.map((step, idx) => [step, [...path, {id: step.id, idx}]] as const);
    else if (isPipelineSelfRef(data)) {
      const next = getPipelineRef(refMap, data.nqName, data.version);
      if (!next)
        throw new Error(`Failed to deref nqName ${data.nqName} version ${data.version} on path ${JSON.stringify(path)}`);

      return [[next, [...path, {id: next.id, idx: 0}]]] as const;
    }
    return [] as const;
  }, new Set<ConfigTraverseItem>());

  const tree = traverse(startNode, (acc, state, path) => {
    const [node, ppath, idx] = makeTreeNode(config, refMap, path, isReadonly, logger);
    if (isFuncCallNode(node)) {
      if (!isPipelineStepConfig(state))
        throw new Error(`Wrong FuncCall node state type ${state.type} on path ${JSON.stringify(path)}`);
      node.initState(state);
    }
    return addTreeNodeOrCreate({acc, config, node, ppath, pos: idx, defaultValidators, mockMode, logger});
  }, startState);
  return tree!;
}

export function fromPipelineInstanceState({
  state,
  config,
  isReadonly,
  mockMode = false,
  defaultValidators = false,
  logger,
}: {
  state: PipelineSerializedState;
  config: PipelineConfigurationProcessed;
  isReadonly: boolean,
  defaultValidators?: boolean,
  mockMode?: boolean,
  logger?: DriverLogger
}): StateTree {
  const refMap = buildRefMap(config);

  const traverse = buildTraverseD([] as Readonly<NodePath>, (item: PipelineSerializedState, path) => {
    if (isFuncCallSerializedState(item))
      return [];
    else
      return item.steps.map((step, idx) => [step, [...path, {id: step.configId, idx}]] as const);
  });

  const tree = traverse(state, (acc, state, path) => {
    const [node, ppath, idx] = makeTreeNode(config, refMap, path, isReadonly || state.isReadonly, logger);
    node.restoreState(state);
    return addTreeNodeOrCreate({acc, config, node, ppath, pos: idx, defaultValidators, mockMode, logger});
  }, undefined as StateTree | undefined);
  return tree!;
}

export function fromPipelineInstanceConfig({
  instanceConfig,
  config,
  isReadonly = false,
  defaultValidators = false,
  mockMode = false,
  logger,
}: {
  instanceConfig: PipelineInstanceConfig,
  config: PipelineConfigurationProcessed,
  isReadonly?: boolean,
  defaultValidators?: boolean,
  mockMode?: boolean,
  logger?: DriverLogger
}): StateTree {
  const refMap = buildRefMap(config);

  // TODO: fix static pipeline missing steps
  const traverse = buildTraverseD([] as Readonly<NodePath>, (data: PipelineInstanceConfig, path, visited) => {
    if (visited!.has(data))
      throw new Error(`Initial config cycle on node ${data.id} path ${JSON.stringify(path)}`);
    visited!.add(data);
    if (data.steps)
      return data.steps.map((step, idx) => [step, [...path, {id: step.id, idx}], visited] as const);
    else
      return [];
  }, new Set<PipelineInstanceConfig>());

  const tree = traverse(instanceConfig, (acc, state, path) => {
    const [node, ppath, idx] = makeTreeNode(config, refMap, path, isReadonly, logger);
    if (isFuncCallNode(node))
      node.initState(state);
    return addTreeNodeOrCreate({acc, config, node, ppath, pos: idx, defaultValidators, mockMode, logger});
  }, undefined as StateTree | undefined);
  return tree!;
}

export function loadOrCreateCalls(stateTree: StateTree, mockMode: boolean) {
  return defer(() => {
    const pendingFuncCallLoads = stateTree.nodeTree.traverse(stateTree.nodeTree.root, (acc, node) => {
      const item = node.getItem();
      if (isFuncCallNode(item)) {
        if (item.instancesWrapper.getInstance())
          return acc;
        else if (mockMode) {
          const obs$ = defer(() => {
            const adapter = new FuncCallMockAdapter(item.config.io!, item.isReadonly);
            item.initAdapter({adapter, restrictions: {}, isOutputOutdated: true, runError: undefined}, true);
            return of(undefined);
          });
          return [...acc, obs$];
        } else if (item.instancesWrapper.id == null && item.pendingId == null) {
          const obs$ = defer(() => {
            return from(makeFuncCall(item.config.nqName, false)).pipe(
              map((data) => item.initAdapter(data, true)),
            );
          });
          return [...acc, obs$];
        } else if (item.instancesWrapper.id == null && item.pendingId) {
          const savedId = item.pendingId;
          const obs$ = defer(() => {
            return from(loadFuncCall(savedId, item.isReadonly)).pipe(
              map((data) => item.initAdapter(data, false)),
            );
          });
          return [...acc, obs$];
        }
      }
      return acc;
    }, [] as Array<Observable<undefined | void>>);
    return merge(...pendingFuncCallLoads, MAX_CONCURRENT_LOADS).pipe(
      toArray(),
      mapTo(stateTree),
    );
  });
}

export function findPipelineNode(stateTree: StateTree, uuid?: string) {
  const data = stateTree.nodeTree.find((item) => uuid == null || item.uuid === uuid);
  if (data == null)
    throw new Error(`Node uuid ${uuid} not found`);
  const [root, path] = data;
  const item = root.getItem();
  if (isFuncCallNode(item))
    throw new Error(`FuncCall node ${uuid}, but pipeline is expected`);

  const nqName = root.getItem().config.nqName;
  return [root, nqName, path] as const;
}

export function getSubConfig(config: PipelineConfigurationProcessed, path: NodePath, id: string) {
  const refMap = buildRefMap(config);
  const subConfPath = [...path.map((s) => s.id), id];
  const subConfig = getConfigByInstancePath(subConfPath, config, refMap);
  return subConfig;
}

function makeTreeNode(
  config: PipelineConfigurationProcessed,
  refMap: PipelineRefMap,
  path: Readonly<NodePath>,
  isReadonly: boolean,
  logger?: DriverLogger,
) {
  const confPath = path.map((p) => p.id);
  const nodeConf = getConfigByInstancePath(confPath, config, refMap);
  const node = makeNode(nodeConf, isReadonly, logger);
  const ppath = path.slice(0, -1);
  const idx = indexFromEnd(path)?.idx ?? 0;
  return [node, ppath, idx] as const;
}

function makeNode(
  nodeConf: PipelineConfigurationProcessed | PipelineStepConfigurationProcessed,
  isReadonly: boolean,
  logger?: DriverLogger,
) {
  if (isPipelineStepConfig(nodeConf))
    return new FuncCallNode(nodeConf, isReadonly, logger);
  else if (isPipelineStaticConfig(nodeConf))
    return new StaticPipelineNode(nodeConf, isReadonly, logger);
  else if (isPipelineDynamicConfig(nodeConf))
    return new DynamicPipelineNode(nodeConf, isReadonly, logger);

  throw new Error(`Wrong node type ${nodeConf}`);
}

function addTreeNodeOrCreate(
  {
    acc,
    config,
    node,
    ppath,
    pos,
    defaultValidators,
    mockMode,
    logger,
  } : {
    acc: StateTree | undefined;
    config: PipelineConfigurationProcessed;
    node: StateTreeNode;
    ppath: NodePath;
    pos: number;
    defaultValidators: boolean;
    mockMode: boolean;
    logger?: DriverLogger
  },
) {
  if (acc)
    acc.nodeTree.addItem(ppath, node, node.config.id, pos);
  else
    return new StateTree(node, config, mockMode, defaultValidators, logger);
  return acc;
}

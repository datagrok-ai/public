import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, defer, identity, of, merge, Subject} from 'rxjs';
import {delay, finalize, map, mapTo, toArray, concatMap} from 'rxjs/operators';
import {NodePath, BaseTree, TreeNode, NodeAddressSegment} from '../data/BaseTree';
import {PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {isFuncCallState, PipelineInstanceConfig, PipelineSerializedState, PipelineState} from '../config/PipelineInstance';
import {buildTraverseD} from '../data/traversable';
import {buildRefMap, ConfigTraverseItem, getConfigByInstancePath, isPipelineParallelConfig, isPipelineSelfRef, isPipelineSequentialConfig, isPipelineStaticConfig, isPipelineStepConfig, PipelineStepConfigurationProcessed} from '../config/config-utils';
import {FuncCallAdapter, FuncCallMockAdapter, IFuncCallAdapter} from './FuncCallAdapters';
import {loadFuncCall, loadInstanceState, makeFuncCall, saveFuncCall, saveInstanceState} from './adapter-utils';
import {FuncCallNode, isFuncCallNode, ParallelPipelineNode, PipelineNodeBase, SequentialPipelineNode, StateTreeNode, StaticPipelineNode} from './StateTreeNodes';
import {indexFromEnd} from '../utils';

export type StateTreeSaveOptions = {
  isSerialized?: boolean,
  disableUUID?: boolean,
}

export class StateTree extends BaseTree<StateTreeNode> {
  public makeStateRequests = new Subject<true>();

  constructor(
    item: StateTreeNode,
    private config: PipelineConfigurationProcessed,
    private mockMode = false,
  ) {
    super(item);
  }

  public toSerializedState(disableUUID = false): PipelineSerializedState {
    return this.toStateRec(this.getRoot(), {isSerialized: true, disableUUID});
  }

  public toState(disableUUID = false): PipelineState {
    return this.toStateRec(this.getRoot(), {isSerialized: false, disableUUID});
  }

  public save(uuid?: string, mockDelay?: number): Observable<string> {
    let rootItem: PipelineNodeBase;
    return defer(() => {
      this.globalStructureLock();
      const [root, nqName] = this.findPipelineRootNode(uuid);
      if (nqName == null)
        throw new Error(`Attempting to save pipeline with no nqName`);
      rootItem = root.getItem() as PipelineNodeBase;
      if (rootItem.config.provider == null)
        throw new Error(`Attempting to save pipeline with no nqName provider`);
      const pendingFuncCallSaves = this.traverse(root, (acc, node) => {
        const item = node.getItem();
        if (isFuncCallNode(item)) {
          if (item.instancesWrapper.isRunning$.value || item.instancesWrapper.isLoading$.value)
            throw new Error(`FuncCall node ${item.uuid} saving during being updated`);
          const bridge = item.instancesWrapper.getInstance();
          if (bridge == null)
            throw new Error(`FuncCall node ${item.uuid} has no instance`);
          item.instancesWrapper.isLoading$.next(true);
          if (this.mockMode) {
            const obs$ = of(bridge).pipe(mockDelay ? delay(mockDelay) : identity).pipe(
              map((adapter) => [item, adapter] as const),
            );
            return [...acc, obs$];
          } else {
            const obs$ = defer(() => saveFuncCall(bridge)).pipe(map((nextCall) => [item, new FuncCallAdapter(nextCall)] as const));
            return [...acc, obs$];
          }
        } else {
          if (item.isUpdating$.value && item !== rootItem)
            throw new Error(`Pipeline node ${item.uuid} saving during being updated`);
        }
        return acc;
      }, [] as Array<Observable<Readonly<[FuncCallNode, IFuncCallAdapter]>>>);
      this.emitNewState();
      return merge(...pendingFuncCallSaves, 5).pipe(
        map(([node, nAdapter]) => {
          node.setAdapter(nAdapter);
          node.instancesWrapper.isLoading$.next(false);
        }),
        toArray(),
        concatMap(() => {
          if (this.mockMode)
            return of('');
          const state = this.toStateRec(root, {isSerialized: true});
          const json = JSON.stringify(state);
          return defer(() => saveInstanceState(nqName, json));
        }),
      );
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
      if (rootItem)
        rootItem.isUpdating$.next(false);
      this.emitNewState();
    }));
  }

  public init() {
    return defer(() => {
      this.globalStructureLock();
      return this.loadOrCreateCalls(this);
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
      this.emitNewState();
    }));
  }

  public moveMoveSubtree(uuid: string, pos: number) {
    return defer(() => {
      this.globalStructureLock();
      const data = this.find((item) => item.uuid === uuid);
      if (data == null)
        throw new Error(`Node uuid ${uuid} not found`);
      const [node, path] = data;
      this.removeBrunch(path);
      const ppath = path.slice(0, -1);
      this.attachBrunch(ppath, node, node.getItem().config.id, pos);
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
      this.emitNewState();
    }));
  }

  public removeSubtree(uuid: string) {
    return defer(() => {
      this.globalStructureLock();
      const data = this.find((item) => item.uuid === uuid);
      if (data == null)
        throw new Error(`Node uuid ${uuid} not found`);
      const [, path] = data;
      this.removeBrunch(path);
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
      this.emitNewState();
    }));
  }

  public addSubTree(puuid: string, id: string, pos: number) {
    return defer(() => {
      this.globalStructureLock();
      const refMap = buildRefMap(this.config);
      const [_root, _nqName, path] = this.findPipelineRootNode(puuid);
      const subConfPath = [...path.map((s) => s.id), id];
      const subConfig = getConfigByInstancePath(subConfPath, this.config, refMap);
      if (isPipelineStepConfig(subConfig))
        this.attachBrunch(path, new TreeNode(new FuncCallNode(subConfig)), id, pos);
      else
        StateTree.fromConfig(this.config, subConfig, [...path, {id, idx: pos}], this, this.mockMode);
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
      this.emitNewState();
    }));
  }

  public loadSubTree(puuid: string, stateId: string, id: string, pos: number) {
    return defer(() => {
      this.globalStructureLock();
      const [_root, _nqName, path] = this.findPipelineRootNode(puuid);
      const subPath = [...path, {id, idx: pos}];
      const tree = StateTree.load(stateId, this.config, subPath, this.mockMode).pipe(
        concatMap((tree) => this.loadOrCreateCalls(tree)),
      ).pipe(
        map((tree) => this.attachBrunch(path, tree.root, id, pos)),
      );
      return tree;
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
      this.emitNewState();
    }));
  }

  static load(id: string, config: PipelineConfigurationProcessed, startPath: Readonly<NodePath> = [], mockMode = false): Observable<StateTree> {
    return defer(async () => {
      const state = await loadInstanceState(id);
      const tree = StateTree.fromState(state, config, startPath, undefined, mockMode);
      return tree;
    });
  }

  static fromConfig(config: PipelineConfigurationProcessed, startNode = config, startPath: Readonly<NodePath> = [], currentState?: StateTree, mockMode = false): StateTree {
    const refMap = buildRefMap(config);

    // TODO: initial infinite cycles detection
    const traverse = buildTraverseD(startPath, (data: ConfigTraverseItem, path) => {
      if (isPipelineParallelConfig(data) || isPipelineSequentialConfig(data)) {
        const items = (data?.initialSteps ?? []).map((step, idx) => {
          const item = data.stepTypes.find((t) => {
            if (isPipelineSelfRef(t)) {
              const derefedStep = refMap.get(t.selfRef);
              return derefedStep?.id === step.id;
            }
            return t.id === step.id;
          });
          if (!item)
            throw new Error(`Node ${step.id} not found on path ${JSON.stringify(path)}`);
          const nextPath = [...path, {id: step.id, idx}];
          return [item, nextPath] as const;
        });
        return items;
      } else if (isPipelineStaticConfig(data))
        return data.steps.map((step, idx) => [step, [...path, {id: step.id, idx}]] as const);
      else if (isPipelineSelfRef(data)) {
        const next = refMap.get(data.selfRef);
        if (!next)
          throw new Error(`Failed to deref ${data.selfRef} on path ${JSON.stringify(path)}`);

        return [[next, [...path, {id: next.id, idx: 0}]]] as const;
      }
      return [] as const;
    }, new Set<ConfigTraverseItem>());

    const tree = traverse(startNode, (acc, state, path) => {
      const [node, ppath, idx] = StateTree.makeTreeNode(config, refMap, path);
      if (isFuncCallNode(node)) {
        if (!isPipelineStepConfig(state))
          throw new Error(`Wrong FuncCall node state type ${state.type} on path ${JSON.stringify(path)}`);
        node.initState(state);
      }
      return StateTree.addTreeNodeOrCreate(acc, config, node, ppath, idx, mockMode);
    }, currentState);
    return tree!;
  }

  static fromState(state: PipelineSerializedState, config: PipelineConfigurationProcessed, startPath: Readonly<NodePath> = [], currentState?: StateTree, mockMode = false): StateTree {
    const refMap = buildRefMap(config);

    const traverse = buildTraverseD(startPath, (item: PipelineState, path) => {
      if (isFuncCallState(item))
        return [];
      else
        return item.steps.map((step, idx) => [step, [...path, {id: step.configId, idx}]] as const);
    });

    const tree = traverse(state, (acc, state, path) => {
      const [node, ppath, idx] = StateTree.makeTreeNode(config, refMap, path);
      if (isFuncCallNode(node)) {
        if (!isFuncCallState(state))
          throw new Error(`Wrong FuncCall state type ${state.type} on path ${JSON.stringify(path)}`);
        node.restoreState(state);
      }
      return StateTree.addTreeNodeOrCreate(acc, config, node, ppath, idx, mockMode);
    }, currentState);
    return tree!;
  }

  static fromInstanceConfig(instanceConfig: PipelineInstanceConfig, config: PipelineConfigurationProcessed, mockMode = false): StateTree {
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
      const [node, ppath, idx] = StateTree.makeTreeNode(config, refMap, path);
      if (isFuncCallNode(node))
        node.initState(state);
      return StateTree.addTreeNodeOrCreate(acc, config, node, ppath, idx, mockMode);
    }, undefined as StateTree | undefined);
    return tree!;
  }

  private static makeTreeNode(config: PipelineConfigurationProcessed, refMap: Map<string, PipelineConfigurationProcessed>, path: Readonly<NodePath>) {
    const confPath = path.map((p) => p.id);
    const nodeConf = getConfigByInstancePath(confPath, config, refMap);
    if (nodeConf == null)
      throw new Error(`Unable to find configuration node on path ${JSON.stringify(path)}`);
    const node = StateTree.makeNode(nodeConf);
    const ppath = path.slice(0, -1);
    const idx = indexFromEnd(path)?.idx ?? 0;
    return [node, ppath, idx] as const;
  }

  private static addTreeNodeOrCreate(acc: StateTree | undefined, config: PipelineConfigurationProcessed, node: StateTreeNode, ppath: NodePath, pos: number, mockMode: boolean) {
    if (acc)
      acc.addItem(ppath, node, node.config.id, pos);
    else
      return new StateTree(node, config, mockMode);
    return acc;
  }

  private findPipelineRootNode(uuid?: string) {
    const data = this.find((item) => uuid == null || item.uuid === uuid);
    if (data == null)
      throw new Error(`Node uuid ${uuid} not found`);
    const [root, path] = data;
    const item = root.getItem();
    if (isFuncCallNode(item))
      throw new Error(`FuncCall node ${uuid} cannot be pipeline root`);

    const nqName = root.getItem().config.nqName;
    return [root, nqName, path] as const;
  }

  private loadOrCreateCalls(tree: StateTree) {
    return defer(() => {
      const pendingFuncCallLoads = tree.traverse(tree.root, (acc, node) => {
        const item = node.getItem();
        if (isFuncCallNode(item)) {
          if (this.mockMode) {
            const obs$ = of([item, new FuncCallMockAdapter(item.config.io!), false] as const);
            return [...acc, obs$];
          } else if (item.instancesWrapper.id == null) {
            const obs$ = defer(() => {
              const savedId = item.pendingId;
              if (savedId == null)
                return makeFuncCall(item.config.nqName);
              else {
                item.pendingId = undefined;
                return loadFuncCall(savedId);
              }
            }).pipe(
              map((nitem) => [item, nitem, true] as const),
            );
            return [...acc, obs$];
          }
        }
        return acc;
      }, [] as Array<Observable<Readonly<[FuncCallNode, IFuncCallAdapter, boolean]>>>);
      return merge(...pendingFuncCallLoads, 5).pipe(
        map(([node, nAdapter, isNew]) => {
          node.setAdapter(nAdapter, isNew);
          node.instancesWrapper.isLoading$.next(false);
        }),
        toArray(),
        mapTo(tree),
      );
    });
  }

  private emitNewState() {
    this.makeStateRequests.next(true);
  }

  private globalStructureLock() {
    const item = this.root.getItem() as PipelineNodeBase;
    if (item.isUpdating$.value)
      throw new Error(`Global structure double lock`);
    item.setLoading(true);
  }

  private globalStructureUnlock() {
    const item = this.root.getItem() as PipelineNodeBase;
    if (!item.isUpdating$.value)
      throw new Error(`Global structure double unlock`);
    item.setLoading(false);
  }

  private static makeNode(nodeConf: PipelineConfigurationProcessed | PipelineStepConfigurationProcessed) {
    if (isPipelineStepConfig(nodeConf))
      return new FuncCallNode(nodeConf);
    else if (isPipelineStaticConfig(nodeConf))
      return new StaticPipelineNode(nodeConf);
    else if (isPipelineParallelConfig(nodeConf))
      return new ParallelPipelineNode(nodeConf);
    else if (isPipelineSequentialConfig(nodeConf))
      return new SequentialPipelineNode(nodeConf);

    throw new Error(`Wrong node type ${nodeConf}`);
  }

  private toStateRec(node: TreeNode<StateTreeNode>, options: StateTreeSaveOptions = {}): PipelineState {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return options.isSerialized ? item.toSerializedState(options.disableUUID) : item.toState(options.disableUUID);

    const state = item.toState(options.disableUUID);
    const steps = node.getChildren().map((node) => {
      const item = this.toStateRec(node.item, options);
      return item;
    });
    return {...state, steps} as PipelineState;
  }
}

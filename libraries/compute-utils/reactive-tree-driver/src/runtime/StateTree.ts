import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, defer, identity, of, merge, Subject, BehaviorSubject} from 'rxjs';
import {delay, finalize, map, mapTo, toArray, concatMap, filter, tap} from 'rxjs/operators';
import {NodePath, BaseTree, TreeNode} from '../data/BaseTree';
import {PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {isFuncCallState, PipelineInstanceConfig, PipelineSerializedState, PipelineState} from '../config/PipelineInstance';
import {buildTraverseD} from '../data/traversable';
import {buildRefMap, ConfigTraverseItem, getConfigByInstancePath, isPipelineParallelConfig, isPipelineSelfRef, isPipelineSequentialConfig, isPipelineStaticConfig, isPipelineStepConfig, PipelineStepConfigurationProcessed} from '../config/config-utils';
import {FuncCallAdapter, FuncCallMockAdapter, IFuncCallAdapter} from './FuncCallAdapters';
import {loadFuncCall, loadInstanceState, makeFuncCall, makeMetaCall, saveFuncCall, saveInstanceState} from './adapter-utils';
import {FuncCallNode, isFuncCallNode, ParallelPipelineNode, PipelineNodeBase, SequentialPipelineNode, StateTreeNode, StateTreeSerializationOptions, StaticPipelineNode} from './StateTreeNodes';
import {indexFromEnd} from '../utils';

const MAX_CONCURENT_SAVES = 5;

export class StateTree extends BaseTree<StateTreeNode> {
  public makeStateRequests = new Subject<true>();
  public metaCall$ = new BehaviorSubject<DG.FuncCall | undefined>(undefined);

  constructor(
    item: StateTreeNode,
    private config: PipelineConfigurationProcessed,
    private mockMode = false,
  ) {
    super(item);
  }

  public toSerializedState(options: StateTreeSerializationOptions = {}): PipelineSerializedState {
    return StateTree.toStateRec(this.getRoot(), true, options);
  }

  public toState(options: StateTreeSerializationOptions = {}): PipelineState {
    return StateTree.toStateRec(this.getRoot(), false, options);
  }

  public save(uuid?: string, mockDelay?: number) {
    return this.updateWithLock(() => {
      const [root, nqName] = StateTree.findPipelineNode(this, uuid);
      if (nqName == null)
        throw new Error(`Attempting to save pipeline with no nqName`);
      const rootItem = root.getItem() as PipelineNodeBase;
      if (rootItem.config.provider == null)
        throw new Error(`Attempting to save pipeline with no nqName provider`);
      const isRoot = root === this.root;
      const pendingFuncCallSaves = this.traverse(root, (acc, node) => {
        const item = node.getItem();
        if (isFuncCallNode(item)) {
          if (item.instancesWrapper.isRunning$.value || item.instancesWrapper.isLoading$.value)
            throw new Error(`FuncCall node ${item.uuid} saving during being updated`);
          const inst = item.instancesWrapper.getInstance();
          if (inst == null)
            throw new Error(`FuncCall node ${item.uuid} has no instance`);
          item.instancesWrapper.isLoading$.next(true);
          if (this.mockMode) {
            const obs$ = of(inst).pipe(mockDelay ? delay(mockDelay) : identity).pipe(
              map((adapter) => [item, adapter] as const),
            );
            return [...acc, obs$];
          } else {
            const obs$ = defer(() => saveFuncCall(inst)).pipe(map((nextCall) => [item, new FuncCallAdapter(nextCall)] as const));
            return [...acc, obs$];
          }
        } else {
          if (item.isUpdating$.value && item !== rootItem)
            throw new Error(`Pipeline node ${item.uuid} saving during being updated`);
        }
        return acc;
      }, [] as Array<Observable<Readonly<[FuncCallNode, IFuncCallAdapter]>>>);
      return merge(...pendingFuncCallSaves, MAX_CONCURENT_SAVES).pipe(
        map(([node, nAdapter]) => {
          node.setAdapter(nAdapter);
          node.instancesWrapper.isLoading$.next(false);
        }),
        toArray(),
        concatMap(() => this.saveMetaCall(root, nqName, isRoot ? this.metaCall$.value : undefined)),
        tap((call) => isRoot ? this.metaCall$.next(call) : void(0)),
      );
    });
  }

  public initFuncCalls() {
    return this.updateWithLock(() => {
      return StateTree.loadOrCreateCalls(this, this.mockMode);
    });
  }

  public initMetaCall() {
    return this.updateWithLock(() => {
      return this.makeMetaCall().pipe(
        tap((call) => this.metaCall$.next(call)),
      );
    });
  }

  public moveSubtree(uuid: string, pos: number) {
    return this.updateWithLock(() => {
      const data = this.find((item) => item.uuid === uuid);
      if (data == null)
        throw new Error(`Node uuid ${uuid} not found`);
      const [node, path] = data;
      this.removeBrunch(path);
      const ppath = path.slice(0, -1);
      this.attachBrunch(ppath, node, node.getItem().config.id, pos);
      return of(this);
    });
  }

  public removeSubtree(uuid: string) {
    return this.updateWithLock(() => {
      const data = this.find((item) => item.uuid === uuid);
      if (data == null)
        throw new Error(`Node uuid ${uuid} not found`);
      const [, path] = data;
      this.removeBrunch(path);
      return of(this);
    });
  }

  public addSubTree(puuid: string, id: string, pos: number, initCalls = true) {
    return this.updateWithLock(() => {
      const [_root, _nqName, path] = StateTree.findPipelineNode(this, puuid);
      const subConfig = this.getSubConfig(path, id);
      if (isPipelineStepConfig(subConfig))
        this.attachBrunch(path, new TreeNode(new FuncCallNode(subConfig)), id, pos);
      else
        StateTree.fromConfig({config: this.config, startNode: subConfig, startPath: [...path, {id, idx: pos}], startState: this, mockMode: this.mockMode});
      return initCalls ? StateTree.loadOrCreateCalls(this, this.mockMode) : of(this);
    });
  }

  public loadSubTree(puuid: string, dbId: string, id: string, pos: number) {
    return this.updateWithLock(() => {
      const [_root, _nqName, path] = StateTree.findPipelineNode(this, puuid);
      const subConfig = this.getSubConfig(path, id);
      if (isPipelineStepConfig(subConfig))
        throw new Error(`FuncCall node ${JSON.stringify(path)}, but pipeline is expected`);
      const tree = StateTree.load(dbId, subConfig, this.mockMode).pipe(
        concatMap((tree) => StateTree.loadOrCreateCalls(tree, this.mockMode)),
      ).pipe(
        map((tree) => this.attachBrunch(path, tree.root, id, pos)),
      );
      return tree.pipe(mapTo(this));
    });
  }

  public runStep(uuid: string, mockResults?: Record<string, any>, mockDelay?: number) {
    if (!this.mockMode && mockResults)
      throw new Error(`Mock results passed to runStep while mock mode is off`);
    const res = this.find((item) => item.uuid === uuid);
    if (!res)
      throw new Error(`Step uuid ${uuid} not found`);
    const [snode] = res;
    const node = snode.getItem();
    if (!isFuncCallNode(node))
      throw new Error(`Step uuid ${uuid} is not FuncCall`);
    return node.instancesWrapper.run(mockResults, mockDelay);
  }

  public close() {
    this.traverse(this.root, (acc, node) => {
      const item = node.getItem();
      if (isFuncCallNode(item))
        item.instancesWrapper.close();
      return acc;
    }, undefined);
  }

  static load(dbId: string, config: PipelineConfigurationProcessed, mockMode = false): Observable<StateTree> {
    return defer(async () => {
      const [metaCall, state] = await loadInstanceState(dbId);
      const tree = StateTree.fromState({state, config, metaCall, mockMode});
      return tree;
    });
  }

  static fromConfig({config, startNode = config, startPath = [], startState, mockMode = false}: { config: PipelineConfigurationProcessed; startNode?: PipelineConfigurationProcessed; startPath?: Readonly<NodePath>; startState?: StateTree; mockMode?: boolean; }): StateTree {
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
    }, startState);
    return tree!;
  }

  static fromState({state, config, metaCall, mockMode = false}: { state: PipelineSerializedState; config: PipelineConfigurationProcessed; metaCall: DG.FuncCall; mockMode?: boolean }): StateTree {
    const refMap = buildRefMap(config);

    const traverse = buildTraverseD([] as Readonly<NodePath>, (item: PipelineState, path) => {
      if (isFuncCallState(item))
        return [];
      else
        return item.steps.map((step, idx) => [step, [...path, {id: step.configId, idx}]] as const);
    });

    const tree = traverse(state, (acc, state, path) => {
      const [node, ppath, idx] = StateTree.makeTreeNode(config, refMap, path);
      node.restoreState(state);
      return StateTree.addTreeNodeOrCreate(acc, config, node, ppath, idx, mockMode);
    }, undefined as StateTree | undefined);
    tree!.metaCall$.next(metaCall);
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

  private static loadOrCreateCalls(tree: StateTree, mockMode: boolean) {
    return defer(() => {
      const pendingFuncCallLoads = tree.traverse(tree.root, (acc, node) => {
        const item = node.getItem();
        if (isFuncCallNode(item)) {
          if (item.instancesWrapper.getInstance()?.getFuncCall())
            return [...acc, of(null)];
          else if (mockMode) {
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
      }, [] as Array<Observable<Readonly<[FuncCallNode, IFuncCallAdapter, boolean]>> | Observable<null>>);
      return merge(...pendingFuncCallLoads, MAX_CONCURENT_SAVES).pipe(
        filter((x) => !!x),
        map(([node, nAdapter, isNew]) => {
          node.setAdapter(nAdapter, isNew);
          node.instancesWrapper.isLoading$.next(false);
        }),
        toArray(),
        mapTo(tree),
      );
    });
  }

  private static findPipelineNode(tree: StateTree, uuid?: string) {
    const data = tree.find((item) => uuid == null || item.uuid === uuid);
    if (data == null)
      throw new Error(`Node uuid ${uuid} not found`);
    const [root, path] = data;
    const item = root.getItem();
    if (isFuncCallNode(item))
      throw new Error(`FuncCall node ${uuid}, but pipeline is expected`);

    const nqName = root.getItem().config.nqName;
    return [root, nqName, path] as const;
  }

  private updateWithLock<R>(fn: () => Observable<R>) {
    return defer(() => {
      this.globalStructureLock();
      return fn();
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
      this.makeStateRequests.next(true);
    }));
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

  private getSubConfig(path: NodePath, id: string) {
    const refMap = buildRefMap(this.config);
    const subConfPath = [...path.map((s) => s.id), id];
    const subConfig = getConfigByInstancePath(subConfPath, this.config, refMap);
    return subConfig;
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

  private makeMetaCall() {
    return defer(() => {
      if (this.mockMode || !this.config.nqName)
        return undefined;
      return makeMetaCall(this.config.nqName);
    });
  }

  private saveMetaCall(root: TreeNode<StateTreeNode>, nqName: string, currentMetaCall?: DG.FuncCall) {
    return defer(() => {
      if (this.mockMode || !nqName)
        return undefined;
      const state = StateTree.toStateRec(root, true, {disableNodesUUID: true});
      const json = JSON.stringify(state);
      return saveInstanceState(nqName, json, currentMetaCall);
    });
  }

  public static toStateRec(node: TreeNode<StateTreeNode>, isSerialized: boolean, options: StateTreeSerializationOptions = {}): PipelineState {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return isSerialized ? item.toSerializedState(options) : item.toState(options);

    const state = item.toState(options);
    const steps = node.getChildren().map((node) => {
      const item = this.toStateRec(node.item, isSerialized, options);
      return item;
    });
    return {...state, steps} as PipelineState;
  }
}

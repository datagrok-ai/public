import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, defer, of, merge, Subject, BehaviorSubject, from} from 'rxjs';
import {finalize, map, mapTo, toArray, concatMap, tap} from 'rxjs/operators';
import {NodePath, BaseTree, TreeNode} from '../data/BaseTree';
import {PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {isFuncCallSerializedState, PipelineInstanceConfig, PipelineSerializedState, PipelineState} from '../config/PipelineInstance';
import {buildTraverseD} from '../data/graph-traverse-utils';
import {buildRefMap, ConfigTraverseItem, getConfigByInstancePath, isPipelineParallelConfig, isPipelineSelfRef, isPipelineSequentialConfig, isPipelineStaticConfig, isPipelineStepConfig, PipelineStepConfigurationProcessed} from '../config/config-utils';
import {FuncCallAdapter, FuncCallMockAdapter} from './FuncCallAdapters';
import {loadFuncCall, loadInstanceState, makeFuncCall, makeMetaCall, saveFuncCall, saveInstanceState} from './funccall-utils';
import {ConsistencyInfo, FuncCallNode, FuncCallStateInfo, isFuncCallNode, ParallelPipelineNode, PipelineNodeBase, SequentialPipelineNode, StateTreeNode, StateTreeSerializationOptions, StaticPipelineNode} from './StateTreeNodes';
import {indexFromEnd} from '../utils';
import {ValidationResultBase} from '../../../shared-utils/validation';
import {LinksState} from './LinksState';

const MAX_CONCURENT_SAVES = 5;

export class StateTree extends BaseTree<StateTreeNode> {
  public makeStateRequests$ = new Subject<true>();
  public metaCall$ = new BehaviorSubject<DG.FuncCall | undefined>(undefined);
  public isLocked$ = new BehaviorSubject(false);
  public linksState = new LinksState();

  constructor(
    item: StateTreeNode,
    private config: PipelineConfigurationProcessed,
    private mockMode = false,
  ) {
    super(item);
  }

  public toSerializedState(options: StateTreeSerializationOptions = {}): PipelineSerializedState {
    return StateTree.toStateRec(this.getRoot(), true, options) as PipelineSerializedState;
  }

  public toState(options: StateTreeSerializationOptions = {}): PipelineState {
    return StateTree.toStateRec(this.getRoot(), false, options) as PipelineState;
  }

  public getValidations() {
    const entries = this.traverse(this.getRoot(), (acc, node) => {
      const item = node.getItem();
      if (isFuncCallNode(item))
        return [...acc, [item.uuid, item.validationInfo$] as const];

      return acc;
    }, [] as (readonly [string, BehaviorSubject<Record<string, ValidationResultBase>>])[]);
    return Object.fromEntries(entries);
  }

  public getConsistency() {
    const entries = this.traverse(this.getRoot(), (acc, node) => {
      const item = node.getItem();
      if (isFuncCallNode(item))
        return [...acc, [item.uuid, item.consistencyInfo$] as const];

      return acc;
    }, [] as (readonly [string, BehaviorSubject<Record<string, ConsistencyInfo>>])[]);
    return Object.fromEntries(entries);
  }

  public getMeta() {
    const entries = this.traverse(this.getRoot(), (acc, node) => {
      const item = node.getItem();
      if (isFuncCallNode(item))
        return [...acc, [item.uuid, item.metaInfo$] as const];
      return acc;
    }, [] as (readonly [string, BehaviorSubject<Record<string, BehaviorSubject<any | undefined>>>])[]);
    return Object.fromEntries(entries);
  }

  public getFuncCallStates() {
    const entries = this.traverse(this.getRoot(), (acc, node) => {
      const item = node.getItem();
      if (isFuncCallNode(item))
        return [...acc, [item.uuid, item.funcCallState$] as const];
      return acc;
    }, [] as (readonly [string, BehaviorSubject<FuncCallStateInfo | undefined>])[]);
    return Object.fromEntries(entries);
  }

  public save(uuid?: string) {
    return this.updateWithTreeLock(() => {
      const [root, nqName] = StateTree.findPipelineNode(this, uuid);
      if (this.mockMode)
        throw new Error(`Cannot save mock tree`);
      if (nqName == null)
        throw new Error(`Attempting to save pipeline with no nqName`);
      const rootItem = root.getItem() as PipelineNodeBase;
      if (rootItem.config.provider == null)
        throw new Error(`Attempting to save pipeline with no nqName provider`);
      const isRoot = root === this.root;
      const pendingFuncCallSaves = this.traverse(root, (acc, node) => {
        const item = node.getItem();
        if (isFuncCallNode(item)) {
          const obs$ = defer(() => saveFuncCall(item.instancesWrapper)).pipe(
            map((nextCall) => {
              const nadapter = new FuncCallAdapter(nextCall, item.isReadonly);
              item.changeAdapter(nadapter);
            }),
          );
          return [...acc, obs$];
        }
        return acc;
      }, [] as Array<Observable<void>>);
      return merge(...pendingFuncCallSaves, MAX_CONCURENT_SAVES).pipe(
        toArray(),
        concatMap(() => this.saveMetaCall(root, nqName, isRoot ? this.metaCall$.value : undefined)),
        tap((call) => isRoot ? this.metaCall$.next(call) : void(0)),
      );
    });
  }

  public initAll() {
    return this.updateWithTreeLock(() => {
      return merge(StateTree.loadOrCreateCalls(this, this.mockMode), this.initMetaCall()).pipe(toArray(), mapTo(this));
    });
  }

  public initFuncCalls() {
    return this.updateWithTreeLock(() => {
      return StateTree.loadOrCreateCalls(this, this.mockMode);
    });
  }

  public moveSubtree(uuid: string, pos: number) {
    return this.updateWithTreeLock(() => {
      const data = this.find((item) => item.uuid === uuid);
      if (data == null)
        throw new Error(`Node uuid ${uuid} not found`);
      const [node, path] = data;
      this.removeBrunch(path);
      const ppath = path.slice(0, -1);
      this.attachBrunch(ppath, node, node.getItem().config.id, pos);
      return of(ppath);
    });
  }

  public removeSubtree(uuid: string) {
    return this.updateWithTreeLock(() => {
      const data = this.find((item) => item.uuid === uuid);
      if (data == null)
        throw new Error(`Node uuid ${uuid} not found`);
      const [, path] = data;
      this.removeBrunch(path);
      return of(path);
    });
  }

  public addSubTree(puuid: string, id: string, pos: number, initCalls = true) {
    return this.updateWithTreeLock(() => {
      const [_root, _nqName, path] = StateTree.findPipelineNode(this, puuid);
      const subConfig = this.getSubConfig(path, id);
      if (isPipelineStepConfig(subConfig))
        this.attachBrunch(path, new TreeNode(new FuncCallNode(subConfig, false)), id, pos);
      else
        StateTree.fromPipelineConfig({config: this.config, startNode: subConfig, startPath: [...path, {id, idx: pos}], startState: this, isReadonly: false, mockMode: this.mockMode});
      return (initCalls ? StateTree.loadOrCreateCalls(this, this.mockMode) : of(this)).pipe(mapTo(path));
    });
  }

  public loadSubTree(puuid: string, dbId: string, id: string, pos: number, isReadonly: boolean) {
    return this.updateWithTreeLock(() => {
      const [_root, _nqName, path] = StateTree.findPipelineNode(this, puuid);
      const subConfig = this.getSubConfig(path, id);
      if (isPipelineStepConfig(subConfig))
        throw new Error(`FuncCall node ${JSON.stringify(path)}, but pipeline is expected`);
      const tree = StateTree.load(dbId, subConfig, {mockMode: this.mockMode, isReadonly}).pipe(
        concatMap((tree) => StateTree.loadOrCreateCalls(tree, this.mockMode)),
      ).pipe(
        map((tree) => this.attachBrunch(path, tree.root, id, pos)),
      );
      return tree.pipe(mapTo(path));
    });
  }

  public runStep(uuid: string, mockResults?: Record<string, any>, mockDelay?: number) {
    return this.updateWithTreeLock(() => {
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
    });
  }

  public close() {
    this.traverse(this.root, (acc, node) => {
      const item = node.getItem();
      if (isFuncCallNode(item))
        item.close();
      return acc;
    }, undefined);
  }

  static load(dbId: string, config: PipelineConfigurationProcessed, {isReadonly = false, mockMode = false}: {isReadonly?: boolean, mockMode?: boolean } = {}): Observable<StateTree> {
    return defer(async () => {
      const [metaCall, state] = await loadInstanceState(dbId);
      const tree = StateTree.fromInstanceState({state, config, metaCall, isReadonly, mockMode});
      return tree;
    });
  }

  static fromPipelineConfig({config, startNode = config, startPath = [], startState, isReadonly = false, mockMode = false}: { config: PipelineConfigurationProcessed; startNode?: PipelineConfigurationProcessed; startPath?: Readonly<NodePath>; startState?: StateTree; isReadonly?: boolean, mockMode?: boolean; }): StateTree {
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
      const [node, ppath, idx] = StateTree.makeTreeNode(config, refMap, path, isReadonly);
      if (isFuncCallNode(node)) {
        if (!isPipelineStepConfig(state))
          throw new Error(`Wrong FuncCall node state type ${state.type} on path ${JSON.stringify(path)}`);
        node.initState(state);
      }
      return StateTree.addTreeNodeOrCreate(acc, config, node, ppath, idx, mockMode);
    }, startState);
    return tree!;
  }

  static fromInstanceState({state, config, metaCall, isReadonly, mockMode = false}: { state: PipelineSerializedState; config: PipelineConfigurationProcessed; metaCall: DG.FuncCall; isReadonly: boolean, mockMode?: boolean }): StateTree {
    const refMap = buildRefMap(config);

    const traverse = buildTraverseD([] as Readonly<NodePath>, (item: PipelineSerializedState, path) => {
      if (isFuncCallSerializedState(item))
        return [];
      else
        return item.steps.map((step, idx) => [step, [...path, {id: step.configId, idx}]] as const);
    });

    const tree = traverse(state, (acc, state, path) => {
      const [node, ppath, idx] = StateTree.makeTreeNode(config, refMap, path, isReadonly || state.isReadonly);
      node.restoreState(state);
      return StateTree.addTreeNodeOrCreate(acc, config, node, ppath, idx, mockMode);
    }, undefined as StateTree | undefined);
    tree!.metaCall$.next(metaCall);
    return tree!;
  }

  static fromInstanceConfig(instanceConfig: PipelineInstanceConfig, config: PipelineConfigurationProcessed, {isReadonly = false, mockMode = false}: {isReadonly?: boolean, mockMode?: boolean } = {}): StateTree {
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
      const [node, ppath, idx] = StateTree.makeTreeNode(config, refMap, path, isReadonly);
      if (isFuncCallNode(node))
        node.initState(state);
      return StateTree.addTreeNodeOrCreate(acc, config, node, ppath, idx, mockMode);
    }, undefined as StateTree | undefined);
    return tree!;
  }

  private static makeTreeNode(config: PipelineConfigurationProcessed, refMap: Map<string, PipelineConfigurationProcessed>, path: Readonly<NodePath>, isReadonly: boolean) {
    const confPath = path.map((p) => p.id);
    const nodeConf = getConfigByInstancePath(confPath, config, refMap);
    const node = StateTree.makeNode(nodeConf, isReadonly);
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
          if (item.instancesWrapper.getInstance())
            return acc;
          else if (mockMode) {
            const obs$ = defer(() => {
              const adapter = new FuncCallMockAdapter(item.config.io!, item.isReadonly);
              item.initAdapter(adapter, {}, true, true);
              return of(undefined);
            });
            return [...acc, obs$];
          } else if (item.instancesWrapper.id == null && item.pendingId == null) {
            const obs$ = defer(() => {
              return from(makeFuncCall(item.config.nqName, false)).pipe(
                map(([adapter, restrictions, outputState]) => {
                  item.initAdapter(adapter, restrictions, outputState, true);
                }),
              );
            });
            return [...acc, obs$];
          } else if (item.instancesWrapper.id == null && item.pendingId) {
            const savedId = item.pendingId;
            const obs$ = defer(() => {
              return from(loadFuncCall(savedId, item.isReadonly)).pipe(
                map(([adapter, restrictions, outputState]) => {
                  item.initAdapter(adapter, restrictions, outputState, false);
                }),
              );
            });
            return [...acc, obs$];
          }
        }
        return acc;
      }, [] as Array<Observable<undefined | void>>);
      return merge(...pendingFuncCallLoads, MAX_CONCURENT_SAVES).pipe(
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

  private updateWithTreeLock<R>(fn: () => Observable<R>) {
    return defer(() => {
      this.treeLock();
      return fn();
    }).pipe(finalize(() => {
      this.treeUnlock();
      this.makeStateRequests$.next(true);
    }));
  }

  private treeLock() {
    if (this.isLocked$.value)
      throw new Error(`Global structure double lock`);
    this.isLocked$.next(true);
  }

  private treeUnlock() {
    if (!this.isLocked$.value)
      throw new Error(`Global structure double unlock`);
    this.isLocked$.next(false);
  }

  private getSubConfig(path: NodePath, id: string) {
    const refMap = buildRefMap(this.config);
    const subConfPath = [...path.map((s) => s.id), id];
    const subConfig = getConfigByInstancePath(subConfPath, this.config, refMap);
    return subConfig;
  }

  private static makeNode(nodeConf: PipelineConfigurationProcessed | PipelineStepConfigurationProcessed, isReadonly: boolean) {
    if (isPipelineStepConfig(nodeConf))
      return new FuncCallNode(nodeConf, isReadonly);
    else if (isPipelineStaticConfig(nodeConf))
      return new StaticPipelineNode(nodeConf, isReadonly);
    else if (isPipelineParallelConfig(nodeConf))
      return new ParallelPipelineNode(nodeConf, isReadonly);
    else if (isPipelineSequentialConfig(nodeConf))
      return new SequentialPipelineNode(nodeConf, isReadonly);

    throw new Error(`Wrong node type ${nodeConf}`);
  }

  private initMetaCall() {
    return this.makeMetaCall().pipe(
      tap((call) => this.metaCall$.next(call)),
    );
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
      if (this.mockMode)
        return of(currentMetaCall?.clone());
      const state = StateTree.toStateRec(root, true, {disableNodesUUID: true});
      return saveInstanceState(nqName, state);
    });
  }

  public static toStateRec(node: TreeNode<StateTreeNode>, isSerialized: boolean, options: StateTreeSerializationOptions = {}): PipelineState | PipelineSerializedState {
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

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, Observable, defer, identity, of, from, merge, Subject} from 'rxjs';
import {v4 as uuidv4} from 'uuid';
import {delay, finalize, map, mapTo, toArray, concatMap} from 'rxjs/operators';
import {NodePath, BaseTree, TreeNode} from '../data/BaseTree';
import {PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed} from '../config/config-processing-utils';
import {isFuncCallState, PipelineInstanceConfig, PipelineSerializedState, PipelineState, PipelineStateParallel, PipelineStateSequential, PipelineStateStatic, StepFunCallInitialConfig, StepFunCallSerializedState, StepFunCallState} from '../config/PipelineInstance';
import {buildTraverseD} from '../data/traversable';
import {buildRefMap, ConfigTraverseItem, getConfigByInstancePath, isPipelineParallelConfig, isPipelineSelfRef, isPipelineSequentialConfig, isPipelineStaticConfig, isPipelineStepConfig, PipelineStepConfigurationProcessed} from '../config/config-utils';
import {FuncCallAdapter, IFuncCallAdapter, isMockAdapter, IStateStore, MemoryStore} from './FuncCallAdapters';
import {FuncCallInstancesBridge} from './FuncCallInstancesBridge';
import {historyUtils} from '../../../history-utils';

export interface IStoreProvider {
  getStateStore(): IStateStore;
}

export class FuncCallNode implements IStoreProvider {
  public readonly uuid = uuidv4();
  public readonly nodeType = 'funccall';

  public instancesWrapper = new FuncCallInstancesBridge();

  constructor(
    public readonly config: PipelineStepConfigurationProcessed,
  ) {}

  setAdapter(fc?: IFuncCallAdapter, initValues = false) {
    this.instancesWrapper.setInstance(fc, initValues);
  }

  getStateStore() {
    return this.instancesWrapper;
  }

  restoreState(state: StepFunCallState) {
    this.instancesWrapper.isOutputOutdated$.next(!!state.isOuputOutdated);
    this.instancesWrapper.isCurrent$.next(!!state.isCurrent);
  }

  initState(initialConfig: StepFunCallInitialConfig) {
    this.instancesWrapper.inputRestrictions$.next(initialConfig.inputRestrictions ?? {});
    this.instancesWrapper.initialValues = initialConfig.values ?? {};
  }

  toState(): StepFunCallState {
    const res: StepFunCallState = {
      type: 'funccall',
      uuid: this.uuid,
      nqName: this.config.nqName,
      configId: this.config.id,
      friendlyName: this.config.friendlyName,
      funcCallId: this.instancesWrapper.getInstance()?.id,
      funcCall: this.instancesWrapper.getInstance()?.getFuncCall(),
      isRunning: false,
      isRunable: this.instancesWrapper.isRunable$.value,
      isOuputOutdated: this.instancesWrapper.isOutputOutdated$.value,
      isCurrent: this.instancesWrapper.isCurrent$.value,
    };
    return res;
  }

  toSerializedState(): StepFunCallSerializedState {
    const res: StepFunCallSerializedState = {
      type: 'funccall',
      uuid: this.uuid,
      nqName: this.config.nqName,
      configId: this.config.id,
      friendlyName: this.config.friendlyName,
      funcCallId: this.instancesWrapper.id,
      isOuputOutdated: this.instancesWrapper.isOutputOutdated$.value,
      isCurrent: this.instancesWrapper.isCurrent$.value,
    };
    return res;
  }
}

export class PipelineNodeBase implements IStoreProvider {
  public readonly uuid = uuidv4();
  private store: MemoryStore;

  public isUpdating$ = new BehaviorSubject(true);
  public isLoadable$ = new BehaviorSubject(false);
  public isSavable$ = new BehaviorSubject(false);
  public isReadonly$ = new BehaviorSubject(false);
  // TODO: implement (?)
  public metaCallId$ = new BehaviorSubject<string | undefined>(undefined);

  constructor(
    public readonly config: PipelineConfigurationProcessed,
  ) {
    this.store = new MemoryStore(config.states ?? []);
  }

  // TODO: count running links as well as loading
  setLoading(val: boolean) {
    this.isUpdating$.next(val);
  }

  getStateStore() {
    return this.store;
  }

  toState() {
    return {
      configId: this.config.id,
      uuid: this.uuid,
    };
  }
}

export class StaticPipelineNode extends PipelineNodeBase {
  public readonly nodeType = 'static';

  constructor(
    public readonly config: PipelineConfigurationStaticProcessed,
  ) {
    super(config);
  }

  toState(): PipelineStateStatic {
    const base = super.toState();
    const res: PipelineStateStatic = {
      ...base,
      nqName: this.config.nqName,
      type: this.nodeType,
      steps: [],
    };
    return res;
  }
}

export class ParallelPipelineNode extends PipelineNodeBase {
  public readonly nodeType = 'parallel';

  constructor(
    public readonly config: PipelineConfigurationParallelProcessed,
  ) {
    super(config);
  }

  toState(): PipelineStateParallel {
    const base = super.toState();
    const res: PipelineStateParallel = {
      ...base,
      type: this.nodeType,
      steps: [],
      stepTypes: this.config.stepTypes.map((s) => {
        const {id: configId, allowAdding} = s;
        return {configId, allowAdding};
      }),
    };
    return res;
  }
}

export class SequentialPipelineNode extends PipelineNodeBase {
  public readonly nodeType = 'sequential';

  constructor(
    public readonly config: PipelineConfigurationSequentialProcessed,
  ) {
    super(config);
  }

  toState() {
    const base = super.toState();
    const res: PipelineStateSequential = {
      ...base,
      type: this.nodeType,
      steps: [],
      stepTypes: this.config.stepTypes.map((s) => {
        const {id: configId, allowAdding} = s;
        return {configId, allowAdding};
      }),
    };
    return res;
  }
}

export type StateTreeNode = FuncCallNode | StaticPipelineNode | ParallelPipelineNode | SequentialPipelineNode;

export function isFuncCallNode(node: StateTreeNode): node is FuncCallNode {
  return node.nodeType === 'funccall';
}

export function isStaticPipelineNode(node: StateTreeNode): node is StaticPipelineNode {
  return node.nodeType === 'static';
}

export function isParallelPipelineNode(node: StateTreeNode): node is ParallelPipelineNode {
  return node.nodeType === 'parallel';
}

export function isSequentialPipelineNode(node: StateTreeNode): node is SequentialPipelineNode {
  return node.nodeType === 'sequential';
}

const RESTRICTIONS_PATH = 'INPUT_RESTRICTIONS';
const CONFIG_PATH = 'PIPELINE_CONFIG';
const PROVIDER_PATH = 'PIPELINE_PROVIDER';

export class StateTree extends BaseTree<StateTreeNode> {
  public makeStateRequests = new Subject<true>();

  constructor(
    item: StateTreeNode,
    private mockMode = false,
  ) {
    super(item);
  }

  public toSerializedState(): PipelineSerializedState {
    return this.toStateRec(this.getRoot(), true);
  }

  public toState(): PipelineState {
    return this.toStateRec(this.getRoot(), false);
  }

  public save(uuid: string, mockDelay?: number): Observable<string> {
    let item: PipelineNodeBase;
    return defer(() => {
      this.globalStructureLock();
      const [root, nqName] = this.findPipelineRootNode(uuid);
      const pendingFuncCallSaves = this.traverse(root, (acc, node) => {
        const item = node.getItem();
        if (isFuncCallNode(item)) {
          if (item.instancesWrapper.isRunning$.value || item.instancesWrapper.isLoading$.value)
            throw new Error(`FuncCall node ${item.uuid} saving during being updated`);
          const bridge = item.instancesWrapper.getInstance();
          if (bridge == null)
            throw new Error(`FuncCall node ${item.uuid} has no instance`);
          item.instancesWrapper.isLoading$.next(true);
          if (isMockAdapter(bridge)) {
            const obs$ = of(bridge).pipe(mockDelay ? delay(mockDelay) : identity).pipe(
              map((adapter) => [item, adapter] as const),
            );
            return [...acc, obs$];
          } else {
            const obs$ = defer(() => this.saveFuncCall(bridge)).pipe(map((nextCall) => [item, new FuncCallAdapter(nextCall)] as const));
            return [...acc, obs$];
          }
        } else {
          if (item.isUpdating$.value)
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
          const state = this.toStateRec(root, true);
          const json = JSON.stringify(state);
          return defer(() => StateTree.saveMetaCall(nqName, json));
        }),
      );
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
      if (item)
        item.isUpdating$.next(false);
      this.emitNewState();
    }));
  }

  private findPipelineRootNode(uuid: string) {
    const data = this.find((item) => item.uuid === uuid);
    if (data == null)
      throw new Error(`Node uuid ${uuid} not found`);
    const [root, path] = data;
    const item = root.getItem();
    if (isFuncCallNode(item))
      throw new Error(`Attempting to save FuncCall node ${uuid}`);

    const nqName = root.getItem().config.nqName;
    if (nqName == null)
      throw new Error(`Attempting to save pipeline with no nqName`);
    return [root, nqName, path] as const;
  }

  private async saveFuncCall(bridge: IFuncCallAdapter) {
    const fc = bridge.getFuncCall();
    // TODO: DF handling (?)
    const restrictions = bridge.inputRestrictions$.value;
    fc.options[RESTRICTIONS_PATH] = JSON.stringify(restrictions);
    return historyUtils.saveRun(fc);
  }

  private static async loadFuncCall(id: string): Promise<IFuncCallAdapter> {
    const fc = await historyUtils.loadRun(id, false);
    const restrictions = JSON.parse(fc.options[RESTRICTIONS_PATH] ?? {});
    const bridge = new FuncCallAdapter(fc);
    bridge.inputRestrictions$.next(restrictions);
    return bridge;
  }

  private static async saveMetaCall(nqName: string, stateJson: string) {
    const metaCall = await this.makeFuncCall(nqName);
    metaCall.options[CONFIG_PATH] = stateJson;
    const fc = await historyUtils.saveRun(metaCall);
    return fc.id;
  }

  private static async loadMetaConfig(id: string) {
    const fc = await historyUtils.loadRun(id, false);
    const config = JSON.parse(fc.options[CONFIG_PATH] ?? {});
    return config as PipelineSerializedState;
  }

  private static async makeFuncCall(nqName: string) {
    const func: DG.Func = await grok.functions.eval(nqName);
    return func.prepare({});
  }

  public load(id: string, uuid: string, config: PipelineConfigurationProcessed, pos?: number) {
    return defer(() => {
      this.globalStructureLock();
      const [root, _nqName, path] = this.findPipelineRootNode(uuid);
      const tree = StateTree.load(id, config).pipe(
        concatMap((tree) => this.loadOrCreateCalls(tree)),
      );
      const branch = root.getItem();
      this.attachBrunch(path, root, branch.config.id, pos);
      return tree;
    }).pipe(finalize(() => {
      this.globalStructureUnlock();
    }));
  }

  private loadOrCreateCalls(tree: StateTree) {
    return defer(() => {
      const pendingFuncCallLoads = tree.traverse(tree.root, (acc, node) => {
        const item = node.getItem();
        if (isFuncCallNode(item)) {
          const obs$ = defer(() => StateTree.loadFuncCall(item.instancesWrapper.id!)).pipe(
            map((nitem) => [item, nitem] as const),
          );
          return [...acc, obs$];
        }
        return acc;
      }, [] as Array<Observable<Readonly<[FuncCallNode, IFuncCallAdapter]>>>);
      return merge(...pendingFuncCallLoads, 5).pipe(
        map(([node, nAdapter]) => {
          node.setAdapter(nAdapter);
          node.instancesWrapper.isLoading$.next(false);
        }),
        toArray(),
        mapTo(tree),
      );
    });
  }

  static load(id: string, config: PipelineConfigurationProcessed): Observable<StateTree> {
    return defer(async () => {
      const state = await StateTree.loadMetaConfig(id);
      const tree = StateTree.fromState(state, config);
      return tree;
    });
  }

  static fromState(state: PipelineSerializedState, config: PipelineConfigurationProcessed): StateTree {
    const refMap = buildRefMap(config);

    const traverse = buildTraverseD([] as Readonly<NodePath>, (item: PipelineState, path) => {
      if (isFuncCallState(item))
        return [];
      else
        return item.steps.map((step, idx) => [step, [...path, {id: step.configId, idx}]] as const);
    });

    const tree = traverse(state, (acc, state, path) => {
      const confPath = path.map((p) => p.id);
      const treePath = path.map((p) => ({idx: p.idx}));
      const nodeConf = getConfigByInstancePath(confPath, config, refMap);
      if (nodeConf == null)
        throw new Error(`Unable to find find configuration node on path ${JSON.stringify(path)}`);
      const node = StateTree.makeNode(nodeConf);
      if (isFuncCallNode(node)) {
        if (!isFuncCallState(state))
          throw new Error(`Wrong FuncCall state type ${state.type} on path ${JSON.stringify(path)}`);
        node.restoreState(state);
      }
      if (acc)
        acc.addItem(treePath, node, node.config.id);
      else
        return new StateTree(node);
      return acc;
    }, undefined as StateTree | undefined);
    return tree!;
  }

  static fromInstanceConfig(instanceConfig: PipelineInstanceConfig, config: PipelineConfigurationProcessed): StateTree {
    const refMap = buildRefMap(config);

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
      const confPath = path.map((p) => p.id);
      const treePath = path.map((p) => ({idx: p.idx}));
      const nodeConf = getConfigByInstancePath(confPath, config, refMap);
      if (nodeConf == null)
        throw new Error(`Unable to find configuration node on path ${JSON.stringify(path)}`);
      const node = StateTree.makeNode(nodeConf);
      if (isFuncCallNode(node))
        node.initState(state);
      if (acc)
        acc.addItem(treePath, node, node.config.id);
      else
        return new StateTree(node);
      return acc;
    }, undefined as StateTree | undefined);
    return tree!;
  }

  static makeInitialTree(config: PipelineConfigurationProcessed): StateTree {
    const refMap = buildRefMap(config);

    const traverse = buildTraverseD([] as Readonly<NodePath>, (data: ConfigTraverseItem, path, visited) => {
      if (visited!.has(data))
        throw new Error(`Initial config cycle on node ${data} on path ${JSON.stringify(path)}`);
      visited!.add(data);
      if (isPipelineParallelConfig(data) || isPipelineSequentialConfig(data)) {
        const items = (data?.initialSteps ?? []).map((step, idx) => {
          const item = data.stepTypes.find((t) => t.id === step.id);
          if (!item)
            throw new Error(`Node ${step.id} not found on path ${JSON.stringify(path)}`);
          const nextPath = [...path, {id: step.id, idx}];
          return [item, nextPath, visited] as const;
        });
        return items;
      } else if (isPipelineStaticConfig(data))
        return data.steps.map((step, idx) => [step, [...path, {id: step.id, idx}], visited] as const);
      else if (isPipelineSelfRef(data)) {
        const next = refMap.get(data.selfRef);
        if (!next)
          throw new Error(`Failed to deref ${data.selfRef} on path ${JSON.stringify(path)}`);

        return [[next, [...path, {id: next.id, idx: 0}], visited!]] as const;
      }
      return [] as const;
    }, new Set<ConfigTraverseItem>());

    const tree = traverse(config, (acc, state, path) => {
      const confPath = path.map((p) => p.id);
      const treePath = path.map((p) => ({idx: p.idx}));
      const nodeConf = getConfigByInstancePath(confPath, config, refMap);
      if (nodeConf == null)
        throw new Error(`Unable to find configuration node on path ${JSON.stringify(path)}`);
      const node = StateTree.makeNode(nodeConf);
      if (isFuncCallNode(node)) {
        if (!isPipelineStepConfig(state))
          throw new Error(`Wrong FuncCall node state type ${state.type} on path ${JSON.stringify(path)}`);
        node.initState(state);
      }
      if (acc)
        acc.addItem(treePath, node, node.config.id);
      else
        return new StateTree(node);
      return acc;
    }, undefined as StateTree | undefined);
    return tree!;
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

  private toStateRec(node: TreeNode<StateTreeNode>, isSerialized: boolean): PipelineState {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return isSerialized ? item.toSerializedState() : item.toState();

    const state = item.toState();
    const steps = node.getChildren().map((node) => {
      const item = this.toStateRec(node.item, isSerialized);
      return item;
    });
    return {...state, steps} as PipelineState;
  }
}

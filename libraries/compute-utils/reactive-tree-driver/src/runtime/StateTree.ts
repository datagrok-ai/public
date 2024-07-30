import {BehaviorSubject, Observable, of, combineLatest, from, defer} from 'rxjs';
import {v4 as uuidv4} from 'uuid';
import {switchMap, map} from 'rxjs/operators';
import {NodePath, NodeTree, TreeNode} from '../data/NodeTree';
import {PipelineStepConfiguration} from '../config/PipelineConfiguration';
import {FuncallStateItem, PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed } from '../config/config-processing-utils';
import {ValidationResultBase} from '../../../shared-utils/validation';
import { isFuncCallState, PipelineInstanceConfig, PipelineSerializedState, PipelineState, PipelineStateParallel, PipelineStateSequential, PipelineStateStatic, StepFunCallInitialConfig, StepFunCallSerializedState, StepFunCallState} from '../config/PipelineInstance';
import {buildTraverseD} from '../data/traversable';
import { buildRefMap, ConfigTraverseItem, getConfigByInstancePath, isPipelineParallelConfig, isPipelineSelfRef, isPipelineSequentialConfig, isPipelineStaticConfig, isPipelineStepConfig} from '../config/config-utils';
import {IFuncCallAdapter, IRunnableWrapper, IStateStore, IValidationStore, MemoryStore, RestrictionType} from './FuncCallAdapters';

export class FuncCallInstancesBridge implements IStateStore, IValidationStore, IRunnableWrapper {
  private instance$ = new BehaviorSubject<IFuncCallAdapter | undefined>(undefined);
  public isRunning$ = new BehaviorSubject(false);
  public isRunable$ = new BehaviorSubject(false);
  public isOutputOutdated$ = new BehaviorSubject(false);
  public isCurrent$ = new BehaviorSubject(false);
  public validations$ = new BehaviorSubject({});
  public inputRestrictions$ = new BehaviorSubject({});

  constructor() {
    this.instance$.pipe(
      switchMap((instance) => instance? instance.isRunning$ : of(false)),
    ).subscribe(this.isRunning$);

    // TODO: add global lock as well
    this.instance$.pipe(
      switchMap((instance) => {
        if (instance == null)
          return of(false);
        return combineLatest([
          instance.isRunning$,
          instance.validations$.pipe(map(validations => this.isRunnable(validations))),
        ]).pipe(map((isRunning, isValid) => !isRunning && isValid))
      }),
    ).subscribe(this.isRunable$);

    this.instance$.pipe(
      switchMap((instance) => instance? instance.isOutputOutdated$ : of(false)),
    ).subscribe(this.isOutputOutdated$);
  }

  get id() {
    return this.instance$.value?.id;
  }

  setInstance(instance: IFuncCallAdapter | undefined) {
    this.instance$.next(instance);
  }

  getInstance() {
    return this.instance$.value;
  }

  getState<T = any>(id: string): T | undefined {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return currentInstance.getState<T>(id);
  }

  getStateChanges<T = any>(id: string): Observable<T | undefined> {
    return this.instance$.pipe(
      switchMap((instance) => instance ? instance.getStateChanges(id) : of(undefined)),
    );
  }

  setState<T = any>(id: string, val: T | undefined, restrictionType: RestrictionType = 'none') {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      currentInstance.setState(id, val, restrictionType);
  }

  setValidation(id: string, validatorId: string, validation: ValidationResultBase | undefined) {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return currentInstance.setValidation(id, validatorId, validation);
  }

  run() {
    const currentInstance = this.instance$.value;
    if (currentInstance) {
      return from(defer(() => currentInstance.run()));
    } else {
      throw new Error(`Attempting to run an empty FuncCallInstancesBridge`);
    }
  }

  private isRunnable(validations: Record<string, Record<string, ValidationResultBase | undefined>>) {
    for (const validatorResults of Object.values(validations)) {
      for (const res of Object.values(validatorResults)) {
        if (res?.errors?.length)
          return false;
      }
    }
    return true;
  }
}

export interface IStoreProvider {
  getStore(): IStateStore;
}

export class FuncCallNode implements IStoreProvider {
  public readonly uuid = uuidv4();
  public readonly nodeType = 'funccall';
  public initialValues: Record<string, any> = {};
  public initialInputRestriction: Record<string, RestrictionType> = {};

  private instancesWrapper = new FuncCallInstancesBridge();

  constructor(
    public readonly config: PipelineStepConfiguration<FuncallStateItem[]>,
  ) {}

  setFuncall(fc?: IFuncCallAdapter) {
    this.instancesWrapper.setInstance(fc);
  }

  getStore() {
    return this.instancesWrapper;
  }

  restoreState(state: StepFunCallState) {
    this.instancesWrapper.isOutputOutdated$.next(!!state.isOuputOutdated);
    this.instancesWrapper.isCurrent$.next(!!state.isCurrent);
  }

  initState(initialConfig: StepFunCallInitialConfig) {
    this.initialValues = initialConfig.values ?? {};
    this.initialInputRestriction = initialConfig.inputRestrictions ?? {};
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
      isRunning: false, //this.adapter.getInstance()?.isRunning$.value,
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

  constructor(
    public readonly config: PipelineConfigurationProcessed,
  ) {
    this.store = new MemoryStore(config.states ?? []);
  }

  getStore() {
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
      stepTypes: this.config.stepTypes.map(s => {
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
      stepTypes: this.config.stepTypes.map(s => {
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


export class StateTree extends NodeTree<StateTreeNode> {
  public toSerializedState(): PipelineState {
    return this.toStateRec(this.getRoot(), true);
  }

  public toState(): PipelineState {
    return this.toStateRec(this.getRoot(), false);
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
      const confPath = path.map(p => p.id);
      const treePath = path.map(p => ({idx: p.idx}));
      const nodeConf = getConfigByInstancePath(confPath, config, refMap);
      if (nodeConf == null)
        throw new Error(`Unable to find conf node for ${path}`);
      const node = StateTree.makeNode(nodeConf);
      if (isFuncCallNode(node)) {
        if (!isFuncCallState(state))
          throw new Error(`Wrong FuncCall node state: ${state}`);
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
        throw new Error(`Initial config cycle on node ${data}`);
      visited!.add(data);
      if (data.steps)
        return data.steps.map((step, idx) => [step, [...path, {id: step.id, idx}], visited] as const);
      else
        return [];
    }, new Set<PipelineInstanceConfig>());

    const tree = traverse(instanceConfig, (acc, state, path) => {
      const confPath = path.map(p => p.id);
      const treePath = path.map(p => ({idx: p.idx}));
      const nodeConf = getConfigByInstancePath(confPath, config, refMap);
      if (nodeConf == null)
        throw new Error(`Unable to find conf node for ${path}`);
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

  static stateMakeInitialTree(config: PipelineConfigurationProcessed): StateTree {
    const refMap = buildRefMap(config);

    const traverse = buildTraverseD([] as Readonly<NodePath>, (data: ConfigTraverseItem, path, visited) => {
      if (visited!.has(data))
        throw new Error(`Initial config cycle on node ${data}`);
      visited!.add(data);
      if (isPipelineParallelConfig(data) || isPipelineSequentialConfig(data)) {
        const items = (data?.initialSteps ?? []).map((step, idx) => {
          const item = data.stepTypes.find(t => t.id === step.id);
          if (!item)
            throw new Error(`Item ${step.id} not found on path ${path}`);
          const nextPath = [...path, {id: step.id, idx}];
          return [item, nextPath, visited] as const;
        });
        return items;
      } else if (isPipelineStaticConfig(data)) {
        return data.steps.map((step, idx) => [step, [...path, {id: step.id, idx}], visited] as const);
      } else if (isPipelineSelfRef(data)) {
        const next = refMap.get(data.selfRef);
        if (!next) {
          throw new Error(`Failed to deref ${data.selfRef} on path ${path}`);
        }
        return [[next, [...path, {id: next.id, idx: 0}], visited!]] as const
      }
      return [] as const;
    }, new Set<ConfigTraverseItem>());

    const tree = traverse(config, (acc, state, path) => {
      const confPath = path.map(p => p.id);
      const treePath = path.map(p => ({idx: p.idx}));
      const nodeConf = getConfigByInstancePath(confPath, config, refMap);
      if (nodeConf == null)
        throw new Error(`Unable to find conf node for ${path}`);
      const node = StateTree.makeNode(nodeConf);
      if (isFuncCallNode(node)) {
        if (!isPipelineStepConfig(state))
          throw new Error(`Wrong FuncCall node state: ${state}`);
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

  private static makeNode(nodeConf: PipelineConfigurationProcessed | PipelineStepConfiguration<FuncallStateItem[]>) {
    if (isPipelineStepConfig(nodeConf)) {
      return new FuncCallNode(nodeConf);
    } else if (isPipelineStaticConfig(nodeConf)) {
      return new StaticPipelineNode(nodeConf);
    } else if (isPipelineParallelConfig(nodeConf)) {
      return new ParallelPipelineNode(nodeConf);
    } else if (isPipelineSequentialConfig(nodeConf)) {
      return new SequentialPipelineNode(nodeConf);
    }
    throw new Error(`Wrong node type ${nodeConf}`);
  }

  private toStateRec(node: TreeNode<StateTreeNode>, isSerialized: boolean): PipelineState {
    const item = node.getItem();
    if (isFuncCallNode(item)) {
      return isSerialized ? item.toSerializedState() : item.toState();
    }
    const state = item.toState();
    const steps = node.getChildren().map((node) => {
      const item = this.toStateRec(node.item , isSerialized);
      return item;
    });
    return {...state, steps} as PipelineState;
  }
}

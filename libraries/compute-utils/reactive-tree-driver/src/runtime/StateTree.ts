import {BehaviorSubject, Observable, of} from 'rxjs';
import {InputState, ItemPathArray} from '../data/common-types';
import {v4 as uuidv4} from 'uuid';
import {switchMap} from 'rxjs/operators';
import {NodeTree, TreeNode} from '../data/NodeTree';
import {PipelineStepConfiguration, StateItem} from '../config/PipelineConfiguration';
import {FuncallStateItem, PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed } from '../config/config-processing-utils';
import {ValidationResultBase} from '../../../shared-utils/validation';
import {PipelineState, PipelineStateParallel, PipelineStateSequential, PipelineStateStatic, StepFunCallState} from '../config/PipelineInstance';

export interface ICallable {
  run(): void
  isRunning$: BehaviorSubject<boolean>;
  isRunable$: BehaviorSubject<boolean>;
  isOuputOutdated$: BehaviorSubject<boolean>;
  isCurrent$: BehaviorSubject<boolean>;
}

export interface IStateStore {
  getStateChanges<T = any>(id: string): Observable<T | undefined>;
  getState<T = any>(id: string): T | undefined;
  setState<T = any >(id: string, value: T | undefined, inputState?: InputState): void;
}

export interface IValidationDisplay {
  setValidation(id: string, validation: ValidationResultBase | undefined): void;
}

export interface IIdWrapper {
  id?: string;
}

export class StoreItem<T = any> {
  state$ = new BehaviorSubject<T | undefined>(undefined);
}

export class MemoryStore implements IStateStore {
  public readonly uuid = uuidv4();
  states: Record<string, StoreItem> = {};

  constructor(private statesDescriptions: StateItem[]) {
    for (const description of this.statesDescriptions)
      this.states[description.id] = new StoreItem();
  }

  getState<T = any>(id: string): T {
    return this.states[id]?.state$?.value;
  }

  getStateChanges<T = any>(id: string): Observable<T | undefined> {
    return this.states[id]?.state$.asObservable();
  }

  setState<T = any>(id: string, val: T | undefined, _inputState: InputState) {
    this.states[id]?.state$.next(val);
  }
}

export type IFuncCallBridge = IStateStore & IValidationDisplay & ICallable & IIdWrapper;

export class FuncCallInstanceAdapter implements IFuncCallBridge {
  private instance$ = new BehaviorSubject<IFuncCallBridge | undefined>(undefined);
  public isRunning$ = new BehaviorSubject(false);
  public isRunable$ = new BehaviorSubject(false);
  public isOuputOutdated$ = new BehaviorSubject(false);
  public isCurrent$ = new BehaviorSubject(false);

  constructor() {
    this.instance$.pipe(
      switchMap((instance) => instance? instance.isRunning$ : of(false)),
    ).subscribe(this.isRunning$);
    this.instance$.pipe(
      switchMap((instance) => instance? instance.isRunable$ : of(false)),
    ).subscribe(this.isRunable$);
    this.instance$.pipe(
      switchMap((instance) => instance? instance.isOuputOutdated$ : of(false)),
    ).subscribe(this.isOuputOutdated$);
  }

  get id() {
    return this.instance$.value?.id;
  }

  setInstance(instance: IFuncCallBridge | undefined) {
    this.instance$.next(instance);
  }

  getInstance() {
    return this.instance$.value;
  }

  getState<T = any>(id: string): T | undefined {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return currentInstance.getState(id);
  }

  getStateChanges<T = any>(id: string): Observable<T | undefined> {
    return this.instance$.pipe(
      switchMap((instance) => instance ? instance.getStateChanges(id) : of(undefined)),
    );
  }

  setState<T = any>(id: string, val: T | undefined, inputState: InputState) {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      currentInstance.setState(id, val, inputState);
  }

  setValidation(id: string, validation: ValidationResultBase | undefined) {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return currentInstance.setValidation(id, validation);
  }

  run() {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      currentInstance.run();
  }
}

export interface IStoreProvider {
  getStore(): IStateStore;
}

export class FuncCallNode implements IStoreProvider {
  public readonly uuid = uuidv4();
  public readonly nodeType = 'funccall';
  private adapter = new FuncCallInstanceAdapter();

  constructor(
    public readonly config: PipelineStepConfiguration<FuncallStateItem>,
  ) {}

  setFuncall(fc?: IFuncCallBridge) {
    this.adapter.setInstance(fc);
  }

  getStore() {
    return this.adapter;
  }

  static fromState() {

  }

  toState(): StepFunCallState {
    const res: StepFunCallState = {
      type: 'funccall',
      uuid: this.uuid,
      nqName: this.config.nqName,
      configId: this.config.id,
      friendlyName: this.config.friendlyName,
      funcCallId: this.adapter.getInstance()?.id,
      isRunning: this.adapter.getInstance()?.isRunning$.value,
      isRunable: this.adapter.getInstance()?.isRunable$.value,
      isOuputOutdated: this.adapter.getInstance()?.isOuputOutdated$.value,
      isCurrent: this.adapter.getInstance()?.isCurrent$.value,
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

  static fromState(state: PipelineStateStatic) {

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

  static fromState(state: PipelineStateParallel) {

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

  static fromState(state: PipelineStateSequential) {

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
  public toState(): PipelineState {
    return this.toStateRec(this.getRoot());
  }

  static fromState(state: PipelineState): StateTree {
    // TODO
  }

  private toStateRec(node: TreeNode<StateTreeNode>): PipelineState {
    const item = node.getItem();
    if (isFuncCallNode(item))
      return item.toState();

    const state = item.toState();
    const steps = node.getChildren().map((node) => {
      const item = this.toStateRec(node.item);
      return item;
    });
    return {...state, steps} as PipelineState;
  }
}

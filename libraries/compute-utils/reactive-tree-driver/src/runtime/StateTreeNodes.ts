import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {v4 as uuidv4} from 'uuid';
import { isFuncCallState, PipelineInstanceConfig, PipelineState, PipelineStateParallel, PipelineStateSequential, PipelineStateStatic, StepFunCallInitialConfig, StepFunCallSerializedState, StepFunCallState} from '../config/PipelineInstance';
import {PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed} from '../config/config-processing-utils';
import {IFuncCallAdapter, IStateStore, MemoryStore} from './FuncCallAdapters';
import {FuncCallInstancesBridge} from './FuncCallInstancesBridge';
import { isPipelineStepConfig, PipelineStepConfigurationProcessed} from '../config/config-utils';


export interface IStoreProvider {
  getStateStore(): IStateStore;
}

export class FuncCallNode implements IStoreProvider {
  public uuid = uuidv4();
  public readonly nodeType = 'funccall';

  public instancesWrapper = new FuncCallInstancesBridge();
  public pendingId?: string;

  constructor(
    public readonly config: PipelineStepConfigurationProcessed,
  ) {}

  setAdapter(fc?: IFuncCallAdapter, initValues = false) {
    this.instancesWrapper.setInstance(fc, initValues);
  }

  getStateStore() {
    return this.instancesWrapper;
  }

  restoreState(state: PipelineState) {
    if (!isFuncCallState(state))
      throw new Error(`Wrong FuncCall node state ${JSON.stringify(state)}`);
    this.instancesWrapper.isOutputOutdated$.next(!!state.isOuputOutdated);
    this.instancesWrapper.isCurrent$.next(!!state.isCurrent);
    this.pendingId = state.funcCallId;
    this.uuid = state.uuid;
  }

  initState(initialConfig: StepFunCallInitialConfig) {
    this.instancesWrapper.inputRestrictions$.next(initialConfig.inputRestrictions ?? {});
    this.instancesWrapper.initialValues = initialConfig.values ?? {};
  }

  toState(disableUUID = false): StepFunCallState {
    const res: StepFunCallState = {
      type: 'funccall',
      uuid: this.uuid,
      nqName: this.config.nqName,
      configId: this.config.id,
      friendlyName: this.config.friendlyName,
      funcCallId: this.instancesWrapper.getInstance()?.getFuncCall()?.id,
      funcCall: this.instancesWrapper.getInstance()?.getFuncCall(),
      isRunning: false,
      isRunable: this.instancesWrapper.isRunable$.value,
      isOuputOutdated: this.instancesWrapper.isOutputOutdated$.value,
      isCurrent: this.instancesWrapper.isCurrent$.value,
    };
    if (disableUUID)
      res.uuid = '';
    return res;
  }

  toSerializedState(disableUUID = false): StepFunCallSerializedState {
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
    if (disableUUID)
      res.uuid = '';
    return res;
  }
}

export class PipelineNodeBase implements IStoreProvider {
  public uuid = uuidv4();
  private store: MemoryStore;

  public isUpdating$ = new BehaviorSubject(false);
  public isReadonly$ = new BehaviorSubject(false);

  constructor(
    public readonly config: PipelineConfigurationProcessed,
  ) {
    this.store = new MemoryStore(config.states ?? []);
  }

  setLoading(val: boolean) {
    this.isUpdating$.next(val);
  }

  restoreState(state: PipelineState) {
    if (isFuncCallState(state))
      throw new Error(`Wrong pipeline node state ${JSON.stringify(state)}`);
    this.uuid = state.uuid;
  }

  getStateStore() {
    return this.store;
  }

  toState(disableUUID = false) {
    const state = {
      configId: this.config.id,
      uuid: this.uuid,
    };
    if (disableUUID)
      state.uuid = '';
    return state;
  }
}

export class StaticPipelineNode extends PipelineNodeBase {
  public readonly nodeType = 'static';

  constructor(
    public readonly config: PipelineConfigurationStaticProcessed,
  ) {
    super(config);
  }

  toState(disableUUID = false): PipelineStateStatic {
    const base = super.toState(disableUUID);
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

  toState(disableUUID = false): PipelineStateParallel {
    const base = super.toState(disableUUID);
    const res: PipelineStateParallel = {
      ...base,
      type: this.nodeType,
      steps: [],
      stepTypes: this.config.stepTypes.map((s) => {
        const {id: configId, disableUIAdding} = s;
        return {configId, disableUIAdding};
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

  toState(disableUUID = false) {
    const base = super.toState(disableUUID);
    const res: PipelineStateSequential = {
      ...base,
      type: this.nodeType,
      steps: [],
      stepTypes: this.config.stepTypes.map((s) => {
        const {id: configId, disableUIAdding} = s;
        return {configId, disableUIAdding};
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

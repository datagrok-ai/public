import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {v4 as uuidv4} from 'uuid';
import {isFuncCallState, PipelineState, PipelineStateParallel, PipelineStateSequential, PipelineStateStatic, StepFunCallInitialConfig, StepFunCallSerializedState, StepFunCallState} from '../config/PipelineInstance';
import {PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed} from '../config/config-processing-utils';
import {IFuncCallAdapter, IStateStore, MemoryStore} from './FuncCallAdapters';
import {FuncCallInstancesBridge} from './FuncCallInstancesBridge';
import {isPipelineConfig, PipelineStepConfigurationProcessed} from '../config/config-utils';

export type StateTreeSerializationOptions = {
  disableNodesUUID?: boolean,
  disableCallsUUID?: boolean,
};

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
    if (state.uuid)
      this.uuid = state.uuid;
  }

  initState(initialConfig: StepFunCallInitialConfig) {
    this.instancesWrapper.inputRestrictions$.next(initialConfig.inputRestrictions ?? {});
    this.instancesWrapper.initialValues = initialConfig.values ?? {};
  }

  toState(options: StateTreeSerializationOptions): StepFunCallState {
    const res: StepFunCallState = {
      type: 'funccall',
      uuid: this.uuid,
      nqName: this.config.nqName,
      configId: this.config.id,
      friendlyName: this.config.friendlyName,
      funcCallId: this.instancesWrapper.getInstance()?.getFuncCall()?.id,
      funcCall: this.instancesWrapper.getInstance()?.getFuncCall(),
      isRunning: this.instancesWrapper.isRunning$.value,
      isRunable: this.instancesWrapper.isRunable$.value,
      isOuputOutdated: this.instancesWrapper.isOutputOutdated$.value,
      isCurrent: this.instancesWrapper.isCurrent$.value,
    };
    if (options.disableNodesUUID)
      res.uuid = '';
    if (options.disableCallsUUID)
      res.funcCallId = '';
    return res;
  }

  toSerializedState(options: StateTreeSerializationOptions): StepFunCallSerializedState {
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
    if (options.disableNodesUUID)
      res.uuid = '';
    if (options.disableCallsUUID)
      res.funcCallId = '';
    return res;
  }
}

export class PipelineNodeBase implements IStoreProvider {
  public uuid = uuidv4();
  private store: MemoryStore;

  public isReadonly$ = new BehaviorSubject(false);

  constructor(
    public readonly config: PipelineConfigurationProcessed,
  ) {
    this.store = new MemoryStore(config.states ?? []);
  }

  restoreState(state: PipelineState) {
    if (isFuncCallState(state))
      throw new Error(`Wrong pipeline node state ${JSON.stringify(state)}`);
    if (state.uuid)
      this.uuid = state.uuid;
  }

  getStateStore() {
    return this.store;
  }

  toState(options: StateTreeSerializationOptions) {
    const res = {
      configId: this.config.id,
      uuid: this.uuid,
      provider: this.config.provider,
      version: this.config.version,
      nqName: this.config.nqName,
    };
    if (options.disableNodesUUID)
      res.uuid = '';

    return res;
  }
}

export class StaticPipelineNode extends PipelineNodeBase {
  public readonly nodeType = 'static';

  constructor(
    public readonly config: PipelineConfigurationStaticProcessed,
  ) {
    super(config);
  }

  toState(options: StateTreeSerializationOptions): PipelineStateStatic {
    const base = super.toState(options);
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

  toState(options: StateTreeSerializationOptions): PipelineStateParallel {
    const base = super.toState(options);
    const res: PipelineStateParallel = {
      ...base,
      type: this.nodeType,
      steps: [],
      stepTypes: this.config.stepTypes.map((s) => {
        if (isPipelineConfig(s)) {
          const {id: configId, disableUIAdding, nqName, friendlyName} = s;
          return {configId, disableUIAdding, nqName, friendlyName};
        } else {
          const {id: configId, disableUIAdding} = s;
          return {configId, disableUIAdding};
        }
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

  toState(options: StateTreeSerializationOptions) {
    const base = super.toState(options);
    const res: PipelineStateSequential = {
      ...base,
      type: this.nodeType,
      steps: [],
      stepTypes: this.config.stepTypes.map((s) => {
        if (isPipelineConfig(s)) {
          const {id: configId, disableUIAdding, nqName, friendlyName} = s;
          return {configId, disableUIAdding, nqName, friendlyName};
        } else {
          const {id: configId, disableUIAdding} = s;
          return {configId, disableUIAdding};
        }
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

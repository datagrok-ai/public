import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, merge, Observable, of, Subject} from 'rxjs';
import {v4 as uuidv4} from 'uuid';
import {PipelineStateParallel, PipelineStateSequential, PipelineStateStatic, ConsistencyInfo, StepFunCallInitialConfig, StepFunCallSerializedState, StepFunCallState, PipelineSerializedState, isFuncCallSerializedState} from '../config/PipelineInstance';
import {PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed} from '../config/config-processing-utils';
import {IFuncCallAdapter, IStateStore, MemoryStore, RestrictionState} from './FuncCallAdapters';
import {FuncCallInstancesBridge} from './FuncCallInstancesBridge';
import {isPipelineConfig, PipelineStepConfigurationProcessed} from '../config/config-utils';
import {mergeValidationResultsBase, ValidationResultBase} from '../../../shared-utils/validation';
import {map, mapTo, scan, skip, switchMap, takeUntil, withLatestFrom} from 'rxjs/operators';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

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

  public instancesWrapper = new FuncCallInstancesBridge(this.config.io!, this.isReadonly);
  public pendingId?: string;
  private pendingState?: StepFunCallSerializedState;

  public consistencyInfo$ = new BehaviorSubject<Record<string, ConsistencyInfo>>({});
  public validationInfo$ = new BehaviorSubject<Record<string, ValidationResultBase>>({});

  private closed$ = new Subject<true>();

  constructor(
    public readonly config: PipelineStepConfigurationProcessed,
    public isReadonly: boolean,
  ) {
    this.instancesWrapper.instance$.pipe(
      switchMap((instance) => this.getConsistencyUpdater(instance)),
      takeUntil(this.closed$),
    ).subscribe(this.consistencyInfo$);

    this.instancesWrapper.instance$.pipe(
      switchMap((instance) => this.getValidationsUpdater(instance)),
      takeUntil(this.closed$),
    ).subscribe(this.validationInfo$);
  }

  setAdapter(adapter?: IFuncCallAdapter, initValues = false) {
    if (this.pendingState) {
      if (adapter) {
        adapter.inputRestrictions$.next(this.pendingState.inputRestrictions);
        adapter.isOutputOutdated$.next(!!this.pendingState.isOuputOutdated);
      }
      this.pendingState = undefined;
    }
    this.instancesWrapper.setInstance(adapter, initValues);
  }

  getStateStore() {
    return this.instancesWrapper;
  }

  restoreState(state: PipelineSerializedState) {
    if (!isFuncCallSerializedState(state))
      throw new Error(`Wrong FuncCall node state ${JSON.stringify(state)}`);
    this.pendingState = state;
    this.pendingId = state.funcCallId;
    if (state.uuid)
      this.uuid = state.uuid;
  }

  initState(initialConfig: StepFunCallInitialConfig) {
    this.instancesWrapper.inputRestrictions$.next(initialConfig.inputRestrictions ?? {});
    this.instancesWrapper.initialValues = initialConfig.values ?? {};
  }

  toState(options: StateTreeSerializationOptions): StepFunCallState {
    const instance = this.instancesWrapper.getInstance();
    const res: StepFunCallState = {
      type: 'funccall',
      uuid: this.uuid,
      nqName: this.config.nqName,
      configId: this.config.id,
      friendlyName: this.config.friendlyName,
      funcCall: instance?.getFuncCall(),
      funcCallId: instance?.getFuncCall()?.id,
      isRunning: this.instancesWrapper.isRunning$.value,
      isRunable: this.instancesWrapper.isRunable$.value,
      isOuputOutdated: this.instancesWrapper.isOutputOutdated$.value,
      inputRestrictions: {}, // TODO: Remove later
      isReadonly: this.isReadonly,
    };
    if (options.disableNodesUUID)
      res.uuid = '';
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
      inputRestrictions: this.instancesWrapper.inputRestrictions$.value,
      isReadonly: this.isReadonly,
    };
    if (options.disableNodesUUID)
      res.uuid = '';
    if (options.disableCallsUUID)
      res.funcCallId = '';
    return res;
  }

  close() {
    this.instancesWrapper.close();
    this.closed$.next(true);
  }

  private getConsistencyUpdater(instance: IFuncCallAdapter | undefined): Observable<Record<string, ConsistencyInfo>> {
    if (!instance)
      return of({});
    const valueUpdates = instance.getStateNames().map((name) => instance.getStateChanges(name).pipe(skip(1), mapTo(name)));
    const restrictionUpdates = instance.inputRestrictionsUpdates$.pipe(map(([name]) => name));
    const inputUpdatesChecker$ = merge(...[...valueUpdates, restrictionUpdates]).pipe(
      withLatestFrom(this.instancesWrapper.inputRestrictions$),
      map(([name, restrictions]) => [name, this.getConsistencyState(instance, name, restrictions)] as const),
    );
    const state$ = inputUpdatesChecker$.pipe(
      scan((acc, [name, info]) => {
        if (!info) {
          const {[name]: omitted, ...rest} = acc;
          return rest;
        } else
          return {...acc, [name]: info};
      }, this.initConsistencyStates(instance, this.instancesWrapper.inputRestrictions$.value)),
    );
    return state$;
  }

  private getValidationsUpdater(instance: IFuncCallAdapter | undefined): Observable<Record<string, ValidationResultBase>> {
    if (!instance)
      return of({});
    return instance.validations$.pipe(map((validations) => this.convertValidations(validations)));
  }

  private convertValidations(validationsIn: Record<string, Record<string, ValidationResultBase | undefined>>) {
    const validationArrays = Object.values(validationsIn).reduce((acc, val) => {
      for (const [k, v] of Object.entries(val)) {
        if (v) {
          if (acc[k])
            acc[k].push(v);
          else
            acc[k] = [v];
        }
      }
      return acc;
    }, {} as Record<string, ValidationResultBase[]>);
    const validationEntries = Object.entries(validationArrays).map(([k, validations]) => [k, mergeValidationResultsBase(...validations)] as const);
    const validations = Object.fromEntries(validationEntries);
    return validations;
  }

  private initConsistencyStates(instance: IFuncCallAdapter, restrictions: Record<string, RestrictionState | undefined>) {
    const res: Record<string, ConsistencyInfo> = {};
    for (const name of instance.getStateNames()) {
      const cinfo = this.getConsistencyState(instance, name, restrictions);
      if (cinfo)
        res[name] = cinfo;
    }
    return res;
  }

  private getConsistencyState(instance: IFuncCallAdapter, inputName: string, restrictions: Record<string, RestrictionState | undefined>): ConsistencyInfo | undefined {
    const restriction = restrictions[inputName];
    if (!restriction)
      return undefined;
    else {
      const {assignedValue, type} = restriction;
      const currentVal = instance.getState(inputName);
      if (!this.deepEq(assignedValue, currentVal)) {
        return {
          restriction: type,
          inconsistent: true,
          assignedValue,
        };
      } else
        return undefined;
    }
  }

  private deepEq(val1: any, val2: any) {
    try {
      expectDeepEqual(val1, val2);
    } catch {
      return false;
    }
    return true;
  }
}

export class PipelineNodeBase implements IStoreProvider {
  public uuid = uuidv4();
  private store: MemoryStore;

  constructor(
    public readonly config: PipelineConfigurationProcessed,
    public readonly isReadonly: boolean,
  ) {
    this.store = new MemoryStore(config.states ?? [], false);
  }

  restoreState(state: PipelineSerializedState) {
    if (isFuncCallSerializedState(state))
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
      isReadonly: this.isReadonly,
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
    public readonly isReadonly: boolean,
  ) {
    super(config, isReadonly);
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
    public readonly isReadonly: boolean,
  ) {
    super(config, isReadonly);
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
    public readonly isReadonly: boolean,
  ) {
    super(config, isReadonly);
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

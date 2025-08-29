import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, combineLatest, merge, Observable, Subject, of} from 'rxjs';
import dayjs from 'dayjs';
import {v4 as uuidv4} from 'uuid';
import {PipelineStateParallel, PipelineStateSequential, PipelineStateStatic, StepFunCallInitialConfig, StepFunCallSerializedState, StepFunCallState, PipelineSerializedState, isFuncCallSerializedState, ViewAction, PipelineInstanceRuntimeData, PipelineOutline} from '../config/PipelineInstance';
import {PipelineConfigurationParallelProcessed, PipelineConfigurationProcessed, PipelineConfigurationSequentialProcessed, PipelineConfigurationStaticProcessed} from '../config/config-processing-utils';
import {IFuncCallAdapter, IStateStore, MemoryStore} from './FuncCallAdapters';
import {FuncCallInstancesBridge, RestrictionState} from './FuncCallInstancesBridge';
import {isPipelineConfig, isPipelineStepConfig, PipelineStepConfigurationProcessed} from '../config/config-utils';
import {map, mapTo, scan, skip, switchMap, takeUntil, withLatestFrom} from 'rxjs/operators';
import {RestrictionType, ValidationResult} from '../data/common-types';
import {customDeepEqual, mergeValidationResults} from '../utils';

export const descriptionOutputs = ['title', 'description', 'tags'] as const;
const descriptionStates = descriptionOutputs.map((id) => ({id}));

export type StateTreeSerializationOptions = {
  disableNodesUUID?: boolean,
  disableCallsUUID?: boolean,
  skipFuncCalls?: boolean,
};

export type ConsistencyInfo = {
  restriction: RestrictionType,
  inconsistent: boolean,
  assignedValue: any,
}

export type FuncCallStateInfo = {
  isRunning: boolean,
  isRunnable: boolean,
  isOutputOutdated: boolean,
  runError: string | undefined,
  pendingDependencies: string[],
}

export type MetaCallInfo = {
  id?: string,
  title?: string,
  description?: string,
  isFavorite?: boolean,
  tags?: string[],
  started?: dayjs.Dayjs
}

export interface AdapterInitData {
  adapter: IFuncCallAdapter,
  restrictions: Record<string, RestrictionState | undefined>,
  isOutputOutdated: boolean,
  runError?: string,
}

export interface IStoreProvider {
  getStateStore(): IStateStore;
}

export class NodeMetaDescription extends MemoryStore {}

export class FuncCallNode implements IStoreProvider {
  private depsData$ = new BehaviorSubject<(readonly [string, Observable<boolean>])[]>([]);

  public uuid = uuidv4();
  public readonly nodeType = 'funccall';

  public instancesWrapper = new FuncCallInstancesBridge(this.config.io!, this.isReadonly);
  public pendingId?: string;

  public consistencyInfo$ = new BehaviorSubject<Record<string, ConsistencyInfo>>({});
  public validationInfo$ = new BehaviorSubject<Record<string, ValidationResult>>({});
  public metaInfo$ = new BehaviorSubject<Record<string, BehaviorSubject<Record<string, any> | undefined>>>({});
  public funcCallState$ = new BehaviorSubject<FuncCallStateInfo | undefined>(undefined);
  public pendingDependencies$ = new BehaviorSubject<string[]>([]);
  public nodeDescription = new NodeMetaDescription(descriptionStates, false);

  private closed$ = new Subject<true>();

  constructor(
    public readonly config: PipelineStepConfigurationProcessed,
    public isReadonly: boolean,
  ) {
    this.getConsistencyChanges().pipe(
      takeUntil(this.closed$),
    ).subscribe(this.consistencyInfo$);

    this.instancesWrapper.validations$.pipe(
      map((validations) => this.convertValidations(validations)),
      takeUntil(this.closed$),
    ).subscribe(this.validationInfo$);

    this.metaInfo$.next(this.instancesWrapper.meta);

    combineLatest([
      this.instancesWrapper.isRunning$,
      this.instancesWrapper.isRunable$,
      this.instancesWrapper.isOutputOutdated$,
      this.instancesWrapper.runError$,
      this.pendingDependencies$,
    ]).pipe(
      map(([isRunning, isRunnable, isOutputOutdated, runError, pendingDependencies]) =>
        ({isRunning, isRunnable, isOutputOutdated, runError, pendingDependencies})),
      takeUntil(this.closed$),
    ).subscribe(this.funcCallState$);

    this.depsData$.pipe(
      switchMap((deps) => this.getPendingDeps(deps)),
      takeUntil(this.closed$),
    ).subscribe(this.pendingDependencies$);
  }

  initAdapter(
    data: AdapterInitData,
    initValues: boolean,
  ) {
    this.instancesWrapper.init({...data, initValues});
  }

  changeAdapter(adapter: IFuncCallAdapter, isNew = false) {
    this.instancesWrapper.change(adapter, isNew);
  }

  setDeps(deps: (readonly [string, Observable<boolean>])[]) {
    this.depsData$.next(deps);
  }

  getStateStore() {
    return this.instancesWrapper;
  }

  restoreState(state: PipelineSerializedState) {
    if (!isFuncCallSerializedState(state))
      throw new Error(`Wrong FuncCall node state ${JSON.stringify(state)}`);
    this.pendingId = state.funcCallId;
    if (state.uuid)
      this.uuid = state.uuid;
  }

  initState(initialConfig: StepFunCallInitialConfig) {
    const initialRestrictions: Record<string, RestrictionState> = {};
    for (const [k, restrictionType] of Object.entries(initialConfig.inputRestrictions ?? {})) {
      const initalVal = (initialConfig.initialValues ?? {})[k];
      initialRestrictions[k] = {
        assignedValue: initalVal,
        type: restrictionType,
      };
    }
    this.instancesWrapper.setPreInitialData({
      initialValues: initialConfig.initialValues ?? {},
      initialRestrictions,
    });
  }

  toState(options: StateTreeSerializationOptions, actions?: ViewAction[]): StepFunCallState {
    const instance = this.instancesWrapper.getInstance();
    const res: StepFunCallState = {
      type: 'funccall',
      uuid: this.uuid,
      configId: this.config.id,
      friendlyName: this.config.friendlyName,
      funcCall: options?.skipFuncCalls ? undefined : instance?.getFuncCall(),
      isReadonly: this.isReadonly,
      viewersHook: this.config.viewersHook,
      actions,
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
      isReadonly: this.isReadonly,
    };
    if (options.disableNodesUUID)
      res.uuid = '';
    if (options.disableCallsUUID)
      res.funcCallId = '';
    return res;
  }

  clearOldMetas(currentIds: Set<string>) {
    for (const meta$ of Object.values(this.instancesWrapper.metaStates ?? {})) {
      if (!meta$?.value)
        continue;
      let needsUpdate = false;
      for (const handlerId of Object.keys(meta$.value)) {
        if (!currentIds.has(handlerId)) {
          delete meta$.value[handlerId];
          needsUpdate = true;
        }
      }
      if (needsUpdate)
        meta$.next(Object.keys(meta$.value).length === 0 ? undefined : meta$.value);
    }
  }

  clearIORestriction(name: string) {
    this.instancesWrapper.removeRestriction(name);
  }

  clearOldValidations(currentIds: Set<string>) {
    const cval = this.instancesWrapper.validations$.value;
    const nval: Record<string, Record<string, ValidationResult | undefined>> = {};
    let needsUpdate = false;
    for (const [k, v] of Object.entries(cval)) {
      if (currentIds.has(k))
        nval[k] = v;
      else
        needsUpdate = true;
    }
    if (needsUpdate)
      this.instancesWrapper.validations$.next(nval);
  }

  close() {
    this.instancesWrapper.close();
    this.closed$.next(true);
  }

  private convertValidations(
    validationsIn: Record<string, Record<string, ValidationResult | undefined>>,
  ) {
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
    }, {} as Record<string, ValidationResult[]>);
    const validationEntries = Object.entries(validationArrays).map(
      ([k, validations]) => [k, mergeValidationResults(...validations)] as const);
    const validations = Object.fromEntries(validationEntries);
    return validations;
  }

  private getConsistencyChanges(): Observable<Record<string, ConsistencyInfo>> {
    const valueUpdates = this.instancesWrapper.getStateNames().map(
      (name) => this.instancesWrapper.getStateChanges(name, true).pipe(skip(1), mapTo(name)));
    const restrictionUpdates = this.instancesWrapper.inputRestrictionsUpdates$.pipe(map(([name]) => name));
    const inputUpdatesChecker$ = merge(...[...valueUpdates, restrictionUpdates]).pipe(
      withLatestFrom(this.instancesWrapper.inputRestrictions$),
      map(([name, restrictions]) => [name, this.getConsistencyState(name, restrictions)] as const),
    );
    const state$ = inputUpdatesChecker$.pipe(
      scan((acc, [name, info]) => {
        if (!info) {
          const {[name]: omitted, ...rest} = acc;
          return rest;
        } else
          return {...acc, [name]: info};
      }, this.initConsistencyStates(this.instancesWrapper.inputRestrictions$.value)),
    );
    return state$;
  }

  private initConsistencyStates(
    restrictions: Record<string, RestrictionState | undefined>,
  ) {
    const res: Record<string, ConsistencyInfo> = {};
    for (const name of this.instancesWrapper.getStateNames()) {
      const cinfo = this.getConsistencyState(name, restrictions);
      if (cinfo)
        res[name] = cinfo;
    }
    return res;
  }

  private getConsistencyState(
    inputName: string,
    restrictions: Record<string, RestrictionState | undefined>,
  ): ConsistencyInfo | undefined {
    const restriction = restrictions[inputName];
    if (!restriction)
      return undefined;
    else {
      const {assignedValue, type} = restriction;
      const currentVal = this.instancesWrapper.getState(inputName);
      const inconsistent = !customDeepEqual(currentVal, assignedValue);
      return {
        restriction: type,
        inconsistent,
        assignedValue,
      };
    }
  }

  private getPendingDeps(deps: (readonly [string, Observable<boolean>])[]) {
    const pendingStates = deps.map(([uuid, state$]) => state$.pipe(map((state) => [uuid, state] as const)));
    if (pendingStates.length === 0)
      return of([]);
    const pending$ = merge(...pendingStates).pipe(
      scan((acc, [uuid, isPending]) => {
        if (isPending)
          acc.add(uuid);
        else
          acc.delete(uuid);

        return acc;
      }, new Set<string>()),
      map((s) => [...s]),
    );
    return pending$;
  }
}

export class PipelineNodeBase implements IStoreProvider {
  public uuid = uuidv4();
  private store: MemoryStore;
  public nodeDescription = new NodeMetaDescription(descriptionStates, false);

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

  clearOldTags(currentIds: Set<string>) {
    const cval = this.nodeDescription.getState<Record<string, string[]>>('tags') ?? [];
    const nval: Record<string, string[]> = {};
    let needsUpdate = false;
    for (const [k, v] of Object.entries(cval)) {
      if (currentIds.has(k))
        nval[k] = v;
      else
        needsUpdate = true;
    }
    if (needsUpdate)
      this.nodeDescription.setState('tags', nval);
  }

  toSerializedState(options: StateTreeSerializationOptions) {
    const res = {
      configId: this.config.id,
      uuid: this.uuid,
      version: this.config.version,
      nqName: this.config.nqName,
      friendlyName: this.config.friendlyName,
      isReadonly: this.isReadonly,
    };
    if (options.disableNodesUUID)
      res.uuid = '';
    return res;
  }

  toState(options: StateTreeSerializationOptions) {
    const state = this.toSerializedState(options);
    const res = {
      ...state,
      nqName: this.config.nqName,
      approversGroup: this.config.approversGroup,
      disableHistory: !!this.config.disableHistory,
      customExports: this.config.customExports,
    };
    return res;
  }

  getStructureCheck(state: PipelineOutline) {
    if (this.config.structureCheck) {
      try {
        return this.config.structureCheck(state);
      } catch (e: any) {
        grok.shell.error(e);
        console.error(e);
      }
    }
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

  toState(options: StateTreeSerializationOptions, actions?: ViewAction[]): PipelineStateStatic<StepFunCallState, PipelineInstanceRuntimeData> {
    const state = super.toState(options);
    return {
      ...state,
      type: this.nodeType,
      steps: [],
      actions,
    };
  }

  toSerializedState(options: StateTreeSerializationOptions): PipelineStateStatic<StepFunCallSerializedState, {}> {
    const base = super.toSerializedState(options);
    const res: PipelineStateStatic<StepFunCallSerializedState, {}> = {
      ...base,
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

  toState(options: StateTreeSerializationOptions, actions?: ViewAction[]): PipelineStateParallel<StepFunCallState, PipelineInstanceRuntimeData> {
    const state = super.toState(options);
    const res: PipelineStateParallel<StepFunCallState, PipelineInstanceRuntimeData> = {
      ...state,
      type: this.nodeType,
      steps: [],
      stepTypes: getStepTypes(this.config),
      actions,
    };
    return res;
  }

  toSerializedState(options: StateTreeSerializationOptions): PipelineStateParallel<StepFunCallSerializedState, {}> {
    const base = super.toSerializedState(options);
    const res: PipelineStateParallel<StepFunCallSerializedState, {}> = {
      ...base,
      type: this.nodeType,
      steps: [],
      stepTypes: getStepTypes(this.config),
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

  toState(options: StateTreeSerializationOptions, actions?: ViewAction[]): PipelineStateSequential<StepFunCallState, PipelineInstanceRuntimeData> {
    const state = super.toState(options);
    const res: PipelineStateSequential<StepFunCallState, PipelineInstanceRuntimeData> = {
      ...state,
      type: this.nodeType,
      steps: [],
      stepTypes: getStepTypes(this.config),
      actions,
    };
    return res;
  }

  toSerializedState(options: StateTreeSerializationOptions): PipelineStateSequential<StepFunCallSerializedState, {}> {
    const base = super.toSerializedState(options);
    const res: PipelineStateSequential<StepFunCallSerializedState, {}> = {
      ...base,
      type: this.nodeType,
      steps: [],
      stepTypes: getStepTypes(this.config),
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

function getStepTypes(conf: PipelineConfigurationParallelProcessed | PipelineConfigurationSequentialProcessed) {
  return conf.stepTypes.map((s) => {
    if (isPipelineConfig(s) || isPipelineStepConfig(s)) {
      const {id: configId, disableUIAdding, disableUIDragging, disableUIRemoving, nqName, friendlyName} = s;
      return {configId, disableUIAdding, disableUIDragging, disableUIRemoving, nqName, friendlyName};
    } else {
      const {id: configId, disableUIAdding, disableUIDragging, disableUIRemoving} = s;
      return {configId, disableUIAdding, disableUIDragging, disableUIRemoving};
    }
  });
}

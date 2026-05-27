import * as DG from 'datagrok-api/dg';
import {DynamicPipelineType, isDynamicType, ItemId, NqName, RestrictionType, ValidationResult} from '../data/common-types';
import {ActionInfoBase, CustomExport, NestedItemContext, ViewersHook} from './PipelineConfiguration';

//
// initial steps config for dynamic pipelines
//

export type StepDynamicInitialConfig = {
  id: ItemId;
  initialValues?: Record<string, any>;
  inputRestrictions?: Record<string, RestrictionType>;
}

export type StepFunCallInitialConfig = {
  id: ItemId;
  initialValues?: Record<string, any>;
  inputRestrictions?: Record<string, RestrictionType>;
}

export type InstanceConfRec<C> = {
  steps?: InstanceConfRec<C>[];
} & C;

export type InstanceConfRecInput<C> = {
  steps?: Array<ItemId | InstanceConfRecInput<C>>;
} & C;

export type PipelineInstanceConfig = InstanceConfRec<StepDynamicInitialConfig | StepFunCallInitialConfig>;
export type PipelineInstanceConfigInput = InstanceConfRecInput<StepDynamicInitialConfig | StepFunCallInitialConfig>;

export function normalizeIdRef<T extends {id: ItemId}>(s: ItemId | T): T {
  return typeof s === 'string' ? ({id: s} as T) : s;
}

export function normalizePipelineInstanceConfig(c: PipelineInstanceConfigInput): PipelineInstanceConfig {
  return {
    ...c,
    steps: c.steps?.map((s) => normalizePipelineInstanceConfig(normalizeIdRef(s))),
  };
}


//
// instance state for view/serialization
//

export interface PipelineInstanceState {
  provider: Function | NqName;
  version?: string;
  state: PipelineState;
};

export type StateTypes = PipelineState['type'];

export type PipelineOutline = PipelineStateRec<StepFunCallStateBase, {}>;
export type PipelineState = PipelineStateRec<StepFunCallState, PipelineInstanceRuntimeData>;
export type PipelineSerializedState = PipelineStateRec<StepFunCallSerializedState, {}>;

export function isFuncCallState(state: PipelineState): state is StepFunCallState {
  return state.type === 'funccall';
}

export function isFuncCallSerializedState(state: PipelineSerializedState): state is StepFunCallSerializedState {
  return state.type === 'funccall';
}

export function isStaticPipelineState(state: PipelineState): state is PipelineStateStatic<StepFunCallState, PipelineInstanceRuntimeData> {
  return state.type === 'static';
}

export function isDynamicPipelineState(state: PipelineState): state is PipelineStateDynamic<StepFunCallState, PipelineInstanceRuntimeData> {
  return isDynamicType(state.type);
}

export function isStaticSerializedPipelineState(state: PipelineSerializedState): state is PipelineStateStatic<StepFunCallSerializedState, {}> {
  return state.type === 'static';
}

export function isDynamicSerializedPipelineState(state: PipelineSerializedState): state is PipelineStateDynamic<StepFunCallSerializedState, {}> {
  return isDynamicType(state.type);
}

export type PipelineStateRec<S, T> = PipelineStateStatic<S, T> | PipelineStateDynamic<S, T> | S;

// funccall

export type ViewAction = ActionInfoBase & {
  uuid: string;
  /** Result of evaluating showWhen / hideWhen against the current tree.
   *  `true` when no condition is set. Compute2 uses this to decide whether
   *  to render the action; RTD itself does not filter or gate on it. */
  visible: boolean;
};

export type StepFunCallStateBase = {
  type: 'funccall';
  configId: string;
  uuid: string;
  friendlyName?: string;
  isReadonly: boolean;
}

export type StepFunCallSerializedState = {
  nqName: string;
  funcCallId?: string;
} & StepFunCallStateBase;


export type StepFunCallState = {
  funcCall?: DG.FuncCall;
  viewersHook?: ViewersHook;
  actions?: ViewAction[];
  enableHistory?: boolean;
} & StepFunCallStateBase;

// pipeline base

export type PipelineInstanceRuntimeData = {
  actions: ViewAction[] | undefined;
  approversGroup: string | undefined;
  disableHistory: boolean;
  customExports: CustomExport[] | undefined;
  disableDefaultExport?: boolean;
  structureCheckResults?: ValidationResult | undefined;
  forceNavigate: boolean;
}

export type PipelineInstanceBase<I, T> = {
  uuid: string;
  configId: string;
  isReadonly: boolean;
  friendlyName: string | undefined;
  description: string | undefined;
  version: string | undefined;
  nqName: string | undefined;
} & I & T;

// static

export type PipelineStateStatic<S, T> = PipelineInstanceBase<{
  type: 'static',
  steps: PipelineStateRec<S, T>[];
  isActionStep?: boolean;
}, T>;

// dynamic (unified type for both parallel and sequential)

export type StepDynamicDescription = {
  configId: string;
  nqName?: string;
  friendlyName?: string;
} & NestedItemContext;

export type StepDynamicState<S, T> = PipelineStateRec<S, T> & StepDynamicDescription;

export type PipelineStateDynamic<S, T> = PipelineInstanceBase<{
  type: DynamicPipelineType;
  steps: StepDynamicState<S, T>[];
  stepTypes: StepDynamicDescription[];
}, T>;

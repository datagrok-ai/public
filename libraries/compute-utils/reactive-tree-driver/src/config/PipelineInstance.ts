import * as DG from 'datagrok-api/dg';
import {ItemId, NqName, RestrictionType, ValidationResult} from '../data/common-types';
import {ActionInfo, CustomExport, NestedItemContext, ViewersHook} from './PipelineConfiguration';

//
// initial steps config for dynamic pipelines
//

export type StepDynamicInitialConfig = {
  id: ItemId;
  initialValues?: Record<string, any>;
  inputRestrictions?: Record<string, RestrictionType>;
}

/** @deprecated Use StepDynamicInitialConfig */
export type StepParallelInitialConfig = StepDynamicInitialConfig;
/** @deprecated Use StepDynamicInitialConfig */
export type StepSequentialInitialConfig = StepDynamicInitialConfig;

export type StepFunCallInitialConfig = {
  id: ItemId;
  initialValues?: Record<string, any>;
  inputRestrictions?: Record<string, RestrictionType>;
}

export type InstanceConfRec<C> = {
  steps?: InstanceConfRec<C>[];
} & C;

export type PipelineInstanceConfig = InstanceConfRec<StepDynamicInitialConfig | StepFunCallInitialConfig>;


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
  return state.type === 'dynamic' || state.type === 'parallel' || state.type === 'sequential';
}

/** @deprecated Use isDynamicPipelineState */
export const isParallelPipelineState = isDynamicPipelineState;
/** @deprecated Use isDynamicPipelineState */
export const isSequentialPipelineState = isDynamicPipelineState;

export function isStaticSerializedPipelineState(state: PipelineSerializedState): state is PipelineStateStatic<StepFunCallSerializedState, {}> {
  return state.type === 'static';
}

export function isDynamicSerializedPipelineState(state: PipelineSerializedState): state is PipelineStateDynamic<StepFunCallSerializedState, {}> {
  return state.type === 'dynamic' || state.type === 'parallel' || state.type === 'sequential';
}

/** @deprecated Use isDynamicSerializedPipelineState */
export const isParallelSerializedPipelineState = isDynamicSerializedPipelineState;
/** @deprecated Use isDynamicSerializedPipelineState */
export const isSequentialSerializedPipelineState = isDynamicSerializedPipelineState;

export type PipelineStateRec<S, T> = PipelineStateStatic<S, T> | PipelineStateDynamic<S, T> | S;

// funccall

export type ViewAction = {
  uuid: string;
} & ActionInfo;

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
  version: string | undefined;
  nqName: string | undefined;
} & I & T;

// static

export type PipelineStateStatic<S, T> = PipelineInstanceBase<{
  type: 'static',
  steps: PipelineStateRec<S, T>[];
}, T>;

// dynamic (unified type for both parallel and sequential)

export type StepDynamicDescription = {
  configId: string;
  nqName?: string;
  friendlyName?: string;
} & NestedItemContext;

export type StepDynamicState<S, T> = PipelineStateRec<S, T> & StepDynamicDescription;

export type PipelineStateDynamic<S, T> = PipelineInstanceBase<{
  type: 'dynamic' | 'parallel' | 'sequential';
  steps: StepDynamicState<S, T>[];
  stepTypes: StepDynamicDescription[];
}, T>;

/** @deprecated Use StepDynamicDescription */
export type StepSequentialDescription = StepDynamicDescription;
/** @deprecated Use StepDynamicDescription */
export type StepParallelDescription = StepDynamicDescription;
/** @deprecated Use StepDynamicState */
export type StepSequentialState<S, T> = StepDynamicState<S, T>;
/** @deprecated Use StepDynamicState */
export type StepParallelState<S, T> = StepDynamicState<S, T>;
/** @deprecated Use PipelineStateDynamic */
export type PipelineStateSequential<S, T> = PipelineStateDynamic<S, T>;
/** @deprecated Use PipelineStateDynamic */
export type PipelineStateParallel<S, T> = PipelineStateDynamic<S, T>;

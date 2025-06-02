import * as DG from 'datagrok-api/dg';
import {ItemId, NqName, RestrictionType, ValidationResult} from '../data/common-types';
import {ActionInfo, CustomExport, NestedItemContext, ViewersHook} from './PipelineConfiguration';

//
// initial steps config for dynamic pipelines
//

export type StepParallelInitialConfig = {
  id: ItemId;
  initialValues?: Record<string, any>;
  inputRestrictions?: Record<string, RestrictionType>;
}

export type StepSequentialInitialConfig = {
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

export type PipelineInstanceConfig = InstanceConfRec<StepParallelInitialConfig | StepSequentialInitialConfig | StepFunCallInitialConfig>;


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

export function isParallelPipelineState(state: PipelineState): state is PipelineStateParallel<StepFunCallState, PipelineInstanceRuntimeData> {
  return state.type === 'parallel';
}

export function isSequentialPipelineState(state: PipelineState): state is PipelineStateSequential<StepFunCallState, PipelineInstanceRuntimeData> {
  return state.type === 'sequential';
}

export function isStaticSerializedPipelineState(state: PipelineSerializedState): state is PipelineStateStatic<StepFunCallSerializedState, {}> {
  return state.type === 'static';
}

export function isParallelSerializedPipelineState(state: PipelineSerializedState): state is PipelineStateParallel<StepFunCallSerializedState, {}> {
  return state.type === 'parallel';
}

export function isSequentialSerializedPipelineState(state: PipelineSerializedState): state is PipelineStateSequential<StepFunCallSerializedState, {}> {
  return state.type === 'sequential';
}

export type PipelineStateRec<S, T> = PipelineStateStatic<S, T> | PipelineStateSequential<S, T> | PipelineStateParallel<S, T> | S;

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
  customExports: CustomExport[] | undefined;
  structureCheckResults?: ValidationResult | undefined;
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

// sequential

export type StepSequentialDescription = {
  configId: string;
  nqName?: string;
  friendlyName?: string;
} & NestedItemContext;

export type StepSequentialState<S, T> = PipelineStateRec<S, T> & StepSequentialDescription;

export type PipelineStateSequential<S, T> = PipelineInstanceBase<{
  type: 'sequential';
  steps: StepSequentialState<S, T>[];
  stepTypes: StepSequentialDescription[];
}, T>;

// parallel

export type StepParallelDescription = {
  configId: string;
  nqName?: string;
  friendlyName?: string;
} & NestedItemContext;

export type StepParallelState<S, T> = PipelineStateRec<S, T> & StepParallelDescription;

export type PipelineStateParallel<S, T> = PipelineInstanceBase<{
  type: 'parallel';
  steps: StepParallelState<S, T>[];
  stepTypes: StepParallelDescription[];
}, T>;

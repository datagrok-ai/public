import * as DG from 'datagrok-api/dg';
import {ItemId, NqName, RestrictionType} from '../data/common-types';
import {Action} from '../runtime/Link';

//
// initial steps config for dynamic pipelines
//

export type StepParallelInitialConfig = {
  id: ItemId;
  disableUIRemoving?: boolean;
}

export type StepSequentialInitialConfig = {
  id: ItemId;
  disableUIRemoving?: boolean;
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

export type PipelineState = PipelineStateRec<StepFunCallState>;
export type PipelineSerializedState = PipelineStateRec<StepFunCallSerializedState>;

export function isFuncCallState(state: PipelineState): state is StepFunCallState {
  return state.type === 'funccall';
}
export function isFuncCallSerializedState(state: PipelineSerializedState): state is StepFunCallSerializedState {
  return state.type === 'funccall';
}


export function isStaticPipelineState(state: PipelineState): state is PipelineStateStatic<StepFunCallState> {
  return state.type === 'static';
}
export function isParallelPipelineState(state: PipelineState): state is PipelineStateParallel<StepFunCallState> {
  return state.type === 'parallel';
}
export function isSequentialPipelineState(state: PipelineState): state is PipelineStateSequential<StepFunCallState> {
  return state.type === 'sequential';
}

export type PipelineStateRec<S> = PipelineStateStatic<S> | PipelineStateSequential<S> | PipelineStateParallel<S> | S;


// funccall

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
  actions?: Action[];
} & StepFunCallStateBase;

// pipeline base

export type PipelineInstanceBase<I> = {
  uuid: string;
  configId: string;
  isReadonly: boolean;
  friendlyName: string | undefined;
  provider: NqName | undefined;
  version: string | undefined;
  nqName: string | undefined;
} & I;

// static

export type PipelineStateStatic<S> = PipelineInstanceBase<{
  type: 'static',
  steps: PipelineStateRec<S>[];
}>;

// sequential

export type StepSequentialDescription = {
  configId: string;
  nqName?: string;
  friendlyName?: string;
  disableUIAdding?: boolean;
};

export type StepSequentialState<S> = PipelineStateRec<S> & StepSequentialDescription;

export type PipelineStateSequential<S> = PipelineInstanceBase<{
  type: 'sequential';
  steps: StepSequentialState<S>[];
  stepTypes: StepSequentialDescription[];
}>;

// parallel

export type StepParallelDescription = {
  configId: string;
  nqName?: string;
  friendlyName?: string;
  disableUIAdding?: boolean;
};

export type StepParallelState<S> = PipelineStateRec<S> & StepParallelDescription;

export type PipelineStateParallel<S> = PipelineInstanceBase<{
  type: 'parallel';
  steps: StepParallelState<S>[];
  stepTypes: StepParallelDescription[];
}>;

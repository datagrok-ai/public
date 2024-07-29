import * as DG from 'datagrok-api/dg';
import {InputState, ItemId, NqName} from '../data/common-types';
import { ValidationResultBase } from '../../../shared-utils/validation';

//
// initial steps config for dynamic pipelines
//

export type StepParallelInitialConfig = {
  id: ItemId;
  allowRemoving: boolean;
}

export type StepSequentialInitialConfig = {
  id: ItemId;
  allowRemoving: boolean;
}

export type StepFunCallInitialConfig = {
  id: ItemId;
  values?: Record<string, any>;
  inputStates?: Record<string, InputState>;
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

export type PipelineState = PipelineStateStatic | PipelineStateSequential | PipelineStateParallel | StepFunCallState;
export type PipelineSerializedState = PipelineStateStatic | PipelineStateSequential | PipelineStateParallel | StepFunCallSerializedState;

export function isFuncCallState(state: PipelineState): state is StepFunCallState {
  return state.type === 'funccall';
}
export function isStaticPipelineState(state: PipelineState): state is PipelineStateStatic {
  return state.type === 'static';
}
export function isParallelPipelineState(state: PipelineState): state is PipelineStateParallel {
  return state.type === 'parallel';
}
export function isSequentialPipelineState(state: PipelineState): state is PipelineStateSequential {
  return state.type === 'static';
}


// funccall

export type StepFunCallSerializedState = {
  type: 'funccall';
  configId: string;
  uuid: string;
  nqName: string;
  friendlyName?: string;
  funcCallId?: string;
  isOuputOutdated?: boolean;
  isCurrent?: boolean;
};

export type StepFunCallState = {
  funcCall?: DG.FuncCall;
  validations?: Record<string, ValidationResultBase>;
  isLoading?: boolean;
  isRunning?: boolean;
  isRunable?: boolean;
} & StepFunCallSerializedState;

// pipeline base

export type PipelineInstanceBase<S> = {
  uuid: string;
  configId: string;
  friendlyName?: string;
  nqName?: string;
} & S;

// static

export type PipelineStateStatic = PipelineInstanceBase<{
  type: 'static',
  steps: PipelineState[];
}>;

// sequential

export type StepSequentialDescription = {
  configId: string;
  allowAdding: boolean;
};

export type StepSequentialState = PipelineState & StepSequentialDescription;

export type PipelineStateSequential = PipelineInstanceBase<{
  type: 'sequential';
  steps: StepSequentialState[];
  stepTypes: StepSequentialDescription[];
}>;

// parallel

export type StepParallelDescription = {
  configId: string;
  allowAdding: boolean;
};

export type StepParallelState = PipelineState & StepParallelDescription;

export type PipelineStateParallel = PipelineInstanceBase<{
  type: 'parallel';
  steps: StepParallelState[];
  stepTypes: StepParallelDescription[];
}>;

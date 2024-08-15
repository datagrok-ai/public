import * as DG from 'datagrok-api/dg';
import {ItemId, NqName} from '../data/common-types';
import {ValidationResultBase} from '../../../shared-utils/validation';
import {RestrictionState} from '../runtime/FuncCallAdapters';
import {ActionPositions} from './PipelineConfiguration';

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
  values?: Record<string, any>;
  inputRestrictions?: Record<string, RestrictionState>;
}

export type InstanceConfRec<C> = {
  steps?: InstanceConfRec<C>[];
} & C;

export type PipelineInstanceConfig = InstanceConfRec<StepParallelInitialConfig | StepSequentialInitialConfig | StepFunCallInitialConfig>;


//
// instance state for view/serialization
//

// TODO: save backref to pipeline provider when possible
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
export function isFuncCallSerializedState(state: PipelineSerializedState): state is StepFunCallSerializedState {
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

export const pipelineMenuCategories = ['export', 'mock', 'template'] as const;
export type PipelineMenuCategories = typeof pipelineMenuCategories[number];


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

export type StepAction = {
  uuid: string,
  friendlyName?: string;
  position: ActionPositions,
}

export type StepFunCallState = {
  funcCall?: DG.FuncCall;
  validations?: Record<string, ValidationResultBase>;
  inputRestrictions?: Record<string, RestrictionState | undefined>;
  isLoading?: boolean;
  isRunning?: boolean;
  isRunable?: boolean;
  actions?: StepAction;
} & StepFunCallSerializedState;

// pipeline base

export type PipelineAction = {
  uuid: string,
  friendlyName?: string;
  position: ActionPositions,
}

export type PipelineInstanceBase<S> = {
  uuid: string;
  configId: string;
  isReadonly?: boolean;
  friendlyName?: string;
  provider: NqName | undefined;
  version: string | undefined;
  nqName: string | undefined;
} & S;

// static

export type PipelineStateStatic = PipelineInstanceBase<{
  type: 'static',
  steps: PipelineState[];
}>;

// sequential

export type StepSequentialDescription = {
  configId: string;
  nqName?: string;
  friendlyName?: string;
  disableUIAdding?: boolean;
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
  nqName?: string;
  friendlyName?: string;
  disableUIAdding?: boolean;
};

export type StepParallelState = PipelineState & StepParallelDescription;

export type PipelineStateParallel = PipelineInstanceBase<{
  type: 'parallel';
  steps: StepParallelState[];
  stepTypes: StepParallelDescription[];
}>;

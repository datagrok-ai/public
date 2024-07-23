import * as DG from 'datagrok-api/dg';
import { ItemId, ItemType, NqName } from '../data/common-types';

//
// initial steps config for dynamic pipelines
//

export type StepParallelInitialConfig = {
  type: ItemType;
  allowRemoving: boolean;
}

export type StepSequentialInitialConfig = {
  type: ItemType;
  allowRemoving: boolean;
}

export type StepFunCallInitialConfig = {
  id: string;
  defaultInputs?: Record<string,any>;
};

export type ConfRec<S> = {
  steps?: ConfRec<S>[];
}

export type PipelineInstanceConfig = ConfRec<StepParallelInitialConfig | StepSequentialInitialConfig | StepFunCallInitialConfig>;


//
// instance state for view/serialization
//

export interface PipelineInstanceState {
  provider: Function | NqName;
  version?: string;
  state: PipelineState;
};

export type PipelineState = PipelineStateStatic | PipelineStateSequential | PipelineStateParallel | StepFunCallState;

// funcall

export type StepFunCallState = {
  instanceId: string;
  type: 'funccall';
  nqName: string;
  funcCallId?: string;
  funcCall?: DG.FuncCall;
  isInconsistent?: boolean;
  ioSynced?: boolean;
  isCurrent?: boolean;
} & StepFunCallInitialConfig;

// pipeline base

type PipelineInstanceBase<S> = {
  id: string;
  instanceId: string;
  nqName?: string;
} & S;

// static

export type PipelineStateStatic = PipelineInstanceBase<{
  type: 'staticPipeline',
  steps: StepFunCallState[];
}>;

// sequential

export type StepSequentialDescription = {
  typeName: string;
  inputTypeId: ItemId;
  outputTypeId: ItemId;
  dynamic: boolean;
} & StepSequentialInitialConfig;

export type StepSequentialState = {
  state: PipelineState;
} & StepSequentialDescription;

export type PipelineStateSequential = PipelineInstanceBase<{
  type: 'sequentialPipeline';
  steps: StepSequentialState[];
  stepTypes: StepSequentialDescription[];
}>;

// parallel

export type StepParallelDescription = {
  typeName: string;
  dynamic: boolean;
} & StepParallelInitialConfig;

export type StepParallelState = {
  state: PipelineState
} & StepParallelDescription;


export type PipelineStateParallel = PipelineInstanceBase<{
  type: 'parallelPipeline';
  steps: StepParallelState[];
  stepTypes: StepParallelDescription[];
}>;

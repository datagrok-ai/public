import * as DG from 'datagrok-api/dg';
import { InputState, ItemId, ItemPathArray, ItemType, NqName} from '../data/common-types';

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

export type PipelineState = PipelineStateStatic | PipelineStateSequential | PipelineStateParallel | StepFunCallState;

// funcall

export type StepFunCallState = {
  type: 'funccall';
  uuid: string;
  nqName: string;
  configPath: ItemPathArray;
  funcCallId?: string;
  funcCall?: DG.FuncCall;
  isLoading?: boolean;
  isRunning?: boolean;
  isRunable?: boolean;
  isOuputOutdated?: boolean;
  isCurrent?: boolean;
};

// pipeline base

export type PipelineInstanceBase<S> = {
  configPath: ItemPathArray;
  uuid: string;
  nqName: string | undefined;
} & S;

// static

export type PipelineStateStatic = PipelineInstanceBase<{
  type: 'static',
  steps: PipelineState[];
}>;

// sequential

export type StepSequentialDescription = {
  type: string;
  inputTypeId: ItemId;
  outputTypeId: ItemId;
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
  type: string;
  allowAdding: boolean;
};

export type StepParallelState = PipelineState & StepParallelDescription;

export type PipelineStateParallel = PipelineInstanceBase<{
  type: 'parallel';
  steps: StepParallelState[];
  stepTypes: StepParallelDescription[];
}>;

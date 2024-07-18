import * as DG from 'datagrok-api/dg';
import {ItemId, NqName, ItemType} from '../config/CommonTypes';
import {ActionPositions} from '../config/PipelineConfiguration';
import {ValidationResult} from '../../../shared-utils/validation';

// view config

export interface ViewConfig {
  provider: NqName;
  pipelineConfig: StepsPipelineConfig;
};


export type PipelineViewConfigVariants = PipelineViewConfigStatic | PipelineViewConfigParallel | PipelineViewConfigSequential;


export interface StepFunCallConfig {
  id: string;
  nqName: string;
  isStep: true;
  funcCallId?: string;
  funcCall?: DG.FuncCall;
  isInconsistent?: boolean;
  ioSynced?: boolean;
  isCurrent?: boolean;
}

export interface StepsPipelineConfig {
  id: string;
  funcCallId?: string;
  nqName?: string
  isStep: false;
  isCurrent?: boolean;
  pipelineConfig: PipelineViewConfigVariants;
}


export interface StepSequentialType {
  type: ItemType;
  typeName: string;
}

export interface StepSequentialConfig extends StepFunCallConfig, StepSequentialType {}

export interface StepsSequentialConfig extends StepsPipelineConfig, StepSequentialType {}


export interface StepParallelType {
  type: ItemType;
  typeName: string;
  inputTypeId: ItemId;
  outputTypeId: ItemId;
}

export interface StepParallelConfig extends StepFunCallConfig, StepParallelType {}

export interface StepsParallelConfig extends StepsPipelineConfig, StepParallelType {}


export interface PipelineViewConfigStatic {
  id: string;
  dynamic: false;
  isStep: false;
  steps: (StepFunCallConfig | StepsPipelineConfig)[];
}

export interface PipelineViewConfigSequential {
  id: string;
  dynamic: 'sequential';
  isStep: false;
  steps: (StepFunCallConfig | StepsSequentialConfig)[];
  stepTypes: StepSequentialType[];
}

export interface PipelineViewConfigParallel {
  id: string;
  dynamic: 'parallel';
  isStep: false;
  steps: (StepFunCallConfig | StepsParallelConfig)[];
  stepTypes: StepParallelType[];
}


// view config updates

export interface AddGroupItem {
  event: 'addGroupItem';
  id: string;
  type: ItemType;
  insertBefore?: string;
}

export interface RemoveGroupItem {
  event: 'removeGroupItem';
  id: string;
}

export interface MoveGroupItem {
  event: 'moveGroupItem';
  id: string;
  insertBefore?: ItemId;
}

export interface IOSyncChange {
  event: 'ioSyncChange';
  id: string;
  ioSynced: boolean;
}

export interface CurrentStepChange {
  event: 'currentStepChange';
  id: string;
}

export interface ConsistencyChange {
  event: 'consistencyChange';
  id: string;
  input: string;
  value: boolean;
}

export interface FuncCallLoaded {
  event: 'funcCallLoaded';
  id: string;
  funcCall: DG.FuncCall;
}

export type ViewConfigChanges = AddGroupItem | RemoveGroupItem | MoveGroupItem | IOSyncChange | CurrentStepChange | ConsistencyChange | FuncCallLoaded;


// additional controls

export const globalCategories = ['export', 'mock', 'template'] as const;
export type GlobalCategories = typeof globalCategories[number];

export interface Actions {
  stepActions: {
    id: string,
    name: string,
    position: ActionPositions,
    stepId: string;
  }[];
  globalActions: {
    id: string,
    name: string,
    category?: GlobalCategories,
  }[];
}

// additional control events

export interface ActionTriggered {
  event: 'actionClicked',
  id: string,
}

// validatation results

export interface ViewValidationResult {
  results: {
    id: string,
    input: string,
    validation: ValidationResult,
  }[];
}


// utils

export function isStepFunCallConfig(config: StepFunCallConfig | StepsPipelineConfig): config is StepFunCallConfig {
  return config.isStep;
}

export function traverseViewConfig<T>(
  config: ViewConfig,
  handler: (acc: T, node: (StepFunCallConfig | StepsPipelineConfig)) => T,
  acc: T,
): T {
  const stk: (StepFunCallConfig | StepsPipelineConfig)[] = [config.pipelineConfig];

  while (stk.length) {
    const conf = stk.pop()!;
    if (isStepFunCallConfig(conf))
      acc = handler(acc, conf);
    else {
      acc = handler(acc, conf);
      const rsteps = [...conf.pipelineConfig.steps].reverse();
      stk.push(...rsteps);
    }
  }
  return acc;
}

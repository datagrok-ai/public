import * as DG from 'datagrok-api/dg';
import {ItemType} from '../data/common-types';
import {ActionPositions} from '../config/PipelineConfiguration';
import {ValidationResult} from '../../../shared-utils/validation';

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
  insertBefore?: string;
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

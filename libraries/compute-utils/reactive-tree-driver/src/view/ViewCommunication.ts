import * as DG from 'datagrok-api/dg';
import {ItemId} from '../data/common-types';
import {ActionPositions} from '../config/PipelineConfiguration';
import {ValidationResult} from '../../../shared-utils/validation';

// view config updates

export interface AddGroupItem {
  event: 'addGroupItem';
  uuid: string;
  id: ItemId;
  insertBefore?: string;
}

export interface RemoveGroupItem {
  event: 'removeGroupItem';
  uuid: string;
}

export interface MoveGroupItem {
  event: 'moveGroupItem';
  uuid: string;
  insertBefore?: string;
}

export interface IOSyncChange {
  event: 'ioSyncChange';
  uuid: string;
  ioSynced: boolean;
}

export interface CurrentStepChange {
  event: 'currentStepChange';
  uuid: string;
}

export interface ConsistencyChange {
  event: 'consistencyChange';
  uuid: string;
  input: string;
  value: boolean;
}

export interface FuncCallLoaded {
  event: 'funcCallLoaded';
  uuid: string;
  funcCall: DG.FuncCall;
}

export type ViewConfigChanges = AddGroupItem | RemoveGroupItem | MoveGroupItem | IOSyncChange | CurrentStepChange | ConsistencyChange | FuncCallLoaded;


// additional controls

export const globalCategories = ['export', 'mock', 'template'] as const;
export type GlobalCategories = typeof globalCategories[number];

export interface Actions {
  stepActions: {
    uuid: string,
    name: string,
    position: ActionPositions,
    step: string;
  }[];
  globalActions: {
    uuid: string,
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
    uuid: string,
    input: string,
    validation: ValidationResult,
  }[];
}

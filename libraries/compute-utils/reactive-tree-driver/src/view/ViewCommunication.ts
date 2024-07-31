import {ItemId} from '../data/common-types';

// view config update requests

export interface AddGroupItem {
  event: 'addGroupItem';
  parentPipelineUuid: string;
  itemId: ItemId;
  postion?: number;
}

export interface RemoveGroupItem {
  event: 'removeGroupItem';
  uuid: string;
}

export interface MoveGroupItem {
  event: 'moveGroupItem';
  uuid: string;
  postion: number;
}

export interface ChangeCurrentStep {
  event: 'changeCurrentStep';
  uuid: string;
}

export interface RunStep {
  event: 'runStep';
  uuid: string;
}

export interface SavePipeline {
  event: 'savePipeline';
  uuid: string;
}

export interface LoadPipeline {
  event: 'loadPipeline';
  funcCallId: string;
  parentPipelineUuid?: string;
  itemId?: ItemId;
  postion?: number;
  isReadonly?: boolean;
}

export interface RunAction {
  event: 'runAction';
  uuid: string;
}

export type ViewConfigChanges = AddGroupItem | RemoveGroupItem | MoveGroupItem | ChangeCurrentStep | RunStep | SavePipeline | LoadPipeline | RunAction;

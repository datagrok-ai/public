import {PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {ItemId, NqName} from '../data/common-types';

// view config update requests

export interface AddDynamicItem {
  event: 'addDynamicItem';
  parentUuid: string;
  itemId: ItemId;
  postion: number;
}

export interface LoadDynamicItem {
  event: 'loadDynamicItem';
  parentUuid: string;
  dbId: string;
  itemId: ItemId;
  postion: number;
  readonly?: boolean; // TODO
}

export interface SaveDynamicItem {
  event: 'saveDynamicItem';
  uuid: string;
}

export interface RemoveDynamicItem {
  event: 'removeDynamicItem';
  uuid: string;
}

export interface MoveDynamicItem {
  event: 'moveDynamicItem';
  uuid: string;
  postion: number;
}

export interface RunStep {
  event: 'runStep';
  uuid: string;
  mockResults?: Record<string, any>;
  mockDelay?: number;
}

export interface SavePipeline {
  event: 'savePipeline';
}

export interface LoadPipeline {
  event: 'loadPipeline';
  funcCallId: string;
  config?: PipelineConfigurationProcessed;
}

export interface InitPipeline {
  event: 'initPipeline';
  provider: NqName;
  version?: number;
}

export type ViewConfigCommands = AddDynamicItem | LoadDynamicItem | SaveDynamicItem | RemoveDynamicItem | MoveDynamicItem | RunStep | SavePipeline | LoadPipeline | InitPipeline;

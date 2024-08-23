import {ActionItem, ValidationResult} from '../../shared-utils/validation';
import {PipelineState} from './config/PipelineInstance';
import {ItemId, RestrictionType} from './data/common-types';


export interface IRuntimeLinkController {
  // state manipulation
  getState<T = any>(name: string): T[] | undefined;
  setState<T = any>(name: string, state: T, restricion?: RestrictionType): void;
  updateRestrictions(name: string, restricion?: RestrictionType): void;
  updateOptions(name: string, options: string[], selected: string): void;
}

export interface IRuntimeValidatorController {
  getState<T = any>(name: string): T[] | undefined;
  getValidationAction(name: string, action: string): ActionItem;
  setValidation(name: string, validation?: ValidationResult | undefined): void;
}

export interface IRuntimePipelineController {
  // pipeline manipulation
  getDynamicItems(name: string): PipelineState;
  addDynamicItem(name: string, type: ItemId, position: number): void;
  moveDynamicItem(name: string, item: ItemId, position: number): void;
  removeDynamicItem(name: string, item: ItemId): void;
}

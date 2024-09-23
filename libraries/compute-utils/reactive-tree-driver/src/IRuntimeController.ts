import {ActionItem, ValidationResult} from '../../shared-utils/validation';
import {PipelineState} from './config/PipelineInstance';
import {ItemId, RestrictionType} from './data/common-types';


export interface IRuntimeController {
  // state manipulation
  getState<T = any>(name: string | number): T[] | undefined;
  setState<T = any>(name: string | number, state: T, restricion?: RestrictionType): void;
  updateRestrictions(name: string | number, restricion?: RestrictionType): void;
  getValidationAction(name: string | number, action: string): ActionItem;
  setValidation(name: string | number, validation?: ValidationResult | undefined): void;
  // pipeline manipulation
  getDynamicItems(name: string): PipelineState;
  addDynamicItem(name: string, type: ItemId, position: number): void;
  moveDynamicItem(name: string, item: ItemId, position: number): void;
  removeDynamicItem(name: string, item: ItemId): void;
}

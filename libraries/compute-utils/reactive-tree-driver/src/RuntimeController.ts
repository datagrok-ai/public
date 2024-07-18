import {ActionItem, ValidationResult} from '../../shared-utils/validation';
import {ItemPath, InputState, ItemType, GroupState, ItemId} from './config/CommonTypes';


export interface RuntimeController {
  // individual links
  enableLink(path: ItemPath): void;
  disableLink(path: ItemPath): void;
  triggerLink(path: ItemPath): void;
  isLinkEnabled(path: ItemPath): void;
  // individual steps
  enableStep(path: ItemPath): void;
  disableStep(path: ItemPath): void;
  isStepEnabled(path: ItemPath): boolean;
  // pipeline level
  enablePipeline(path: ItemPath): void;
  disablePipeline(path: ItemPath): void;
  loadNestedPipeline(path: ItemPath, runId: string): void;
  // individual state manipulation
  getState<T = any>(path: ItemPath): T | undefined;
  setState<T = any>(path: ItemPath, state: T, inputState?: InputState): void;
  updateStateInputState(path: ItemPath, inputState?: InputState): void;
  getValidationAction(path: ItemPath, name?: string): ActionItem;
  setValidation(path: ItemPath, validation?: ValidationResult | undefined): void;
  // dynamic groups state manipulation
  getGroupStates<T = any>(path: ItemPath, type?: ItemType): GroupState<T>;
  setGroupStates<T = any>(path: ItemPath, type: ItemType, state: T): void;
  // dynamic groups manipulation
  addGroupItem(path: ItemPath, type: ItemType, insertBefore?: ItemId): ItemId;
  removeGroupItem(path: ItemPath, item: ItemId): void;
  moveGroupItem(path: ItemPath, item: ItemId, insertBefore?: ItemId): void;
  getGroupConfig(path: ItemPath): any;
}

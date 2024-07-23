import {ActionItem, ValidationResult} from '../../shared-utils/validation';
import {ItemPath, InputState, ItemType, ItemId} from './data/common-types';

export interface IRuntimeController {
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
  loadPipelineRun(path: ItemPath, runId: string): void;
  // individual state manipulation
  getState<T = any>(path: ItemPath): T | undefined;
  setState<T = any>(path: ItemPath, state: T, inputState?: InputState): void;
  updateConsistency(path: ItemPath, inputState?: InputState): void;
  getValidationAction(path: ItemPath, name?: string): ActionItem;
  setValidation(path: ItemPath, validation?: ValidationResult | undefined): void;
  // dynamic pipeline manipulation
  getGroupSteps(path: ItemPath): any;
  addGroupStep(path: ItemPath, type: ItemType, insertBefore?: ItemId): ItemId;
  moveGroupStep(path: ItemPath, item: ItemId, insertBefore?: ItemId): void;
  deleteGroupStep(path: ItemPath, item: ItemId): void;
}

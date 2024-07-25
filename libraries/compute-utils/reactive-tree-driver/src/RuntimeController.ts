import {ActionItem, ValidationResult} from '../../shared-utils/validation';
import {ItemPathArray, InputState, ItemId} from './data/common-types';

export interface IRuntimeController {
  // individual links
  enableLink(path: ItemPathArray): void;
  disableLink(path: ItemPathArray): void;
  triggerLink(path: ItemPathArray): void;
  isLinkEnabled(path: ItemPathArray): void;
  // individual steps
  enableStep(path: ItemPathArray): void;
  disableStep(path: ItemPathArray): void;
  isStepEnabled(path: ItemPathArray): boolean;
  // pipeline level
  enablePipeline(path: ItemPathArray): void;
  disablePipeline(path: ItemPathArray): void;
  loadPipelineRun(path: ItemPathArray, runId: string): void;
  // individual state manipulation
  getState<T = any>(path: ItemPathArray): T | undefined;
  setState<T = any>(path: ItemPathArray, state: T, inputState?: InputState): void;
  updateConsistency(path: ItemPathArray, inputState?: InputState): void;
  getValidationAction(path: ItemPathArray, name?: string): ActionItem;
  setValidation(path: ItemPathArray, validation?: ValidationResult | undefined): void;
  // dynamic pipeline manipulation
  getGroupSteps(path: ItemPathArray): any;
  addGroupStep(path: ItemPathArray, type: ItemId, insertBefore?: ItemId): ItemId;
  moveGroupStep(path: ItemPathArray, item: ItemId, insertBefore?: ItemId): void;
  deleteGroupStep(path: ItemPathArray, item: ItemId): void;
}

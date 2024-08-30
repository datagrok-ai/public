import {ActionItem, ValidationResult} from '../../shared-utils/validation';
import {PipelineState} from './config/PipelineInstance';
import {ItemId, RestrictionType} from './data/common-types';


export interface IRuntimeLinkController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setAll<T = any>(name: string, state: T, restriction?: RestrictionType): void;
  // updateOptions(name: string, options: string[], selected: string): void;
}

export interface IRuntimeValidatorController {
  // getValidationAction(name: string, action: string): ActionItem | undefined;
  setValidation(name: string, validation?: ValidationResult | undefined): void;
}

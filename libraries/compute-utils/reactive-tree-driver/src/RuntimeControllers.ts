import {ActionItem, ValidationResult} from '../../shared-utils/validation';
import {PipelineState} from './config/PipelineInstance';
import {ItemId, RestrictionType} from './data/common-types';


export interface IRuntimeLinkController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setAll<T = any>(name: string, state: T, restriction?: RestrictionType): void;
}

export interface IRuntimeValidatorController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setValidation(name: string, validation?: ValidationResult | undefined): void;
  // getValidationAction(name: string, action: string): ActionItem | undefined;
}

export interface IRuntimeMetaController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setViewMeta(name: string, meta: any): void;
}

import {PipelineInstanceConfig} from './config/PipelineInstance';
import {RestrictionType, ValidationResult} from './data/common-types';

export interface IRuntimeLinkController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setAll<T = any>(name: string, state: T, restriction?: RestrictionType): void;
}

export interface IRuntimeValidatorController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setValidation(name: string, validation?: ValidationResult | undefined): void;
  getValidationAction(name: string, actionId: string): string | undefined;
}

export interface IRuntimeMetaController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setViewMeta(name: string, meta: any): void;
}

export interface IRuntimePipelineMutationController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setPipelineState(name: any, state: PipelineInstanceConfig): void;
}

export interface INameSelectorController {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  setDescriptionItem(name: string, description: any): void;
}

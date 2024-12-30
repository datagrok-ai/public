import {PipelineInstanceConfig} from './config/PipelineInstance';
import {RestrictionType, ValidationResult} from './data/common-types';
import * as DG from 'datagrok-api/dg';

export interface IControllerBase {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  getAdditionalParam(name: string): any | undefined;
}

export interface IRuntimeLinkController extends IControllerBase {
  setAll<T = any>(name: string, state: T, restriction?: RestrictionType): void;
}

export interface IRuntimeValidatorController extends IControllerBase {
  setValidation(name: string, validation?: ValidationResult | undefined): void;
  getValidationAction(name: string, actionId: string): string | undefined;
}

export interface IRuntimeMetaController extends IControllerBase {
  setViewMeta(name: string, meta: any): void;
}

export interface IRuntimePipelineMutationController extends IControllerBase {
  setPipelineState(name: any, state: PipelineInstanceConfig): void;
}

export interface INameSelectorController extends IControllerBase {
  setDescriptionItem(name: string, description: any): void;
}

export interface IFuncallActionController extends IControllerBase {
  setFuncCall(name: string, state: DG.FuncCall): void;
}

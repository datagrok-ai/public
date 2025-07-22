import {PipelineInstanceConfig} from './config/PipelineInstance';
import {RestrictionType, ValidationResult} from './data/common-types';
import * as DG from 'datagrok-api/dg';

export interface IControllerBase {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  getMatchedInputs(): Readonly<Set<string>>;
  getMatchedOutputs(): Readonly<Set<string>>;
  getAdditionalParam(name: string): any | undefined;
}

export interface IRuntimeLinkController extends IControllerBase {
  setAll<T = any>(name: string, state: T, restriction?: RestrictionType): void;
}

export interface IRuntimeReturnController extends IControllerBase {
  returnResult<T = any>(data: T): void;
}

export interface IRuntimeValidatorController extends IControllerBase {
  setValidation(name: string, validation?: ValidationResult | undefined): void;
  getValidationAction(id: string, actionId: string): string | undefined;
}

export interface IRuntimeMetaController extends IControllerBase {
  setViewMeta(name: string, meta: Record<string, any>): void;
}

export interface INameSelectorController extends IControllerBase {
  setDescriptionItem(name: string, description: any): void;
}

export interface IRuntimePipelineMutationController extends IControllerBase {
  setPipelineState(name: any, state: PipelineInstanceConfig): void;
}

export interface IFuncallActionController extends IControllerBase {
  setFuncCall(name: string, state: DG.FuncCall): void;
}

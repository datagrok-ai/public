import {PipelineInstanceConfigInput, PipelineOutline} from './config/PipelineInstance';
import {RestrictionType, StepHandle, ValidationResult} from './data/common-types';
import {NodePath} from './data/BaseTree';
import * as DG from 'datagrok-api/dg';

export type MatchedNodeInfo = {
  /** Node path relative to the node the link is defined on. */
  path: Readonly<NodePath>;
  /** Index of the matched node among its parent's children,
   *  -1 when the match is the node the link is defined on. */
  position: number;
  ioName?: string;
};

export interface IControllerBase {
  getAll<T = any>(name: string): T[] | undefined;
  getFirst<T = any>(name: string): T | undefined;
  getMatchedInputs(): Readonly<Set<string>>;
  getMatchedOutputs(): Readonly<Set<string>>;
  /** Positions of the matched nodes for a `from`/`to` query name,
   *  in the same order as `getAll` values. */
  getMatchedPositions(name: string): MatchedNodeInfo[];
  /** Position of the base node for base-instantiated links. */
  getBasePosition(): MatchedNodeInfo | undefined;
  getAdditionalParam(name: string): any | undefined;
  /** Static per-link parameter from the link `params` config. */
  getParam(name: string): any | undefined;
  hasCall(name: string): boolean;
}

export type TemplateId = string | number;

export type TemplateInfo = {
  name: TemplateId;
  ios: {ioName: string, scriptIoId: string}[];
};

export interface IRuntimeLinkController extends IControllerBase {
  setAll<T = any>(name: string, state: T, restriction?: RestrictionType): void;
  clearRestriction(name: string): void;
  getInputTemplates(): TemplateInfo[];
  getOutputTemplates(): TemplateInfo[];
  propagateTemplatePair(
    inputTemplate: TemplateId,
    outputTemplate: TemplateId,
    defaultRestrictions?: Record<string, RestrictionType> | RestrictionType,
  ): void;
}

export interface IRuntimeReturnController extends IControllerBase {
  returnResult<T = any>(data: T): void;
}

export interface IRuntimeValidatorController extends IControllerBase {
  setValidation(name: string, validation?: ValidationResult | undefined): void;
  getValidationAction(id: string, actionId: string): string | undefined;
  /** True when the action has no `showWhen`/`hideWhen` or its condition currently matches.
   *  Returns false if the action is not visible OR the name/actionId pair is not known.
   *  This is a pure query — `getValidationAction` is not gated by visibility, so a validator
   *  can still surface hidden actions in its result if it chooses to. */
  isActionVisible(name: string, actionId: string): boolean;
}

export interface IRuntimePipelineValidatorController extends IControllerBase {
  setValidation(name: string, validation?: ValidationResult): void;
  /** Outline of the pipeline where this link is defined (not the `to` target). */
  getOutline(): PipelineOutline;
}

export interface IRuntimeMetaController extends IControllerBase {
  setViewMeta(name: string, meta: Record<string, any>): void;
}

export interface INameSelectorController extends IControllerBase {
  setDescriptionItem(name: string, description: any): void;
}

export interface IRuntimePipelineMutationController extends IControllerBase {
  setPipelineState(name: string, state: PipelineInstanceConfigInput): void;
  getSteps(name: string): StepHandle[];
  addStep(name: string, configId: string, position?: number): void;
  removeStep(name: string, step: StepHandle): void;
  moveStep(name: string, step: StepHandle, position: number): void;
}

export interface IFuncallActionController extends IControllerBase {
  setFuncCall(name: string, state: DG.FuncCall): void;
}

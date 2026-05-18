import * as DG from 'datagrok-api/dg';
import {TreeNode} from '../data/BaseTree';
import {IRuntimeLinkController, IRuntimeMetaController, IRuntimePipelineMutationController, INameSelectorController, IRuntimeValidatorController, IFuncallActionController, IRuntimeReturnController, IRuntimePipelineValidatorController} from '../RuntimeControllers';
import {GranularMutationOp, RestrictionType, StepHandle, ValidationResult} from '../data/common-types';
import {StateTreeNode} from './StateTreeNodes';
import {ScopeInfo} from './Link';
import {PipelineInstanceConfig, PipelineOutline} from '../config/PipelineInstance';
import {NodePath} from '../data/BaseTree';

export class ControllerCancelled extends Error { };

export interface ControllerBaseArgs {
  inputs: Record<string, any[]>;
  inputsSet: Set<string>;
  outputsSet: Set<string>;
  callInputs: Set<string>;
  id: string;
  scopeInfo?: ScopeInfo;
}

export interface ValidatorControllerArgs extends ControllerBaseArgs {
  actions: Record<string, Map<string, string>>;
  baseNode?: TreeNode<StateTreeNode>;
}

export interface PipelineValidatorControllerArgs extends ControllerBaseArgs {
  outline: PipelineOutline;
}

export interface MutationControllerArgs extends ControllerBaseArgs {
  outputNodes: Record<string, {node: TreeNode<StateTreeNode>, path: NodePath}[]>;
}

export class ControllerBase<T> {
  private isActive = true;

  public outputs: Record<string, T> = {};

  public inputs: Record<string, any[]>;
  public inputsSet: Set<string>;
  public outputsSet: Set<string>;
  public callInputs: Set<string>;
  public id: string;
  public scopeInfo?: ScopeInfo;

  constructor(args: ControllerBaseArgs) {
    this.inputs = args.inputs;
    this.inputsSet = args.inputsSet;
    this.outputsSet = args.outputsSet;
    this.callInputs = args.callInputs;
    this.id = args.id;
    this.scopeInfo = args.scopeInfo;
  }

  getAll<T = any>(name: string): T[] {
    this.checkIsClosed();
    this.checkInput(name);
    return this.inputs[name];
  }

  getFirst<T = any>(name: string): T {
    return this.getAll<T>(name)?.[0];
  }

  hasCall(name: string): boolean {
    this.checkIsClosed();
    if (!this.callInputs.has(name))
      throw new Error(`Handler for Link ${this.id} is trying to check an unknown call input ${name}`);
    return this.inputsSet.has(name);
  }

  protected checkInput(name: string) {
    if (!this.inputsSet.has(name))
      throw new Error(`Handler for Link ${this.id} is trying to get an unknown input ${name}`);
  }

  protected checkOutput(name: string) {
    if (!this.outputsSet.has(name))
      throw new Error(`Handler for Link ${this.id} is trying to set an unknown output ${name}`);
  }

  protected checkIsClosed() {
    if (!this.isActive)
      throw new ControllerCancelled();
  }

  getAdditionalParam(name: string) {
    return this.scopeInfo?.additionalParams?.[name];
  }

  getMatchedInputs() {
    return this.inputsSet;
  }

  getMatchedOutputs() {
    return this.outputsSet;
  }

  close() {
    this.isActive = false;
  }
}

export class LinkController extends ControllerBase<[any, RestrictionType]> implements IRuntimeLinkController {
  setAll<T = any>(name: string, state: T, restriction: RestrictionType = 'restricted') {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = [state, restriction] as const;
  }
}

export class ValidatorController extends ControllerBase<ValidationResult | undefined> implements IRuntimeValidatorController {
  public actions: Record<string, Map<string, string>>;
  public baseNode?: TreeNode<StateTreeNode>;

  constructor(args: ValidatorControllerArgs) {
    super(args);
    this.actions = args.actions;
    this.baseNode = args.baseNode;
  }

  getValidationAction(name: string, actionId: string): string | undefined {
    this.checkIsClosed();
    const actions = this.actions[name];
    const actionUUID = actions?.get(actionId);
    return actionUUID;
  }

  setValidation(name: string, validation?: ValidationResult | undefined) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = validation;
  }
}

export class PipelineValidatorController extends ControllerBase<ValidationResult | undefined> implements IRuntimePipelineValidatorController {
  public output: ValidationResult | undefined;
  public outline: PipelineOutline;

  constructor(args: PipelineValidatorControllerArgs) {
    super(args);
    this.outline = args.outline;
  }

  setValidation(validation?: ValidationResult) {
    this.checkIsClosed();
    this.output = validation;
  }

  getOutline(): PipelineOutline {
    this.checkIsClosed();
    return this.outline;
  }
}

export class MetaController extends ControllerBase<any | undefined> implements IRuntimeMetaController {
  setViewMeta(name: string, meta?: any | undefined) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = meta;
  }
}

export class MutationController extends ControllerBase<PipelineInstanceConfig | undefined> implements IRuntimePipelineMutationController {
  public granularOps: Record<string, GranularMutationOp[]> = {};
  private removedUuids = new Set<string>();
  private usedMode: Record<string, 'replace' | 'granular'> = {};
  public outputNodes: Record<string, {node: TreeNode<StateTreeNode>, path: NodePath}[]>;

  constructor(args: MutationControllerArgs) {
    super(args);
    this.outputNodes = args.outputNodes;
  }

  private checkExclusivity(name: string, mode: 'replace' | 'granular') {
    const current = this.usedMode[name];
    if (current && current !== mode) {
      throw new Error(
        `Handler for action ${this.id}: cannot mix setPipelineState and granular ops (addStep/removeStep/moveStep) on the same output "${name}"`,
      );
    }
    this.usedMode[name] = mode;
  }

  setPipelineState(name: string, state?: PipelineInstanceConfig) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.checkExclusivity(name, 'replace');
    this.outputs[name] = state;
  }

  getSteps(name: string): StepHandle[] {
    this.checkIsClosed();
    this.checkOutput(name);
    const nodes = this.outputNodes[name];
    if (!nodes?.length)
      return [];
    const result: StepHandle[] = [];
    for (const {node} of nodes) {
      const children = node.getChildren();
      for (let i = 0; i < children.length; i++) {
        const child = children[i];
        result.push({
          configId: child.id,
          position: i,
          _uuid: child.item.getItem().uuid,
        });
      }
    }
    return result;
  }

  addStep(name: string, configId: string, position?: number) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.checkExclusivity(name, 'granular');
    (this.granularOps[name] ??= []).push({op: 'add', configId, position});
  }

  removeStep(name: string, step: StepHandle) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.checkExclusivity(name, 'granular');
    if (this.removedUuids.has(step._uuid))
      throw new Error(`Handler for action ${this.id}: step handle (configId="${step.configId}") was already removed — stale handle`);
    this.removedUuids.add(step._uuid);
    (this.granularOps[name] ??= []).push({op: 'remove', _uuid: step._uuid});
  }

  moveStep(name: string, step: StepHandle, position: number) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.checkExclusivity(name, 'granular');
    if (this.removedUuids.has(step._uuid))
      throw new Error(`Handler for action ${this.id}: step handle (configId="${step.configId}") was already removed — stale handle`);
    (this.granularOps[name] ??= []).push({op: 'move', _uuid: step._uuid, position});
  }
}

export class NodeMetaController extends ControllerBase<any | undefined> implements INameSelectorController {
  setDescriptionItem(name: string, val: string) {
    this.checkIsClosed();
    this.outputs[name] = val;
  }
}

export class RuntimeReturnController extends ControllerBase<any | undefined> implements IRuntimeReturnController {
  public result: any;

  returnResult(result: any) {
    this.checkIsClosed();
    this.result = result;
  }
}

export class FuncallActionController extends ControllerBase<any | undefined> implements IFuncallActionController {
  setFuncCall(name: string, state: DG.FuncCall) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = state;
  }
}

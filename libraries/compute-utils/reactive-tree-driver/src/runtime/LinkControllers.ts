import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TreeNode} from '../data/BaseTree';
import {IRuntimeLinkController, IRuntimeMetaController, IRuntimeValidatorController} from '../RuntimeControllers';
import {RestrictionType, ValidationResult} from '../data/common-types';
import {StateTreeNode} from './StateTreeNodes';
import {Action, ScopeInfo} from './Link';

export class ControllerCancelled extends Error { };

export class ControllerBase<T> {
  private isClosed = true;

  public outputs: Record<string, T> = {};

  constructor(
    public inputs: Record<string, any[]>,
    public inputsSet: Set<string>,
    public outputsSet: Set<string>,
    public id: string,
    public scopeInfo?: ScopeInfo,
  ) {}

  protected checkInput(name: string) {
    if (!this.inputsSet.has(name)) {
      const e = new Error(`Handler for Link ${this.id} is trying to set an unknown input ${name}`);
      grok.shell.error(e.message);
      throw e;
    }
  }

  protected checkOutput(name: string) {
    if (!this.outputsSet.has(name)) {
      const e = new Error(`Handler for Link ${this.id} is trying to set an unknown output ${name}`);
      grok.shell.error(e.message);
      throw e;
    }
  }

  protected checkIsClosed() {
    if (!this.isClosed)
      throw new ControllerCancelled();
  }

  close() {
    this.isClosed = false;
  }
}

export class LinkController extends ControllerBase<[any, RestrictionType]> implements IRuntimeLinkController {
  constructor(
    public inputs: Record<string, any[]>,
    public inputsSet: Set<string>,
    public outputsSet: Set<string>,
    public id: string,
    public scopeInfo?: ScopeInfo,
  ) {
    super(inputs, inputsSet, outputsSet, id, scopeInfo);
  }

  getAll<T = any>(name: string): T[] {
    this.checkIsClosed();
    this.checkInput(name);
    return this.inputs[name];
  }

  getFirst<T = any>(name: string): T {
    return this.getAll<T>(name)?.[0];
  }

  setAll<T = any>(name: string, state: T, restriction: RestrictionType = 'none') {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = [state, restriction] as const;
  }
}

export class ValidatorController extends ControllerBase<ValidationResult | undefined> implements IRuntimeValidatorController {
  constructor(
    public inputs: Record<string, any[]>,
    public inputsSet: Set<string>,
    public outputsSet: Set<string>,
    public id: string,
    public actions: Action[],
    public baseNode?: TreeNode<StateTreeNode>,
    public scopeInfo?: ScopeInfo,
  ) {
    super(inputs, inputsSet, outputsSet, id, scopeInfo);
  }

  getAll<T = any>(name: string): T[] {
    this.checkIsClosed();
    this.checkInput(name);
    return this.inputs[name];
  }

  getFirst<T = any>(name: string): T {
    return this.getAll<T>(name)?.[0];
  }

  getValidationAction(actionId: string): string | undefined {
    this.checkIsClosed();
    if (!this.baseNode)
      return undefined;
    const action = this.actions.find((a) => a.matchInfo.spec.id === actionId);
    return action?.uuid;
  }

  setValidation(name: string, validation?: ValidationResult | undefined) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = validation;
  }
}

export class MetaController extends ControllerBase<any | undefined> implements IRuntimeMetaController {
  constructor(
    public inputs: Record<string, any[]>,
    public inputsSet: Set<string>,
    public outputsSet: Set<string>,
    public id: string,
    public scopeInfo?: ScopeInfo,
  ) {
    super(inputs, inputsSet, outputsSet, id, scopeInfo);
  }

  getAll<T = any>(name: string): T[] {
    this.checkIsClosed();
    this.checkInput(name);
    return this.inputs[name];
  }

  getFirst<T = any>(name: string): T {
    return this.getAll<T>(name)?.[0];
  }

  setViewMeta(name: string, meta?: any | undefined) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = meta;
  }
}

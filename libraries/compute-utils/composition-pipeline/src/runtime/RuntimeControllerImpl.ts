import * as grok from 'datagrok-api/grok';
import {RichFunctionView} from '../../../function-views';
import {ActionItem, ValidationResult} from '../../../shared-utils/validation';
import {ItemPath, RuntimeController} from '../PipelineConfiguration';
import {keyToPath, pathToKey, pathJoin} from '../config-processing-utils';
import {Aborted} from './NodeConf';
import {ControllerConfig} from './ControllerConfig';
import {PipelineRuntime} from './PipelineRuntime';
import {debuglog} from '../utils';

export class RuntimeControllerImpl implements RuntimeController {
  constructor(private handlerId: string, private config: ControllerConfig, private rt: PipelineRuntime, private signal?: AbortSignal) {
  }

  public enableLink(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setLinkState(fullPath, true);
  }

  public disableLink(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setLinkState(fullPath, false);
  }

  public isLinkEnabled(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setLinkState(fullPath, false);
  }

  public enableStep(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath))
      return;

    debuglog(`handler: ${this.handlerId}, try enable step: ${pathToKey(fullPath)}`);

    return this.rt.setStepState(fullPath, true);
  }

  public disableStep(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath))
      return;

    debuglog(`handler: ${this.handlerId}, try enable step: ${pathToKey(fullPath)}`);

    return this.rt.setStepState(fullPath, false);
  }

  public isStepEnabled(path: ItemPath): boolean {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.getStepState(fullPath);
  }

  public enablePipeline(path: ItemPath) {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath))
      return;

    debuglog(`handler: ${this.handlerId}, try enable pipeline: ${pathToKey(fullPath)}`);

    return this.rt.setPipelineState(fullPath, true);
  }

  public disablePipeline(path: ItemPath) {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath))
      return;

    debuglog(`handler: ${this.handlerId}, try disable pipeline: ${pathToKey(fullPath)}`);

    return this.rt.setPipelineState(fullPath, false);
  }

  public triggerLink(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath))
      return;


    return this.rt.triggerLink(fullPath);
  }

  public getState<T = any>(path: ItemPath): T | void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkInput(fullPath))
      return;

    return this.rt.getState(fullPath);
  }

  public setState<T>(path: ItemPath, value: T, inputState?: 'disabled' | 'restricted' | 'user input'): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath))
      return;

    debuglog(`handler: ${this.handlerId}, try update state: ${pathToKey(fullPath)}, value: ${value}`);
    return this.rt.setState(fullPath, value, inputState);
  }

  public updateStateState(path: ItemPath, inputState: 'disabled' | 'restricted' | 'user input'): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath))
      return;

    this.rt.setState(fullPath, this.rt.getState(fullPath), inputState);
  }

  public getView(path: ItemPath): RichFunctionView | void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.getView(fullPath);
  }

  public goToStep(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.goToStep(fullPath);
  }

  public setValidation(path: ItemPath, validation?: ValidationResult): void {
    this.checkAborted();
    const fullPathTarget = pathJoin(this.config.pipelinePath, path);
    const fullPathLinkPath = keyToPath(this.handlerId);
    if (!this.checkValidation(fullPathTarget))
      return;

    debuglog(`handler: ${this.handlerId}: update validation: ${fullPathTarget}, value: ${validation}`);
    this.rt.setValidation(fullPathTarget, fullPathLinkPath, validation);
  }

  public getValidationAction(path: ItemPath, name?: string): ActionItem {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.getAction(fullPath, name);
  }

  public loadNestedPipeline(path: ItemPath, runId: string) {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    this.rt.loadNestedPipeline(fullPath, runId);
  }

  private checkAborted() {
    if (this.signal?.aborted)
      throw new Aborted();
  }

  private checkInput(path: ItemPath) {
    const key = pathToKey(path);
    if (!this.config.fromKeys.has(key))
      grok.shell.warning(`Handler ${this.handlerId} has not declared access to input ${key}`);

    return true;
  }

  private checkOutput(path: ItemPath) {
    const key = pathToKey(path);
    if (!this.config.toKeys.has(key))
      grok.shell.warning(`Handler ${this.handlerId} has not declared access to output ${key}`);

    if (this.config.fromKeys.has(key)) {
      grok.shell.error(`Handler ${this.handlerId} trying to set output value which is input ${key}`);
      return false;
    }

    return true;
  }

  private checkValidation(path: ItemPath) {
    const key = pathToKey(path);
    if (!this.config.toKeys.has(key) || !this.config.fromKeys.has(key))
      grok.shell.warning(`Handler ${this.handlerId} should declare both input and output access to set validation ${key}`);

    return true;
  }
}

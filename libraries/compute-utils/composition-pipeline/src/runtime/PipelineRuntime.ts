import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, Observable, of, from, Subject} from 'rxjs';
import {filter, switchMap, debounceTime, takeUntil, take, finalize, mapTo, concatMap, startWith} from 'rxjs/operators';
import {ActionItem, ValidationResult, mergeValidationResults} from '../../../shared-utils/validation';
import {historyUtils} from '../../../history-utils';
import {PipelineStepConfiguration, PipelineActionConfiguraion, Handler} from '../config/PipelineConfiguration';
import {ItemPath, InputState} from '../config/CommonTypes';
import {keyToPath, pathToKey, PathKey, getSuffix, pathJoin} from '../config/config-processing-utils';
import {NodeConfTypes, Aborted} from './NodeConf';
import {NodeState} from './NodeState';
import {LinkState} from './LinkState';
import {debuglog, callHandler} from '../utils';
import {RuntimeControllerImpl} from './RuntimeControllerImpl';
import {PipelineDriverState} from './PipelineDriverState';
import {IViewBridge} from '../view/IViewBridge';

export class PipelineRuntime {
  private runningLinks = new Set<string>();
  private disabledSettersPrefixes = new Set<string>();
  private validationState = new Map<string, Map<string, ValidationResult | undefined>>();
  private pipelineLoads = new Subject<{ path: ItemPath; runId: string; }>();
  public isUpdating = new BehaviorSubject<boolean>(false);

  constructor(
    private nodes: Map<PathKey, NodeState>,
    private links: Map<PathKey, LinkState>,
    private viewBridge: IViewBridge,
    public pipelineState: PipelineDriverState,
  ) {
    this.pipelineLoads.pipe(
      concatMap(({path, runId}) => from(this.pipelineLoader(path, runId))),
      takeUntil(this.pipelineState.closed),
    ).subscribe();
  }

  public setLinkState(path: ItemPath, enabled: boolean): void {
    const k = pathToKey(path);
    const link = this.links.get(k);
    if (link)
      link.enabled.next(enabled);
  }

  public getLinkState(path: ItemPath): boolean {
    const k = pathToKey(path);
    const link = this.links.get(k);
    return !!link?.enabled.value;
  }

  public loadNestedPipeline(path: ItemPath, runId: string) {
    this.pipelineLoads.next({path, runId});
  }

  public triggerLink(path: ItemPath): void {
    const k = pathToKey(path);
    const link = this.links.get(k);
    if (link)
      link.trigger();
  }

  public getState<T = any>(path: ItemPath): T | undefined {
    const {state} = this.getNodeState(path)!;
    const val = state?.getValue();
    return val;
  }

  public setState<T>(path: ItemPath, value: T, inputState?: InputState): void {
    const {state, node} = this.getNodeState(path)!;
    if (state && !this.isSetterDisabled(path, node.type))
      state.setValue(value, inputState);
  }

  public setValidation(targetPath: ItemPath, linkPath: ItemPath, validation: ValidationResult | undefined): void {
    const targetId = pathToKey(targetPath);
    const linkId = pathToKey(linkPath);
    const targetState = this.validationState.get(targetId) ?? new Map<string, ValidationResult | undefined>();
    targetState.set(linkId, validation);
    this.validationState.set(targetId, targetState);
    this.updateViewValidation(targetPath);
  }

  public getAction(path: ItemPath, name?: string): ActionItem {
    const k = pathToKey(path);
    const node = this.nodes.get(k);
    if (node?.type !== 'action') {
      const msg = `Node ${path.join(',')} is not an action`;
      grok.shell.error(msg);
      throw new Error(msg);
    }

    const conf = node.conf as PipelineActionConfiguraion;
    const handler = conf.handler;
    const ctrlConf = node.controllerConfig!;
    return {
      actionName: name ?? conf.friendlyName,
      action: async () => {
        try {
          const controller = new RuntimeControllerImpl(node.conf.id, ctrlConf, this!);
          await callHandler(handler, {controller}).toPromise();
        } catch (e) {
          grok.shell.error(String(e));
          throw (e);
        }
      },
    };
  }

  public wireViews() {
    for (const [, node] of this.nodes) {
      for (const [, state] of node.states) {
        if ((node.type === 'step') && (state.conf.stateType === 'input' || state.conf.stateType === 'output')) {
          const stateId = state.conf.id;
          const {changes, setter} = this.viewBridge.getStateBindings(node.conf.id, stateId);
          state.linkState(changes, setter);
        }
      }
    }
  }

  public wireLinks() {
    for (const [, link] of this.links) {
      const changes = link.controllerConfig.from.map((path) => this.getValueChanges(path, link));
      const handler = link.conf.handler ?? this.getDefaultHandler(link);
      link.getValuesChanges().pipe(
        debounceTime(0),
        switchMap(() => this.getHandlerObservable(link, handler)),
        takeUntil(this.pipelineState.closed),
      ).subscribe();
      link.setSources(changes);
    }
  }

  private getHandlerObservable(link: LinkState, handler: Handler) {
    const obs$ = new Observable((observer) => {
      this.addRunningHandler(link.conf.id);
      const abortController = new AbortController();
      const signal = abortController.signal;
      const controller = new RuntimeControllerImpl(link.conf.id, link.controllerConfig, this, signal);
      let done = false;
      const sub = from(callHandler(handler, {controller})).pipe(
        finalize(() => this.removeRunningHandler(link.conf.id)),
      ).subscribe(
        (val) => {
          done = true;
          observer.next(val);
        },
        (error: any) => {
          if (!(error instanceof Aborted)) {
            grok.shell.error(error);
            console.error(error);
          }
        });
      return () => {
        if (!done)
          abortController.abort();
        sub.unsubscribe();
      };
    });
    return obs$;
  }

  private getValueChanges(path: ItemPath, link: LinkState) {
    const {state} = this.getNodeState(path)!;
    const valueChanges = link.conf.includeDataFrameMutations ?
      state.value.pipe(switchMap((param) => {
        if (param?.onDataChanged && param?.onDataChanged?.pipe)
          return param.onDataChanged.pipe(startWith(null), mapTo(param)) as Observable<any>;

        return of(param);
      })) : state.value;
    return valueChanges;
  }

  private getDefaultHandler(link: LinkState) {
    return (async () => {
      const nLinks = Math.min(link.controllerConfig.from.length, link.controllerConfig.to.length);
      for (let idx = 0; idx < nLinks; idx++) {
        let state = this.getState(link.controllerConfig.from[idx]);
        if (state instanceof DG.DataFrame)
          state = state.clone();
        this.setState(link.controllerConfig.to[idx], state, link.conf.inputState);
      }
    });
  }

  public getRunningUpdates() {
    return [...this.runningLinks];
  }

  public disableSubStepsIOSetters(prefixPath: ItemPath) {
    const k = pathToKey(prefixPath);
    debuglog(`disable substeps setters: ${k}`);
    this.disabledSettersPrefixes.add(k);
  }

  public enableSubStepsIOSetters(prefixPath: ItemPath) {
    const k = pathToKey(prefixPath);
    debuglog(`enable substeps setters: ${k}`);
    this.disabledSettersPrefixes.delete(k);
  }

  public async awaitForAllUpdatesDone() {
    await this.isUpdating.pipe(
      debounceTime(1), // wait more than state debounce time
      filter((x) => !x),
      take(1),
    ).toPromise();
  }

  private async pipelineLoader(path: ItemPath, runId: string) {
    this.disableSubStepsIOSetters(path);
    try {
      const calls = await this.loadPipelineFuncCalls(runId);
      const nestedSteps = this.getNestedSteps(path);
      await this.loadNestedSteps(calls, nestedSteps);
      await this.awaitForAllUpdatesDone();
    } finally {
      this.enableSubStepsIOSetters(path);
    }
  }

  private getNestedSteps(pipelinePath: ItemPath) {
    const pipelineKey = pathToKey(pipelinePath);
    const pipelineName = pipelinePath[pipelinePath.length - 1];
    const nestedSteps = ([...this.nodes.entries()])
      .filter(([key, node]) => node.type === 'step' && getSuffix(key, pipelineKey))
      .map(([key, node]) => {
        const pathKey = key;
        const suffix = getSuffix(key, pipelineKey)!;
        const nestedStepId = pathToKey(pathJoin([pipelineName], keyToPath(suffix)));
        const nqName = (node.conf as PipelineStepConfiguration).nqName;
        return {pathKey, nestedStepId, nqName};
      });
    return nestedSteps;
  }

  private async loadPipelineFuncCalls(funcCallId: string) {
    const runs = await historyUtils.loadChildRuns(funcCallId);
    return runs.childRuns;
  }

  private addRunningHandler(id: string) {
    this.runningLinks.add(id);
    this.isUpdating.next(true);
  }

  private removeRunningHandler(id: string) {
    this.runningLinks.delete(id);
    this.isUpdating.next(this.runningLinks.size !== 0);
  }

  private isSetterDisabled(path: ItemPath, type: NodeConfTypes) {
    if (type === 'step') {
      const k = pathToKey(path);
      for (const prefix of this.disabledSettersPrefixes) {
        const suffix = getSuffix(k, prefix);
        if (suffix)
          return true;
      }
    }
    return false;
  }

  private getNodeState(path: ItemPath) {
    const nodePath = path.slice(0, path.length - 1);
    const stateName = path[path.length - 1];
    const k = pathToKey(nodePath);
    const node = this.nodes.get(k);
    if (node && node.states) {
      const state = node.states.get(stateName);
      if (state)
        return {node, state};
    }
  }

  private updateViewValidation(targetPath: ItemPath) {
    const targetId = pathToKey(targetPath);
    const targetState = this.validationState.get(targetId)!;
    const keys = [...targetState.keys()].sort();
    const validations = keys.map((k) => targetState.get(k));
    const {state, node} = this.getNodeState(targetPath)!;
    this.viewBridge.setExternalValidationResults(node.conf.id, state.conf.id, mergeValidationResults(...validations));
  }
}

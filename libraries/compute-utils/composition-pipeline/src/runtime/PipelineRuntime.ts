import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, Observable, of, from, Subject} from 'rxjs';
import {RichFunctionView} from '../../../function-views';
import {filter, switchMap, debounceTime, takeUntil, take, sample, finalize, mapTo, concatMap, startWith} from 'rxjs/operators';
import {ActionItem, ValidationResult, mergeValidationResults} from '../../../shared-utils/validation';
import {historyUtils} from '../../../history-utils';
import {ABILITY_STATE, VISIBILITY_STATE} from '../../../shared-utils/consts';
import {isIncomplete} from '../../../shared-utils/utils';
import {ItemPath, PipelineStepConfiguration, PipelinePopupConfiguration, PipelineActionConfiguraion, InputState, Handler} from '../PipelineConfiguration';
import {keyToPath, pathToKey, PathKey, getSuffix, pathJoin} from '../config-processing-utils';
import {NodeConfTypes, Aborted} from './NodeConf';
import {NodeState} from './NodeState';
import {LinkState} from './LinkState';
import {debuglog, callHandler} from '../utils';
import {ICompositionView} from '../view/CompositionPipelineView';
import {RuntimeControllerImpl} from './RuntimeControllerImpl';
import {PipelineGlobalState} from './PipelineGlobalState';

interface NestedStepData {
  pathKey: string;
  nestedStepId: string;
  nqName: string;
}

export class PipelineRuntime {
  private runningLinks = new Set<string>();
  private disabledSettersPrefixes = new Set<string>();
  private validationState = new Map<string, Map<string, ValidationResult | undefined>>();
  private pipelineLoads = new Subject<{ path: ItemPath; runId: string; }>();
  public isUpdating = new BehaviorSubject<boolean>(false);

  constructor(
    private nodes: Map<PathKey, NodeState>,
    private links: Map<PathKey, LinkState>,
    private view: ICompositionView,
    public pipelineState: PipelineGlobalState,
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

  public setStepState(path: ItemPath, enabled: boolean): void {
    const k = pathToKey(path);

    if (this.isSetterDisabled(path, 'step'))
      return;

    debuglog(`step visibility updated: ${k}, new value: ${enabled}`);

    if (enabled)
      this.view.showSteps(k);

    else
      this.view.hideSteps(k);
  }

  public getStepState(path: ItemPath): boolean {
    const k = pathToKey(path);
    const value = this.view.getStepState(k);
    return value === VISIBILITY_STATE.VISIBLE;
  }

  public setPipelineState(path: ItemPath, enabled: boolean): void {
    const steps = this.getNestedSteps(path);
    for (const data of steps)
      this.setStepState(keyToPath(data.pathKey), enabled);
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

  public getState<T = any>(path: ItemPath): T | void {
    const {state} = this.getNodeState(path)!;
    const val = state?.getValue();
    return val;
  }

  public setState<T>(path: ItemPath, value: T, inputState?: InputState): void {
    const {state, node} = this.getNodeState(path)!;
    if (state && !this.isSetterDisabled(path, node.type))
      state.setValue(value, inputState);
  }

  public getView(path: ItemPath): RichFunctionView | void {
    const k = pathToKey(path);
    return this.view.getStepView<RichFunctionView>(k);
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
    if (node?.type !== 'action' && node?.type !== 'popup') {
      const msg = `Node ${path.join(',')} is not an action/popup`;
      grok.shell.error(msg);
      throw new Error(msg);
    }

    if (node.type === 'action') {
      const conf = node.conf as PipelineActionConfiguraion;
      const handler = conf.handler;
      const ctrlConf = node.controllerConfig!;
      return {
        actionName: name ?? conf.friendlyName,
        action: async () => {
          try {
            const controller = new RuntimeControllerImpl(node.conf.id, ctrlConf, this!);
            await callHandler(handler, {controller});
          } catch (e) {
            grok.shell.error(String(e));
            throw (e);
          }
        },
      };
    } else {
      const conf = node.conf as PipelinePopupConfiguration;
      const view = this.view.getStepView(conf.id);
      if (!view) {
        const msg = `Popup view ${conf.id} not found`;
        grok.shell.error(msg);
        throw new Error(msg);
      }
      return {
        actionName: name ?? conf.friendlyName,
        action: async () => {
          try {
            await new Promise((resolve, _reject) => {
              ui.dialog(conf.friendlyName)
                .add(view.root)
                .onOK(() => {
                  node.notifier.next(true);
                  resolve(true);
                })
                .onCancel(() => {
                  resolve(false);
                })
                .showModal(true);
            });
          } catch (e) {
            grok.shell.error(String(e));
            throw (e);
          }
        },
      };
    }
  }

  public wireViews() {
    for (const [, node] of this.nodes) {
      for (const [, state] of node.states) {
        if ((node.type === 'popup' || node.type === 'step') && (state.conf.stateType === 'input' || state.conf.stateType === 'output')) {
          const stateId = state.conf.id;
          const {changes, setter} = this.view.getStateBindings(node.conf.id, stateId);
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
    const notifier = state.notifier;
    if (link.conf.ignoreNotifier || !notifier)
      return valueChanges;
    return valueChanges.pipe(sample(notifier));
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

  public goToStep(_path: ItemPath): void {
    console.log('TODO: not implemented');
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
    this.setPipelineState(path, false);
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

  private async loadNestedSteps(calls: DG.FuncCall[], stepsData: NestedStepData[]) {
    const stepsMappings = new Map<string, NestedStepData>();
    for (const data of stepsData)
      stepsMappings.set(data.nestedStepId, data);

    for (const call of calls) {
      const customId = call.options['customId'];
      const data = stepsMappings.get(customId);
      if (data) {
        const rfv = this.view.getStepView(data.pathKey);
        if (rfv) {
          debuglog(`loaded step: ${data.pathKey}, nqName: ${data.nqName}, id: ${call.id}`);
          const ncall = await rfv.loadRun(call.id);
          this.view.showSteps(data.pathKey);
          if (!isIncomplete(ncall))
            this.view.enableStep(data.pathKey);
        }
      }
    }

    const nextStep = this.view.getFollowingStep();
    if (nextStep)
      nextStep.ability.next(ABILITY_STATE.ENABLED);
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
    this.view.setExternalValidationResults(node.conf.id, state.conf.id, mergeValidationResults(...validations));
  }
}

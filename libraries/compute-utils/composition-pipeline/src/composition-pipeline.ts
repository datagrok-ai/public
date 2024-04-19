import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import cloneDeepWith from 'lodash.clonedeepwith';
import {BehaviorSubject, Observable, merge, of, from, Subject} from 'rxjs';
import {PipelineView, RichFunctionView} from '../../function-views';
import {withLatestFrom, filter, switchMap, debounceTime, takeUntil, take, sample, finalize, mapTo, concatMap} from 'rxjs/operators';
import {ActionItem, ValidationResult, mergeValidationResults} from '../../shared-utils/validation';
import {historyUtils} from '../../history-utils';

//
// State values

const primitiveTypes = ['bool', 'int', 'number', 'string'] as const;
type SupportedPrimitives = typeof primitiveTypes[number];
const containerTypes = ['dataframe', 'object'] as const;
type SupportedContainers = typeof primitiveTypes[number];
const supportedTypes = [...primitiveTypes, ...containerTypes] as const;
type SupportedTypes = typeof supportedTypes[number];
export type TypeSpec = {
  [key: string]: SupportedTypes | { type: SupportedContainers, spec?: TypeSpec };
}

//
// Common interfaces

export type ItemName = string;
export type ItemPath = ItemName[];
export type NqName = string;

export interface RuntimeController {
  enableLink(path: ItemPath): void;
  disableLink(path: ItemPath): void;
  triggerLink(path: ItemPath): void;
  enableStep(path: ItemPath): void;
  disableStep(path: ItemPath): void;
  loadNestedPipeline(path: ItemPath, runId: string): void;
  getState<T = any>(path: ItemPath): T | void;
  setState<T>(path: ItemPath, state: T, inputState?: 'disabled' | 'restricted' | 'user input'): void;
  getValidationAction(path: ItemPath, name?: string): ActionItem
  setValidation(path: ItemPath, validation: ValidationResult | undefined): void
  getView(path: ItemPath): RichFunctionView | void;
  goToStep(path: ItemPath): void;
}

//
// CompositionPipeline individual RFV configurations

export type ExportConfig = {
  export?: ((format: string) => Promise<Blob>);
  filename?: ((format: string) => string);
  supportedFormats?: string[];
  supportedExtensions?: Record<string, string>;
}

export type StateType = 'input' | 'output' | 'state';

export type StateItemConfiguration = {
  id: ItemName;
  type?: TypeSpec;
  stateType?: StateType; // system field
}

export type Handler = ((params: { controller: RuntimeController }) => Promise<void>) | NqName;

export type PipelineLinkConfiguration = {
  id: ItemName;
  from: ItemPath | ItemPath[];
  to: ItemPath | ItemPath[];
  handler?: Handler;
  ignoreNotifier?: boolean;
}

export type PipelineHookConfiguration = {
  id: ItemName;
  from?: ItemPath | ItemPath[];
  to?: ItemPath | ItemPath[];
  handler: Handler;
}

export type PipelineActionConfiguraion = PipelineHookConfiguration & {
  position: 'buttons' | 'menu' | 'none';
  friendlyName: string;
}

export type PipelinePopupConfiguration = {
  id: ItemName;
  nqName: NqName;
  position: 'buttons' | 'menu' | 'none';
  friendlyName: string;
  helpUrl?: string;
  states?: StateItemConfiguration[];
  actions?: PipelineActionConfiguraion[];
}

export type PipelineStepConfiguration = {
  id: ItemName;
  nqName: NqName;
  friendlyName?: string;
  helpUrl?: string;
  states?: StateItemConfiguration[];
  actions?: PipelineActionConfiguraion[];
  popups?: PipelinePopupConfiguration[];
}

export type PipelineHooks = {
  beforeInit?: PipelineHookConfiguration[];
  afterInit?: PipelineHookConfiguration[];
  beforeFuncCallReady?: PipelineHookConfiguration[];
  afterFuncCallReady?: PipelineHookConfiguration[];
  beforeLoadRun?: PipelineHookConfiguration[];
  afterLoadRun?: PipelineHookConfiguration[];
  onViewReady?: PipelineHookConfiguration[];
}

export type PipelineConfiguration = {
  id: ItemName;
  nqName: NqName;
  // composable
  steps: PipelineStepConfiguration[];
  hooks?: PipelineHooks;
  links?: PipelineLinkConfiguration[];
  actions?: PipelineActionConfiguraion[];
  // non-composable
  states?: StateItemConfiguration[];
  exportConfig?: ExportConfig;
  helpUrl?: string;
}

//
// Composition Pipeline configuration

export type ItemsToAdd = {
  popupsToAdd?: [PipelinePopupConfiguration, ItemPath][];
  actionsToAdd?: [PipelineActionConfiguraion, ItemPath][];
}

export type ItemsToRemove = {
  itemsToRemove?: ItemPath[];
}

export type NestedPipelineConfig = {
  insertBeforeStep?: ItemName,
  namePrefix?: string;
}

export type ComposedPipelinesConfig = {
  nestedPipelinesConfig?: {
    [id: string]: NestedPipelineConfig;
  }
}

export type PipelineCompositionConfiguration = ItemsToAdd & PipelineConfiguration & ComposedPipelinesConfig & ItemsToRemove;

//
// Internal config

export type PipelineConfigVariants = PipelineConfiguration | CompositionGraphConfig;

export type GraphNestedPipelines = {
  nestedPipelines?: {
    [key: ItemName]: PipelineConfiguration | CompositionGraphConfig;
  };
}

export type CompositionGraphConfig = (PipelineConfiguration | PipelineCompositionConfiguration) & GraphNestedPipelines;

export function isCompositionConfig(config: CompositionGraphConfig | PipelineCompositionConfiguration | PipelineConfiguration): config is PipelineCompositionConfiguration {
  return (config as any)?.nestedPipelines;
}

export function isPipelineConfig(config: CompositionGraphConfig | PipelineCompositionConfiguration | PipelineConfiguration): config is PipelineConfiguration {
  return !(config as any)?.nestedPipelines;
}

export function cloneConfig<T>(config: T): T {
  return cloneDeepWith(config, (val) => {
    if (typeof val === 'function')
      return val;
  });
}

export type PathKey = string;

export function pathJoin(path: ItemPath, ...restPaths: ItemPath[]): ItemPath {
  return path.concat(...restPaths);
}

export function pathToKey(path: ItemPath): PathKey {
  return path.filter((x) => x).join('/');
}

export function keyToPath(key: PathKey) {
  return key.split('/');
}

export function getParentKey(key: PathKey): PathKey {
  const path = keyToPath(key);
  const ppath = path.slice(0, path.length - 1);
  return pathToKey(ppath);
}

export function getSuffix(key: PathKey, prefix: PathKey) {
  if (key.startsWith(prefix))
    return key.substring(prefix.length);
}

export function normalizePaths(paths: ItemPath | ItemPath[]): ItemPath[] {
  if (Array.isArray(paths[0]))
    return paths as ItemPath[];
  else
    return [paths as ItemPath];
}

export function traverseConfigPipelines<T>(
  graph: CompositionGraphConfig,
  nodeHandler: (
    acc: T,
    node: CompositionGraphConfig,
    path: ItemPath) => T,
  acc: T,
) {
  const stk: [{
    node: CompositionGraphConfig,
    path: ItemPath,
  }] = [{
    node: graph,
    path: [] as ItemPath,
  }];

  const queuedNested = new Set<object>();

  while (stk.length) {
    const {node, path} = stk.pop()!;
    if (isCompositionConfig(node)) {
      if (queuedNested.has(node))
        acc = nodeHandler(acc, node, path);
      else {
        stk.push({node: node, path});
        for (const nextNode of Object.values(node.nestedPipelines ?? {}))
          stk.push({node: nextNode, path: [...path, node.id]});
        queuedNested.add(node);
      }
    } else if (isPipelineConfig(node))
      acc = nodeHandler(acc, node, path);
  }

  return acc;
}

//
// Internal runtime state

export class Aborted extends Error {}

export type ItemsToMerge = {
  step?: PipelineStepConfiguration[];
  popup?: PipelinePopupConfiguration[];
  action?: PipelineActionConfiguraion[];
}

export type SubNodeConfTypes = 'action' | 'popup' | 'step';
export type SubNodeConf = PipelineActionConfiguraion | PipelinePopupConfiguration | PipelineStepConfiguration;
export type NodeConfTypes = SubNodeConfTypes | 'pipeline';
export type NodeConf = SubNodeConf | PipelineConfiguration;

export type InputState = 'disabled' | 'restricted' | 'user input';

export class NodeItemState<T = any> {
  public currentSource = new BehaviorSubject<Observable<T>>(of());
  private valueChanges = this.currentSource.pipe(
    switchMap((source) => source),
  );

  public value = new BehaviorSubject<T | undefined>(undefined);

  private setter = (x: T, _inputState?: InputState) => {
    this.currentSource.next(of(x));
  };

  constructor(
    public conf: StateItemConfiguration,
    public pipelineState: PipelineGlobalState,
    public notifier?: Observable<true>,
  ) {
    this.valueChanges.pipe(
      takeUntil(this.pipelineState.closed),
    ).subscribe((x) => {
      // console.log('state', this.conf.id, x);
      this.value.next(x);
    });
  }

  linkState(source: Observable<T>, setter?: (x: T, inputState?: InputState) => void) {
    this.currentSource.next(source);
    if (setter)
      this.setter = setter;
  }

  getValue() {
    return this.value.value;
  }

  setValue(val: T, inputState?: InputState) {
    this.setter(val, inputState);
  }
}

export class NodeState {
  public controllerConfig?: ControllerConfig;
  public states = new Map<ItemName, NodeItemState>();
  public notifier = new Subject<true>();

  constructor(
    public conf: NodeConf,
    public type: NodeConfTypes,
    public pipelinePath: ItemPath,
    public pipelineState: PipelineGlobalState,
  ) {
    if (type === 'pipeline') {
      const states = (conf as PipelineStepConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeItemState(state, this.pipelineState));
    }

    if (type === 'step') {
      const states = (conf as PipelineStepConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeItemState(state, this.pipelineState));
    }

    if (type === 'popup') {
      const states = (conf as PipelinePopupConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeItemState(state, this.pipelineState, this.notifier));
    }

    if (type === 'action') {
      const link = conf as PipelineActionConfiguraion;
      this.controllerConfig = new ControllerConfig(pipelinePath, link.from, link.to);
    }
  }
}


export class LinkState {
  public controllerConfig: ControllerConfig;
  public enabled = new BehaviorSubject(true);
  private currentSource = new BehaviorSubject<Observable<any>>(of());
  private externalTrigger = new Subject<true>();
  public valueChanges = this.currentSource.pipe(
    switchMap((source) => source),
  );

  constructor(
    public conf: PipelineLinkConfiguration,
    public pipelinePath: ItemPath,
    public pipelineState: PipelineGlobalState,
  ) {
    this.controllerConfig = new ControllerConfig(pipelinePath, conf.from, conf.to);
  }

  public setSources(sources: Observable<any>[]) {
    this.currentSource.next(merge(...sources, this.externalTrigger));
  }

  public trigger() {
    this.externalTrigger.next(true);
  }

  public getValuesChanges() {
    return this.valueChanges.pipe(
      withLatestFrom(this.enabled),
      filter(([, enabled]) => {
        if (!enabled)
          return false;

        return true;
      }),
      // tap(() => console.log('link', this.conf.id)),
      mapTo(true),
    );
  }
}

export class PipelineGlobalState {
  closed = new Subject<true>();
}

export interface StepData {
  path: string;
  subPath: string;
  nqName: string;
}

export class ControllerConfig {
  public from: ItemPath[] = [];
  public to: ItemPath[] = [];
  public fromKeys = new Set<string>();
  public toKeys = new Set<string>();

  constructor(
    public pipelinePath: ItemPath,
    from?: ItemPath | ItemPath[],
    to?: ItemPath | ItemPath[],
  ) {
    if (from) {
      this.from = normalizePaths(from);
      this.fromKeys = new Set(this.from.map((path) => pathToKey(path)));
    }

    if (to) {
      this.to = normalizePaths(to);
      this.toKeys = new Set(this.to.map((path) => pathToKey(path)));
    }
  }
}

export class PipelineRuntime {
  private runningLinks = new Set<string>();
  private disabledStepsPaths = new Set<string>();
  private validationState = new Map<string, Map<string, ValidationResult | undefined>>();
  private pipelineLoads = new Subject<{path: ItemPath, runId: string}>();
  public isUpdating = new BehaviorSubject<boolean>(false);

  constructor(
    private nodes: Map<PathKey, NodeState>,
    private links: Map<PathKey, LinkState>,
    private view: ICompositionView,
    public pipelineState: PipelineGlobalState,
  ) {
    this.pipelineLoads.pipe(
      concatMap(({path, runId}) => from(this.loadNestedPipeline(path, runId))),
      takeUntil(this.pipelineState.closed),
    ).subscribe();
  }

  setLinkState(path: ItemPath, enabled: boolean): void {
    const k = pathToKey(path);
    const link = this.links.get(k);
    if (link)
      link.enabled.next(enabled);
  }

  setStepState(path: ItemPath, enabled: boolean): void {
    const k = pathToKey(path);

    if (enabled)
      this.view.showSteps(k);
    else
      this.view.hideSteps(k);
  }

  async loadNestedPipeline(path: ItemPath, runId: string) {
    this.disableSubStepsIOSetters(path);
    try {
      const calls = await this.loadPipelineFuncCalls(runId);
      const nestedSteps = this.getNestedSteps(path);
      this.insertLoadedPipelineSteps(calls, nestedSteps);
      await this.awaitForAllUpdatesDone();
    } finally {
      this.enableSubStepsIOSetters(path);
    }
  }

  triggerLink(path: ItemPath): void {
    const k = pathToKey(path);
    const link = this.links.get(k);
    if (link)
      link.trigger();
  }

  getState<T = any>(path: ItemPath): T | void {
    const {state} = this.getNodeState(path)!;
    const val = state?.getValue();
    return val;
  }

  updateState<T>(path: ItemPath, value: T, inputState: InputState = 'user input'): void {
    const {state, node} = this.getNodeState(path)!;
    if (state && !this.isSetterDisabled(path, node.type))
      state.setValue(value, inputState);
  }

  getView(path: ItemPath): RichFunctionView | void {
    const k = pathToKey(path);
    return this.view.getStepView<RichFunctionView>(k);
  }

  setValidation(targetPath: ItemPath, linkPath: ItemPath, validation: ValidationResult | undefined): void {
    const targetId = pathToKey(targetPath);
    const linkId = pathToKey(linkPath);
    const targetState = this.validationState.get(targetId) ?? new Map<string, ValidationResult | undefined>();
    targetState.set(linkId, validation);
    this.validationState.set(targetId, targetState);
    this.updateViewValidation(targetPath);
  }


  getAction(path: ItemPath, name?: string): ActionItem {
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

  wireViews() {
    // wiring by nqName
    for (const [, node] of this.nodes) {
      for (const [, state] of node.states) {
        if ((node.type === 'popup' || node.type === 'step') && (state.conf.stateType === 'input' || state.conf.stateType === 'output')) {
          const stateId = state.conf.id;
          const {changes, setter} = this.view.getStateBindings(node.conf.id, stateId);
          if (node.type === 'step')
            state.linkState(changes, setter);
        }
      }
    }
  }

  wireLinks() {
    for (const [, link] of this.links) {
      const changes = link.controllerConfig.from.map((path) => {
        const {state} = this.getNodeState(path)!;
        const value = state.value;
        const notifier = state.notifier;
        if (link.conf.ignoreNotifier || !notifier)
          return value;
        return value.pipe(sample(notifier));
      });
      const handler = link.conf.handler ?? (async () => {
        // console.log('default handler', link.conf.id);
        const nLinks = Math.min(link.controllerConfig.from.length, link.controllerConfig.to.length);
        for (let idx = 0; idx < nLinks; idx++) {
          let state = this.getState(link.controllerConfig.from[idx]);
          if (state instanceof DG.DataFrame)
            state = state.clone();
          this.updateState(link.controllerConfig.to[idx], state);
        }
      });
      link.getValuesChanges().pipe(
        debounceTime(0),
        switchMap(() => {
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
        }),
        takeUntil(this.pipelineState.closed),
      ).subscribe();
      link.setSources(changes);
    }
  }

  getRunningUpdates() {
    return [...this.runningLinks];
  }

  goToStep(_path: ItemPath): void {
    console.log('TODO: goToStep not implemented');
  }

  disableSubStepsIOSetters(prefixPath: ItemPath) {
    const k = pathToKey(prefixPath);
    this.disabledStepsPaths.add(k);
  }

  enableSubStepsIOSetters(prefixPath: ItemPath) {
    const k = pathToKey(prefixPath);
    this.disabledStepsPaths.delete(k);
  }

  private insertLoadedPipelineSteps(calls: DG.FuncCall[], stepsData: StepData[]) {
    const stepsMappings = new Map<string, StepData>();
    for (const data of stepsData)
      stepsMappings.set(data.subPath, data);

    for (const call of calls) {
      const id = call.options['customId'];
      const data = stepsMappings.get(id);
      if (data) {
        const rfv = this.view.getStepView(id);
        if (rfv)
          rfv.linkFunccall(call);
      }
    }
  }

  private async awaitForAllUpdatesDone() {
    await this.isUpdating.pipe(
      debounceTime(0),
      filter((x) => !x),
      take(1),
    ).toPromise();
  }

  private getNestedSteps(pipelinePath: ItemPath) {
    const k = pathToKey(pipelinePath);
    const nestedSteps = ([...this.nodes.entries()])
      .filter(([key, node]) => node.type === 'step' && getSuffix(k, key))
      .map(([key, node]) => {
        const path = key;
        const subPath = getSuffix(k, key)!;
        const nqName = (node.conf as PipelineStepConfiguration).nqName;
        return {path, subPath, nqName};
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
      for (const prefix of this.disabledStepsPaths) {
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

export interface StepSpec {
  id: string;
  funcName: string;
  friendlyName?: string;
  helpUrl?: string | HTMLElement;
}

export interface HookSpec {
  hooks: PipelineHooks;
  pipelinePath: ItemPath;
}


export class RuntimeControllerImpl implements RuntimeController {
  public disabledSteps = new Set<string>();

  constructor(private handlerId: string, private config: ControllerConfig, private rt: PipelineRuntime, private signal?: AbortSignal) {
  }

  enableLink(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setLinkState(fullPath, true);
  }

  disableLink(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setLinkState(fullPath, false);
  }

  enableStep(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    this.disabledSteps.delete(pathToKey(fullPath));
    return this.rt.setStepState(fullPath, true);
  }

  disableStep(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    this.disabledSteps.add(pathToKey(fullPath));
    return this.rt.setStepState(fullPath, false);
  }

  triggerLink(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.triggerLink(fullPath);
  }

  getState<T = any>(path: ItemPath): T | void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkInput(fullPath))
      return;

    return this.rt.getState(fullPath);
  }

  setState<T>(path: ItemPath, value: T, inputState?: 'disabled' | 'restricted' | 'user input'): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath))
      return;

    return this.rt.updateState(fullPath, value, inputState);
  }

  getView(path: ItemPath): RichFunctionView | void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.getView(fullPath);
  }

  goToStep(path: ItemPath): void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.goToStep(fullPath);
  }

  setValidation(path: ItemPath, validation: ValidationResult | undefined): void {
    this.checkAborted();
    const fullPathTarget = pathJoin(this.config.pipelinePath, path);
    const fullPathLinkPath = keyToPath(this.handlerId);
    if (!this.checkValidation(fullPathTarget))
      return;
    this.rt.setValidation(fullPathTarget, fullPathLinkPath, validation);
  }

  getValidationAction(path: ItemPath, name?: string): ActionItem {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.getAction(fullPath, name);
  }

  async loadNestedPipeline(path: ItemPath, runId: string) {
    this.checkAborted();
    await this.rt.loadNestedPipeline(path, runId);
    this.checkAborted();
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


//
// Pipeline classes

export interface ICompositionView {
  injectConfiguration(steps: StepSpec[], hooks: HookSpec[], rt: PipelineRuntime, exportConfig?: ExportConfig): void;
  getStateBindings<T = any>(viewId: PathKey, stateId: string): {changes: Observable<T>, setter: (x: any) => void};
  showSteps(...id: NqName[]): void;
  hideSteps(...id: NqName[]): void;
  getStepView<T = RichFunctionView>(viewId?: PathKey): T;
  setExternalValidationResults(viewId: PathKey, stateId: string, results: ValidationResult): void;
  isUpdating: BehaviorSubject<boolean>;
  getRunningUpdates(): string[];
}

export class CompositionPipelineView extends PipelineView implements ICompositionView {
  private hooks: HookSpec[] = [];
  private rt?: PipelineRuntime;
  public customViews = new Map<string, RFVPopup>();

  constructor(funcName: string) {
    super(funcName, [], {historyEnabled: true, isTabbed: false, skipInit: true});
  }

  public injectConfiguration(steps: StepSpec[], hooks: HookSpec[], rt: PipelineRuntime, exportConfig?: ExportConfig) {
    if (exportConfig)
      this.exportConfig = {...this.exportConfig, ...exportConfig};

    this.initialConfig = steps.map((step) => ({...step, customId: step.id}));
    this.hooks = hooks;
    this.rt = rt;
    this.rt.isUpdating.pipe(takeUntil(this.rt.pipelineState.closed)).subscribe((val) => this.isUpdating.next(val));
  }

  public getStateBindings(viewId: PathKey, stateId: string) {
    const stepView = this.getStepView<RichFunctionView>(viewId);
    const changes = stepView.getParamChanges(stateId);
    const setter = (x: any, inputState?: 'disabled' | 'restricted' | 'user input') => stepView.setInput(stateId, x, inputState);
    return {changes, setter};
  }

  public setExternalValidationResults(viewId: PathKey, inputName: string, results: ValidationResult) {
    const view = this.getStepView(viewId);
    if (view)
      view.setExternalValidationResults(inputName, results);
  }

  override getRunningUpdates() {
    return this.rt!.getRunningUpdates();
  }

  override getStepView<T = RichFunctionView>(viewId: PathKey) {
    const view = super.getStepView<RichFunctionView>(viewId);
    if (!view)
      return this.customViews.get(viewId)! as T;

    return view as T;
  }

  override async init() {
    this.rt!.disableSubStepsIOSetters([]);
    await this.execHooks('beforeInit');
    await super.init();
    await this.execHooks('afterInit');
  }

  override async onBeforeStepFuncCallApply(nqName: string, scriptCall: DG.FuncCall, editorFunc: DG.Func) {
    await this.execHooks('beforeFuncCallReady', {nqName, scriptCall, editorFunc});
  }

  override async onAfterStepFuncCallApply(nqName: string, scriptCall: DG.FuncCall, view: RichFunctionView) {
    await this.execHooks('afterFuncCallReady', {nqName, scriptCall, view});
  }

  override async onBeforeLoadRun() {
    this.rt!.disableSubStepsIOSetters([]);
    await super.onBeforeLoadRun();
    await this.execHooks('beforeLoadRun');
  }

  override async onAfterLoadRun(run: DG.FuncCall) {
    await super.onAfterLoadRun(run);
    this.rt!.enableSubStepsIOSetters([]);
    await this.execHooks('afterLoadRun', {run});
  }

  override build() {
    super.build();
    this.rt!.wireViews();
    this.rt!.wireLinks();
    this.rt!.enableSubStepsIOSetters([]);
    this.execHooks('onViewReady', {view: this});
  }

  override close() {
    this.rt!.pipelineState.closed.next(true);
    for (const v of this.customViews.values())
      v.customClose();

    super.close();
  }

  private async execHooks(category: keyof PipelineHooks, additionalParams: Record<string, any> = {}) {
    for (const {hooks, pipelinePath} of this.hooks!) {
      const items = hooks[category];
      for (const item of items ?? []) {
        const handler = item.handler!;
        const ctrlConf = new ControllerConfig(pipelinePath, item.from, item.to);
        const controller = new RuntimeControllerImpl(item.id, ctrlConf, this.rt!);
        const params = {...additionalParams, controller};
        await callHandler(handler, params);
      }
    }
  }
}

export class RFVPopup extends RichFunctionView {
  override detach() {}

  customClose() {
    super.detach();
    super.close();
  }
}

export class CompositionPipeline {
  private id?: ItemName;
  private nqName?: NqName;
  private exportConfig?: ExportConfig;

  private pipelineState = new PipelineGlobalState();

  private config?: CompositionGraphConfig;
  private ioInfo = new Map<NqName, StateItemConfiguration[]>();
  private nodes = new Map<PathKey, NodeState>();
  private links = new Map<PathKey, LinkState>();
  private hooks: HookSpec[] = [];
  private steps: { id: string, funcName: string, friendlyName?: string, helpUrl?: string | HTMLElement }[] = [];
  private rt?: PipelineRuntime;

  private viewInst?: CompositionPipelineView;
  private isInit = false;

  static compose(
    compositionConfig: PipelineCompositionConfiguration,
    nestedConfigs: PipelineConfigVariants[],
  ): CompositionGraphConfig {
    const nestedPipelines: {
      [key: PathKey]: PipelineConfigVariants;
    } = {};
    for (const conf of nestedConfigs) {
      const {id} = conf;
      nestedPipelines[id] = cloneConfig(conf);
    }
    return {...compositionConfig, nestedPipelines};
  }

  constructor(conf?: PipelineConfigVariants) {
    if (conf)
      this.setConfig(conf);
  }

  public setConfig(conf: PipelineConfigVariants) {
    if (this.config)
      throw new Error('Pipeline already has config');
    this.config = cloneConfig(conf);
    this.id = this.config.id;
    if (this.nqName && this.nqName !== this.config.nqName)
      throw new Error(`Config different wrapper nqName ${this.config.nqName}, already set to ${this.nqName}`);
    this.nqName = this.config.nqName;
    this.exportConfig = this.config.exportConfig;
  }

  public makePipelineView(nqName = this.nqName) {
    if (this.viewInst)
      throw new Error(`View has been already created for pipeline ${nqName}`);

    if (!nqName)
      throw new Error('No nqName for pipeline');

    if (this.nqName && this.nqName !== nqName)
      throw new Error(`View different wrapper nqName ${nqName}, already set to ${this.nqName}`);

    this.viewInst = new CompositionPipelineView(nqName);
    this.nqName = nqName;
    return this.viewInst;
  }

  public async init() {
    if (this.isInit)
      throw new Error(`Double init for pipeline ${this.nqName} ${this.id}`);

    this.isInit = true;
    if (!this.viewInst)
      throw new Error(`No view for pipeline ${this.nqName} ${this.id}`);

    await this.loadFuncCallsIO();
    this.addScriptStatesToConfig();
    this.processConfig();
    this.addSystemHooks();
    this.rt = new PipelineRuntime(this.nodes, this.links, this.viewInst, this.pipelineState);
    this.viewInst.injectConfiguration(this.steps, this.hooks, this.rt, this.exportConfig);
    this.isInit = true;
    await this.viewInst.init();
  }

  private async loadFuncCallsIO() {
    const names = await this.gatherIONqNames();
    for (const name of names) {
      const spec = await this.getFuncIOSpec(name);
      this.ioInfo.set(name, spec);
    }
  }

  private addSystemHooks() {
    const stepIdsToBtnNodes = new Map<PathKey, NodeState[]>();
    const additionalViewNames = new Map<PathKey, NqName>();
    const additionalViews = this.viewInst!.customViews;

    for (const node of this.nodes.values()) {
      if (node.type === 'popup' || node.type === 'action') {
        if (node.type === 'popup') {
          const conf = node.conf as PipelinePopupConfiguration;
          additionalViewNames.set(conf.id, conf.nqName);
        }

        const conf = node.conf as (PipelineActionConfiguraion | PipelinePopupConfiguration);

        if (conf.position === 'none')
          continue;

        const parentPathKey = getParentKey(node.conf.id);

        // TODO: menu handling
        if (conf.position === 'buttons') {
          const nodes = stepIdsToBtnNodes.get(parentPathKey) ?? [];
          nodes.push(node);
          stepIdsToBtnNodes.set(parentPathKey, nodes);
        }
      }
    }

    const beforeInit = [{
      id: '_SystemViewsAdd_',
      handler: async () => {
        for (const [k, nqName] of additionalViewNames.entries()) {
          const view = new RFVPopup(nqName, {historyEnabled: false, isTabbed: true});
          await view.isReady.pipe(filter((x) => x), take(1)).toPromise();
          (view.root).style.width = '100%';
          (view.root).style.height = '100%';
          (view.root).style.overflow = 'hidden';
          additionalViews.set(k, view);
        }
      },
    }];

    const onViewReady = [{
      id: '_SystemButtonsAdd_',
      handler: async ({controller}: { controller: RuntimeController }) => {
        for (const [viewKey, nodes] of stepIdsToBtnNodes.entries()) {
          const btns = nodes.map((node) => {
            const conf = node.conf as (PipelineActionConfiguraion | PipelinePopupConfiguration);
            const action = this.rt!.getAction(keyToPath(conf.id));
            return ui.button(action.actionName, action.action);
          });
          const view = controller.getView(keyToPath(viewKey))!;
          view.setAdditionalButtons(btns);
        }
      },
    }];

    const hooks = {
      beforeInit,
      onViewReady,
    };
    this.hooks.unshift({
      pipelinePath: [],
      hooks,
    });
  }

  private processConfig() {
    if (!this.config)
      throw new Error('No pipeline config');

    // process merge config
    const {toRemove, toAdd, stepConfig} = traverseConfigPipelines(
      this.config,
      (acc, node, path) => {
        if (isCompositionConfig(node))
          this.processMergeConfig(node, path, acc.toRemove, acc.toAdd, acc.stepConfig);
        return acc;
      },
      {
        toRemove: new Set<string>(),
        toAdd: new Map<string, ItemsToMerge>(),
        stepConfig: new Map<string, NestedPipelineConfig>(),
      },
    );

    // process hoooks
    this.hooks = traverseConfigPipelines(
      this.config,
      (acc, node, path) => {
        const hooks = this.getPipelineHooks(node, toRemove, path);
        acc.unshift({hooks, pipelinePath: pathJoin(path, [node.id])});
        return acc;
      },
      [] as HookSpec[],
    );

    // process node items
    traverseConfigPipelines(
      this.config,
      (acc, node, path) => {
        this.processPipelineConfig(node, path, toRemove, toAdd);
        return acc;
      },
      undefined,
    );

    // get steps sequence
    const {pipelineSteps, pipelineSeq} = traverseConfigPipelines(
      this.config,
      (acc, node) => {
        const id = node.id;
        const specs = node.steps.map((s) => ({id: s.id, funcName: s.nqName, friendlyName: s.friendlyName, helpUrl: s.helpUrl}));
        acc.pipelineSteps.set(id, specs);
        acc.pipelineSeq.unshift(id);
        return acc;
      },
      {
        pipelineSteps: new Map<ItemName, StepSpec[]>(),
        pipelineSeq: [] as ItemName[],
      },
    );

    this.steps = pipelineSeq.reduce((acc, id) => {
      const steps = pipelineSteps.get(id)!;
      const pos = stepConfig.get(id)?.insertBeforeStep;
      if (pos) {
        const idx = acc.findIndex((spec) => spec.id === pos);
        if (idx > 0)
          acc.splice(idx, 0, ...steps);

        return acc;
      } else
        return [...acc, ...steps ?? []];
    }, [] as StepSpec[]);
  }

  // deal with hooks

  private getPipelineHooks(node: PipelineConfiguration, toRemove: Set<string>, path: ItemPath) {
    node.hooks = node.hooks ?? {};
    for (const [type, hooks] of Object.entries(node.hooks ?? {})) {
      const filteredHooks = [];
      for (const hook of hooks ?? []) {
        const hookPath = this.updateFullPathLink(hook, pathJoin(path, [node.id]));
        const id = pathToKey(hookPath);
        if (toRemove.has(id))
          continue;

        filteredHooks.push(hook);
      }
      node.hooks[type as keyof PipelineHooks] = filteredHooks;
    }
    return node.hooks;
  }

  // pipeline config processing

  private processPipelineConfig(pipelineConf: PipelineConfiguration, pipelinePath: ItemPath, toRemove: Set<string>, toAdd: Map<string, ItemsToMerge>) {
    const subPath = this.proccessPipelineNodeConfig(pipelineConf, pipelinePath);
    const steps: PipelineStepConfiguration[] = [];
    for (const stepConf of pipelineConf.steps) {
      const sKey = stepConf.id;
      const popups: PipelinePopupConfiguration[] = [];
      for (const popupConf of this.getAllItems(stepConf.popups ?? [], sKey, toAdd, 'popup')) {
        const pPopupConf = this.processNodeConfig(popupConf, keyToPath(sKey), subPath, toRemove, 'popup');
        if (!pPopupConf)
          continue;

        const actions: PipelineActionConfiguraion[] = [];
        for (const actionConf of this.getAllItems(popupConf.actions ?? [], pPopupConf.id, toAdd, 'action')) {
          const pActionConf = this.processActionNodeConfig(actionConf, keyToPath(pPopupConf.id), subPath, toRemove);
          if (!pActionConf)
            continue;

          actions.push(pActionConf);
        }
        popupConf.actions = actions;
        popups.push(popupConf);
      }
      const actions: PipelineActionConfiguraion[] = [];
      for (const actionConf of this.getAllItems(stepConf.actions ?? [], sKey, toAdd, 'action')) {
        const pActionConf = this.processActionNodeConfig(actionConf, keyToPath(sKey), subPath, toRemove);
        if (!pActionConf)
          continue;

        actions.push(pActionConf);
      }
      stepConf.popups = popups;
      stepConf.actions = actions;
      steps.push(stepConf);
    }
    const links: PipelineLinkConfiguration[] = [];
    for (const link of pipelineConf.links ?? []) {
      const plink = this.processLinkConfig(link, subPath, toRemove);
      if (!plink)
        continue;

      links.push(plink);
    }
    const actions: PipelineActionConfiguraion[] = [];
    for (const actionConf of this.getAllItems(pipelineConf.actions ?? [], pipelineConf.id, toAdd, 'action')) {
      const pActionConf = this.processActionNodeConfig(actionConf, subPath, subPath, toRemove);
      if (!pActionConf)
        continue;

      actions.push(pActionConf);
    }

    pipelineConf.steps = steps;
    pipelineConf.links = links;
    pipelineConf.actions = actions;
  }

  private getAllItems<T extends SubNodeConf>(data: T[], prefix: string, toAdd: Map<string, ItemsToMerge>, type: SubNodeConfTypes) {
    const moreItems = toAdd.get(prefix)?.[type] ?? [];
    return [...data, ...moreItems] as T[];
  }

  private proccessPipelineNodeConfig(conf: PipelineConfiguration, pipelinePath: ItemPath) {
    const nodePath = this.updateFullPathNode(conf, pipelinePath);
    const nodeKey = pathToKey(nodePath);

    this.nodes.set(nodeKey, new NodeState(conf, 'pipeline', pipelinePath, this.pipelineState));
    return nodePath;
  }

  private processActionNodeConfig(conf: PipelineActionConfiguraion, path: ItemPath, pipelinePath: ItemPath, toRemove: Set<string>) {
    this.updateFullPathLink(conf, pipelinePath, false);
    const nodePath = this.updateFullPathNode(conf, path);
    const nodeKey = pathToKey(nodePath);
    if (toRemove.has(nodeKey))
      return;

    this.nodes.set(nodeKey, new NodeState(conf, 'action', pipelinePath, this.pipelineState));
    return conf;
  }

  private processNodeConfig<T extends NodeConf>(conf: T, path: ItemPath, pipelinePath: ItemPath, toRemove: Set<string>, type: NodeConfTypes) {
    const nodePath = this.updateFullPathNode(conf, path);
    const nodeKey = pathToKey(nodePath);

    if (toRemove.has(nodeKey))
      return;

    this.nodes.set(nodeKey, new NodeState(conf, type, pipelinePath, this.pipelineState));
    return conf;
  }

  private processLinkConfig(conf: PipelineLinkConfiguration, pipelinePath: ItemPath, toRemove: Set<string>) {
    const linkPath = this.updateFullPathLink(conf, pipelinePath);
    const linkKey = pathToKey(linkPath);
    if (toRemove.has(linkKey))
      return;

    this.links.set(linkKey, new LinkState(conf, pipelinePath, this.pipelineState));
    return conf;
  }

  private updateFullPathLink(target: { id: string, from?: ItemPath | ItemPath[], to?: ItemPath | ItemPath[] }, currentPath: ItemPath, updateId = true) {
    if (target.from) {
      if (Array.isArray(target.from[0]))
        target.from = target.from.map((path) => [...currentPath, ...path]);
      else
        target.from = [...currentPath, ...(target.from as ItemPath)];
    }

    if (target.to) {
      if (Array.isArray(target.to[0]))
        target.to = target.to.map((path) => [...currentPath, ...path]);
      else
        target.to = [...currentPath, ...(target.to as ItemPath)];
    }

    if (updateId)
      return this.updateFullPathNode(target, currentPath);

    return [target.id];
  }

  private updateFullPathNode(target: { id: string }, path: ItemPath) {
    const nPath = pathJoin(path, [target.id]);
    target.id = pathToKey(nPath);
    return nPath;
  }

  // merge config processing

  private processMergeConfig(mergeConf: PipelineCompositionConfiguration,
    path: ItemPath, toRemove: Set<string>, toAdd: Map<string, ItemsToMerge>, stepConfig: Map<string, NestedPipelineConfig>,
  ) {
    const subPath = pathJoin(path, [mergeConf.id]);
    for (const item of mergeConf.itemsToRemove ?? []) {
      const itemPath = pathJoin(subPath, item);
      const itemKey = pathToKey(itemPath);
      toRemove.add(itemKey);
    }
    delete mergeConf.itemsToRemove;
    for (const [pipelineId, conf] of Object.entries(mergeConf.nestedPipelinesConfig ?? {})) {
      const pipelineIdFull = pathToKey(pathJoin(subPath, [pipelineId]));
      if (conf.insertBeforeStep) {
        const insertBeforeStep = pathToKey(pathJoin(subPath, [conf.insertBeforeStep]));
        stepConfig.set(pipelineIdFull, { ...conf, insertBeforeStep });
      } else {
        stepConfig.set(pipelineIdFull, { ...conf });
      }
    }
    delete mergeConf.nestedPipelinesConfig;
    const {popupsToAdd, actionsToAdd} = mergeConf;
    this.processMergeNodeList(popupsToAdd ?? [], subPath, toAdd, 'popup');
    this.processMergeNodeList(actionsToAdd ?? [], subPath, toAdd, 'action');
    delete mergeConf.popupsToAdd;
    delete mergeConf.actionsToAdd;
  }

  private processMergeNodeList(itemsToAdd: [SubNodeConf, ItemPath?][], pipelinePath: ItemPath,
    toAdd: Map<string, ItemsToMerge>, type: SubNodeConfTypes,
  ) {
    for (const [item, path] of itemsToAdd) {
      const itemPath = pathJoin(pipelinePath, path ?? []);
      const tagetKey = pathToKey(itemPath);
      const addConf = toAdd.get(tagetKey) ?? {};
      const addArray = addConf[type] ?? [];
      addArray.unshift(item as any);
      addConf[type] = addArray as any;
      toAdd.set(tagetKey, addConf);
    }
  }

  // funcall states

  private addScriptStatesToConfig() {
    if (!this.config)
      throw new Error('No pipeline config');

    traverseConfigPipelines(
      this.config,
      (acc, node) => {
        for (const conf of node.steps) {
          this.addNodeIOInfo(conf);
          for (const pconf of conf.popups ?? [])
            this.addNodeIOInfo(pconf);
        }
        if (isCompositionConfig(node)) {
          for (const [conf] of node.popupsToAdd ?? [])
            this.addNodeIOInfo(conf);
        }
        return acc;
      },
      undefined,
    );
  }

  private addNodeIOInfo(conf: PipelineStepConfiguration | PipelinePopupConfiguration) {
    const info = cloneConfig(this.ioInfo.get(conf.nqName)) ?? [];
    if (!conf.states)
      conf.states = info;
    else
      conf.states.push(...info);
  }

  private async getFuncIOSpec(nqName: NqName): Promise<StateItemConfiguration[]> {
    const func: DG.Func = await grok.functions.eval(nqName);
    const inputs = func.inputs.map((input) => ({id: input.name, type: input.propertyType as any, stateType: 'input' as StateType}));
    const outputs = func.outputs.map((output) => ({id: output.name, type: output.propertyType as any, stateType: 'output' as StateType}));
    const io = [...inputs, ...outputs];
    return io;
  }

  private async gatherIONqNames() {
    if (!this.config)
      throw new Error('No pipeline config');

    const names = new Set<string>();
    traverseConfigPipelines(
      this.config,
      (acc, node) => {
        for (const conf of node.steps) {
          names.add(conf.nqName);
          for (const pconf of conf.popups ?? [])
            names.add(pconf.nqName);
        }
        if (isCompositionConfig(node)) {
          for (const [conf] of node.popupsToAdd ?? [])
            names.add(conf.nqName);
        }
        return acc;
      },
      names,
    );
    return names;
  }
}

async function callHandler(fn: string | Function, params: Record<string, any>) {
  if (typeof fn === 'string') {
    const f: DG.Func = await grok.functions.eval(fn);
    const call = f.prepare({params});
    await call.call();
    return call.outputs;
  } else {
    const res = await fn(params);
    return res;
  }
}

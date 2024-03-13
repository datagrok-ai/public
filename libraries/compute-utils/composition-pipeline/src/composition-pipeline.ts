import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import cloneDeepWith from 'lodash.clonedeepwith';
import {BehaviorSubject, Observable, merge, of, from, Subject} from 'rxjs';
import {FunctionView, PipelineView, RichFunctionView} from '../../function-views';
import {withLatestFrom, filter, map, switchMap, debounceTime, takeUntil, take, sample} from 'rxjs/operators';

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
  enableStep(path: ItemPath): void;
  disableStep(path: ItemPath): void;
  getState<T = any>(path: ItemPath): T | void;
  updateState<T>(path: ItemPath, state: T, inputState?: 'disabled' | 'restricted' | 'user input'): void;
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
  enableOnStart?: boolean;
  enableOnLoadRun?: boolean;
}

export type PipelineHookConfiguration = {
  id: ItemName;
  from?: ItemPath | ItemPath[];
  to?: ItemPath | ItemPath[];
  handler: Handler;
}

export type ActionConfiguraion = PipelineHookConfiguration & {
  position: 'buttons' | 'menu';
  friendlyName: string;
  runStepOnComplete?: boolean;
}

export type PipelinePopupConfiguration = {
  id: ItemName;
  nqName: NqName;
  position: 'buttons' | 'menu';
  friendlyName: string;
  helpUrl?: string;
  states?: StateItemConfiguration[];
  actions?: ActionConfiguraion[];
}

export type PipelineStepConfiguration = {
  id: ItemName;
  nqName: NqName;
  friendlyName?: string;
  helpUrl?: string;
  states?: StateItemConfiguration[];
  actions?: ActionConfiguraion[];
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
  // non-composable
  states?: StateItemConfiguration[];
  exportConfig?: ExportConfig;
  helpUrl?: string;
}

//
// Composition Pipeline configuration

export type ItemsToAdd = {
  stepsToAdd?: [PipelineStepConfiguration, ItemPath][];
  popupsToAdd?: [PipelinePopupConfiguration, ItemPath][];
  actionsToAdd?: [ActionConfiguraion, ItemPath][];
}

export type ItemsToRemove = {
  itemsToRemove?: ItemPath[];
}

export type PipelineCompositionConfiguration = ItemsToAdd & PipelineConfiguration & ItemsToRemove;

//
// Internal config

export type PipelineConfigVariants = PipelineConfiguration | CompositionGraphConfig;

export type CompositionGraphConfig = (PipelineConfiguration | PipelineCompositionConfiguration) & {
  nestedPipelines?: {
    [key: ItemName]: PipelineConfiguration | CompositionGraphConfig;
  };
}

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

export function normalizePaths(paths: ItemPath | ItemPath[]): ItemPath[] {
  if (Array.isArray(paths[0]))
    return paths as ItemPath[];
  else
    return [paths as ItemPath];
}

function getStepId(conf: PipelinePopupConfiguration | PipelineStepConfiguration) {
  return conf.nqName;
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
  action?: ActionConfiguraion[];
}
export type SubNodeConf = ActionConfiguraion | PipelinePopupConfiguration | PipelineStepConfiguration;
export type SubNodeConfTypes = 'action' | 'popup' | 'step';
export type NodeConfTypes = SubNodeConfTypes | 'pipeline';
export type NodeConf = SubNodeConf | PipelineConfiguration;

export type InputState = 'disabled' | 'restricted' | 'user input';

export class NodeItemState<T = any> {
  public currentSource = new BehaviorSubject<Observable<T>>(of());
  private valueChanges = this.currentSource.pipe(
    switchMap((source) => source),
  );

  public value = new BehaviorSubject<T | undefined>(undefined);

  private setter = (x: T, inputState?: InputState) => {
    this.currentSource.next(of(x));
  };

  constructor(public conf: StateItemConfiguration, public pipelineState: PipelineState, public notifier?: Observable<true>) {
    this.valueChanges.pipe(
      this.notifier ? sample(this.notifier) : map((x) => x),
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
    public type: SubNodeConfTypes,
    public pipelinePath: ItemPath,
    public parrentNodePath: ItemPath,
    public pipelineState: PipelineState,
  ) {
    if (type === 'step' || type === 'popup') {
      const states = (conf as PipelineStepConfiguration | PipelinePopupConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeItemState(state, this.pipelineState, type === 'popup' ? this.notifier : undefined));
    }
    if (type === 'action') {
      const link = conf as ActionConfiguraion;
      this.controllerConfig = new ControllerConfig(pipelinePath, link.from, link.to);
    }
  }
}

export class PipelineState {
  initialLoading = new BehaviorSubject<boolean>(true);
  runChanging = new BehaviorSubject<boolean>(false);
  closed = new Subject<true>();
}

export class LinkState {
  public controllerConfig: ControllerConfig;
  public enabled = new BehaviorSubject(true);
  private currentSource = new BehaviorSubject<Observable<any>>(of());
  public valueChanges = this.currentSource.pipe(
    switchMap((source) => source),
  );

  constructor(
    public conf: PipelineLinkConfiguration,
    public pipelinePath: ItemPath,
    public pipelineState: PipelineState,
  ) {
    this.controllerConfig = new ControllerConfig(pipelinePath, conf.from, conf.to);
  }

  public setSources(sources: Observable<any>[]) {
    this.currentSource.next(merge(...sources));
  }

  public getValuesChanges() {
    return this.valueChanges.pipe(
      withLatestFrom(this.enabled, this.pipelineState.initialLoading, this.pipelineState.runChanging),
      filter(([, enabled, initialLoading, runChanging]) => {
        if (!enabled)
          return false;

        if ((initialLoading && !this.conf.enableOnStart) || (runChanging && !this.conf.enableOnLoadRun))
          return false;

        return true;
      }),
      // tap(() => console.log('link', this.conf.id)),
      map(() => true),
    );
  }
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
  constructor(
    private nodes: Map<PathKey, NodeState>,
    private links: Map<PathKey, LinkState>,
    private view: ICompositionView,
    public pipelineState: PipelineState,
  ) {}

  setLinkState(path: ItemPath, enabled: boolean): void {
    const k = pathToKey(path);
    const link = this.links.get(k);
    if (link)
      link.enabled.next(enabled);
  }

  setStepState(path: ItemPath, enabled: boolean): void {
    const k = pathToKey(path);
    const conf = this.nodes.get(k)?.conf;
    if (!conf)
      return;

    if (enabled)
      this.view.showSteps((conf as PipelineStepConfiguration).nqName);
    else
      this.view.hideSteps((conf as PipelineCompositionConfiguration).nqName);
  }

  getState<T = any>(path: ItemPath): T | void {
    const {state} = this.getNodeState(path)!;
    const val = state?.getValue();
    return val;
  }

  updateState<T>(path: ItemPath, value: T, inputState: InputState = 'user input'): void {
    const {state} = this.getNodeState(path)!;
    if (state)
      state.setValue(value, inputState);
  }

  getView(path: ItemPath): RichFunctionView | void {
    const k = pathToKey(path);
    const conf = this.nodes.get(k)?.conf;
    if (conf)
      return this.view.getStepView<RichFunctionView>((conf as PipelineStepConfiguration).nqName, k);
  }

  wireView() {
    // wiring by nqName
    for (const [, node] of this.nodes) {
      for (const [, state] of node.states) {
        if ((node.type === 'popup' || node.type === 'step') && (state.conf.stateType === 'input' || state.conf.stateType === 'output')) {
          const conf = node.conf as PipelinePopupConfiguration | PipelineStepConfiguration;
          const stepId = getStepId(conf);
          const stateId = state.conf.id;
          const {changes, setter} = this.view.getStateBindings(stepId, stateId, node.conf.id);
          state.linkState(changes, setter);
        }
      }
    }
  }

  wireLinks() {
    for (const [, link] of this.links) {
      const changes = link.controllerConfig.from.map((input) => {
        return this.getNodeState(input)!.state.value;
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
          const abortController = new AbortController();
          const signal = abortController.signal;
          let done = false;
          const obs$ = new Observable((observer) => {
            const controller = new RuntimeControllerImpl(link.conf.id, link.controllerConfig, this, signal);
            const sub = from(callHandler(handler, {controller})).subscribe((val) => {
              done = true;
              observer.next(val);
            }, (error: any) => {
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

  goToStep(path: ItemPath): void {
    console.log('TODO: goToStep not implemented');
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
}

export interface StepSpec {
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

  getState<T = any>(path: ItemPath): T | void {
    this.checkAborted();
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkInput(fullPath))
      return;

    return this.rt.getState(fullPath);
  }

  updateState<T>(path: ItemPath, value: T, inputState?: 'disabled' | 'restricted' | 'user input'): void {
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

    return true;
  }
}


//
// Pipeline classes

export interface ICompositionView {
  showSteps(...id: string[]): void;
  hideSteps(...id: string[]): void;
  getStepView<T extends RichFunctionView>(id: string): T;
  injectConfiguration(steps: StepSpec[], hooks: HookSpec[], rt: PipelineRuntime, exportConfig?: ExportConfig): void;
  getStateBindings<T = any>(stepId: string, stateId: string, k: string): {changes: Observable<T>, setter: (x: any) => void};
  getStepView<T = RichFunctionView>(name: string, k?: string): T;
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

    this.initialConfig = steps;
    this.hooks = hooks;
    this.rt = rt;
  }

  public getStateBindings(stepId: string, stateId: string, k: string) {
    const stepView = this.getStepView<RichFunctionView>(stepId, k);
    const changes = stepView.getParamChanges(stateId);
    const setter = (x: any, inputState?: 'disabled' | 'restricted' | 'user input') => stepView.setInput(stateId, x, inputState);
    return {changes, setter};
  }

  override getStepView<T = RichFunctionView>(name: string, k?: string) {
    const view = super.getStepView<RichFunctionView>(name);
    if (!view && k)
      return this.customViews.get(k)! as T;

    return view as T;
  }

  override async init() {
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
    this.rt!.pipelineState.runChanging.next(true);
    await super.onBeforeLoadRun();
    await this.execHooks('beforeLoadRun');
  }

  override async onAfterLoadRun(run: DG.FuncCall) {
    await super.onAfterLoadRun(run);
    this.rt!.pipelineState.runChanging.next(false);
    this.rt!.wireView();
    await this.execHooks('afterLoadRun', {run});
  }

  override build() {
    super.build();
    this.rt!.wireView();
    this.rt!.wireLinks();
    this.rt!.pipelineState.initialLoading.next(false);
    this.execHooks('onViewReady');
  }

  override close() {
    this.rt!.pipelineState.closed.next(true);
    for (const v of this.customViews.values()) {
      v.customClose();
    }
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
  override detach() {

  }

  customClose() {
    super.detach();
    super.close();
  }
}

export class CompositionPipeline {
  private id: ItemName;
  private nqName: NqName;
  private exportConfig?: ExportConfig;

  private pipelineState = new PipelineState();

  private config: CompositionGraphConfig;
  private ioInfo = new Map<NqName, StateItemConfiguration[]>();
  private nodes = new Map<PathKey, NodeState>();
  private links = new Map<PathKey, LinkState>();
  private hooks: HookSpec[] = [];
  private steps: { funcName: string, friendlyName?: string, helpUrl?: string | HTMLElement }[] = [];
  private rt?: PipelineRuntime;

  private viewInst?: CompositionPipelineView;
  private isInit = false;

  static compose(
    config: PipelineCompositionConfiguration,
    items: PipelineConfigVariants[],
  ): CompositionGraphConfig {
    const nestedPipelines: {
      [key: PathKey]: PipelineConfigVariants;
    } = {};
    for (const conf of items) {
      const {id} = conf;
      nestedPipelines[id] = cloneConfig(conf);
    }
    return {...config, nestedPipelines};
  }

  constructor(conf: PipelineConfigVariants) {
    this.config = cloneConfig(conf);
    this.id = this.config.id;
    this.nqName = this.config.nqName;
    this.exportConfig = this.config.exportConfig;
  }

  public makePipelineView() {
    if (this.viewInst)
      throw new Error(`View has been already created for pipeline ${this.nqName} ${this.id}`);

    this.viewInst = new CompositionPipelineView(this.nqName);
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
    // TODO: set menu as buttons for now;
    const stepIdsToNodes = new Map<string, NodeState[]>();
    const additionalViewNames = new Map<string, NqName>();
    const additionalViews = this.viewInst!.customViews;

    for (const node of this.nodes.values()) {
      if (node.type === 'popup' || node.type === 'action') {
        const parrentPathKey = pathToKey(node.parrentNodePath);
        const nodes = stepIdsToNodes.get(parrentPathKey) ?? [];
        nodes.push(node);
        stepIdsToNodes.set(parrentPathKey, nodes);
        if (node.type === 'popup') {
          const conf = node.conf as PipelinePopupConfiguration;
          additionalViewNames.set(conf.id, conf.nqName);
        }
      }
    }

    const beforeInit = [{
      id: '_SystemViewsAdd_',
      handler: async () => {
        for (const [k, nqName] of additionalViewNames.entries()) {
          const view = new RFVPopup(nqName, {historyEnabled: false, isTabbed: true});
          await view.isReady.pipe(filter((x) => x), take(1)).toPromise();
          additionalViews.set(k, view);
        }
      },
    }];

    const onViewReady = [{
      id: '_SystemButtonsAdd_',
      handler: async ({controller}: { controller: RuntimeController }) => {
        for (const [viewKey, nodes] of stepIdsToNodes.entries()) {
          const btns = nodes.map((node) => {
            if (node.type === 'action') {
              const conf = node.conf as ActionConfiguraion;
              const handler = conf.handler;
              const ctrlConf = node.controllerConfig!;
              return ui.button(conf.friendlyName, async () => {
                const controller = new RuntimeControllerImpl(node.conf.id, ctrlConf, this.rt!);
                await callHandler(handler, {controller});
              });
            } else if (node.type === 'popup') {
              const conf = node.conf as PipelinePopupConfiguration;
              const view = additionalViews.get(conf.id)!;
              (view.root).style.width = '100%';
              (view.root).style.height = '100%';
              (view.root).style.overflow = 'hidden';
              const buttonName = conf.friendlyName;
              return ui.button(buttonName, async () => {
                await new Promise((resolve, _reject) => {
                  ui.dialog(buttonName)
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
              });
            }
            throw new Error(`No element for ${node.conf.id} ${node.type}`);
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
    // process merge config
    const {toRemove, toAdd} = traverseConfigPipelines(
      this.config,
      (acc, node, path) => {
        this.processMergeConfig(node, path, acc.toRemove, acc.toAdd);
        return acc;
      },
      {
        toRemove: new Set<string>(),
        toAdd: new Map<string, ItemsToMerge>(),
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
    this.steps = traverseConfigPipelines(
      this.config,
      (acc, node) => {
        acc.unshift(...node.steps.map((s) => ({funcName: s.nqName, friendlyName: s.friendlyName, helpUrl: s.helpUrl})));
        return acc;
      },
      [] as StepSpec[],
    );
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
    const subPath = this.updateFullPathNode(pipelineConf, pipelinePath);
    const steps: PipelineStepConfiguration[] = [];
    for (const stepConf of this.getAllSteps(pipelineConf.steps, subPath, toRemove, toAdd)) {
      const sKey = stepConf.id;
      const popups: PipelinePopupConfiguration[] = [];
      for (const popupConf of this.getAllConfItems(stepConf.popups ?? [], sKey, toAdd, 'popup')) {
        const pPopupConf = this.processNodeConfig(popupConf, keyToPath(sKey), subPath, toRemove, 'popup');
        if (!pPopupConf)
          continue;

        const actions: ActionConfiguraion[] = [];
        for (const actionConf of this.getAllConfItems(popupConf.actions ?? [], pPopupConf.id, toAdd, 'action')) {
          const pActionConf = this.processActionNodeConfig(actionConf, keyToPath(pPopupConf.id), subPath, toRemove);
          if (!pActionConf)
            continue;

          actions.push(pActionConf);
        }
        popupConf.actions = actions;
        popups.push(popupConf);
      }
      const actions: ActionConfiguraion[] = [];
      for (const actionConf of this.getAllConfItems(stepConf.actions ?? [], sKey, toAdd, 'action')) {
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
    pipelineConf.steps = steps;
    pipelineConf.links = links;
  }

  private getAllSteps(data: PipelineStepConfiguration[], pipelinePath: ItemPath, toRemove: Set<string>, toAdd: Map<string, ItemsToMerge>) {
    const nData: PipelineStepConfiguration[] = data.flatMap((conf) => {
      const pconf = this.processNodeConfig<PipelineStepConfiguration>(conf, pipelinePath, pipelinePath, toRemove, 'step');
      const addConfs = (toAdd.get(conf.id)?.step ?? [])
        .map((step) => this.processNodeConfig(step, pipelinePath, pipelinePath, toRemove, 'step'))
        .filter((x) => x) as PipelineStepConfiguration[];
      return [...(pconf ? [pconf] : []), ...addConfs];
    });
    return nData;
  }

  private getAllConfItems<T extends SubNodeConf>(data: T[], prefix: string,
    toAdd: Map<string, ItemsToMerge>, type: SubNodeConfTypes,
  ) {
    const moreItems = toAdd.get(prefix)?.[type] ?? [];
    return [...data, ...moreItems] as T[];
  }


  private processActionNodeConfig(conf: ActionConfiguraion, path: ItemPath, pipelinePath: ItemPath, toRemove: Set<string>) {
    this.updateFullPathLink(conf, pipelinePath, false);
    const nodePath = this.updateFullPathNode(conf, path);
    const parrentPath = nodePath.slice(0, nodePath.length - 1);
    const nodeKey = pathToKey(nodePath);
    if (toRemove.has(nodeKey))
      return;

    this.nodes.set(nodeKey, new NodeState(conf, 'action', pipelinePath, parrentPath, this.pipelineState));
    return conf;
  }

  private processNodeConfig<T extends NodeConf>(conf: T, path: ItemPath, pipelinePath: ItemPath, toRemove: Set<string>, type: SubNodeConfTypes) {
    const nodePath = this.updateFullPathNode(conf, path);
    const nodeKey = pathToKey(nodePath);
    const parrentPath = nodePath.slice(0, nodePath.length - 1);

    if (toRemove.has(nodeKey))
      return;

    this.nodes.set(nodeKey, new NodeState(conf, type, pipelinePath, parrentPath, this.pipelineState));
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
    path: ItemPath, toRemove: Set<string>, toAdd: Map<string, ItemsToMerge>,
  ) {
    const subPath = pathJoin(path, [mergeConf.id]);
    for (const item of mergeConf.itemsToRemove ?? []) {
      const itemPath = pathJoin(subPath, item);
      const itemKey = pathToKey(itemPath);
      toRemove.add(itemKey);
    }
    delete mergeConf.itemsToRemove;
    const {stepsToAdd, popupsToAdd, actionsToAdd} = mergeConf;
    this.processMergeNodeList(stepsToAdd ?? [], subPath, toAdd, 'step');
    this.processMergeNodeList(popupsToAdd ?? [], subPath, toAdd, 'popup');
    this.processMergeNodeList(actionsToAdd ?? [], subPath, toAdd, 'action');
    delete mergeConf.stepsToAdd;
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
    traverseConfigPipelines(
      this.config,
      (acc, node) => {
        for (const conf of node.steps) {
          this.addNodeIOInfo(conf);
          for (const pconf of conf.popups ?? [])
            this.addNodeIOInfo(pconf);
        }
        if (isCompositionConfig(node)) {
          for (const [conf] of node.stepsToAdd ?? [])
            this.addNodeIOInfo(conf);

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
          for (const [conf] of node.stepsToAdd ?? [])
            names.add(conf.nqName);

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

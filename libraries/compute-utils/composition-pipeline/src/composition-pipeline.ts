import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import cloneDeepWith from 'lodash.clonedeepwith';
import {BehaviorSubject, Observable, merge} from 'rxjs';
import {PipelineView, RichFunctionView} from '../../function-views';
import {withLatestFrom, filter, map} from 'rxjs/operators';

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

export type ActionConfiguraion = PipelineLinkConfiguration & {
  position: 'buttons' | 'menu';
  runStepOnComplete?: boolean;
}

export type PipelinePopupConfiguration = {
  id: ItemName;
  nqName: NqName;
  position: 'buttons' | 'menu';
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
  id: ItemName;
  nqName: NqName;
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

export function pathsToKeys(paths: ItemPath | ItemPath[]): PathKey[] {
  if (Array.isArray(paths[0]))
    return paths.map((path) => pathToKey(path as ItemPath));
  else
    return [pathToKey(paths as ItemPath)];
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

export type ItemsToMerge = {
  step?: PipelineStepConfiguration[];
  popup?: PipelinePopupConfiguration[];
  action?: ActionConfiguraion[];
}
export type AddNodeConf = ActionConfiguraion | PipelinePopupConfiguration | PipelineStepConfiguration;
export type NodeConf = AddNodeConf | PipelineConfiguration;
export type AddNodeConfTypes = 'action' | 'popup' | 'step';
export type NodeConfTypes = AddNodeConfTypes | 'pipeline';

export type NodeItemState<T = any> = StateItemConfiguration & {
  value: BehaviorSubject<T>,
  inputState?: 'disabled' | 'restricted' | 'user input',
}

export class NodeState {
  public controllerConfig?: ControllerConfig;
  public states?: Map<ItemName, NodeItemState>;
  public enabled = new BehaviorSubject(true);

  constructor(
    public conf: NodeConf,
    public type: NodeConfTypes,
    public pipelinePath: ItemPath,
  ) {
    if (type === 'step' || type === 'popup') {
      const states = (conf as PipelineStepConfiguration | PipelinePopupConfiguration).states;
      for (const state of states ?? [])
        this.states?.set(state.id, {...state, value: new BehaviorSubject<any>(undefined)});
    }
    if (type = 'action') {
      const link = conf as ActionConfiguraion;
      this.controllerConfig = new ControllerConfig(pipelinePath, link.from, link.to);
    }
  }
}

export class PipelineState {
  initialLoading = new BehaviorSubject<boolean>(true);
  runChanging = new BehaviorSubject<boolean>(false);
}

export class LinkState {
  public controllerConfig: ControllerConfig;
  public enabled = new BehaviorSubject(true);
  public valueChanges?: Observable<any>;

  constructor(
    public globalState: PipelineState,
    public config: PipelineLinkConfiguration,
    public pipelinePath: ItemPath,
  ) {
    this.controllerConfig = new ControllerConfig(pipelinePath, config.from, config.to);
  }

  public linkValues(sources: Observable<any>[]) {
    this.valueChanges = merge(sources);
  }

  public getValuesChanges() {
    return this.valueChanges!.pipe(
      withLatestFrom(this.enabled, this.globalState.initialLoading, this.globalState.runChanging),
      filter(([_val, enabled, initialLoading, runChanging]) => {
        if (!enabled)
          return false;

        if ((initialLoading && !this.config.enableOnStart) || (runChanging && !this.config.enableOnLoadRun))
          return false;

        return true;
      }),
      map(() => true),
    );
  }
}

export class ControllerConfig {
  public from?: Set<PathKey>;
  public to?: Set<PathKey>;

  constructor(
    public pipelinePath: ItemPath,
    from?: ItemPath | ItemPath[],
    to?: ItemPath | ItemPath[],
  ) {
    if (from)
      this.from = new Set(pathsToKeys(from));

    if (to)
      this.to = new Set(pathsToKeys(to));
  }
}

export class PipelineRuntime {
  constructor(
    private nodes: Map<PathKey, NodeState>,
    private links: Map<PathKey, LinkState>,
    private compositionPipelineView: CompositionPipelineView,
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
      this.compositionPipelineView.showSteps((conf as PipelineStepConfiguration).nqName);
    else
      this.compositionPipelineView.hideSteps((conf as PipelineCompositionConfiguration).nqName);
  }

  getState<T = any>(path: ItemPath): T | void {
    const nodePath = path.slice(0, path.length - 2);
    const stateName = path[path.length - 1];
    const k = pathToKey(nodePath);
    const node = this.nodes.get(k);
    if (node && node.states) {
      const state = node.states.get(stateName);
      if (state)
        return state.value.value;
    }
  }

  updateState<T>(path: ItemPath, value: T, inputState: 'disabled' | 'restricted' | 'user input' = 'user input'): void {
    const nodePath = path.slice(0, path.length - 2);
    const stateName = path[path.length - 1];
    const k = pathToKey(nodePath);
    const node = this.nodes.get(k);
    if (node && node.states) {
      const state = node.states.get(stateName);
      if (state) {
        state.inputState = inputState;
        state.value.next(value);
      }
    }
  }

  getView(path: ItemPath): RichFunctionView | void {
    const k = pathToKey(path);
    const conf = this.nodes.get(k)?.conf;
    if (conf)
      return this.compositionPipelineView.getStepView<RichFunctionView>((conf as PipelineStepConfiguration).nqName);
  }

  goToStep(path: ItemPath): void {
    console.log('TODO: goToStep not implemented');
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
  constructor(private config: ControllerConfig, private rt: PipelineRuntime) {
  }

  enableLink(path: ItemPath): void {
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setLinkState(fullPath, true);
  }

  disableLink(path: ItemPath): void {
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setLinkState(fullPath, false);
  }

  enableStep(path: ItemPath): void {
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setStepState(fullPath, true);
  }

  disableStep(path: ItemPath): void {
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.setStepState(fullPath, false);
  }

  getState<T = any>(path: ItemPath): T | void {
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkInput(fullPath)) {
      grok.shell.warning(`No input link ${fullPath}, pipeline: ${this.config.pipelinePath}`);
      return;
    }
    return this.rt.getState(fullPath);
  }

  updateState<T>(path: ItemPath, value: T, inputState?: 'disabled' | 'restricted' | 'user input'): void {
    const fullPath = pathJoin(this.config.pipelinePath, path);
    if (!this.checkOutput(fullPath)) {
      grok.shell.warning(`No input link ${fullPath}, pipeline: ${this.config.pipelinePath}`);
      return;
    }
    return this.rt.updateState(fullPath, value, inputState);
  }

  getView(path: ItemPath): RichFunctionView | void {
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.getView(fullPath);
  }

  goToStep(path: ItemPath): void {
    const fullPath = pathJoin(this.config.pipelinePath, path);
    return this.rt.goToStep(fullPath);
  }

  private checkInput(path: ItemPath) {
    // TODO: check if has an active link
    return true;
  }

  private checkOutput(path: ItemPath) {
    // TODO: check if has an active link
    return true;
  }
}


//
// Pipeline classes

export class CompositionPipelineView extends PipelineView {
  private hooks: HookSpec[] = [];
  private rt?: PipelineRuntime;

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
    await super.onBeforeLoadRun();
    await this.execHooks('beforeLoadRun');
  }

  override async onAfterLoadRun(run: DG.FuncCall) {
    await super.onAfterLoadRun(run);
    await this.execHooks('afterLoadRun', {run});
  }

  override build() {
    super.build();
    this.execHooks('onViewReady');
  }

  private async execHooks(category: keyof PipelineHooks, additionalParams: Record<string, any> = {}) {
    for (const {hooks, pipelinePath} of this.hooks!) {
      const items = hooks[category];
      for (const item of items ?? []) {
        const handler = item.handler!;
        const ctrlConf = new ControllerConfig(pipelinePath, item.from, item.to);
        const controller = new RuntimeControllerImpl(ctrlConf, this.rt!);
        const params = {...additionalParams, controller};
        await callHandler(handler, params);
      }
    }
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
    items: PipelineConfigVariants[],
    config: PipelineCompositionConfiguration,
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
    this.rt = new PipelineRuntime(this.nodes, this.links, this.viewInst);
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
    for (const stepConf of this.processSteps(pipelineConf.steps, subPath, toRemove, toAdd)) {
      const sKey = stepConf.id;
      const popups: PipelinePopupConfiguration[] = [];
      for (const popupConf of this.getAllConfItems(stepConf.popups ?? [], sKey, toAdd, 'popup')) {
        const pPopupConf = this.processNodeConfig(popupConf, keyToPath(sKey), subPath, toRemove, 'action');
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

  private processSteps(data: PipelineStepConfiguration[], pipelinePath: ItemPath, toRemove: Set<string>, toAdd: Map<string, ItemsToMerge>) {
    const nData: PipelineStepConfiguration[] = data.flatMap((conf) => {
      const pconf = this.processNodeConfig<PipelineStepConfiguration>(conf, pipelinePath, pipelinePath, toRemove, 'step');
      // TODO: fix prev step suffix
      const addConfs = (toAdd.get(conf.id)?.step ?? [])
        .map((step) => this.processNodeConfig(step, pipelinePath, pipelinePath, toRemove, 'step'))
        .filter((x) => x) as PipelineStepConfiguration[];
      return [...(pconf ? [pconf] : []), ...addConfs];
    });
    return nData;
  }

  private getAllConfItems<T extends AddNodeConf>(data: T[], prefix: string,
    toAdd: Map<string, ItemsToMerge>, type: AddNodeConfTypes,
  ) {
    const moreItems = toAdd.get(prefix)?.[type] ?? [];
    return [...data, ...moreItems] as T[];
  }


  private processActionNodeConfig(conf: ActionConfiguraion, path: ItemPath, pipelinePath: ItemPath, toRemove: Set<string>) {
    const nodePath = this.updateFullPathLink(conf, path);
    const nodeKey = pathToKey(nodePath);
    if (toRemove.has(nodeKey))
      return;

    this.nodes.set(nodeKey, new NodeState(conf, 'action', pipelinePath));
    return conf;
  }

  private processNodeConfig<T extends NodeConf>(conf: T, path: ItemPath, pipelinePath: ItemPath, toRemove: Set<string>, type: NodeConfTypes) {
    const nodePath = this.updateFullPathNode(conf, path);
    const nodeKey = pathToKey(nodePath);
    if (toRemove.has(nodeKey))
      return;

    this.nodes.set(nodeKey, new NodeState(conf, type, pipelinePath));
    return conf;
  }

  private processLinkConfig(conf: PipelineLinkConfiguration, pipelinePath: ItemPath, toRemove: Set<string>) {
    const linkPath= this.updateFullPathLink(conf, pipelinePath);
    const linkKey = pathToKey(linkPath);
    if (toRemove.has(linkKey))
      return;

    this.links.set(linkKey, new LinkState(this.pipelineState, conf, pipelinePath));
    return conf;
  }

  private updateFullPathLink(target: { id: string, from?: ItemPath | ItemPath[], to?: ItemPath | ItemPath[] }, currentPath: ItemPath) {
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

    return this.updateFullPathNode(target, currentPath);
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
    this.processMergeNodeList(stepsToAdd ?? [], path, toAdd, 'step');
    this.processMergeNodeList(popupsToAdd ?? [], path, toAdd, 'popup');
    this.processMergeNodeList(actionsToAdd ?? [], path, toAdd, 'action');
    delete mergeConf.stepsToAdd;
    delete mergeConf.popupsToAdd;
    delete mergeConf.actionsToAdd;
  }

  private processMergeNodeList(itemsToAdd: [AddNodeConf, ItemPath?][], pipelinePath: ItemPath,
    toAdd: Map<string, ItemsToMerge>, type: AddNodeConfTypes,
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

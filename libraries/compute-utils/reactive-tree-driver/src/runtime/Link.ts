import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {HandlerBase} from '../config/PipelineConfiguration';
import {BaseTree, NodePath, NodePathSegment, TreeNode} from '../data/BaseTree';
import {descriptionOutputs, isFuncCallNode, StateTreeNode} from './StateTreeNodes';
import {ActionSpec, MatchedIO, MatchedNodePaths, MatchInfo} from './link-matching';
import {BehaviorSubject, combineLatest, defer, EMPTY, merge, Subject, of, asapScheduler} from 'rxjs';
import {map, filter, takeUntil, withLatestFrom, switchMap, catchError, mapTo, finalize, debounceTime, timestamp, distinctUntilChanged, take} from 'rxjs/operators';
import {callHandler, indexFromEnd} from '../utils';
import {defaultLinkHandler, Slot} from './default-handler';
import {ControllerCancelled, FuncallActionController, LinkController, MetaController, MutationController, NodeMetaController, PipelineValidatorController, RuntimeReturnController, ValidatorController} from './LinkControllers';
import {TemplateInfo} from '../RuntimeControllers';
import {LinkIOParsed, LinkSelectorSegment} from '../config/LinkSpec';
import {FuncCallAdapter, FuncCallMockAdapter, MemoryStore} from './FuncCallAdapters';
import {LinksState} from './LinksState';
import {PipelineInstanceConfig} from '../config/PipelineInstance';
import {GranularMutationOp} from '../data/common-types';
import {DriverLogger, reportError} from '../data/Logger';
import {FuncCallInstancesBridge} from './FuncCallInstancesBridge';
import {toStateRec} from './StateTreeSerializer';
import {PipelineNodeBase} from './StateTreeNodes';
import {PipelineOutline} from '../config/PipelineInstance';

const VALIDATOR_DEBOUNCE_TIME = 250;

export type ScopeInfo = {
  additionalParams?: Record<string, any>
}

export class Link {
  protected destroyed$ = new Subject<true>();
  protected isActive$ = new BehaviorSubject(false);
  protected trigger$ = new Subject<ScopeInfo>();

  public uuid = uuidv4();
  public readonly isValidator = this.matchInfo.spec.type === 'validator';
  public readonly isPipelineValidator = this.matchInfo.spec.type === 'pipelineValidator';
  public readonly isMeta = this.matchInfo.spec.type === 'meta';
  public readonly isMutation = this.matchInfo.spec.type === 'pipeline';
  public readonly isNodeMeta = this.matchInfo.spec.type === 'nodemeta' || this.matchInfo.spec.type === 'selector';
  public readonly isFuncallAction = this.matchInfo.spec.type === 'funccall';
  public readonly isReturn = this.matchInfo.spec.type === 'return';
  public readonly isBatchable: boolean;

  // probably a better api
  public lastPipelineMutations?: {
    path: NodePath,
    initConfig: PipelineInstanceConfig,
  }[];
  public lastGranularMutations?: {
    path: NodePath,
    ops: GranularMutationOp[],
  }[];
  public returnResult: any;

  private nextScheduled$ = new BehaviorSubject(-1);
  private lastFinished$ = new BehaviorSubject(-1);
  public isRunning$ = combineLatest([this.nextScheduled$, this.lastFinished$]).pipe(
    map(([next, last]) => next > last),
    distinctUntilChanged(),
  );

  private isWired = false;

  constructor(
    public prefix: NodePath,
    public matchInfo: MatchInfo,
    private customDebounceTime?: number,
    private logger?: DriverLogger,
  ) {
    const spec = this.matchInfo.spec;
    if (spec.type === 'validator' || spec.type === 'pipelineValidator') {
      const effectiveDebounce = this.customDebounceTime ?? VALIDATOR_DEBOUNCE_TIME;
      this.isBatchable = effectiveDebounce === 0 && !spec.sequential;
    } else if (spec.type === 'meta' || spec.type === 'nodemeta' || spec.type === 'selector')
      this.isBatchable = !spec.sequential;
    else
      this.isBatchable = false;
  }

  wire(state: BaseTree<StateTreeNode>, linksState?: LinksState) {
    if (this.isWired)
      return;

    const {slots: inputSlots, templates: inputTemplates} =
      this.buildSlotsAndTemplates(this.matchInfo.spec.from ?? [], this.matchInfo.inputs);
    const {slots: outputSlots, templates: outputTemplates} =
      this.buildSlotsAndTemplates(this.matchInfo.spec.to ?? [], this.matchInfo.outputs);
    const inputSet = new Set(this.getOrderedIO(this.matchInfo.inputs));
    const outputSet = new Set(this.getOrderedIO(this.matchInfo.outputs));
    const callInputs = new Set(
      (this.matchInfo.spec.from ?? []).filter((io) => io.flags?.includes('call')).map((io) => io.name),
    );

    const inputsChanges$ = this.makeInputsChanges(state, linksState);
    const baseNode = this.matchInfo.basePath ?
      state.getNode(([...this.prefix, ...this.matchInfo.basePath])) :
      undefined;

    const actions: Record<string, Map<string, string>> = {};

    if (linksState) {
      for (const [name, minfos] of Object.entries(this.matchInfo.actions)) {
        if (minfos.length > 1)
          reportError('warning', `link:${this.matchInfo.spec.id}`, `Multiple action nodes with the same name ${name}`, this.logger, [this.matchInfo.spec.id]);
        const nodeActions = minfos.map((minfo) => {
          const node = state.getNode([...this.prefix, ...minfo.path]);
          const actions = linksState.nodesActions.get(node.getItem().uuid) ?? [];
          return new Map(actions.map((action) => [action.matchInfo.spec.id, action.uuid]));
        })[0];
        actions[name] = nodeActions;
      }
    }

    inputsChanges$.pipe(
      timestamp(),
      map(({timestamp}) => timestamp),
      takeUntil(this.destroyed$),
    ).subscribe(this.nextScheduled$);

    const actionsVisibility = linksState?.actionsVisibility ?? new Map<string, boolean>();

    inputsChanges$.pipe(
      switchMap(
        ([scope, inputs]) =>
          this.runHandler(inputs, inputSet, outputSet, callInputs, inputSlots, outputSlots, inputTemplates, outputTemplates, actions, actionsVisibility, baseNode, scope, state).pipe(
            map((controller) => this.setHandlerResults(controller, state)),
            catchError((error) => {
              reportError('recoverable', `link:${this.matchInfo.spec.id}`, error, this.logger, [this.matchInfo.spec.id]);
              return EMPTY;
            }),
          ),
      ),
      timestamp(),
      map(({timestamp}) => timestamp),
      takeUntil(this.destroyed$),
    ).subscribe(this.lastFinished$);

    this.isWired = true;
  }

  setActive() {
    this.isActive$.next(true);
  }

  setInactive() {
    this.isActive$.next(false);
  }

  trigger(scope?: ScopeInfo) {
    this.trigger$.next(scope);
  }

  destroy() {
    this.destroyed$.next(true);
  }

  private makeInputsChanges(state: BaseTree<StateTreeNode>, linksState?: LinksState) {
    const inputs = Object.entries(this.matchInfo.inputs).map(([inputAlias, inputItems]) => {
      const nodes = inputItems.map(
        (input) => [
          input.ioName!,
          state.getNode([...this.prefix, ...input.path]),
        ] as const,
      );
      const inputStates = nodes.map(([ioName, node]) => {
        const item = node.getItem();
        const {dataFrameMutations} = this.matchInfo.spec;
        const includeDFMutations = Array.isArray(dataFrameMutations) ?
          dataFrameMutations.includes(inputAlias) :
          !!dataFrameMutations;
        const store = item.getStateStore();
        const state$ = ioName ? store.getStateChanges(ioName, includeDFMutations) :
          (store instanceof FuncCallInstancesBridge ? store.instance$.pipe(map((x) => {
            const adapter = x?.adapter;
            return (adapter && !(adapter instanceof FuncCallMockAdapter)) ? adapter.getFuncCall() : undefined;
          })) : of(undefined));
        return state$;
      });
      return [inputAlias, combineLatest(inputStates)] as const;
    });

    const inputEntries = inputs.map(
      ([name, values$]) => values$.pipe(map((val) => [name, val] as const)));

    const inputsEntries$ = combineLatest(inputEntries);

    const inputsTriggered$ = inputEntries.length ?
      this.trigger$.pipe(withLatestFrom(inputsEntries$)) :
      this.trigger$.pipe(map((scope) => [scope, [] as any[]] as const));

    const toInputsChanges = ([scope, entries]: readonly [any, any]) =>
      [scope as ScopeInfo | undefined, Object.fromEntries(entries) as Record<string, any>] as const;

    if (this.isBatchable && linksState?.batchLinks) {
      inputsEntries$.pipe(
        filter(() => this.isActive$.value),
        takeUntil(this.destroyed$),
      ).subscribe(() => linksState.scheduleBatch(this.uuid));

      return inputsTriggered$.pipe(map(toInputsChanges));
    }

    const debounceVal = this.customDebounceTime ?? ((this.isValidator || this.isPipelineValidator) ? VALIDATOR_DEBOUNCE_TIME : 0);

    const activeInputs$ = inputsEntries$.pipe(
      filter(() => this.isActive$.value),
      debounceTime(debounceVal, debounceVal === 0 ? asapScheduler : undefined),
      map((obs) => [undefined, obs] as const),
    );

    return merge(activeInputs$, inputsTriggered$).pipe(map(toInputsChanges));
  }

  private runHandler(
    inputs: Record<string, any>,
    inputSet: Set<string>,
    outputSet: Set<string>,
    callInputs: Set<string>,
    inputSlots: Slot[],
    outputSlots: Slot[],
    inputTemplates: TemplateInfo[],
    outputTemplates: TemplateInfo[],
    actions: Record<string, Map<string, string>>,
    actionsVisibility: ReadonlyMap<string, boolean>,
    baseNode?: TreeNode<StateTreeNode>,
    scope?: ScopeInfo,
    state?: BaseTree<StateTreeNode>,
  ) {
    const controller = this.getControllerInstance(inputs, inputSet, outputSet, callInputs, inputTemplates, outputTemplates, actions, actionsVisibility, baseNode, scope, state);
    if (this.logger) {
      this.logger.logLink('linkRunStarted', {
        prefix: this.prefix,
        id: this.matchInfo.spec.id,
        linkUUID: this.uuid,
        basePath: this.matchInfo.basePath,
        isDefaultValidator: this.matchInfo.isDefaultValidator,
      });
    }
    if (this.matchInfo.spec.handler) {
      return callHandler(this.matchInfo.spec.handler as HandlerBase<any, void>, {controller}).pipe(
        catchError((e) => {
          if (e instanceof ControllerCancelled)
            return EMPTY;
          throw e;
        }),
        mapTo(controller),
        finalize(() => controller.close()),
      );
    } else if (!this.isValidator && !this.isMeta) {
      return defer(() => of(
        defaultLinkHandler(controller as LinkController, inputSlots, outputSlots, this.matchInfo.spec.defaultRestrictions),
      ).pipe(mapTo(controller)));
    }

    throw Error(`Unable to run handler for link ${this.matchInfo.spec.id}`);
  }

  private resolveOutputNodes(state: BaseTree<StateTreeNode>): Record<string, {node: TreeNode<StateTreeNode>, path: NodePath}[]> {
    const result: Record<string, {node: TreeNode<StateTreeNode>, path: NodePath}[]> = {};
    for (const [outputAlias, outputItems] of Object.entries(this.matchInfo.outputs)) {
      result[outputAlias] = outputItems.map((output) => {
        const path = [...this.prefix, ...output.path];
        return {node: state.getNode(path), path};
      });
    }
    return result;
  }

  private getControllerInstance(
    inputs: Record<string, any>,
    inputSet: Set<string>,
    outputSet: Set<string>,
    callInputs: Set<string>,
    inputTemplates: TemplateInfo[],
    outputTemplates: TemplateInfo[],
    actions: Record<string, Map<string, string>>,
    actionsVisibility: ReadonlyMap<string, boolean>,
    baseNode?: TreeNode<StateTreeNode>,
    scope?: ScopeInfo,
    state?: BaseTree<StateTreeNode>,
  ) {
    const baseArgs = {inputs, inputsSet: inputSet, outputsSet: outputSet, callInputs, id: this.matchInfo.spec.id, scopeInfo: scope, inputTemplates, outputTemplates};

    if (this.isValidator)
      return new ValidatorController({...baseArgs, actions, actionsVisibility, baseNode});

    if (this.isPipelineValidator) {
      const outline = this.resolveOutline(state);
      return new PipelineValidatorController({...baseArgs, outline});
    }

    if (this.isMeta)
      return new MetaController(baseArgs);

    if (this.isMutation) {
      const outputNodes = state ? this.resolveOutputNodes(state) : {};
      return new MutationController({...baseArgs, outputNodes});
    }

    if (this.isNodeMeta)
      return new NodeMetaController({...baseArgs, outputsSet: new Set(descriptionOutputs)});

    if (this.isFuncallAction)
      return new FuncallActionController(baseArgs);

    if (this.isReturn)
      return new RuntimeReturnController(baseArgs);

    return new LinkController(baseArgs);
  }

  private resolveOutline(state?: BaseTree<StateTreeNode>): PipelineOutline {
    if (!state)
      throw new Error(`Link ${this.matchInfo.spec.id}: pipeline validator has no tree state`);
    return toStateRec(state.getNode(this.prefix), {skipFuncCalls: true}) as PipelineOutline;
  }

  private buildSlotsAndTemplates(
    specSide: LinkIOParsed[],
    matchSide: Record<string, MatchedNodePaths>,
  ): {slots: Slot[], templates: TemplateInfo[]} {
    const slots: Slot[] = [];
    const templates: TemplateInfo[] = [];
    const templateById = new Map<string | number, TemplateInfo>();
    for (const io of specSide) {
      if (matchSide[io.name] == null) continue;
      const tplId = io.templateName;
      if (tplId != null) {
        let tpl = templateById.get(tplId);
        if (!tpl) {
          tpl = {name: tplId, ios: []};
          templateById.set(tplId, tpl);
          templates.push(tpl);
          slots.push({kind: 'template', name: tplId});
        }
        const lastSeg = indexFromEnd(io.segments) as LinkSelectorSegment;
        tpl.ios.push({ioName: io.name, scriptIoId: lastSeg.ids[0]});
      } else {
        slots.push({kind: 'bare', name: io.name});
      }
    }
    return {slots, templates};
  }

  private getOrderedIO(ioData: Record<string, MatchedNodePaths>) {
    return Object.entries(ioData).sort(([, v1], [, v2]) => {
      const p1 = this.getFirstMatch(v1);
      const p2 = this.getFirstMatch(v2);
      return BaseTree.compareAddresses(p1, p2);
    }).map(([k]) => k);
  }

  private getFirstMatch(matchIO: readonly MatchedIO[]) {
    const p0 = matchIO.reduce((acc, val) => {
      const d = BaseTree.compareAddresses(acc, val.path);
      return d < 0 ? val.path : acc;
    }, [{idx: Infinity, id: ''}] as readonly NodePathSegment[] );
    return p0;
  }

  private setHandlerResults(controller: LinkController | ValidatorController | MetaController | MutationController | NodeMetaController | FuncallActionController | RuntimeReturnController | PipelineValidatorController, state: BaseTree<StateTreeNode>) {
    this.lastPipelineMutations = [];
    this.lastGranularMutations = [];
    if (this.logger) {
      this.logger.logLink('linkRunFinished', {
        prefix: this.prefix,
        id: this.matchInfo.spec.id,
        linkUUID: this.uuid,
        basePath: this.matchInfo.basePath,
        isDefaultValidator: this.matchInfo.isDefaultValidator,
      });
    }
    const outputsEntries = Object.entries(this.matchInfo.outputs).map(([outputAlias, outputItems]) => {
      const nodes = outputItems.map((output) => {
        const path = [...this.prefix, ...output.path];
        return [output.ioName!, path, state.getNode(path)] as const;
      });
      return [outputAlias, nodes] as const;
    });
    if (controller instanceof RuntimeReturnController) {
      this.returnResult = controller.result;
      return;
    }
    if (controller instanceof PipelineValidatorController) {
      for (const matched of Object.values(this.matchInfo.outputs)) {
        for (const target of matched) {
          const item = state.getNode([...this.prefix, ...target.path]).getItem();
          if (item instanceof PipelineNodeBase)
            item.setPipelineValidation(this.uuid, controller.output);
          else
            reportError('warning', `link:${this.matchInfo.spec.id}`, `pipelineValidator \`to\` target ${item.uuid} is not a pipeline node — skipped`, this.logger, [this.matchInfo.spec.id]);
        }
      }
      return;
    }
    for (const [outputAlias, nodesData] of outputsEntries) {
      for (const [ioName, nodePath, node] of nodesData) {
        if (controller instanceof ValidatorController) {
          const store = node.getItem().getStateStore();
          if (store instanceof MemoryStore)
            throw new Error(`Unable to set validations to a raw memory store ${node.getItem().uuid}`);
          store.setValidation(ioName, this.uuid, controller.outputs[outputAlias]);
        } else if (controller instanceof MetaController) {
          const store = node.getItem().getStateStore();
          if (store instanceof MemoryStore)
            throw new Error(`Unable to set meta to a raw memory store ${node.getItem().uuid}`);
          store.setMeta(ioName, this.uuid, controller.outputs[outputAlias]);
        } else if (controller instanceof MutationController) {
          const initConfig = controller.outputs[outputAlias];
          if (initConfig)
            this.lastPipelineMutations.push({path: nodePath, initConfig});
          const ops = controller.granularOps[outputAlias];
          if (ops?.length)
            this.lastGranularMutations!.push({path: nodePath, ops});
        } else if (controller instanceof NodeMetaController) {
          const data = controller.outputs[outputAlias];
          const descrStore = node.getItem().nodeDescription;
          if (ioName === 'tags') {
            const ndata = {[this.uuid]: data};
            const odata = descrStore.getState('tags') ?? {};
            descrStore.setState(ioName, {...odata, ...ndata});
          } else
            descrStore.setState(ioName, data);
        } else if (controller instanceof FuncallActionController) {
          const data = controller.outputs[outputAlias];
          const stateStore = node.getItem().getStateStore();
          if (stateStore instanceof FuncCallInstancesBridge && !stateStore.isReadonly && data != null) {
            const adapter = new FuncCallAdapter(data, false);
            stateStore.change(adapter, true);
          }
        } else {
          const data = controller.outputs[outputAlias];
          if (data) {
            const [state, restriction] = data;
            const item = node.getItem();
            const isDf = state instanceof DG.DataFrame;
            const nextValue = isDf ? state.clone() : state;
            if (isDf)
              (nextValue as DG.DataFrame).id = uuidv4();
            if (isFuncCallNode(item) && !item.instancesWrapper.isOutputOutdated$.value)
              item.instancesWrapper.setRestriction(ioName, nextValue, restriction);
            else
              node.getItem().getStateStore().setState(ioName, nextValue, restriction);
          }
          if (controller instanceof LinkController && controller.consistencyResets.has(outputAlias)) {
            const item = node.getItem();
            if (isFuncCallNode(item))
              item.instancesWrapper.removeRestriction(ioName);
          }
        }
      }
    }
    if (controller instanceof MutationController)
      checkMutationOverlap(this.matchInfo.spec.id, this.lastPipelineMutations, this.lastGranularMutations);
  }
}

function checkMutationOverlap(
  actionId: string,
  pipelineMutations?: Array<{path: NodePath}>,
  granularMutations?: Array<{path: NodePath}>,
) {
  const all: NodePath[] = [...(pipelineMutations ?? []).map((m) => m.path), ...(granularMutations ?? []).map((m) => m.path)];
  for (let i = 0; i < all.length; i++) {
    for (let j = i + 1; j < all.length; j++) {
      if (BaseTree.isNodeChildOrEq(all[i], all[j]) || BaseTree.isNodeChildOrEq(all[j], all[i]))
        throw new Error(`Handler for action ${actionId}: mutations on overlapping subtrees (${JSON.stringify(all[i])} and ${JSON.stringify(all[j])}). Combine into a single setPipelineState at the common ancestor, or operate only at the inner level.`);
    }
  }
}

export class Action extends Link {
  constructor(
    public prefix: NodePath,
    public matchInfo: MatchInfo,
    public spec: ActionSpec,
    logger?: DriverLogger,
  ) {
    super(prefix, matchInfo, undefined, logger);
  }

  exec(additionalParams: Record<string, any>) {
    this.trigger$.next({additionalParams});
    return this.isRunning$.pipe(
      filter((isRunning) => !isRunning),
      take(1),
      mapTo(undefined),
    );
  }

  execPipelineMutations(additionalParams: Record<string, any>) {
    this.trigger$.next({additionalParams});
    return this.isRunning$.pipe(
      filter((isRunning) => !isRunning),
      take(1),
      map(() => ({
        pipelineMutations: this.lastPipelineMutations,
        granularMutations: this.lastGranularMutations,
      })),
    );
  }
}

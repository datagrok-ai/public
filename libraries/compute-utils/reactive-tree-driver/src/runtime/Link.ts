import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {HandlerBase} from '../config/PipelineConfiguration';
import {BaseTree, NodePath, TreeNode} from '../data/BaseTree';
import {descriptionOutputs, isFuncCallNode, StateTreeNode} from './StateTreeNodes';
import {ActionSpec, MatchInfo} from './link-matching';
import {BehaviorSubject, combineLatest, defer, EMPTY, merge, Subject, of} from 'rxjs';
import {map, filter, takeUntil, withLatestFrom, switchMap, catchError, mapTo, finalize, debounceTime, timestamp, distinctUntilChanged, take} from 'rxjs/operators';
import {callHandler, pathToUUID} from '../utils';
import {defaultLinkHandler} from './default-handler';
import {ControllerCancelled, FuncallActionController, LinkController, MetaController, MutationController, NameSelectorController, ValidatorController} from './LinkControllers';
import {FuncCallAdapter, MemoryStore} from './FuncCallAdapters';
import {LinksState} from './LinksState';
import {PipelineInstanceConfig} from '../config/PipelineInstance';
import {DriverLogger} from '../data/Logger';
import {FuncCallInstancesBridge} from './FuncCallInstancesBridge';

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
  public readonly isMeta = this.matchInfo.spec.type === 'meta';
  public readonly isMutation = this.matchInfo.spec.type === 'pipeline';
  public readonly isSelector = this.matchInfo.spec.type === 'selector';
  public readonly isFuncallAction = this.matchInfo.spec.type === 'funccall';

  // probably a better api
  public lastPipelineMutations?: {
    path: NodePath,
    initConfig: PipelineInstanceConfig,
  }[];

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
  ) {}

  wire(state: BaseTree<StateTreeNode>, linksState?: LinksState) {
    if (this.isWired)
      return;

    const inputNames = Object.keys(this.matchInfo.inputs);
    const outputNames = Object.keys(this.matchInfo.outputs);
    const inputSet = new Set(inputNames);
    const outputSet = new Set(outputNames);

    const inputsChanges$ = this.makeInputsChanges(state);
    const baseNode = this.matchInfo.basePath ?
      state.getNode(([...this.prefix, ...this.matchInfo.basePath])) :
      undefined;

    const actions: Record<string, Map<string, string>> = {};

    if (linksState) {
      for (const [name, minfos] of Object.entries(this.matchInfo.actions)) {
        if (minfos.length > 1)
          grok.shell.warning(`Node ${this.matchInfo.spec.id} prefix ${this.prefix} multiple action nodes with the same name ${name}`);
        const nodeActions = minfos.map((minfo) => {
          const node = state.getNode(minfo.path);
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

    inputsChanges$.pipe(
      switchMap(
        ([scope, inputs]) =>
          this.runHandler(inputs, inputSet, outputSet, inputNames, outputNames, actions, baseNode, scope),
      ),
      map((controller) => this.setHandlerResults(controller, state)),
      catchError((error) => {
        console.error(error);
        grok.shell.error(error?.message);
        return EMPTY;
      }),
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

  private makeInputsChanges(state: BaseTree<StateTreeNode>) {
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
          (store instanceof FuncCallInstancesBridge ? store.instance$.pipe(map((x) => x?.adapter.getFuncCall())) : of(undefined));
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

    const debounceVal = this.customDebounceTime ?? (this.isValidator ? VALIDATOR_DEBOUNCE_TIME : 0);

    const activeInputs$ = inputsEntries$.pipe(
      filter(() => this.isActive$.value),
      debounceTime(debounceVal),
      map((obs) => [undefined, obs] as const),
    );

    const inputsChanges$ = merge(activeInputs$, inputsTriggered$).pipe(
      map(([scope, entries]) => [scope, Object.fromEntries(entries)] as const),
    );
    return inputsChanges$;
  }

  private runHandler(
    inputs: Record<string, any>,
    inputSet: Set<string>,
    outputSet: Set<string>,
    inputNames: string[],
    outputNames: string[],
    actions: Record<string, Map<string, string>>,
    baseNode?: TreeNode<StateTreeNode>,
    scope?: ScopeInfo,
  ) {
    const controller = this.getControllerInstance(inputs, inputSet, outputSet, actions, baseNode, scope);
    if (this.logger) {
      this.logger.logLink('linkRunStarted', {
        prefix: this.prefix,
        id: this.matchInfo.spec.id,
        linkUUID: this.uuid,
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
        defaultLinkHandler(controller as LinkController, inputNames, outputNames, this.matchInfo.spec.defaultRestrictions),
      ).pipe(mapTo(controller)));
    }

    throw Error(`Unable to run handler for link ${this.matchInfo.spec.id}`);
  }

  private getControllerInstance(
    inputs: Record<string, any>,
    inputSet: Set<string>,
    outputSet: Set<string>,
    actions: Record<string, Map<string, string>>,
    baseNode?: TreeNode<StateTreeNode>,
    scope?: ScopeInfo,
  ) {
    if (this.isValidator)
      return new ValidatorController(inputs, inputSet, outputSet, this.matchInfo.spec.id, actions, baseNode, scope);

    if (this.isMeta)
      return new MetaController(inputs, inputSet, outputSet, this.matchInfo.spec.id, scope);

    if (this.isMutation)
      return new MutationController(inputs, inputSet, outputSet, this.matchInfo.spec.id, scope);

    if (this.isSelector)
      return new NameSelectorController(inputs, inputSet, new Set(descriptionOutputs), this.matchInfo.spec.id, scope);

    if (this.isFuncallAction)
      return new FuncallActionController(inputs, inputSet, outputSet, this.matchInfo.spec.id, scope);

    return new LinkController(inputs, inputSet, outputSet, this.matchInfo.spec.id, scope);
  }

  private setHandlerResults(controller: LinkController | ValidatorController | MetaController | MutationController | NameSelectorController | FuncallActionController, state: BaseTree<StateTreeNode>) {
    this.lastPipelineMutations = [];
    if (this.logger) {
      this.logger.logLink('linkRunFinished', {
        prefix: this.prefix,
        id: this.matchInfo.spec.id,
        linkUUID: this.uuid,
      });
    }
    const outputsEntries = Object.entries(this.matchInfo.outputs).map(([outputAlias, outputItems]) => {
      const nodes = outputItems.map((output) => {
        const path = [...this.prefix, ...output.path];
        return [output.ioName!, path, state.getNode(path)] as const;
      });
      return [outputAlias, nodes] as const;
    });
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
          store.setMeta(ioName, controller.outputs[outputAlias]);
        } else if (controller instanceof MutationController) {
          const initConfig = controller.outputs[outputAlias];
          if (initConfig)
            this.lastPipelineMutations.push({path: nodePath, initConfig});
        } else if (controller instanceof NameSelectorController) {
          const data = controller.outputs[outputAlias];
          const descrStore = node.getItem().nodeDescription;
          if (ioName === 'tags') {
            const ndata = { [this.uuid]: data };
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
            const nextValue = state instanceof DG.DataFrame ? state.clone() : state;
            if (isFuncCallNode(item) && !item.instancesWrapper.isOutputOutdated$.value)
              item.instancesWrapper.setRestriction(ioName, nextValue, restriction);
            else
              node.getItem().getStateStore().setState(ioName, nextValue, restriction);
          }
        }
      }
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
      map(() => this.lastPipelineMutations),
    );
  }
}

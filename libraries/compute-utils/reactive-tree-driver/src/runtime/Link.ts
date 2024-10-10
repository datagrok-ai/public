import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {ActionPositions, HandlerBase} from '../config/PipelineConfiguration';
import {BaseTree, NodeAddress, NodePath, TreeNode} from '../data/BaseTree';
import {StateTreeNode} from './StateTreeNodes';
import {MatchInfo} from './link-matching';
import {BehaviorSubject, combineLatest, defer, EMPTY, merge, Subject, of} from 'rxjs';
import {map, filter, takeUntil, withLatestFrom, switchMap, catchError, mapTo, finalize, debounceTime, timestamp, distinctUntilChanged} from 'rxjs/operators';
import {callHandler} from '../utils';
import {defaultLinkHandler} from './default-handler';
import {ControllerCancelled, LinkController, MetaController, MutationController, ValidatorController} from './LinkControllers';
import {MemoryStore} from './FuncCallAdapters';
import {LinksState} from './LinksState';
import {PipelineInstanceConfig} from '../config/PipelineInstance';

const VALIDATOR_DEBOUNCE_TIME = 250;

export type ScopeInfo = {
  scope?: NodeAddress, childOffset?: number
}

export class Link {
  private destroyed$ = new Subject<true>();
  private isActive$ = new BehaviorSubject(false);
  private trigger$ = new Subject<ScopeInfo>();

  public uuid = uuidv4();
  public readonly isValidator = !!this.matchInfo.spec.isValidator;
  public readonly isMeta = !!this.matchInfo.spec.isMeta;
  public readonly isMutation = !!this.matchInfo.spec.isPipeline;

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

  constructor(public prefix: NodePath, public matchInfo: MatchInfo) {}

  wire(state: BaseTree<StateTreeNode>, linksState?: LinksState) {
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
  }

  setActive() {
    this.isActive$.next(true);
  }

  trigger(scope?: NodeAddress, childOffset?: number) {
    this.trigger$.next({scope, childOffset});
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
        const state$ = item.getStateStore().getStateChanges(ioName, includeDFMutations);
        return state$;
      });
      return [inputAlias, combineLatest(inputStates)] as const;
    });

    const inputEntries = inputs.map(
      ([name, values$]) => values$.pipe(map((val) => [name, val] as const)));

    const inputsEntries$ = combineLatest(inputEntries);

    const inputsTriggered$ = inputEntries.length ?
      this.trigger$.pipe(withLatestFrom(inputsEntries$)) :
      this.trigger$.pipe(map((scope) => [scope, [] as any[]] as const)) ;

    const activeInputs$ = inputsEntries$.pipe(
      filter(() => this.isActive$.value),
      debounceTime(this.isValidator ? VALIDATOR_DEBOUNCE_TIME : 0),
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

    return new LinkController(inputs, inputSet, outputSet, this.matchInfo.spec.id, scope);
  }

  private setHandlerResults(controller: LinkController | ValidatorController | MetaController | MutationController, state: BaseTree<StateTreeNode>) {
    this.lastPipelineMutations = [];
    const outputsEntries = Object.entries(this.matchInfo.outputs).map(([outputAlias, outputItems]) => {
      const nodes = outputItems.map((output) => {
        const path = [...this.prefix, ...output.path];
        return [output.ioName!, path, state.getNode(path)] as const;
      });
      return [outputAlias, nodes] as const;
    });
    for (const [outputAlias, nodesData] of outputsEntries) {
      const nextData = controller.outputs[outputAlias];
      if (nextData) {
        for (const [ioName, nodePath, node] of nodesData) {
          if (controller.scopeInfo?.scope &&
            !BaseTree.isNodeChildOffseted(controller.scopeInfo.scope, nodePath, controller.scopeInfo.childOffset)
          )
            continue;

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
          } else {
            const [state, restriction] = controller.outputs[outputAlias];
            const nextValue = state instanceof DG.DataFrame ? state.clone() : state;
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
    public position: ActionPositions,
    public isPipelineMutation: boolean,
    public friendlyName?: string,
    public menuCategory?: string,
  ) {
    super(prefix, matchInfo);
  }
}

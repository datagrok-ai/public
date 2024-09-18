import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {ActionPositions, HandlerBase} from '../config/PipelineConfiguration';
import {BaseTree, NodePath, TreeNode} from '../data/BaseTree';
import {StateTree} from './StateTree';
import {StateTreeNode} from './StateTreeNodes';
import {MatchInfo} from './link-matching';
import {BehaviorSubject, combineLatest, defer, EMPTY, merge, Subject, of} from 'rxjs';
import {map, filter, takeUntil, withLatestFrom, switchMap, catchError, mapTo, finalize, debounceTime, timestamp, distinctUntilChanged} from 'rxjs/operators';
import {callHandler} from '../utils';
import {defaultLinkHandler} from './default-handler';
import {ControllerCancelled, LinkController, MetaController, ValidatorController} from './LinkControllers';
import {MemoryStore} from './FuncCallAdapters';
import {LinksState} from './LinksState';

const VALIDATOR_DEBOUNCE_TIME = 250;

export type ScopeInfo = {
  scope?: NodePath, childOffset?: number
}

export class Link {
  private destroyed$ = new Subject<true>();
  private isActive$ = new BehaviorSubject(false);
  private trigger$ = new Subject<ScopeInfo>();

  public uuid = uuidv4();
  public readonly isValidator = !!this.matchInfo.spec.isValidator;
  public readonly isMeta = !!this.matchInfo.spec.isMeta;
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

    inputsChanges$.pipe(
      timestamp(),
      map(({timestamp}) => timestamp),
      takeUntil(this.destroyed$),
    ).subscribe(this.nextScheduled$);

    const actions = ((this.isValidator && baseNode && linksState) &&
      linksState.baseNodeActions.get(baseNode.getItem().uuid)) ||
      [];

    inputsChanges$.pipe(
      switchMap(
        ([scope, inputs]) =>
          this.runHandler(actions, inputs, inputSet, outputSet, inputNames, outputNames, baseNode, scope),
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

  trigger(scope?: NodePath, childOffset?: number) {
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

    const inputsTriggered$ = this.trigger$.pipe(
      withLatestFrom(inputsEntries$),
    );

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
    actions: Action[],
    inputs: Record<string, any>,
    inputSet: Set<string>,
    outputSet: Set<string>,
    inputNames: string[],
    outputNames: string[],
    baseNode?: TreeNode<StateTreeNode>,
    scope?: ScopeInfo,
  ) {
    const controller = this.getControllerInstance(actions, inputs, inputSet, outputSet, baseNode, scope);

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
    actions: Action[],
    inputs: Record<string, any>,
    inputSet: Set<string>,
    outputSet: Set<string>,
    baseNode?: TreeNode<StateTreeNode>,
    scope?: ScopeInfo,
  ) {
    if (this.isValidator)
      return new ValidatorController(inputs, inputSet, outputSet, this.matchInfo.spec.id, actions, baseNode, scope);

    if (this.isMeta)
      return new MetaController(inputs, inputSet, outputSet, this.matchInfo.spec.id, scope);
    return new LinkController(inputs, inputSet, outputSet, this.matchInfo.spec.id, scope);
  }

  private setHandlerResults(controller: LinkController | ValidatorController | MetaController, state: BaseTree<StateTreeNode>) {
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
              throw new Error(`Unable to set validations to raw memory store ${node.getItem().uuid}`);
            store.setValidation(ioName, this.uuid, controller.outputs[outputAlias]);
          } else if (controller instanceof MetaController) {
            const store = node.getItem().getStateStore();
            if (store instanceof MemoryStore)
              throw new Error(`Unable to set meta to raw memory store ${node.getItem().uuid}`);
            store.setMeta(ioName, controller.outputs[outputAlias]);
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
    public friendlyName?: string,
    public menuCategory?: string,
  ) {
    super(prefix, matchInfo);
  }
}

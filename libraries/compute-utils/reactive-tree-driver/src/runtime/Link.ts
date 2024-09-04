import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {HandlerBase} from '../config/PipelineConfiguration';
import {NodePath, TreeNode} from '../data/BaseTree';
import {StateTree} from './StateTree';
import {StateTreeNode} from './StateTreeNodes';
import {MatchInfo} from './link-matching';
import {BehaviorSubject, combineLatest, defer, EMPTY, merge, Subject, of} from 'rxjs';
import {map, filter, takeUntil, withLatestFrom, switchMap, catchError, mapTo, finalize, debounceTime, timestamp, distinctUntilChanged} from 'rxjs/operators';
import {callHandler} from '../utils';
import {defaultLinkHandler} from './default-handler';
import {ControllerCancelled, LinkController, ValidatorController} from './LinkControllers';

const VALIDATOR_DEBOUNCE_TIME = 250;

export class Link {
  private destroyed$ = new Subject<true>();
  private isActive$ = new BehaviorSubject(false);
  private trigger$ = new Subject<true>();

  public uuid = uuidv4();
  public readonly isValidator = !!this.matchInfo.spec.isValidator;
  private nextScheduled$ = new BehaviorSubject(-1);
  private lastFinished$ = new BehaviorSubject(-1);
  public isRunning$ = combineLatest([this.nextScheduled$, this.lastFinished$]).pipe(
    map(([next, last]) => next > last),
    distinctUntilChanged(),
  );

  constructor(public prefix: NodePath, public matchInfo: MatchInfo) {}

  wire(state: StateTree) {
    const inputNames = Object.keys(this.matchInfo.inputs);
    const outputNames = Object.keys(this.matchInfo.outputs);
    const inputSet = new Set(inputNames);
    const outputSet = new Set(outputNames);

    const inputsChanges$ = this.makeInputsChanges(state);
    const baseNode = this.matchInfo.basePath ? state.getNode(([...this.prefix, ...this.matchInfo.basePath])) : undefined;

    inputsChanges$.pipe(
      timestamp(),
      map(({timestamp}) => timestamp),
      takeUntil(this.destroyed$),
    ).subscribe(this.nextScheduled$);

    inputsChanges$.pipe(
      switchMap((inputs) => this.runHandler(inputs, inputSet, outputSet, inputNames, outputNames, baseNode)),
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

  enable() {
    this.isActive$.next(true);
  }

  disable() {
    this.isActive$.next(false);
  }

  trigger() {
    this.trigger$.next(true);
  }

  destroy() {
    this.destroyed$.next(true);
  }

  private makeInputsChanges(state: StateTree) {
    const inputs = Object.entries(this.matchInfo.inputs).map(([inputAlias, inputItems]) => {
      const nodes = inputItems.map((input) => [input.ioName!, state.getNode([...this.prefix, ...input.path])] as const);
      const inputStates = nodes.map(([ioName, node]) => {
        const item = node.getItem();
        const {dataFrameMutations} = this.matchInfo.spec;
        const includeDFMutations = Array.isArray(dataFrameMutations) ? dataFrameMutations.includes(inputAlias) : !!dataFrameMutations;
        const state$ = item.getStateStore().getStateChanges(ioName, includeDFMutations);
        return state$;
      });
      return [inputAlias, combineLatest(inputStates)] as const;
    });

    const inputEntries = inputs.map(([name, values$]) => values$.pipe(map((val) => [name, val] as const)));

    const inputsEntries$ = combineLatest(inputEntries);

    const inputsTriggered$ = this.trigger$.pipe(
      withLatestFrom(inputsEntries$),
      map(([, obs]) => obs),
    );

    const activeInputs$ = inputsEntries$.pipe(
      filter(() => this.isActive$.value),
      debounceTime(this.isValidator ? VALIDATOR_DEBOUNCE_TIME : 0),
    );

    const inputsChanges$ = merge(activeInputs$, inputsTriggered$).pipe(
      map((entries) => Object.fromEntries(entries)),
    );
    return inputsChanges$;
  }

  private runHandler(inputs: Record<string, any>, inputSet: Set<string>, outputSet: Set<string>, inputNames: string[], outputNames: string[], baseNode?: TreeNode<StateTreeNode>) {
    const controller = this.isValidator ?
      new ValidatorController(inputs, inputSet, outputSet, this.matchInfo.spec.id, baseNode) :
      new LinkController(inputs, inputSet, outputSet, this.matchInfo.spec.id);

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
    } else if (!this.isValidator)
      return defer(() => of(defaultLinkHandler(controller as LinkController, inputNames, outputNames, this.matchInfo.spec.defaultRestrictions)).pipe(mapTo(controller)));

    throw Error(`Unable to run handler for link ${this.matchInfo.spec.id}`);
  }

  private setHandlerResults(controller: LinkController | ValidatorController, state: StateTree) {
    const outputsEntries = Object.entries(this.matchInfo.outputs).map(([outputAlias, outputItems]) => {
      const nodes = outputItems.map((output) => [output.ioName!, state.getNode([...this.prefix, ...output.path])] as const);
      return [outputAlias, nodes] as const;
    });
    for (const [outputAlias, nodesData] of outputsEntries) {
      const nextData = controller.outputs[outputAlias];
      if (nextData) {
        for (const [ioName, node] of nodesData) {
          if (controller instanceof ValidatorController)
            node.getItem().getStateStore().setValidation(ioName, this.uuid, controller.outputs[outputAlias]);
          else {
            const [state, restriction] = controller.outputs[outputAlias];
            const nextValue = state instanceof DG.DataFrame ? state.clone() : state;
            node.getItem().getStateStore().setState(ioName, nextValue, restriction);
          }
        }
      }
    }
  }
}

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {PipelineLinkConfigurationBase} from '../config/PipelineConfiguration';
import {NodePath, TreeNode} from '../data/BaseTree';
import {PipelineState} from '../config/PipelineInstance';
import {PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {StateTree} from './StateTree';
import {isSequentialPipelineNode, isStaticPipelineNode, StateTreeNode} from './StateTreeNodes';
import {MatchInfo, matchNodeLink} from './link-matching';
import {BehaviorSubject, combineLatest, defer, EMPTY, merge, Subject} from 'rxjs';
import {IRuntimeLinkController, IRuntimeValidatorController} from '../RuntimeControllers';
import {RestrictionType} from '../data/common-types';
import {map, filter, takeUntil, withLatestFrom, switchMap, catchError, mapTo, tap, finalize, debounceTime} from 'rxjs/operators';
import {callHandler} from '../utils';
import {defaultLinkHandler} from './default-handler';
import {ActionItem, ValidationResultBase} from '../../../shared-utils/validation';


class ControllerCancelled extends Error { };

export class ControllerBase<T> {
  private isClosed = true;

  public outputs: Record<string, T> = {};

  constructor(public inputs: Record<string, any[]>, public inputsSet: Set<string>, public outputsSet: Set<string>, public id?: string) {}

  protected checkInput(name: string) {
    if (!this.inputsSet.has(name)) {
      const e = new Error(`Handler for Link ${this.id} is trying to set an unknown input ${name}`);
      grok.shell.error(e.message);
      throw e;
    }
  }

  protected checkOutput(name: string) {
    if (!this.outputsSet.has(name)) {
      const e = new Error(`Handler for Link ${this.id} is trying to set an unknown output ${name}`);
      grok.shell.error(e.message);
      throw e;
    }
  }

  protected checkIsClosed() {
    if (!this.isClosed)
      throw new ControllerCancelled();
  }

  close() {
    this.isClosed = false;
  }
}

export class LinkController extends ControllerBase<[any, RestrictionType]> implements IRuntimeLinkController {
  constructor(public inputs: Record<string, any[]>, public inputsSet: Set<string>, public outputsSet: Set<string>, public id?: string) {
    super(inputs, inputsSet, outputsSet, id);
  }

  getAll<T = any>(name: string): T[] {
    this.checkIsClosed();
    this.checkInput(name);
    return this.inputs[name];
  }

  getFirst<T = any>(name: string): T {
    return this.getAll<T>(name)?.[0];
  }

  setAll<T = any>(name: string, state: T, restriction: RestrictionType = 'none') {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = [state, restriction] as const;
  }
}

export class ValidatorController extends ControllerBase<ValidationResultBase | undefined> implements IRuntimeValidatorController {
  constructor(public inputs: Record<string, any[]>, public inputsSet: Set<string>, public outputsSet: Set<string>, public id?: string, public baseNode?: TreeNode<StateTreeNode>) {
    super(inputs, inputsSet, outputsSet, id);
  }

  getAll<T = any>(name: string): T[] {
    this.checkIsClosed();
    this.checkInput(name);
    return this.inputs[name];
  }

  getFirst<T = any>(name: string): T {
    return this.getAll<T>(name)?.[0];
  }

  // TODO: needs links structure
  // getValidationAction(id: string, action: string): ActionItem | undefined {
  //   this.checkIsClosed();
  // }

  setValidation(name: string, validation?: ValidationResultBase | undefined) {
    this.checkIsClosed();
    this.checkOutput(name);
    this.outputs[name] = validation;
  }
}

export class Link {
  public uuid = uuidv4();

  private destroyed$ = new Subject<true>();
  private isActive$ = new BehaviorSubject(false);
  private trigger$ = new Subject<true>();

  constructor(public prefix: NodePath, public matchInfo: MatchInfo) {}

  wire(state: StateTree) {
    const inputNames = Object.keys(this.matchInfo.inputs);
    const outputNames = Object.keys(this.matchInfo.outputs);
    const inputSet = new Set(inputNames);
    const outputSet = new Set(outputNames);

    const inputs$ = this.makeInputsChanges(state);
    const baseNode = this.matchInfo.basePath ? state.getNode(([...this.prefix, ...this.matchInfo.basePath])) : undefined;

    inputs$.pipe(
      debounceTime(0),
      switchMap((inputs) => this.runHandler(inputs, inputSet, outputSet, inputNames, outputNames, baseNode)),
      tap((controller) => this.setHandlerResults(controller, state)),
      takeUntil(this.destroyed$),
    ).subscribe();
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

    const activeInputs$ = inputsEntries$.pipe(filter(() => this.isActive$.value));

    const inputsChanges$ = merge(activeInputs$, inputsTriggered$).pipe(
      map((entries) => Object.fromEntries(entries)),
      filter(() => this.isActive$.value),
    );
    return inputsChanges$;
  }

  private runHandler(inputs: Record<string, any>, inputSet: Set<string>, outputSet: Set<string>, inputNames: string[], outputNames: string[], baseNode?: TreeNode<StateTreeNode>) {
    const controller = this.matchInfo.spec.isValidator ?
      new ValidatorController(inputs, inputSet, outputSet, this.matchInfo.spec.id, baseNode) :
      new LinkController(inputs, inputSet, outputSet, this.matchInfo.spec.id);

    if (this.matchInfo.spec.handler) {
      return callHandler(this.matchInfo.spec.handler, {controller}).pipe(
        catchError((e) => {
          if (e instanceof ControllerCancelled)
            return EMPTY;
          throw e;
        }),
        mapTo(controller),
        finalize(() => controller.close()),
      );
    } else if (!this.matchInfo.spec.isValidator)
      return defer(() => defaultLinkHandler(controller as LinkController, inputNames, outputNames)).pipe(mapTo(controller));

    throw Error(`Unable to run handler for link ${this.matchInfo.spec.id}`);
  }

  private setHandlerResults(controller: ControllerBase<any>, state: StateTree) {
    const outputsEntries = Object.entries(this.matchInfo.outputs).map(([outputAlias, outputItems]) => {
      const nodes = outputItems.map((output) => [output.ioName!, state.getNode([...this.prefix, ...output.path])] as const);
      return [outputAlias, nodes] as const;
    });
    for (const [outputAlias, nodesData] of outputsEntries) {
      const nextData = controller.outputs[outputAlias];
      if (nextData) {
        for (const [ioName, node] of nodesData) {
          if (this.matchInfo.spec.isValidator)
            node.getItem().getStateStore().setValidation(ioName, this.uuid, nextData[0]);
          else {
            const nextValue = nextData[0] instanceof DG.DataFrame ? nextData[0].clone() : nextData[0];
            node.getItem().getStateStore().setState(ioName, nextValue, nextData[1]);
          }
        }
      }
    }
  }
}

export class LinksState {
  constructor() {}

  public update(state: StateTree, mutationPath: NodePath) {
    const nlinks = this.createLinks(state);
    // TODO
  }

  public createLinks(state: StateTree) {
    const links = state.traverse(state.getRoot(), (acc, node, path) => {
      const item = node.getItem();
      if (isStaticPipelineNode(item) || isSequentialPipelineNode(item)) {
        const {config} = item;
        const matchedLinks = (config.links ?? [])
          .map((link) => matchNodeLink(node, link))
          .filter((x) => !!x)
          .flat();
        const links = matchedLinks.map((minfo) => new Link(path, minfo));
        return [...acc, ...links];
      }
      return acc;
    }, [] as Link[]);
    return links;
  }
}

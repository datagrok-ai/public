import {BehaviorSubject, Observable, Subject, of} from 'rxjs';
import {InputState} from '../data/common-types';
import {v4 as uuidv4} from 'uuid';
import {switchMap} from 'rxjs/operators';
import {NodeTree} from '../data/NodeTree';
import {PipelineConfigurationProvided, PipelineStepConfiguration, StateItem} from '../config/PipelineConfiguration';
import {FuncallStateItem} from '../config/config-processing-utils';
import {ValidationResultBase} from '../../../shared-utils/validation';

export interface ICallable {
  run(): void
  isRunning$: Observable<boolean>;
  isRunable$: Observable<boolean>;
  isOuputOutdated$: Observable<boolean>;
}

export interface IStateStore {
  getStateChanges<T = any>(id: string): Observable<T | undefined>;
  getState<T = any>(id: string): T | undefined;
  setState<T = any >(id: string, value: T | undefined, inputState?: InputState): void;
}

export interface IValidationDisplay {
  setValidation(id: string, validation: ValidationResultBase | undefined): void;
}

export class StoreItem<T = any> {
  state$ = new BehaviorSubject<T | undefined>(undefined);
}

export class MemoryStore implements IStateStore {
  public readonly uuid = uuidv4();
  states: Record<string, StoreItem> = {};

  constructor(private statesDescriptions: StateItem[]) {
    for (const description of this.statesDescriptions)
      this.states[description.id] = new StoreItem();
  }

  getState<T = any>(id: string): T {
    return this.states[id]?.state$?.value;
  }

  getStateChanges<T = any>(id: string): Observable<T | undefined> {
    return this.states[id]?.state$.asObservable();
  }

  setState<T = any>(id: string, val: T | undefined, _inputState: InputState) {
    this.states[id]?.state$.next(val);
  }
}

export type IFuncCallBridge = IStateStore & IValidationDisplay & ICallable;

export class FuncCallInstanceAdapter implements IFuncCallBridge {
  private instance$ = new BehaviorSubject<IFuncCallBridge | undefined>(undefined);
  public isRunning$ = this.instance$.pipe(
    switchMap((instance) => instance? instance.isRunning$ : of(false)),
  );
  public isRunable$ = this.instance$.pipe(
    switchMap((instance) => instance? instance.isRunable$ : of(false)),
  );
  public isOuputOutdated$ = this.instance$.pipe(
    switchMap((instance) => instance? instance.isOuputOutdated$ : of(false)),
  );

  constructor() {}

  setInstance(instance: IFuncCallBridge | undefined) {
    this.instance$.next(instance);
  }

  getState<T = any>(id: string): T | undefined {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return currentInstance.getState(id);
  }

  getStateChanges<T = any>(id: string): Observable<T | undefined> {
    return this.instance$.pipe(
      switchMap((instance) => instance ? instance.getStateChanges(id) : of(undefined)),
    );
  }

  setState<T = any>(id: string, val: T | undefined, inputState: InputState) {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      currentInstance.setState(id, val, inputState);
  }

  setValidation(id: string, validation: ValidationResultBase | undefined) {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return currentInstance.setValidation(id, validation);
  }

  run() {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      currentInstance.run();
  }
}

export interface IStoreProvider {
  getStore(): IStateStore;
}

export class FuncCallNode implements IStoreProvider {
  public readonly uuid = uuidv4();
  public readonly nodeType = 'funccall';
  private adapter = new FuncCallInstanceAdapter();

  constructor(
    public readonly config: PipelineStepConfiguration<FuncallStateItem>[],
  ) {}

  setFuncall(fc?: IFuncCallBridge) {
    this.adapter.setInstance(fc);
  }

  getStore() {
    return this.adapter;
  }
}


export class PipelineNode implements IStoreProvider {
  public readonly uuid = uuidv4();
  public readonly nodeType = 'pipeline';
  private store: MemoryStore;

  constructor(
    public readonly config: PipelineConfigurationProvided,
  ) {
    this.store = new MemoryStore(config.states ?? []);
  }

  getStore() {
    return this.store;
  }
}

export type StateTree = NodeTree<FuncCallNode | PipelineNode>;

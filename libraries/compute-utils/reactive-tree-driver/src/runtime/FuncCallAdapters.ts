import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {BehaviorSubject, Observable, from, defer, Subject, merge, identity} from 'rxjs';
import {ValidationResultBase} from '../../../shared-utils/validation';
import {StateItem} from '../config/PipelineConfiguration';
import {delay, map, mapTo, skip, startWith, switchMap, takeUntil} from 'rxjs/operators';
import {RestrictionType} from '../data/common-types';


export interface RestrictionState {
  type: RestrictionType,
  assignedValue: any,
}

export interface IRunnableWrapper {
  id?: string;
  instance?: DG.FuncCall;
  run(mockResults?: Record<string, any>, mockDelay?: number): Observable<any>;
  close(): void;
  isRunning$: BehaviorSubject<boolean>;
  isOutputOutdated$: BehaviorSubject<boolean>;
}

export interface IFuncCallWrapper extends IRunnableWrapper {
  getFuncCall(): DG.FuncCall;
}

export interface IStateStore {
  getStateChanges<T = any>(name: string): Observable<T | undefined>;
  getState<T = any>(name: string): T | undefined;
  setState<T = any>(name: string, value: T | undefined, restrictionType?: RestrictionType): void;
  getStateNames(): string[];
}

export interface IValidationStore {
  setValidation(name: string, validatorId: string, validation: ValidationResultBase | undefined): void;
  validations$: BehaviorSubject<Record<string, Record<string, ValidationResultBase | undefined>>>;
  inputRestrictions$: BehaviorSubject<Record<string, RestrictionState | undefined>>;
}

export type IFuncCallAdapter = IStateStore & IValidationStore & IFuncCallWrapper;


// real implementation

export class FuncCallAdapter implements IFuncCallAdapter {
  id = this.instance.id;
  isRunning$ = new BehaviorSubject(false);
  isOutputOutdated$ = new BehaviorSubject(true);
  validations$ = new BehaviorSubject<Record<string, Record<string, ValidationResultBase | undefined>>>({});
  inputRestrictions$ = new BehaviorSubject<Record<string, RestrictionState | undefined>>({});

  private closed$ = new Subject<true>();

  constructor(public instance: DG.FuncCall) {
    const allParamsChanges = Object.keys(instance.inputs).map((inputName) => this.getStateChanges(inputName).pipe(skip(1)));
    merge(...allParamsChanges).pipe(
      mapTo(true),
      takeUntil(this.closed$),
    ).subscribe(this.isOutputOutdated$);
  }

  run() {
    return from(defer(async () => {
      try {
        this.isRunning$.next(true);
        await this.instance.call();
      } finally {
        this.isRunning$.next(false);
        this.isOutputOutdated$.next(false);
      }
    }));
  }

  getStateChanges(name: string) {
    const ptype = this.getPtype(name);
    const param = this.instance[ptype][name];
    return param.onChanged.pipe(
      startWith(null),
      map(() => param.value),
    );
  }

  getState<T = any>(name: string) {
    const ptype = this.getPtype(name);
    return this.instance[ptype][name] as T;
  }

  setState<T = any>(name: string, value: T | undefined, restrictionType: RestrictionType = 'none') {
    this.instance.inputs[name] = value;
    const currentRestrictions = this.inputRestrictions$.value;
    const restrictionState = restrictionType === 'none' ? undefined : {type: restrictionType, assignedValue: value};
    this.inputRestrictions$.next({...currentRestrictions, [name]: restrictionState});
  }

  getStateNames() {
    return [...Object.keys(this.instance.inputs), ...Object.keys(this.instance.outputs)];
  }

  setValidation(name: string, validatorId: string, validation: ValidationResultBase | undefined) {
    const allValidations = this.validations$.value;
    const validatorResults = allValidations[validatorId] ?? {};
    this.validations$.next({...allValidations, [validatorId]: {...validatorResults, [name]: validation}});
  }

  getFuncCall() {
    return this.instance;
  }

  close() {
    this.closed$.next(true);
  }

  private getPtype(name: string) {
    return this.instance['inputParams'][name] ? 'inputParams' : 'outputParams';
  }
}

// non-funcall backed states

export class StoreItem<T = any> {
  state$ = new BehaviorSubject<T | undefined>(undefined);
}

export class MemoryStore implements IStateStore, IValidationStore {
  public readonly uuid = uuidv4();
  states: Record<string, StoreItem> = {};
  inputRestrictions$ = new BehaviorSubject<Record<string, RestrictionState | undefined>>({});
  validations$ = new BehaviorSubject<Record<string, Record<string, ValidationResultBase | undefined>>>({});

  constructor(private statesDescriptions: StateItem[]) {
    for (const description of this.statesDescriptions)
      this.states[description.id] = new StoreItem();
  }

  getState<T = any>(id: string): T {
    return this.states[id]?.state$?.value;
  }

  getStateChanges<T = any>(id: string, includeDataFrameMutations = false): Observable<T | undefined> {
    return this.states[id]?.state$.pipe(
      includeDataFrameMutations ? switchMap((x) => x instanceof DG.DataFrame ? x.onDataChanged.pipe(mapTo(x)) : x) : identity,
    );
  }

  getStateNames() {
    return Object.keys(this.states);
  }

  setState<T = any>(name: string, value: T | undefined, restrictionType: RestrictionType = 'none') {
    const currentRestrictions = this.inputRestrictions$.value;
    const restrictionState = restrictionType === 'none' ? undefined : {type: restrictionType, assignedValue: value};
    this.inputRestrictions$.next({...currentRestrictions, [name]: restrictionState});
    this.states[name]?.state$.next(value);
  }

  setValidation(name: string, validatorId: string, validation: ValidationResultBase | undefined) {
    const allValidations = this.validations$.value;
    const validatorResults = allValidations[validatorId] ?? {};
    this.validations$.next({...allValidations, [validatorId]: {...validatorResults, [name]: validation}});
  }
}

// mock implementation for rxjs testing

export class FuncCallMockAdapter extends MemoryStore implements IFuncCallAdapter {
  id = uuidv4();
  isRunning$ = new BehaviorSubject(false);
  isOutputOutdated$ = new BehaviorSubject(true);
  instance = undefined;

  constructor(statesDescriptions: StateItem[]) {
    super(statesDescriptions);
  }

  run(outputs?: Record<string, any>, delayTime = 0) {
    return from(defer(async () => {
      try {
        this.isRunning$.next(true);
        for (const [k, v] of Object.entries(outputs ?? {}))
          this.setState(k, v);
      } finally {
        this.isRunning$.next(false);
        this.isOutputOutdated$.next(false);
      }
    })).pipe(delay(delayTime));
  }

  getFuncCall(): DG.FuncCall {
    throw new Error(`Not implemented for mocks`);
  }

  close() {}
}

export function isMockAdapter(adapter: IRunnableWrapper): adapter is FuncCallMockAdapter {
  return adapter.instance == null;
}

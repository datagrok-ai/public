import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {BehaviorSubject, Observable, from, defer, Subject, merge} from 'rxjs';
import {ValidationResultBase} from '../../../shared-utils/validation';
import {StateItem} from '../config/PipelineConfiguration';
import {delay, map, mapTo, skip, startWith, takeUntil, tap} from 'rxjs/operators';
import {RestrictionType} from '../data/common-types';


export interface RestrictionState {
  type: RestrictionType,
  assignedValue: any,
}

export interface IRunnableWrapper {
  id?: string;
  run(mockResults?: Record<string, any>, mockDelay?: number): Observable<any>;
  close(): void;
  isRunning$: BehaviorSubject<boolean>;
  isOutputOutdated$: BehaviorSubject<boolean>;
}

export interface IFuncCallWrapper extends IRunnableWrapper {
  getFuncCall(): DG.FuncCall;
}

export interface IStateStore {
  isReadonly: boolean;
  getStateChanges<T = any>(name: string): Observable<T | undefined>;
  getState<T = any>(name: string): T | undefined;
  setState<T = any>(name: string, value: T | undefined, restrictionType?: RestrictionType): void;
  editState<T = any>(name: string, value: T | undefined): void; // for tests
  getStateNames(): string[];
}

export interface IValidationStore {
  setValidation(name: string, validatorId: string, validation: ValidationResultBase | undefined): void;
  validations$: BehaviorSubject<Record<string, Record<string, ValidationResultBase | undefined>>>;
  inputRestrictions$: BehaviorSubject<Record<string, RestrictionState | undefined>>;
  inputRestrictionsUpdates$: Subject<[string, RestrictionState | undefined]>;
}

export type IFuncCallAdapter = IStateStore & IValidationStore & IFuncCallWrapper & {
  transferAdditionalState(newFC: DG.FuncCall) : IFuncCallAdapter
};


// real implementation

export class FuncCallAdapter implements IFuncCallAdapter {
  id = this.instance.id;
  isRunning$ = new BehaviorSubject(false);
  isOutputOutdated$ = new BehaviorSubject(true);
  validations$ = new BehaviorSubject<Record<string, Record<string, ValidationResultBase | undefined>>>({});
  inputRestrictions$ = new BehaviorSubject<Record<string, RestrictionState | undefined>>({});
  inputRestrictionsUpdates$ = new Subject<[string, RestrictionState | undefined]>();

  private closed$ = new Subject<true>();

  constructor(public instance: DG.FuncCall, public readonly isReadonly: boolean) {
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
    return this.instance[ptype][name].value as T;
  }

  setState<T = any>(name: string, value: T | undefined, restrictionType: RestrictionType = 'none') {
    const currentRestrictions = this.inputRestrictions$.value;
    const restrictionState = restrictionType === 'none' ? undefined : {type: restrictionType, assignedValue: value};
    this.inputRestrictions$.next({...currentRestrictions, [name]: restrictionState});

    if (this.isReadonly)
      this.inputRestrictionsUpdates$.next([name, restrictionState]);
    else
      this.instance.inputs[name] = value;
  }

  editState<T = any>(name: string, value: T | undefined) {
    if (!this.isReadonly)
      this.instance.inputs[name] = value;
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

  transferAdditionalState(newFC: DG.FuncCall) : IFuncCallAdapter {
    const nAdapter = new FuncCallAdapter(newFC, this.isReadonly);
    nAdapter.isOutputOutdated$.next(this.isOutputOutdated$.value);
    nAdapter.validations$.next(this.validations$.value);
    nAdapter.inputRestrictions$.next(this.inputRestrictions$.value);
    return nAdapter;
  }

  private getPtype(name: string) {
    return this.instance['inputParams'][name] ? 'inputParams' : 'outputParams';
  }
}

// non-funcall backed states


export class MemoryStore implements IStateStore, IValidationStore {
  public readonly uuid = uuidv4();
  states: Record<string, BehaviorSubject<any | undefined>> = {};
  validations$ = new BehaviorSubject<Record<string, Record<string, ValidationResultBase | undefined>>>({});
  inputRestrictions$ = new BehaviorSubject<Record<string, RestrictionState | undefined>>({});
  inputRestrictionsUpdates$ = new Subject<[string, RestrictionState | undefined]>();

  constructor(private statesDescriptions: StateItem[], public readonly isReadonly: boolean) {
    for (const description of this.statesDescriptions)
      this.states[description.id] = new BehaviorSubject(undefined);
  }

  getState<T = any>(id: string): T {
    return this.states[id]?.value;
  }

  getStateChanges<T = any>(id: string, includeDataFrameMutations = false): Observable<T | undefined> {
    return this.states[id]?.pipe(
      (x) => (includeDataFrameMutations && x instanceof DG.DataFrame) ?
        x.onDataChanged.pipe(startWith(null), mapTo(x)) :
        x,
    );
  }

  getStateNames() {
    return Object.keys(this.states);
  }

  setState<T = any>(name: string, value: T | undefined, restrictionType: RestrictionType = 'none') {
    const currentRestrictions = this.inputRestrictions$.value;
    const restrictionState = restrictionType === 'none' ? undefined : {type: restrictionType, assignedValue: value};
    this.inputRestrictions$.next({...currentRestrictions, [name]: restrictionState});

    if (this.isReadonly)
      this.inputRestrictionsUpdates$.next([name, restrictionState]);
    else
      this.states[name]?.next(value);
  }

  editState<T = any>(name: string, value: T | undefined) {
    if (!this.isReadonly)
      this.states[name]?.next(value);
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

  constructor(statesDescriptions: StateItem[], isReadonly: boolean) {
    super(statesDescriptions, isReadonly);
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
    throw new Error(`getFuncCall is not implemented for mocks`);
  }

  transferAdditionalState( _newFC: DG.FuncCall) : IFuncCallAdapter {
    throw new Error(`transferAdditionalState is not implemented for mocks`);
  }

  close() {}
}

import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {BehaviorSubject, Observable, from, defer} from 'rxjs';
import {ValidationResultBase} from '../../../shared-utils/validation';
import {StateItem} from '../config/PipelineConfiguration';
import {delay, map, startWith} from 'rxjs/operators';

export type RestrictionType = 'disabled' | 'restricted' | 'info' | 'none';

export interface RestrictionState {
  type: RestrictionType,
  assignedValue: any,
}

export interface IRunnableWrapper {
  id?: string;
  instance?: DG.FuncCall;
  run(): Observable<any>;
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

  constructor(public instance: DG.FuncCall) {}

  run() {
    return from(defer(async () => {
      try {
        this.isRunning$.next(true);
        this.instance.call();
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
    this.instance.inputs[name].value = value;
    const currentRestrictions = this.inputRestrictions$.value;
    const restrictionState = restrictionType === 'none' ? undefined : {type: restrictionType, assignedValue: value};
    this.inputRestrictions$.next({...currentRestrictions, [name]: restrictionState});
  }

  setValidation(name: string, validatorId: string, validation: ValidationResultBase | undefined) {
    const allValidations = this.validations$.value;
    const validatorResults = allValidations[validatorId] ?? {};
    this.validations$.next({...allValidations, [validatorId]: {...validatorResults, [name]: validation}});
  }

  getFuncCall() {
    return this.instance;
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

  getStateChanges<T = any>(id: string): Observable<T | undefined> {
    return this.states[id]?.state$.asObservable();
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
}

export function isMockAdapter(adapter: IRunnableWrapper): adapter is FuncCallMockAdapter {
  return adapter.instance == null;
}

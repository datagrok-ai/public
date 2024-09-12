import {BehaviorSubject, of, combineLatest, Observable, defer, Subject, merge} from 'rxjs';
import {switchMap, map, takeUntil, finalize, mapTo, skip, distinctUntilChanged, withLatestFrom, filter} from 'rxjs/operators';
import {ValidationResultBase} from '../../../shared-utils/validation';
import {IFuncCallAdapter, IRunnableWrapper, IStateStore} from './FuncCallAdapters';
import {RestrictionType} from '../data/common-types';
import {FuncallStateItem} from '../config/config-processing-utils';

export interface RestrictionState {
  type: RestrictionType,
  assignedValue: any,
}

export interface BridgePreInitData {
  initialValues: Record<string, any>,
  initialRestrictions: Record<string, RestrictionState | undefined>,
}

export interface BridgeInitData {
  adapter: IFuncCallAdapter,
  restrictions: Record<string, RestrictionState | undefined>,
  isOutputOutdated: boolean,
  initValues: boolean,
}

interface InstanceData {
  adapter: IFuncCallAdapter,
  isNew: boolean
}

export class FuncCallInstancesBridge implements IStateStore, IRunnableWrapper {
  public instance$ = new BehaviorSubject<InstanceData | undefined>(undefined);

  public isRunning$ = new BehaviorSubject(false);
  public isRunable$ = new BehaviorSubject(false);
  public isOutputOutdated$ = new BehaviorSubject(true);

  public validations$ = new BehaviorSubject<Record<string, Record<string, ValidationResultBase | undefined>>>({});
  public meta$ = new BehaviorSubject<Record<string, BehaviorSubject<any | undefined>>>({});

  public inputRestrictions$ = new BehaviorSubject<Record<string, RestrictionState | undefined>>({});
  public inputRestrictionsUpdates$ = new Subject<[string, RestrictionState | undefined]>();

  public initialValues: Record<string, any> = {};

  private closed$ = new Subject<true>();
  public outdatedChanged$ = new BehaviorSubject(true);

  constructor(private io: FuncallStateItem[], public readonly isReadonly: boolean) {}

  get id() {
    return this.instance$.value?.adapter?.id;
  }

  setPreInitialData(data: BridgePreInitData) {
    this.initialValues = data.initialValues;
    this.inputRestrictions$.next(data.initialRestrictions);
  }

  init(data: BridgeInitData) {
    if (this.instance$.value)
      throw new Error(`Double funcCall bridge instance init`);
    for (const [key, val] of Object.entries(this.initialValues))
      data.adapter.setState(key, val);
    if (data.initValues && data.adapter) {
      for (const [key, val] of Object.entries(this.initialValues))
        data.adapter.setState(key, val);
    }

    this.inputRestrictions$.next({...this.inputRestrictions$.value, ...data.restrictions});
    this.instance$.next({adapter: data.adapter, isNew: true});
    this.outdatedChanged$.next(data.isOutputOutdated);

    this.setupStateWatcher();
  }

  change(adapter: IFuncCallAdapter) {
    this.instance$.next({adapter, isNew: false});
  }

  getInstance() {
    return this.instance$.value?.adapter;
  }

  getState<T = any>(id: string): T | undefined {
    const currentInstance = this.instance$.value?.adapter;
    if (currentInstance)
      return currentInstance.getState<T>(id);
  }

  getStateChanges<T = any>(id: string, includeDataFrameMutations = false): Observable<T | undefined> {
    return this.instance$.pipe(
      switchMap((data) => {
        if (data) {
          const changes$ = data.adapter.getStateChanges(id, includeDataFrameMutations);
          if (!data.isNew)
            return changes$.pipe(skip(1));
          return changes$;
        }
        return of(undefined);
      }),
    );
  }

  setState<T = any>(id: string, val: T | undefined, restrictionType: RestrictionType = 'none') {
    const currentInstance = this.instance$.value?.adapter;
    if (currentInstance == null)
      throw new Error(`Attempting to set an empty FuncCallInstancesBridge`);
    this.inputRestrictions$.next({...this.inputRestrictions$.value, [id]: {assignedValue: val, type: restrictionType}});
    if (!this.isReadonly)
      currentInstance.setState(id, val, restrictionType);
    else
      this.inputRestrictionsUpdates$.next([id, {assignedValue: val, type: restrictionType}] as const);
  }

  editState<T = any>(id: string, val: T | undefined) {
    const currentInstance = this.instance$.value?.adapter;
    if (currentInstance)
      currentInstance.editState(id, val);
  }

  setValidation(id: string, validatorId: string, validation: ValidationResultBase | undefined) {
    const allValidations = this.validations$.value;
    const validatorResults = allValidations[validatorId] ?? {};
    this.validations$.next({...allValidations, [validatorId]: {...validatorResults, [id]: validation}});
  }

  setMeta(id: string, meta: any | undefined) {
    const allMeta = this.meta$.value;
    allMeta[id].next(meta);
  }

  run(mockResults?: Record<string, any>, mockDelay?: number) {
    return defer(() => {
      const currentInstance = this.instance$.value?.adapter;
      if (!currentInstance)
        throw new Error(`Attempting to run an empty FuncCallInstancesBridge`);
      if (this.isRunning$.value)
        throw new Error(`Attempting to run a running FuncCallInstancesBridge`);
      if (!this.isRunable$.value)
        throw new Error(`Attempting to run FuncCallInstancesBridge with validation errors`);
      this.isRunning$.next(true);
      return currentInstance.run(mockResults, mockDelay).pipe(
        finalize(() => {
          this.isRunning$.next(false);
          this.outdatedChanged$.next(false);
        }),
      );
    });
  }

  getStateNames() {
    return this.io.map((item) => item.id);
  }

  close() {
    this.closed$.next(true);
  }

  private setupStateWatcher() {
    this.instance$.pipe(
      switchMap((instance) => {
        if (instance == null)
          return of(false);
        return combineLatest([
          this.isRunning$,
          this.validations$.pipe(map((validations) => this.isRunnable(validations))),
        ]).pipe(map(([isRunning, isValid]) => !isRunning && isValid));
      }),
      distinctUntilChanged(),
      takeUntil(this.closed$),
    ).subscribe(this.isRunable$);

    const inputs = this.io.filter((item) => item.direction === 'input').map((item) => item.id);

    const inputsChanges = inputs.map((inputName) => this.getStateChanges(inputName).pipe(skip(1)));
    merge(...inputsChanges).pipe(
      mapTo(true),
      takeUntil(this.closed$),
    ).subscribe(this.outdatedChanged$);

    this.outdatedChanged$.pipe(
      withLatestFrom(this.isOutputOutdated$),
      filter(([next, current]) => next !== current),
      map(([next]) => next),
      takeUntil(this.closed$),
    ).subscribe(this.isOutputOutdated$);

    this.instance$.pipe(
      map((instance) => {
        if (instance == null)
          return {};
        const res: Record<string, BehaviorSubject<any | undefined>> = {};
        for (const item of this.io) {
          res[item.id] = new BehaviorSubject<any>(undefined);
        }
        return res;
      }),
      takeUntil(this.closed$),
    ).subscribe(this.meta$);
  }

  private isRunnable(validations: Record<string, Record<string, ValidationResultBase | undefined>>) {
    for (const validatorResults of Object.values(validations)) {
      for (const res of Object.values(validatorResults)) {
        if (res?.errors?.length)
          return false;
      }
    }
    return true;
  }
}

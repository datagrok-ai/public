import {BehaviorSubject, of, combineLatest, Observable, from, defer, Subject} from 'rxjs';
import {switchMap, map, takeUntil, take, withLatestFrom} from 'rxjs/operators';
import {ValidationResultBase} from '../../../shared-utils/validation';
import {IStateStore, IValidationStore, IRunnableWrapper, IFuncCallAdapter, RestrictionType} from './FuncCallAdapters';

export class FuncCallInstancesBridge implements IStateStore, IValidationStore, IRunnableWrapper {
  private instance$ = new BehaviorSubject<IFuncCallAdapter | undefined>(undefined);
  public isRunning$ = new BehaviorSubject(false);
  public isLoading$ = new BehaviorSubject(true);
  public isRunable$ = new BehaviorSubject(false);
  public isOutputOutdated$ = new BehaviorSubject(false);
  public isCurrent$ = new BehaviorSubject(false);
  public validations$ = new BehaviorSubject({});
  public inputRestrictions$ = new BehaviorSubject({});
  public initialValues: Record<string, any> = {};

  private closed$ = new Subject<true>();

  constructor() {
    this.instance$.pipe(
      switchMap((instance) => instance ? instance.isRunning$ : of(false)),
      takeUntil(this.closed$),
    ).subscribe(this.isRunning$);

    this.instance$.pipe(
      switchMap((instance) => {
        if (instance == null)
          return of(false);
        return combineLatest([
          instance.isRunning$,
          instance.validations$.pipe(map((validations) => this.isRunnable(validations))),
        ]).pipe(map((isRunning, isValid) => !isRunning && isValid));
      }),
      takeUntil(this.closed$),
    ).subscribe(this.isRunable$);

    this.instance$.pipe(
      switchMap((instance) => instance ? instance.isOutputOutdated$ : of(false)),
      takeUntil(this.closed$),
    ).subscribe(this.isOutputOutdated$);

    this.closed$.pipe(
      take(1),
      withLatestFrom(this.instance$),
    ).subscribe(([, inst]) => {
      if (inst)
        inst.close();
    });
  }

  get id() {
    return this.instance$.value?.id;
  }

  setInstance(instance: IFuncCallAdapter | undefined, initValues = false) {
    if (initValues && instance) {
      for (const [key, val] of Object.entries(this.initialValues))
        instance.setState(key, val);
    }
    if (this.instance$.value)
      this.instance$.value.close();

    this.instance$.next(instance);
  }

  getInstance() {
    return this.instance$.value;
  }

  getState<T = any>(id: string): T | undefined {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return currentInstance.getState<T>(id);
  }

  getStateChanges<T = any>(id: string): Observable<T | undefined> {
    return this.instance$.pipe(
      switchMap((instance) => instance ? instance.getStateChanges(id) : of(undefined)),
    );
  }

  setState<T = any>(id: string, val: T | undefined, restrictionType: RestrictionType = 'none') {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      currentInstance.setState(id, val, restrictionType);
  }

  setInitialValues(values: Record<string, any>) {
    this.initialValues = values;
  }

  getInitialValues(): Record<string, any> {
    return this.initialValues;
  }

  setValidation(id: string, validatorId: string, validation: ValidationResultBase | undefined) {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return currentInstance.setValidation(id, validatorId, validation);
  }

  run(mockResults?: Record<string, any>, mockDelay?: number) {
    const currentInstance = this.instance$.value;
    if (currentInstance)
      return from(defer(() => currentInstance.run(mockResults, mockDelay)));

    else
      throw new Error(`Attempting to run an empty FuncCallInstancesBridge`);
  }

  close() {
    this.closed$.next(true);
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

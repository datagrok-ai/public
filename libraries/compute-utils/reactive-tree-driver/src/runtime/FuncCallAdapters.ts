import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {BehaviorSubject, Observable, defer, of} from 'rxjs';
import {StateItem} from '../config/PipelineConfiguration';
import {delay, map, mapTo, startWith, switchMap, tap} from 'rxjs/operators';
import {RestrictionType} from '../data/common-types';


export interface IRunnableWrapper {
  id?: string;
  run(mockResults?: Record<string, any>, mockDelay?: number): Observable<any>;
}

export interface IFuncCallWrapper {
  getFuncCall(): DG.FuncCall;
}

export interface IStateStore {
  isReadonly: boolean;
  getStateChanges<T = any>(name: string, includeDataFrameMutations: boolean): Observable<T | undefined>;
  getState<T = any>(name: string): T | undefined;
  setState<T = any>(name: string, value: T | undefined, restrictionType?: RestrictionType): void;
  editState<T = any>(name: string, value: T | undefined): void; // for tests
  removeRestriction(name: string): void;
  getStateNames(): string[];
}

export type IFuncCallAdapter = IStateStore & IRunnableWrapper & IFuncCallWrapper;

// funcall based implementation

export class FuncCallAdapter implements IFuncCallAdapter {
  id = this.instance.id;

  constructor(public instance: DG.FuncCall, public readonly isReadonly: boolean) {}

  run() {
    return defer(() => this.instance.call());
  }

  getStateChanges<T = any>(
    name: string,
    includeDataFrameMutations: boolean,
  ): Observable<T | undefined> {
    const ptype = this.getPtype(name);
    const param = this.instance[ptype][name];
    const changes$ = param.onChanged.pipe(
      startWith(null),
      map(() => param.value),
    );
    if (includeDataFrameMutations) {
      return changes$.pipe(
        switchMap((x) => x instanceof DG.DataFrame ?
          x.onDataChanged.pipe(startWith(null), mapTo(x)) :
          of(x)),
      );
    }
    return changes$;
  }

  getState<T = any>(name: string) {
    const ptype = this.getPtype(name);
    return this.instance[ptype][name].value as T;
  }

  setState<T = any>(name: string, value: T | undefined, _restrictionType?: RestrictionType) {
    if (!this.isReadonly)
      this.instance.inputs[name] = value;
  }

  editState<T = any>(name: string, value: T | undefined) {
    this.instance.inputs[name] = value;
  }

  removeRestriction(_name: string) {
    return;
  }

  getStateNames() {
    return [...Object.keys(this.instance.inputs), ...Object.keys(this.instance.outputs)];
  }

  getFuncCall() {
    return this.instance;
  }

  private getPtype(name: string) {
    return this.instance['inputParams'][name] ? 'inputParams' : 'outputParams';
  }
}

// non-funcall backed states

export class MemoryStore implements IStateStore {
  public readonly uuid = uuidv4();
  states: Record<string, BehaviorSubject<any | undefined>> = {};

  constructor(private statesDescriptions: StateItem[], public readonly isReadonly: boolean) {
    for (const description of this.statesDescriptions)
      this.states[description.id] = new BehaviorSubject(undefined);
  }

  getState<T = any>(id: string): T {
    return this.states[id]?.value;
  }

  getStateChanges<T = any>(id: string, includeDataFrameMutations = false): Observable<T | undefined> {
    const changes$ = this.states[id];
    if (includeDataFrameMutations) {
      return changes$.pipe(
        switchMap((x) => x instanceof DG.DataFrame ?
          x.onDataChanged.pipe(startWith(null), mapTo(x)) :
          of(x)),
      );
    }
    return changes$;
  }

  getStateNames() {
    return Object.keys(this.states);
  }

  setState<T = any>(name: string, value: T | undefined, _restrictionType?: RestrictionType) {
    if (!this.isReadonly)
      this.states[name]?.next(value);
  }

  editState<T = any>(name: string, value: T | undefined) {
    this.states[name]?.next(value);
  }

  removeRestriction(_name: string) {
    return;
  }
}

// mock implementation for rxjs testing

export class FuncCallMockAdapter extends MemoryStore implements IFuncCallAdapter {
  id = uuidv4();

  constructor(statesDescriptions: StateItem[], isReadonly: boolean) {
    super(statesDescriptions, isReadonly);
  }

  run(outputs?: Record<string, any>, delayTime = 0) {
    return of(null).pipe(
      delay(delayTime),
      tap(() => {
        for (const [k, v] of Object.entries(outputs ?? {}))
          this.setState(k, v);
      }),
    );
  }

  getFuncCall(): DG.FuncCall {
    throw new Error(`getFuncCall is not implemented for mocks`);
  }
}

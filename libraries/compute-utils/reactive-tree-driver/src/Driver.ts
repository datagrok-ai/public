import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, Observable, Subject, EMPTY, of, from} from 'rxjs';
import {isFuncCallSerializedState, PipelineState} from './config/PipelineInstance';
import {AddDynamicItem, InitPipeline, LoadDynamicItem, LoadPipeline, MoveDynamicItem, RemoveDynamicItem, RunAction, RunSequence, RunStep, SaveDynamicItem, SavePipeline, ViewConfigCommands} from './view/ViewCommunication';
import {pairwise, takeUntil, concatMap, catchError, switchMap, map, mapTo, startWith, withLatestFrom, tap} from 'rxjs/operators';
import {StateTree} from './runtime/StateTree';
import {loadInstanceState} from './runtime/funccall-utils';
import {callHandler} from './utils';
import {PipelineConfiguration} from './config/PipelineConfiguration';
import {getProcessedConfig} from './config/config-processing-utils';
import {ConsistencyInfo, FuncCallStateInfo} from './runtime/StateTreeNodes';
import {ValidationResult} from './data/common-types';

export class Driver {
  public currentState$ = new BehaviorSubject<PipelineState | undefined>(undefined);
  public currentValidations$ = new BehaviorSubject<Record<
    string, BehaviorSubject<Record<string, ValidationResult>>
  >>({});
  public currentConsistency$ = new BehaviorSubject<Record<
    string, BehaviorSubject<Record<string, ConsistencyInfo>>
  >>({});
  public currentMeta$ = new BehaviorSubject<Record<string, BehaviorSubject<any | undefined>>>({});
  public currentCallsState$ = new BehaviorSubject<Record<string, BehaviorSubject<FuncCallStateInfo | undefined>>>({});

  public globalROLocked$ = new BehaviorSubject(false);
  public treeMutationsLocked$ = new BehaviorSubject(false);

  private states$ = new BehaviorSubject<StateTree | undefined>(undefined);
  private commands$ = new Subject<ViewConfigCommands>();
  private closed$ = new Subject<true>();

  constructor(private mockMode = false) {
    this.commands$.pipe(
      withLatestFrom(this.states$),
      concatMap(([msg, state]) => this.executeCommand(msg, state)),
      catchError((error) => {
        console.error(error);
        grok.shell.error(error?.message);
        return EMPTY;
      }),
      takeUntil(this.closed$),
    ).subscribe();

    this.states$.pipe(
      pairwise(),
      takeUntil(this.closed$),
    ).subscribe(([oldval]) => {
      if (oldval)
        oldval.close();
    });

    this.states$.pipe(
      switchMap((state) => state ?
        state.makeStateRequests$.pipe(startWith(null), mapTo(state)) :
        of(undefined)),
      map((state) => state ? state.toState() : undefined),
      takeUntil(this.closed$),
    ).subscribe(this.currentState$);

    this.states$.pipe(
      map((state) => state ? state.getConsistency() : {}),
      takeUntil(this.closed$),
    ).subscribe(this.currentConsistency$);

    this.states$.pipe(
      map((state) => state ? state.getMeta() : {}),
      takeUntil(this.closed$),
    ).subscribe(this.currentMeta$);

    this.states$.pipe(
      map((state) => state ? state.getValidations() : {}),
      takeUntil(this.closed$),
    ).subscribe(this.currentValidations$);

    this.states$.pipe(
      map((state) => state ? state.getFuncCallStates() : {}),
      takeUntil(this.closed$),
    ).subscribe(this.currentCallsState$);

    this.states$.pipe(
      switchMap((state) => state ? state.globalROLocked$ : of(false)),
      takeUntil(this.closed$),
    ).subscribe(this.globalROLocked$);

    this.states$.pipe(
      switchMap((state) => state ? state.treeMutationsLocked$ : of(false)),
      takeUntil(this.closed$),
    ).subscribe(this.treeMutationsLocked$);
  }

  public sendCommand(msg: ViewConfigCommands) {
    this.commands$.next(msg);
  }

  public close() {
    this.closed$.next(true);
    this.states$.value?.close();
    this.currentState$.next(undefined);
  }

  private executeCommand(msg: ViewConfigCommands, state?: StateTree): Observable<any> {
    switch (msg.event) {
    case 'addDynamicItem':
      return this.addDynamicItem(msg, state);
    case 'loadDynamicItem':
      return this.loadDynamicItem(msg, state);
    case 'saveDynamicItem':
      return this.saveDynamicItem(msg, state);
    case 'removeDynamicItem':
      return this.removeDynamicItem(msg, state);
    case 'moveDynamicItem':
      return this.moveDynamicItem(msg, state);
    case 'runStep':
      return this.runStep(msg, state);
    case 'runAction':
      return this.runAction(msg, state);
    case 'runSequence':
      return this.runSequence(msg, state);
    case 'savePipeline':
      return this.savePipeline(msg, state);
    case 'loadPipeline':
      return this.loadPipeline(msg);
    case 'initPipeline':
      return this.initPipeline(msg);
    }
    throw new Error(`Unknow tree driver command ${(msg as any).event}`);
  }

  private addDynamicItem(msg: AddDynamicItem, state?: StateTree) {
    this.checkState(msg, state);
    return state.addSubTree(msg.parentUuid, msg.itemId, msg.position);
  }

  private loadDynamicItem(msg: LoadDynamicItem, state?: StateTree) {
    this.checkState(msg, state);
    return state.loadSubTree(msg.parentUuid, msg.dbId, msg.itemId, msg.position, !!msg.readonly);
  }

  public saveDynamicItem(msg: SaveDynamicItem, state?: StateTree) {
    this.checkState(msg, state);
    return state.save(msg.uuid);
  }

  private removeDynamicItem(msg: RemoveDynamicItem, state?: StateTree) {
    this.checkState(msg, state);
    return state.removeSubtree(msg.uuid);
  }

  private moveDynamicItem(msg: MoveDynamicItem, state?: StateTree) {
    this.checkState(msg, state);
    return state.moveSubtree(msg.uuid, msg.position);
  }

  private runStep(msg: RunStep, state?: StateTree) {
    this.checkState(msg, state);
    return state.runStep(msg.uuid, msg.mockResults, msg.mockDelay);
  }

  private runAction(msg: RunAction, state?: StateTree) {
    this.checkState(msg, state);
    return state.runAction(msg.uuid);
  }

  private runSequence(msg: RunSequence, state?: StateTree) {
    this.checkState(msg, state);
    return state.runAction(msg.startUuid);
  }

  public savePipeline(msg: SavePipeline, state?: StateTree) {
    this.checkState(msg, state);
    return state.save();
  }

  public loadPipeline(msg: LoadPipeline) {
    return from(loadInstanceState(msg.funcCallId)).pipe(
      concatMap((stateLoaded) => {
        if (isFuncCallSerializedState(stateLoaded))
          throw new Error(`Wrong pipeline config in wrapper FuncCall ${msg.funcCallId}`);
        if (!stateLoaded.provider)
          throw new Error(`Pipeline config in wrapper FuncCall ${msg.funcCallId} missing provider`);
        if (msg.config)
          return of([stateLoaded, msg.config] as const);
        return callHandler<PipelineConfiguration>(stateLoaded.provider, {version: stateLoaded.version}).pipe(
          concatMap((conf) => from(getProcessedConfig(conf))),
          map((config) => [stateLoaded, config] as const),
        );
      }),
      map(([state, config]) =>
        StateTree.fromInstanceState({state, config, isReadonly: !!msg.readonly, defaultValidators: false, mockMode: this.mockMode})),
      concatMap((state) => state.init()),
      tap((state) => this.states$.next(state)),
    );
  }

  private initPipeline(msg: InitPipeline) {
    return callHandler<PipelineConfiguration>(msg.provider, {version: msg.version}).pipe(
      concatMap((conf) => from(getProcessedConfig(conf))),
      map((config) => StateTree.fromPipelineConfig({config, isReadonly: false, defaultValidators: false, mockMode: this.mockMode})),
      concatMap((state) => state.init()),
      tap((state) => this.states$.next(state)),
    );
  }

  private checkState(msg: ViewConfigCommands, state?: StateTree): asserts state is StateTree {
    if (!state)
      throw new Error(`Command ${msg.event} requires state tree`);
  }
}

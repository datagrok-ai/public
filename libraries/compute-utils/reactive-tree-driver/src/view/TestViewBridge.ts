import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, EMPTY, Subject, merge, of} from 'rxjs';
import {HooksRunner} from './HooksRunner';
import {ActionTriggered, Actions, FuncCallLoaded, ViewConfig, ViewConfigChanges, ViewValidationResult, isStepFunCallConfig, traverseViewConfig} from './ViewCommunication';
import {IViewBridge} from './IViewBridge';
import {delay, switchMap, takeUntil, map} from 'rxjs/operators';

export class TestViewBridge implements IViewBridge {
  private destroyed$ = new Subject<true>();

  actionsInput = new BehaviorSubject<Actions | undefined>(undefined);
  actionsEvents = new Subject<ActionTriggered>();

  configInput = new BehaviorSubject<ViewConfig | undefined>(undefined);
  configChanges = new Subject<ViewConfigChanges>();

  blockRunInput = new BehaviorSubject(false);
  validationInput = new BehaviorSubject<ViewValidationResult | undefined>(undefined);

  constructor(public hooksRunner: HooksRunner, public mockCalls: Record<string, DG.FuncCall>) {
    this.configInput.pipe(
      takeUntil(this.destroyed$),
      switchMap((config) => {
        if (!config)
          return EMPTY;
        return this.loadAllCalls(config);
      }),
    ).subscribe((payload) => this.configChanges.next({event: 'funcCallLoaded', ...payload} as FuncCallLoaded));
  }

  private loadAllCalls(config: ViewConfig) {
    const callsIds = traverseViewConfig(config, (acc, node) => {
      if (isStepFunCallConfig(node))
        return [...acc, {id: node.id, funcCallId: node.funcCallId!}];

      return acc;
    }, [] as {id: string, funcCallId: string}[]);
    return merge(...callsIds.map((data) => this.loadCall(data.id, data.funcCallId)));
  }

  private loadCall(id: string, funcCallId: string) {
    return of(this.mockCalls[funcCallId]).pipe(
      delay(100),
      map((funcCall) => ({id, funcCall})),
    );
  }

  destroy() {
    this.destroyed$.next(true);
    this.destroyed$.complete();
  }
}

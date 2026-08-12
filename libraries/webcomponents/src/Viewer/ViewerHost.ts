import {BehaviorSubject, Subject, Observable, SchedulerLike, EMPTY, asapScheduler, combineLatest} from 'rxjs';
import {concatMap, debounceTime, takeUntil, tap} from 'rxjs/operators';

export interface ViewerLike {
  type: string;
  dataFrame: any;
  setOptions(options: Record<string, any>): void;
}

export type ViewerFactory<V extends ViewerLike> = (type: string, dataFrame: any) => Observable<V>;

/**
 * Framework-agnostic state core of the dg-viewer element. Collects type/dataFrame/options
 * set within one tick into a single transition (microtask-scheduled), so the order in which
 * the host framework patches the properties is irrelevant. A transition applies the frame
 * first and the final options second in one synchronous block: the frame swap's Dart look
 * refresh is deferred, so the corrected look always lands before it runs, while options
 * apply eagerly against the frame they were derived from. Stale options are never applied —
 * re-asserting them against a new frame is what crashed the Dart legend refresh.
 */
export class ViewerHost<V extends ViewerLike> {
  readonly typeSetted$ = new BehaviorSubject<string | undefined>(undefined);
  readonly dfSetted$ = new BehaviorSubject<any>(undefined);
  readonly optionsSetted$ = new BehaviorSubject<Record<string, any>>({});
  readonly viewer$ = new BehaviorSubject<V | undefined>(undefined);
  // fires after an in-place frame swap is applied to the attached viewer: application is
  // asynchronous (batched per tick), so consumers doing imperative post-swap work need a
  // signal that the viewer actually holds the new frame (creation is signaled by viewer$)
  readonly frameApplied$ = new Subject<V>();

  private appliedOptionsJson = '{}';
  private destroyed$ = new Subject<boolean>();

  constructor(
    private readonly createViewer: ViewerFactory<V>,
    scheduler: SchedulerLike = asapScheduler,
  ) {
    combineLatest([this.typeSetted$, this.dfSetted$, this.optionsSetted$]).pipe(
      debounceTime(0, scheduler),
      concatMap(() => this.applyTransition()),
      takeUntil(this.destroyed$),
    ).subscribe();
  }

  destroy() {
    this.destroyed$.next(true);
  }

  private applyTransition(): Observable<unknown> {
    const type = this.typeSetted$.value;
    const df = this.dfSetted$.value;
    const options = this.optionsSetted$.value;
    const optionsJson = JSON.stringify(options);
    const optionsChanged = optionsJson !== this.appliedOptionsJson;
    this.appliedOptionsJson = optionsJson;
    const viewer = this.viewer$.value;

    if (!type || !df) {
      this.viewer$.next(undefined);
      return EMPTY;
    }
    if (!viewer || viewer.type !== type) {
      return this.createViewer(type, df).pipe(
        tap((next: V) => {
          // options may have moved on while the viewer was being created; the latest
          // value is applied here and the snapshot updated so it is not re-applied
          this.appliedOptionsJson = JSON.stringify(this.optionsSetted$.value);
          next.setOptions(this.optionsSetted$.value);
          this.viewer$.next(next);
        }),
      );
    }
    if (viewer.dataFrame !== df) {
      viewer.dataFrame = df;
      viewer.setOptions(options);
      this.frameApplied$.next(viewer);
    } else if (optionsChanged)
      viewer.setOptions(options);
    return EMPTY;
  }
}

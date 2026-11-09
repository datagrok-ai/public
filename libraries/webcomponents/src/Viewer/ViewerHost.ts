import {BehaviorSubject, Subject, Observable, SchedulerLike, EMPTY, asapScheduler, merge, defer} from 'rxjs';
import {catchError, concatMap, debounceTime, takeUntil, tap} from 'rxjs/operators';

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
 * apply eagerly against the frame they were derived from. Stale options are never applied to a
 * new frame — re-asserting them against it is what crashed the Dart legend refresh.
 * Transitions read the latest state when they run, so ones queued behind a slow creation collapse.
 */
export class ViewerHost<V extends ViewerLike> {
  readonly typeSetted$ = new BehaviorSubject<string | undefined>(undefined);
  readonly dfSetted$ = new BehaviorSubject<any>(undefined);
  readonly optionsSetted$ = new BehaviorSubject<Record<string, any>>({});
  readonly viewer$ = new BehaviorSubject<V | undefined>(undefined);

  // compared instead of viewer.dataFrame: Dart may return a different wrapper for the same
  // frame, and a frame set on the viewer directly must survive unrelated transitions
  private appliedFrame: any;
  private appliedOptionsJson = '{}';
  private destroyed$ = new Subject<boolean>();

  constructor(
    private readonly createViewer: ViewerFactory<V>,
    scheduler: SchedulerLike = asapScheduler,
  ) {
    merge(this.typeSetted$, this.dfSetted$, this.optionsSetted$).pipe(
      debounceTime(0, scheduler),
      // a failed transition (unknown type, a throwing frame or options apply) must not end the stream
      concatMap(() => defer(() => this.applyTransition()).pipe(catchError((e) => {
        console.error(e);
        return EMPTY;
      }))),
      takeUntil(this.destroyed$),
    ).subscribe();
  }

  destroy() {
    this.destroyed$.next(true);
  }

  private applyTransition(): Observable<unknown> {
    const type = this.typeSetted$.value;
    const df = this.dfSetted$.value;
    const viewer = this.viewer$.value;

    if (!type || !df) {
      if (viewer)
        this.viewer$.next(undefined);
      return EMPTY;
    }
    if (!viewer || viewer.type !== type) {
      return this.createViewer(type, df).pipe(
        tap((next: V) => {
          // options may have moved on while the viewer was being created; the latest
          // value is applied here so the queued transition does not re-apply it
          this.applyOptions(next, true);
          this.appliedFrame = df;
          this.viewer$.next(next);
        }),
      );
    }
    const frameChanged = this.appliedFrame !== df;
    if (frameChanged) {
      viewer.dataFrame = df;
      this.appliedFrame = df;
    }
    this.applyOptions(viewer, frameChanged);
    return EMPTY;
  }

  private applyOptions(viewer: V, force: boolean) {
    const options = this.optionsSetted$.value;
    const optionsJson = JSON.stringify(options);
    if (!force && optionsJson === this.appliedOptionsJson)
      return;
    viewer.setOptions(options);
    this.appliedOptionsJson = optionsJson;
  }
}

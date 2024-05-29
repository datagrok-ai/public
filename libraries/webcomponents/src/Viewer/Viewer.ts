/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject, BehaviorSubject, from, merge, of, combineLatest, Observable, identity, EMPTY} from 'rxjs';
import {
  distinctUntilChanged, filter, switchMap, takeUntil, withLatestFrom, map, debounceTime, take,
} from 'rxjs/operators';

export class Viewer<T = any> extends HTMLElement {
  private valueSetted$ = new BehaviorSubject<DG.DataFrame | undefined>(undefined);
  private viewerSetted$ = new BehaviorSubject<DG.Viewer<T> | undefined>(undefined);
  private nameSetted$ = new BehaviorSubject<string | undefined>(undefined);

  private viewer$ = new BehaviorSubject<DG.Viewer<T> | undefined>(undefined);

  private destroyed$ = new Subject<boolean>();

  constructor() {
    super();

    const latestViewer$ = this.viewer$.pipe(distinctUntilChanged(), filter((v) => !!v));

    const latestName$ = merge(
      latestViewer$.pipe(map((viewer) => viewer?.type)),
      this.nameSetted$,
    ).pipe(distinctUntilChanged());

    const latestValue$ = merge(
      latestViewer$.pipe(map((viewer) => viewer?.dataFrame)),
      this.getViewerEventObservable('d4-data-frame-changed').pipe(map((ev) => ev.data.args.newValue as DG.DataFrame)),
      this.valueSetted$,
    ).pipe(distinctUntilChanged());

    const latestParams$ = combineLatest([
      latestName$,
      latestValue$,
    ] as const);

    merge(
      latestParams$,
      this.viewerSetted$,
    ).pipe(
      switchMap((payload) => {
        if (Array.isArray(payload)) {
          const [name, value] = payload;
          if (name && value) {
            if (this.viewer?.type !== name || !this.viewer)
              return from(this.createViewer(name, value));
            else
              return EMPTY;
          } else
            return of(undefined);
        } else {
          const viewer = payload;
          return of(viewer);
        }
      }),
      takeUntil(this.destroyed$),
    ).subscribe((viewer) => {
      this.changeAttachedViewer(viewer as DG.Viewer<T>);
    });

    this.valueSetted$.pipe(
      withLatestFrom(this.viewer$),
      takeUntil(this.destroyed$),
    ).subscribe(([value, viewer]) => {
      if (viewer && value)
        viewer.dataFrame = value;
    });

    this.viewer$.pipe(
      distinctUntilChanged(),
      takeUntil(this.destroyed$),
    ).subscribe((viewer) => {
      this.dispatchEvent(new CustomEvent('viewer-changed', {detail: viewer}));
    });

    // try to make sure that detached is the last event
    this.getViewerEventObservable().pipe(
      debounceTime(1000),
      filter((ev) => ev.type === 'd4-viewer-detached'),
      take(1),
    ).subscribe(() => {
      this.destroyed$.next(true);
    });
  }

  connectedCallback() {
  }

  disconnectedCallback() {
  }

  get viewer() {
    return this.viewer$.value;
  }

  set viewer(viewer: DG.Viewer<T> | undefined) {
    this.viewerSetted$.next(viewer);
  }

  get value() {
    return this.viewer$.value?.dataFrame ?? this.valueSetted$.value;
  }

  set value(df: DG.DataFrame | undefined) {
    this.valueSetted$.next(df);
  }

  get name() {
    return this.viewer$.value?.type ?? this.nameSetted$.value;
  }

  set name(name: string | undefined) {
    this.nameSetted$.next(name);
  }

  public getViewerEventObservable<P = any>(eventType?: string): Observable<P> {
    return this.viewer$.pipe(
      switchMap((viewer) => viewer ? viewer.onEvent().pipe(
        eventType ? filter((ev) => ev.type === eventType) : identity,
      ) : of()),
      filter((x) => x),
    );
  }

  private async createViewer(name: string, df: DG.DataFrame) {
    const viewer = await df.plot.fromType(name) as DG.Viewer<T>;
    return viewer;
  }

  private changeAttachedViewer(viewer?: DG.Viewer<T>) {
    this.innerHTML = '';
    this.viewer$.next(viewer);
    if (viewer)
      this.appendChild(viewer.root);
  }
}

export interface ViewerT extends Viewer {};

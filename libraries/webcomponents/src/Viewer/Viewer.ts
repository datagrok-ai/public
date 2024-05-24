/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject, BehaviorSubject, from, merge, of, combineLatest} from 'rxjs';
import {distinctUntilChanged, filter, switchMap, takeUntil, withLatestFrom, map} from 'rxjs/operators';

export class Viewer<T = any> extends HTMLElement {
  private valueSetted$ = new BehaviorSubject<DG.DataFrame | undefined>(undefined);
  private viewerSetted$ = new BehaviorSubject<DG.Viewer<T> | undefined>(undefined);
  private nameSetted$ = new BehaviorSubject<string | undefined>(undefined);

  private viewer$ = new BehaviorSubject<DG.Viewer<T> | undefined>(undefined);

  private isConnected$ = new BehaviorSubject<boolean>(false);
  private destroyed$ = new Subject<boolean>();

  constructor() {
    super();

    const providedViewer$ = this.viewerSetted$.pipe(distinctUntilChanged());

    const newValueForViewer$ = this.valueSetted$.pipe(
      distinctUntilChanged(),
      withLatestFrom(this.viewer$),
    );

    const newValueNoViewer$ = newValueForViewer$.pipe(
      filter(([, viewer]) => !viewer),
      map(([val]) => val));

    const newName$ = this.nameSetted$.pipe(distinctUntilChanged());

    const newViewer$ = combineLatest([
      newName$,
      newValueNoViewer$,
    ]).pipe(
      filter((vals) => vals.every((x) => x)),
    );

    // creating a new or replacing an existing one
    merge(
      newViewer$,
      providedViewer$,
    ).pipe(
      switchMap((payload) => {
        if (Array.isArray(payload)) {
          const [name, value] = payload;
          return from(this.createViewer(name!, value!));
        } else
          return of(payload);
      }),
      takeUntil(this.destroyed$),
    ).subscribe((viewer) => {
      this.viewer$.next(viewer);
      this.changeAttachedViewer(viewer);
    });

    let changingDf = false;

    // replacing viewers DataFrame
    newValueForViewer$.pipe(
      filter((vals) => vals.every((x) => x)),
      takeUntil(this.destroyed$),
    ).subscribe(([value, viewer]) => {
      changingDf = true;
      try {
        viewer!.dataFrame = value!;
      } finally {
        changingDf = false;
      }
    });

    this.viewer$.pipe(
      switchMap((viewer) => viewer ? viewer.onEvent().pipe(filter((ev) => ev.type === 'd4-viewer-detached')) : of()),
      filter((val) => val && !changingDf),
    ).subscribe(() => {
      this.destroyed$.next(true);
      console.log('Viewer webcomponent destroyed');
    });
  }


  connectedCallback() {
    this.isConnected$.next(true);
  }

  disconnectedCallback() {
    this.isConnected$.next(false);
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

  private async createViewer(name: string, df: DG.DataFrame) {
    const viewer = await df.plot.fromType(name) as DG.Viewer<T>;
    return viewer;
  }

  private changeAttachedViewer(viewer?: DG.Viewer<T>) {
    this.innerHTML = '';
    if (viewer)
      this.appendChild(viewer.root);
  }
}

export interface ViewerT extends Viewer {};

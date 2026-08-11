/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject, from, merge, of, combineLatest, Observable, identity} from 'rxjs';
import {distinctUntilChanged, filter, switchMap, takeUntil, map} from 'rxjs/operators';
import {ViewerHost} from './ViewerHost';

export class Viewer<T = any> extends HTMLElement {
  private host = new ViewerHost<DG.Viewer<T>>(
    (type, df) => from((df as DG.DataFrame).plot.fromType(type) as Promise<DG.Viewer<T>>));

  private destroyed$ = new Subject<boolean>();

  constructor() {
    super();

    const latestViewer$ = this.host.viewer$.pipe(distinctUntilChanged());

    const latestType$ = merge(
      latestViewer$.pipe(map((viewer) => viewer?.type)),
      this.host.typeSetted$,
    ).pipe(distinctUntilChanged());

    const latestDf$ = merge(
      latestViewer$.pipe(map((viewer) => viewer?.dataFrame)),
      this.getViewerEventObservable('d4-data-frame-changed').pipe(map((ev) => ev.data.args.newValue as DG.DataFrame)),
      this.host.dfSetted$,
    ).pipe(distinctUntilChanged());

    combineLatest([latestType$, latestDf$]).pipe(
      takeUntil(this.destroyed$),
    ).subscribe(([type, df]) => {
      if (type !== this.host.typeSetted$.value)
        this.dispatchEvent(new CustomEvent('viewer-type-changed', {detail: type}));
      // compare dart handles: Dart events deliver a fresh JS wrapper for the same frame,
      // and a spurious dispatch here runs inside the Dart event fire, where consumers
      // mutating the grid hit "Controller is already firing an event"
      if (df?.dart !== this.host.dfSetted$.value?.dart)
        this.dispatchEvent(new CustomEvent('viewer-data-frame-changed', {detail: df}));
    });

    this.host.viewer$.pipe(
      takeUntil(this.destroyed$),
    ).subscribe((viewer) => {
      ui.empty(this);
      if (viewer)
        this.appendChild(viewer.root);
    });

    // frame application is asynchronous (batched per tick), so an applied in-place swap
    // is announced — consumers doing imperative post-swap work listen to this event
    this.host.frameApplied$.pipe(
      takeUntil(this.destroyed$),
    ).subscribe((viewer) => {
      this.dispatchEvent(new CustomEvent('viewer-data-frame-changed', {detail: viewer.dataFrame}));
    });

    latestViewer$.pipe(
      takeUntil(this.destroyed$),
    ).subscribe((viewer) => {
      this.dispatchEvent(new CustomEvent('viewer-changed', {detail: viewer}));
    });
  }

  public destroy() {
    this.host.destroy();
    this.destroyed$.next(true);
    ui.empty(this);
  }

  connectedCallback() {
  }

  disconnectedCallback() {
  }

  get viewer() {
    return this.host.viewer$.value;
  }

  get dataFrame() {
    return this.host.viewer$.value?.dataFrame ?? this.host.dfSetted$.value;
  }

  set dataFrame(df: DG.DataFrame | undefined) {
    this.host.dfSetted$.next(df);
  }

  get type() {
    return this.host.viewer$.value?.type ?? this.host.typeSetted$.value;
  }

  set type(type: string | undefined) {
    this.host.typeSetted$.next(type);
  }

  set options(options: Record<string, string | boolean>) {
    this.host.optionsSetted$.next(options);
  }

  public getViewerEventObservable<P = any>(eventType?: string): Observable<P> {
    return this.host.viewer$.pipe(
      switchMap((viewer) => viewer ? viewer.onEvent().pipe(
        eventType ? filter((ev) => ev.type === eventType) : identity,
      ) : of()),
      filter((x) => x),
    );
  }
}

export interface ViewerT extends Viewer {};

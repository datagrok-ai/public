/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject, from, EMPTY, asapScheduler} from 'rxjs';
import {distinctUntilChanged, switchMap, takeUntil, map, observeOn, skip, startWith} from 'rxjs/operators';
import {ViewerHost} from './ViewerHost';

export class Viewer<T = any> extends HTMLElement {
  private host = new ViewerHost<DG.Viewer<T>>(
    (type, df) => from((df as DG.DataFrame).plot.fromType(type) as Promise<DG.Viewer<T>>));

  private destroyed$ = new Subject<boolean>();

  constructor() {
    super();

    this.host.viewer$.pipe(
      takeUntil(this.destroyed$),
    ).subscribe((viewer) => {
      ui.empty(this);
      if (viewer)
        this.appendChild(viewer.root);
      this.dispatchEvent(new CustomEvent('viewer-changed', {detail: viewer}));
    });

    // the Dart event fires inside the frame setter, where consumers mutating the grid hit
    // "Controller is already firing an event"; frames are compared by dart handle because
    // Dart can return a different wrapper for the same frame; seeded with the creation frame:
    // creation and replacement are announced by viewer-changed only
    this.host.viewer$.pipe(
      switchMap((viewer) => viewer ? viewer.onDataFrameChanged.pipe(
        observeOn(asapScheduler),
        map(() => viewer.dataFrame),
        startWith(viewer.dataFrame),
        distinctUntilChanged((a, b) => a?.dart === b?.dart),
        skip(1),
      ) : EMPTY),
      takeUntil(this.destroyed$),
    ).subscribe((df) => {
      this.dispatchEvent(new CustomEvent('viewer-data-frame-changed', {detail: df}));
    });
  }

  public destroy() {
    this.host.destroy();
    this.destroyed$.next(true);
    ui.empty(this);
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

  set options(options: Record<string, any>) {
    this.host.optionsSetted$.next(options);
  }
}

export interface ViewerT extends Viewer {};

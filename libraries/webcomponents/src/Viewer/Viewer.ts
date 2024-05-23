/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {filter, switchMap, take} from 'rxjs/operators';

export class Viewer<T = any> extends HTMLElement {
  private df: DG.DataFrame = DG.DataFrame.fromColumns([DG.Column.float('empty', 0)]);
  private viewerInst?: DG.Viewer<T>;

  public viewer$ = new BehaviorSubject<DG.Viewer<T> | undefined>(undefined);

  constructor() {
    super();
    // TODO? closed event?
    // this.attached$.pipe(
    //   filter(val => val),
    //   switchMap(
    // 	() => this.viewerInst!.onEvent().pipe(filter((ev) => ev.type === 'd4-viewer-detached'))
    //   ),
    //   take(1)
    // );
  }

  connectedCallback() {
    this.init();
  }

  disconnectedCallback() {
    // ignoring, using d4-viewer-detached event instead
  }

  attributeChangedCallback(name: string, oldValue: string, newValue: string) {
    if (name === 'name' && oldValue !== newValue) {
      this.viewerInst = undefined;
      this.init();
    }
  }

  get viewer() {
    return this.viewerInst;
  }

  set viewer(viewer: DG.Viewer<T> | undefined) {
    this.viewerInst = viewer;
    if (this.viewerInst) {
      if (this.df)
        this.viewerInst.dataFrame = this.df
      this.attachViewer();
    } else {
      this.clear();
    }
  }

  get value() {
    return this.df;
  }

  set value(df: DG.DataFrame) {
    this.df = df;
    if (this.viewerInst)
      this.viewerInst.dataFrame = this.df;
  }

  private get viewerAttrName() {
    return this.getAttribute('name');
  }

  private async init() {
    if (this.viewerInst || !this.viewerAttrName)
      return;
    await this.createViewer();
    this.attachViewer();
  }

  private async createViewer() {
    // checking if a value has been updated during Ciewer initialization
    const oldVal = this.df;

    this.viewerInst = await this.df.plot.fromType(this.viewerAttrName!) as DG.Viewer;

    if (oldVal !== this.df)
      this.viewerInst.dataFrame = this.df;
  }

  private attachViewer() {
    this.clear();
    this.appendChild(this.viewerInst!.root);
    this.viewer$.next(this.viewerInst!);
  }

  private clear() {
    this.innerHTML = '';
  }
}

export interface ViewerT extends Viewer {};

import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {isLeaf, NodeType} from '@datagrok-libraries/bio';
import * as rxjs from 'rxjs';
import {Unsubscribable} from 'rxjs';
import {HoverType} from './markup';

export type TreeRendererEventArgsType<TNode extends NodeType, THover extends HoverType<TNode>> = {
  target: TreeRendererBase<TNode, THover>,
  context: CanvasRenderingContext2D,
  lengthRatio: number,
};


export abstract class TreeRendererBase<TNode extends NodeType, THover extends HoverType<TNode>> {
  public view?: HTMLElement;

  private _treeRoot: TNode;
  get treeRoot(): TNode { return this._treeRoot; }

  set treeRoot(value: TNode) {
    this._treeRoot = value;

    this._current = null;
    this._mouseOver = null;
    this._selections = [];

    this.render('set treeRoot()');
  }

  //#region current

  public get currentNode(): TNode | null { return this._current ? this._current.node : null; }

  protected _onCurrentChanged: rxjs.Subject<void> = new rxjs.Subject<void>();
  get onCurrentChanged(): rxjs.Observable<void> { return this._onCurrentChanged; }

  protected _current: THover | null = null;
  public get current(): THover | null { return this._current; }

  public set current(value: THover | null) {
    if (
      (value != null && this._current == null) ||
      (value == null && this._current != null) ||
      (value != null && this._current != null && value.node != this._current.node)
    ) {
      this._current = value;
      this.render('set current()');
      this._onCurrentChanged.next();
    }
  }

  //#endregion current

  //#region mouseOver

  public get mouseOverNode(): TNode | null { return this._mouseOver ? this._mouseOver.node : null; }

  protected _onMouseOverChanged: rxjs.Subject<void> = new rxjs.Subject<void>();
  get onMouseOverChanged(): rxjs.Observable<void> { return this._onMouseOverChanged; }

  protected _mouseOver: THover | null = null;
  public get mouseOver(): THover | null { return this._mouseOver; }

  public set mouseOver(value: THover | null) {
    if (
      (value != null && this._mouseOver == null) ||
      (value == null && this._mouseOver != null) ||
      (value != null && this._mouseOver != null && value.node != this._mouseOver.node)
    ) {
      const msg = `this._mouseOver.node.name = '${this._mouseOver ? this._mouseOver.node.name : '<null>'}', ` +
        `value.node.name = '${value ? value.node.name : '<null>'}'`;
      this._mouseOver = value;
      this.render('set mouseOver() ' + msg);
      this._onMouseOverChanged.next();
    }
  }

  //#endregion mouseOver

  //#region selections
  protected _selectedNodes: TNode[] = [];
  get selectedNodes(): TNode[] { return this._selectedNodes; }

  protected _onSelectionChanged: rxjs.Subject<void> = new rxjs.Subject<void>();
  get onSelectionChanged(): rxjs.Observable<void> { return this._onSelectionChanged; }

  protected _selections: THover[] = [];
  public get selections(): THover[] { return this._selections; }

  public set selections(value: THover[]) {
    if (
      value.length != this.selections.length ||
      value.some((v, i) => v.node != this.selections[i].node)
    ) {
      this._selections = value;
      this._selectedNodes = this._selections.map((h) => h.node);
      this.render('set selections()');
      this._onSelectionChanged.next();
    }
  }

  //#endregion selections

  protected subs: Unsubscribable[] = [];

  protected constructor(treeRoot: TNode) {
    this._treeRoot = treeRoot;

    this._onAfterRender = new rxjs.Subject<TreeRendererEventArgsType<TNode, THover>>();
  }

  // -- View --

  public attach(view: HTMLElement): void {
    this.view = view;
    this.subs.push(ui.onSizeChanged(this.view).subscribe(this.viewOnSizeChanged.bind(this)));

    // Postponed call to initial render after attach completed (e.g. canvas created)
    window.setTimeout(() => {
      this.viewOnSizeChanged(); // Calls this.render()
      // this.render('attach()');
    }, 0 /* next event cycle */);
  }

  public detach(): void {
    for (const sub of this.subs) sub.unsubscribe();
    delete this.view;
  }

  public abstract render(purpose: string): void;

  protected readonly _onAfterRender: rxjs.Subject<TreeRendererEventArgsType<TNode, THover>>;

  get onAfterRender(): rxjs.Observable<TreeRendererEventArgsType<TNode, THover>> { return this._onAfterRender; }

  // -- Handle events --

  protected viewOnSizeChanged(): void {
    console.debug('Dendrogram: TreeRendererBase.viewOnSizeChanged()');

    this.render('viewOnSizeChanged');
  }
}


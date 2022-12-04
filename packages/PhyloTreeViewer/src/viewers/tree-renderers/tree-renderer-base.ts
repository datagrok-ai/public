import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {isLeaf, NodeType} from '@datagrok-libraries/bio';
import * as rxjs from 'rxjs';
import {GridTreeRendererBase} from '../grid-tree-renderer';
import {Unsubscribable} from 'rxjs';

export type TreeRendererEventArgsType<TNode extends NodeType> = {
  target: TreeRendererBase<TNode>,
  context: CanvasRenderingContext2D,
  lengthRatio: number,
};


export abstract class TreeRendererBase<TNode extends NodeType> {

  private _treeRoot: TNode;
  private readonly _view: HTMLElement;

  /** Tree hierarchical structure */
  protected _totalLength: number;

  get treeRoot(): TNode { return this._treeRoot; }

  set treeRoot(value: TNode) {
    this._treeRoot = value;
    this.render();
  }

  get totalLength(): number { return this._totalLength; }

  get view(): HTMLElement { return this._view; }

  protected subs: Unsubscribable[] = [];

  protected constructor(treeRoot: TNode, totalLength: number, view: HTMLElement) {
    this._treeRoot = treeRoot;
    this._view = view;
    this._totalLength = totalLength;

    this.subs.push(
      ui.onSizeChanged(this.view).subscribe(this.viewOnSizeChanged.bind(this)));

    this._onAfterRender = new rxjs.Subject<TreeRendererEventArgsType<TNode>>();
  }

  public abstract render(): void;

  protected readonly _onAfterRender: rxjs.Subject<TreeRendererEventArgsType<TNode>>;

  get onAfterRender(): rxjs.Observable<TreeRendererEventArgsType<TNode>> { return this._onAfterRender; }

  // -- Handle events --

  private viewOnSizeChanged(): void {
    console.debug('PhyloTreeViewer: TreeRendererBase.viewOnSizeChanged()');

    this.render();
  }
}


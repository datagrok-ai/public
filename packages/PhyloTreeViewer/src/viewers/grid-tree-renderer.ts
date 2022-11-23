import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable, Subject, Subscribable, Unsubscribable} from 'rxjs';

import {NodeType, isLeaf, ITreeHelper} from '@datagrok-libraries/bio';

import {TreeHelper} from '../utils/tree-helper';

export type GridTreeRendererEventArgsType<TNode extends NodeType> = {
  target: GridTreeRendererBase<TNode>,
  context: CanvasRenderingContext2D,
  lengthRatio: number,
};

export abstract class GridTreeRendererBase<TNode extends NodeType> {
  protected _leftPadding: number = 6;
  protected _rightPadding: number = 6;

  private readonly th: ITreeHelper;

  /** Tree hierarchical structure */
  protected _tree: TNode;
  protected _totalLength: number;

  protected readonly _view: HTMLElement;
  protected readonly _grid: DG.Grid;

  protected subs: Unsubscribable[] = [];

  get tree(): TNode { return this._tree; }

  set tree(value: TNode) {
    this._tree = value;
    this.render();
  }

  get totalLength(): number { return this._totalLength; }

  get view(): HTMLElement { return this._view; }

  get grid(): DG.Grid { return this._grid; }

  get leftPadding(): number { return this._leftPadding; }

  get rightPadding(): number { return this._rightPadding; }

  protected readonly _onAfterRender = new Subject<GridTreeRendererEventArgsType<TNode>>();

  get onAfterRender(): Observable<GridTreeRendererEventArgsType<TNode>> { return this._onAfterRender; }

  protected constructor(tree: TNode, totalLength: number, treeDiv: HTMLElement, grid: DG.Grid) {
    this.th = new TreeHelper();
    this._tree = tree;
    this._view = treeDiv;
    this._grid = grid;
    this._totalLength = totalLength;
    this.subs.push(ui.onSizeChanged(this.view).subscribe(this.render.bind(this)));
    this.subs.push(this.grid.onBeforeDrawContent.subscribe(this.render.bind(this)));
    this.subs.push(ui.onSizeChanged(this.grid.root).subscribe(this.render.bind(this)));
  }

  public abstract render(): void;
}

export type MarkupNodeType = NodeType & {
  children: MarkupNodeType[],
  index: number,
  minIndex: number,
  maxIndex: number,
  /** node's branch_length with max subtreeLength of children */
  subtreeLength?: number,
};

/** Draws only nodes/leaves visible in leaf range */
export class LeafRangeGridTreeRenderer extends GridTreeRendererBase<MarkupNodeType> {
  private readonly canvas: HTMLCanvasElement;

  constructor(tree: MarkupNodeType, totalLength: number, treeDiv: HTMLElement, grid: DG.Grid) {
    super(tree, totalLength, treeDiv, grid);

    this.view.style.setProperty('overflow-y', 'hidden', 'important');

    this.canvas = ui.canvas();
    this.canvas.style.position = 'absolute';
    this.view.appendChild(this.canvas);
  }

  public render(): void {
    const firstRowIndex: number = Math.floor(this._grid.vertScroll.min);
    const rowsGridHeight: number = this._grid.root.clientHeight - this._grid.colHeaderHeight;
    const lastRowIndex: number = firstRowIndex + Math.ceil(rowsGridHeight / this._grid.props.rowHeight);

    console.debug('PhyloTreeViewer: LeafRangeTreeRenderer.render() ' +
      `firstRowIndex = ${firstRowIndex}, lastRowIndex = ${lastRowIndex}`);

    const cw: number = this.view.clientWidth;

    // this.treeDiv.clientHeight is incorrect on start without resize, workaround using this.grid.root.clientHeight
    this.view.style.height = `${this.grid.root.clientHeight}px`;
    const ch: number = this.view.clientHeight - this.grid.colHeaderHeight;
    //const ch: number = this.grid.root.clientHeight - this.grid.colHeaderHeight;

    this.canvas.width = cw * window.devicePixelRatio;
    this.canvas.height = ch * window.devicePixelRatio;
    this.canvas.style.left = `${0}px`;
    this.canvas.style.top = `${this.grid.colHeaderHeight}px`;
    this.canvas.style.width = `${cw}px`;
    this.canvas.style.height = `${ch}px`;

    // TODO:
    const ctx = this.canvas.getContext('2d')!;

    // Here we will render range of leaves
    const lengthRatio: number = window.devicePixelRatio * (cw - this.leftPadding - this.rightPadding) / this.totalLength; // px/[length unit]
    const stepRatio: number = window.devicePixelRatio * this.grid.props.rowHeight; // px/[step unit, row]

    ctx.save();
    try {
      ctx.fillStyle = '#FFFFFF';
      ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);

      // ctx.beginPath();
      // ctx.strokeStyle = '#800000';
      // ctx.moveTo(0, 0);
      // ctx.lineTo(ctx.canvas.width, ctx.canvas.height);
      // ctx.moveTo(0, ctx.canvas.height);
      // ctx.lineTo(ctx.canvas.width, 0);
      // ctx.stroke();

      LeafRangeGridTreeRenderer.renderNode(ctx, this.tree as MarkupNodeType,
        firstRowIndex, lastRowIndex, this.leftPadding, lengthRatio, stepRatio,
        this.tree.subtreeLength!, 0);
    } finally {
      ctx.restore();
      this._onAfterRender.next({target: this, context: ctx, lengthRatio,});
    }
  }

  /**
   * @param ctx
   * @param node
   * @param firstRowIndex
   * @param lastRowIndex
   * @param leftPadding
   * @param lengthRatio
   * @param stepRatio
   * @param {number} totalLength Total (whole) tree length (height)
   * @param currentLength
   * @private
   */
  private static renderNode(ctx: CanvasRenderingContext2D, node: MarkupNodeType,
    firstRowIndex: number, lastRowIndex: number,
    leftPadding: number, lengthRatio: number, stepRatio: number,
    totalLength: number, currentLength: number = 0
  ) {
    const r: number = window.devicePixelRatio;

    if (isLeaf(node)) {
      if (firstRowIndex <= node.index && node.index <= lastRowIndex) {
        const minX = currentLength * lengthRatio + leftPadding * r;
        const maxX = (currentLength + node.branch_length!) * lengthRatio + leftPadding * r;

        // plot leaf grid
        const posY = (node.index - firstRowIndex + 0.5) * stepRatio;
        ctx.beginPath();
        ctx.strokeStyle = '#C0C0C0';
        ctx.lineWidth = 1;
        ctx.moveTo(maxX, posY);
        ctx.lineTo(ctx.canvas.width, posY);
        ctx.stroke();

        // plot branch
        ctx.beginPath();
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 1;
        ctx.moveTo(minX, posY);
        ctx.lineTo(maxX, posY);
        ctx.stroke();


        // plot leaf (marker?)
        ctx.beginPath();
        ctx.fillStyle = 'black';
        ctx.ellipse(maxX, posY, 1.5, 1.5, 0, 0, 2 * Math.PI);
        ctx.fill();
      }
    } else {
      if (firstRowIndex <= node.maxIndex && node.minIndex <= lastRowIndex) {
        for (const childNode of node.children) {
          LeafRangeGridTreeRenderer.renderNode(ctx, childNode,
            firstRowIndex, lastRowIndex,
            leftPadding, lengthRatio, stepRatio,
            totalLength, currentLength + node.branch_length!);
        }

        // plot join
        const joinMinIndex = node.children[0].index;
        const joinMaxIndex = node.children[node.children.length - 1].index;
        const posX = (currentLength + node.branch_length!) * lengthRatio + leftPadding * r;
        const minY = Math.max((joinMinIndex - firstRowIndex + 0.5) * stepRatio, 0);
        const maxY = Math.min((joinMaxIndex - firstRowIndex + 0.5) * stepRatio, ctx.canvas.height);
        //
        ctx.beginPath();
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 1;
        ctx.moveTo(posX, minY);
        ctx.lineTo(posX, maxY);
        ctx.stroke();

        const minX = currentLength * lengthRatio + leftPadding * r;
        const maxX = (currentLength + node.branch_length!) * lengthRatio + leftPadding * r;
        const posY = (node.index - firstRowIndex + 0.5) * stepRatio;

        ctx.beginPath();
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 1;
        ctx.moveTo(minX, posY);
        ctx.lineTo(maxX, posY);
        ctx.stroke();
      }
    }
  }

  /**
   *
   * @param node
   * @param currentLeafIndex
   * @return {number} Index pointing to the next leaf
   */
  public static markupNode(node: MarkupNodeType | NodeType, currentLeafIndex: number = 0): void {
    function markupNodeInt(node: MarkupNodeType, currentLeafIndex: number) {
      if (isLeaf(node)) {
        node.index = currentLeafIndex;
        node.subtreeLength = node.branch_length!;

        return currentLeafIndex + 1;
      } else {
        let maxSubtreeLength = 0;
        let leafIndex: number = currentLeafIndex;
        node.minIndex = leafIndex;
        for (const childNode of node.children!) {
          leafIndex = markupNodeInt(childNode, leafIndex);

          if (maxSubtreeLength < childNode.subtreeLength!) maxSubtreeLength = childNode.subtreeLength!;
        }
        node.maxIndex = leafIndex - 1; // leafIndex points to the next leaf already
        node.index = node.children.map((n) => (n as MarkupNodeType).index)
          .reduce((a, b) => a + b) / node.children.length;
        node.subtreeLength = maxSubtreeLength + node.branch_length!;

        return leafIndex;
      }
    }

    const t1: number = Date.now();
    markupNodeInt(node as MarkupNodeType, currentLeafIndex);
    const t2: number = Date.now();
    console.debug('PhyloTreeViewer: LeafRangeTreeRenderer.markupNode() ' + `ET: ${((t2 - t1) / 1000).toString()} s`);
  }

  public static create(
    tree: NodeType, treeDiv: HTMLElement, grid: DG.Grid
  ): GridTreeRendererBase<MarkupNodeType> {
    // TODO: adapt tree: bio.NodeType to MarkupNodeType
    LeafRangeGridTreeRenderer.markupNode(tree);
    const totalLength: number = (tree as MarkupNodeType).subtreeLength!;
    return new LeafRangeGridTreeRenderer(tree as MarkupNodeType, totalLength, treeDiv, grid);
  }
}

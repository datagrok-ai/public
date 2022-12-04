import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {NodeType} from '@datagrok-libraries/bio';
import {GridTreeRendererBase} from '../grid-tree-renderer';
import {markupNode, MarkupNodeType, renderNode} from './markup';

/** Draws only nodes/leaves visible in leaf range */
export class LeafRangeGridTreeRenderer extends GridTreeRendererBase<MarkupNodeType> {

  constructor(tree: MarkupNodeType, totalLength: number, treeDiv: HTMLElement, grid: DG.Grid) {
    super(tree, totalLength, treeDiv, grid);

    this.view.style.setProperty('overflow-y', 'hidden', 'important');

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

      renderNode(ctx, this.treeRoot as MarkupNodeType,
        firstRowIndex, lastRowIndex, this.leftPadding, lengthRatio, stepRatio,
        this.treeRoot.subtreeLength!, 0);
    } finally {
      ctx.restore();
      this._onAfterRender.next({target: this, context: ctx, lengthRatio,});
    }
  }

  public static create(
    tree: NodeType, treeDiv: HTMLElement, grid: DG.Grid
  ): GridTreeRendererBase<MarkupNodeType> {
    // TODO: adapt tree: bio.NodeType to MarkupNodeType
    markupNode(tree);
    const totalLength: number = (tree as MarkupNodeType).subtreeLength!;
    return new LeafRangeGridTreeRenderer(tree as MarkupNodeType, totalLength, treeDiv, grid);
  }
}
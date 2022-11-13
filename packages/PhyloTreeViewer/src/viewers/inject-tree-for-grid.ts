import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import $ from 'cash-dom';
import {TreeHelper} from '../utils/tree-helper';
import {attachDivToGrid} from './inject-tree-to-grid';
import {Unsubscribable} from 'rxjs';
import {maxIndex, sum} from 'd3';
import {first} from 'rxjs/operators';

export abstract class GridTreeRendererBase<TNode extends bio.NodeType> {
  private readonly th: bio.ITreeHelper;

  /** Tree hierarchical structure */
  protected readonly tree: TNode;
  protected readonly treeDiv: HTMLElement;
  protected readonly grid: DG.Grid;

  protected subs: Unsubscribable[] = [];

  protected constructor(tree: TNode, treeDiv: HTMLElement, grid: DG.Grid,) {
    this.th = new TreeHelper();
    this.tree = tree;
    this.treeDiv = treeDiv;
    this.grid = grid;
    this.subs.push(ui.onSizeChanged(this.treeDiv).subscribe(this.render.bind(this)));
    this.subs.push(this.grid.onBeforeDrawContent.subscribe(this.render.bind(this)));
    this.subs.push(ui.onSizeChanged(this.grid.root).subscribe(this.render.bind(this)));
  }

  public abstract render(): void;
}

export type MarkupNodeType = bio.NodeType & {
  children: MarkupNodeType[],
  index: number,
  minIndex: number,
  maxIndex: number,
  /** node's branch_length with max subtreeLength of children */
  subtreeLength?: number,
};

/** Draws only nodes/leaves visible in leaf range */
export class TreeForGridRenderer extends GridTreeRendererBase<MarkupNodeType> {
  private leftPadding: number = 4;
  private readonly canvas: HTMLCanvasElement;


  constructor(tree: MarkupNodeType, treeDiv: HTMLElement, grid: DG.Grid) {
    super(tree, treeDiv, grid);

    this.treeDiv.style.setProperty('overflow-y', 'hidden', 'important');

    this.canvas = ui.canvas();
    this.canvas.style.position = 'absolute';
    this.treeDiv.appendChild(this.canvas);
  }

  public render(): void {
    const firstRowIndex: number = Math.floor(this.grid.vertScroll.min);
    const rowsGridHeight: number = this.grid.root.clientHeight - this.grid.colHeaderHeight;
    const lastRowIndex: number = firstRowIndex + Math.ceil(rowsGridHeight / this.grid.props.rowHeight);

    console.debug('PhyloTreeViewer: LeafRangeTreeRenderer.render() ' +
      `firstRowIndex = ${firstRowIndex}, lastRowIndex = ${lastRowIndex}`);

    const cw: number = this.treeDiv.clientWidth;

    // this.treeDiv.clientHeight is incorrect on start without resize, workaround using this.grid.root.clientHeight
    this.treeDiv.style.height = `${this.grid.root.clientHeight}px`;
    const ch: number = this.treeDiv.clientHeight - this.grid.colHeaderHeight;
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
    const lengthRatio: number = (this.canvas.width - this.leftPadding) / this.tree.subtreeLength!; // px/[length unit]
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

      TreeForGridRenderer.renderNode(ctx, this.tree as MarkupNodeType,
        firstRowIndex, lastRowIndex, lengthRatio, stepRatio,
        this.tree.subtreeLength!, 0);
    } finally {
      ctx.restore();
    }
  }

  /**
   * @param ctx
   * @param node
   * @param firstRowIndex
   * @param lastRowIndex
   * @param {number} totalLength Total (whole) tree length (height)
   * @param currentLength
   * @private
   */
  private static renderNode(ctx: CanvasRenderingContext2D, node: MarkupNodeType,
    firstRowIndex: number, lastRowIndex: number, lengthRatio: number, stepRatio: number,
    totalLength: number, currentLength: number = 0
  ) {

    if (bio.isLeaf(node)) {
      if (firstRowIndex <= node.index && node.index <= lastRowIndex) {
        const minX = currentLength * lengthRatio;
        const maxX = (currentLength + node.branch_length!) * lengthRatio;

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
          TreeForGridRenderer.renderNode(ctx, childNode, firstRowIndex, lastRowIndex, lengthRatio, stepRatio,
            totalLength, currentLength + node.branch_length!);
        }

        // plot join
        const joinMinIndex = node.children[0].index;
        const joinMaxIndex = node.children[node.children.length - 1].index;
        const posX = (currentLength + node.branch_length!) * lengthRatio;
        const minY = Math.max((joinMinIndex - firstRowIndex + 0.5) * stepRatio, 0);
        const maxY = Math.min((joinMaxIndex - firstRowIndex + 0.5) * stepRatio, ctx.canvas.height);
        //
        ctx.beginPath();
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 1;
        ctx.moveTo(posX, minY);
        ctx.lineTo(posX, maxY);
        ctx.stroke();

        const minX = currentLength * lengthRatio;
        const maxX = (currentLength + node.branch_length!) * lengthRatio;
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
  public static markupNode(node: MarkupNodeType, currentLeafIndex: number = 0): void {
    function markupNodeInt(node: MarkupNodeType, currentLeafIndex: number) {
      if (bio.isLeaf(node)) {
        node.index = currentLeafIndex;
        node.subtreeLength = node.branch_length!;

        return currentLeafIndex + 1;
      } else {
        let maxSubtreeLength = 0;
        let leafIndex: number = currentLeafIndex;
        node.minIndex = leafIndex;
        for (const childNode of node.children!) {
          leafIndex = markupNodeInt(childNode as MarkupNodeType, leafIndex);

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
    markupNodeInt(node, currentLeafIndex);
    const t2: number = Date.now();
    console.debug('PhyloTreeViewer: LeafRangeTreeRenderer.markupNode() ' + `ET: ${((t2 - t1) / 1000).toString()} s`);
  }
}

export function gridTreeRenderer(
  tree: bio.NodeType, treeDiv: HTMLElement, grid: DG.Grid
): GridTreeRendererBase<MarkupNodeType> {
  // TODO: adapt tree: bio.NodeType to MarkupNodeType
  TreeForGridRenderer.markupNode(tree as MarkupNodeType);

  return new TreeForGridRenderer(tree as MarkupNodeType, treeDiv, grid);
}

export function injectTreeForGridUI(
  grid: DG.Grid, newickStr: string, leafColName: string, neighborWidth: number = 100
): GridNeighbor {
  const th: bio.ITreeHelper = new TreeHelper();

  const treeN = attachDivToGrid(grid, neighborWidth);
  const treeRoot = treeN.root!;

  // const treeDiv = ui.div();
  // treeRoot.appendChild(treeDiv);
  // treeRoot.style.backgroundColor = '#FFF0F0';
  // treeRoot.style.setProperty('overflow-y', 'hidden', 'important');

  const newickRoot: bio.NodeType = bio.Newick.parse_newick(newickStr);

  const treeRender: GridTreeRendererBase<MarkupNodeType> = gridTreeRenderer(newickRoot, treeRoot, grid);

  // grid.dataFrame.onFilterChanged.subscribe((args) => {
  //   console.debug('PhyloTreeViewer: injectTreeForGrid() grid.dataFrame.onFilterChanged()');
  //
  //   window.setTimeout(() => {
  //     const [viewedRoot] = th.setGridOrder(newickRoot, grid, leafColName);
  //
  //
  //   }, 0);
  // });

  return treeN;
}


import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

import {HoverType, ITreePlacer, MarkupNodeType} from './markup';
import {Subject} from 'rxjs';
import {isLeaf, NodeType} from '@datagrok-libraries/bio';

export type RectangleTreeHoverType<TNode extends NodeType> = HoverType<TNode> & {
  /** Node position along height axis */
  nodeHeight: number
};


// eslint-disable-next-line max-len
export class RectangleTreePlacer<TNode extends MarkupNodeType> implements ITreePlacer<TNode, RectangleTreeHoverType<TNode>> {
  private _top: number;
  get top(): number { return this._top; }

  private _bottom: number;
  get bottom() { return this._bottom; }

  get height() { return this._bottom - this._top; }

  /** Tree hierarchical structure */
  private _totalLength: number;
  get totalLength() { return this._totalLength; }

  get padding() { return {left: 8, right: 8}; }

  private readonly _onChanged: rxjs.Subject<void>;

  get onPlacingChanged(): rxjs.Observable<void> { return this._onChanged; }

  constructor(top: number, bottom: number, totalLength: number) {
    this._top = top;
    this._bottom = bottom;
    this._totalLength = totalLength;

    this._onChanged = new Subject<void>();
  }

  update(params: { top?: number, bottom?: number, totalLength?: number }): void {
    let changed: boolean = false;

    if (params.top && params.top != this.top) {
      this._top = params.top;
      changed = true;
    }

    if (params.bottom && params.bottom != this.bottom) {
      this._bottom = params.bottom;
      changed = true;
    }

    if (params.totalLength && params.totalLength != this.totalLength) {
      this._totalLength = params.totalLength;
      changed = true;
    }

    if (changed)
      this._onChanged.next();
  }

  /**
   * @param {NodeType} treeRoot
   * @param {DG.Point} point
   * @param {number} nodeSize Size of node in pixels (of styler)
   */
  getNode(
    treeRoot: TNode, canvasPoint: DG.Point, lineWidth: number, nodeSize: number,
    treeToCanvas: (treeP: DG.Point) => DG.Point
  ): RectangleTreeHoverType<TNode> | null {
    function getNodeInt(
      node: MarkupNodeType, canvasPoint: DG.Point, currentHeight: number
    ): RectangleTreeHoverType<TNode> | null {
      const dpr: number = window.devicePixelRatio;
      // console.debug('DendrogramTreePlacer.getNode() ' +
      //   `point = ${JSON.stringify(point)}, currentHeight = ${currentHeight}, ` +
      //   `node = ${JSON.stringify({
      //     index: node.index,
      //     minIndex: node.minIndex,
      //     maxIndex: node.maxIndex,
      //     branch_length: node.branch_length
      //   })}`);
      // tree to canvas is linear transform, so searching node in canvas coords is correct
      const minIndex: number = (node.minIndex ?? node.index) - 0.25;
      const maxIndex: number = (node.maxIndex ?? node.index) + 0.25;
      const leftTop: DG.Point = treeToCanvas(new DG.Point(currentHeight, minIndex));
      const rightBottomP: DG.Point = treeToCanvas(new DG.Point(currentHeight + node.branch_length!, maxIndex));
      if (leftTop.y <= canvasPoint.y && canvasPoint.y <= rightBottomP.y) {
        let res: RectangleTreeHoverType<TNode> | null = null;
        const beginP: DG.Point = treeToCanvas(new DG.Point(currentHeight, node.index));
        const endP: DG.Point = treeToCanvas(new DG.Point(currentHeight + node.branch_length!, node.index));
        const nodeR: number = nodeSize * dpr / 2; // in canvas pixels
        const lineR: number = lineWidth * 1.5 * dpr / 2;
        if (
          // node stick
          (Math.abs(canvasPoint.y - beginP.y) < lineR &&
            beginP.x - lineR <= canvasPoint.x && canvasPoint.x <= endP.x + lineR) ||
          // leaf tip
          (isLeaf(node) &&
            (Math.pow((endP.x - canvasPoint.x) / nodeR, 2) + Math.pow((endP.y - canvasPoint.y) / nodeR, 2)) < 1)
        )
          res = {nodeHeight: currentHeight, node: node as TNode};

        if (!res) {
          for (const childNode of (node.children ?? [])) {
            res = getNodeInt(childNode, canvasPoint, currentHeight + node.branch_length!);
            if (res) break;
          }
        }

        return res;
      }

      return null;
    }

    const res = getNodeInt(treeRoot, canvasPoint, 0);
    return res;
  }

  getNodeHeight(treeRoot: MarkupNodeType, node: MarkupNodeType): number | undefined {
    function getNodeHeightInt(
      currentNode: MarkupNodeType, node: MarkupNodeType, currentHeight: number
    ): number | undefined {
      let res: number | undefined = undefined;
      if (currentNode == node) {
        res = currentHeight;
      } else {
        for (const childNode of (currentNode.children ?? [])) {
          res = getNodeHeightInt(childNode, node, currentHeight + currentNode.branch_length!);
          if (res !== undefined) break;
        }
      }

      return res;
    }

    const res: number | undefined = getNodeHeightInt(treeRoot, node, 0);
    return res;
  }
}

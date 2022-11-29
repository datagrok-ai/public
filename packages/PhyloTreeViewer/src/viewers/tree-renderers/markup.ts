import {isLeaf, NodeType} from '@datagrok-libraries/bio';

/** Markup node for node and its subtree place on axis of leaves. */
export type MarkupNodeType = NodeType & {
  children: MarkupNodeType[],
  index: number,
  minIndex: number,
  maxIndex: number,
  /** node's branch_length with max subtreeLength of children */
  subtreeLength?: number,
};

/**
 *
 * @param node
 * @param currentLeafIndex
 * @return {number} Index pointing to the next leaf
 */
export function markupNode(
  node: MarkupNodeType | NodeType, currentLeafIndex: number = 0
): void {
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
export function renderNode<TNode extends MarkupNodeType>(
  ctx: CanvasRenderingContext2D, node: TNode,
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
        renderNode(ctx, childNode,
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
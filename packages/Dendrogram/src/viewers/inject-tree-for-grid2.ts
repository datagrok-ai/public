import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import {TreeHelper} from '../../src/utils/tree-helper';
import {GridTreeRendererBase} from './tree-renderers/grid-tree-renderer-base';
import {LeafRangeGridTreeRenderer} from '../../src/viewers/tree-renderers/grid-tree-renderer';
import {NodeCuttedType, NodeType} from '@datagrok-libraries/bio/src/trees';
import {TreeCutOptions} from '@datagrok-libraries/bio/src/trees/dendrogram';
import {markupNode, MarkupNodeType} from './tree-renderers/markup';
import {attachDivToGrid} from '../utils';
import {
  PROPS as _D_PROPS,
  PROPS_CATS as _D_PROPS_CATS,
} from './dendrogram';
import {RectangleTreeHoverType} from './tree-renderers/rectangle-tree-placer';
import {GridTreePlacer} from './tree-renderers/grid-tree-placer';
import {Unsubscribable} from 'rxjs';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';

import '../css/injected-dendrogram.css';

export function injectTreeForGridUI2(
  grid: DG.Grid, treeRoot: NodeType | null, leafColName?: string, neighborWidth: number = 100, cut?: TreeCutOptions,
): GridNeighbor {
  const th: ITreeHelper = new TreeHelper();
  const treeNb: GridNeighbor = attachDivToGrid(grid, neighborWidth);
  if (treeRoot) markupNode(treeRoot);
  const totalLength: number = treeRoot ? (treeRoot as MarkupNodeType).subtreeLength! : 1;
  if (Number.isNaN(totalLength))
    throw new Error('Can not calculate totalLength for the tree.');

  const placer: GridTreePlacer<MarkupNodeType> = new GridTreePlacer<MarkupNodeType>(grid, totalLength);

  const renderer: GridTreeRendererBase<MarkupNodeType> =
    LeafRangeGridTreeRenderer.create(grid, treeRoot, placer);
  renderer.attach(treeNb.root!);

  renderer.onAfterRender.subscribe(({target, context, lengthRatio}) => {
    if (cutSlider) {
      const tgt = target as GridTreeRendererBase<MarkupNodeType>;

      cutSlider.root.style.left = `${0}px`;
      cutSlider.root.style.width = `${tgt.view!.clientWidth}px`;
      cutSlider.root.style.height = `${grid.colHeaderHeight}px`;

      const posX = cutSlider.value! * lengthRatio + tgt.leftPadding * window.devicePixelRatio;

      context.strokeStyle = '#A00000';
      context.moveTo(posX, 0);
      context.lineTo(posX, context.canvas.height);
      context.stroke();
    }
  });

  let cutSlider: DG.InputBase<number | null> | null = null;
  if (cut) {
    // TODO: Get max from tree height
    //@ts-ignore
    const treeHeight: number = (treeRoot as MarkupNodeType).subtreeLength;
    cutSlider = ui.input.slider('', {value: 0, min: 0, max: treeHeight});
    $(cutSlider.root).find('input').each((_, el) => {
      el.setAttribute('step', '0.01');
      el.style.width = '100%';
      el.style.height = `100%`;
    });
    cutSlider.root.style.position = 'absolute';
    cutSlider.root.style.top = `${0}px`;
    cutSlider.root.style.backgroundColor = '#FFF0F0';

    treeNb.root!.appendChild(cutSlider.root);

    cutSlider.onChanged.subscribe((value) => {
      // console.debug('Dendrogram: injectTreeToGrid() cutSlider.onChanged() ' + `${cutSlider!.value}`);
      // th.cutTreeToGrid(newickRoot, cutSlider!.value!, grid.dataFrame, leafColName, 'Cluster');

      const t1 = Date.now();
      th.treeCutAsTree(treeRoot, value!, true);
      const t2 = Date.now();
      console.debug('Dendrogram: injectTreeForGrid() cutSlider.onChanged() treeCutAsTree() ' +
        `ET: ${((t2 - t1) / 1000).toString()}`);

      const newickRootCopy = JSON.parse(JSON.stringify(treeRoot));
      const newickRootCutted = th.treeCutAsTree(newickRootCopy, value!);
      th.markClusters(newickRootCutted as NodeCuttedType, cut.dataDf, leafColName ?? null, cut.clusterColName);
      th.buildClusters(newickRootCutted as NodeCuttedType, cut.clusterDf, cut.clusterColName, leafColName);

      markupNode(newickRootCutted!);
      renderer.treeRoot = newickRootCutted as MarkupNodeType;
    });
  }

  function alignGridWithTree(): void {
    const [viewedRoot] = th.setGridOrder(treeRoot, grid, leafColName);
    if (viewedRoot) markupNode(viewedRoot);

    renderer.treeRoot = viewedRoot as MarkupNodeType;
  }

  // initial alignment tree with grid
  alignGridWithTree();

  // -- Handling events --

  // When a node is clicked, all of the leaf nodes come from it should become selected
  function rendererOnCurrentChanged() {
    window.setTimeout(() => {
      if (!renderer || !placer) return;

      // Reset the current row index
      // grid.dataFrame.currentRowIdx = -1;
      const selectionIndexes = new Set<number>();
      // only one current becomes first of sub leaves
      const selectionLeafs: NodeType[] | null = renderer.current ? th.getLeafList(renderer.current.node) : null;
      if (selectionLeafs && selectionLeafs.length > 0) {
        const selectionLeafNames = new Set<string>(selectionLeafs.map((leaf) => leaf.name));
        if (leafColName) {
          const leafCol: DG.Column = grid.dataFrame.getCol(leafColName);
          const rowCount = grid.dataFrame.rowCount;
          for (let rowI: number = 0; rowI < rowCount; rowI++) {
            const rowLeafName: string = leafCol.get(rowI);
            if (selectionLeafNames.has(rowLeafName)) {
              selectionIndexes.add(rowI);
              break;
            }
          }
        } else {
          selectionLeafs.forEach((leaf) => {
            selectionIndexes.add(parseInt(leaf.name));
          });
        }
      }
      if (selectionIndexes.size === 1) {
        grid.dataFrame.currentRowIdx = selectionIndexes.values().next().value;
      } else {
        grid.dataFrame.selection.init((i) => selectionIndexes.has(i));
        grid.dataFrame.selection.fireChanged();
      }
      grid.invalidate(); // fixing stall current on changed currentRowIdx to -1
    });
  }

  function rendererOnMouseOverChanged() {
    window.setTimeout(() => {
      if (!renderer || !placer) return;

      const oldMouseOverRowIdx: number = grid.dataFrame.mouseOverRowIdx;
      // only one mouseOver becomes first of sub leaves
      const mouseOverLeaf: NodeType | null = renderer.mouseOver ? th.getLeafList(renderer.mouseOver.node)[0] : null;
      let newMouseOverRowIdx: number = -1;
      if (mouseOverLeaf) {
        if (leafColName) {
          const leafCol: DG.Column = grid.dataFrame.getCol(leafColName);
          const rowCount = grid.dataFrame.rowCount;
          for (let rowI = 0; rowI < rowCount; rowI++) {
            const rowLeafName: string = leafCol.get(rowI);
            if (rowLeafName == mouseOverLeaf.name) {
              newMouseOverRowIdx = rowI;
              break;
            }
          }
        } else {
          newMouseOverRowIdx = parseInt(mouseOverLeaf.name);
        }
      }
      if (newMouseOverRowIdx != oldMouseOverRowIdx)
        grid.dataFrame.mouseOverRowIdx = newMouseOverRowIdx;
    });
  }

  function rendererOnSelectionChanged() {
    window.setTimeout(() => {
      if (!renderer || !placer) return;

      /* Here we get selected rows from dataFrame (leaves only).
       * Some of selected nodes can be in subtree of renderer.selections nodes.
       * We need to merge nodes to form selections object.
       * Nodes subs can be selected or deselected.
       */

      const oldSelection: DG.BitSet = grid.dataFrame.selection.clone();
      if (renderer.selections.length == 0) {
        grid.dataFrame.selection.init((_) => { return false; }, false);
      } else {
        const leafCol: DG.Column | null = !!leafColName ? grid.dataFrame.getCol(leafColName) : null;
        const nodeNameSet = new Set(
          renderer.selectedNodes
            .map((sn) => th.getNodeList(sn).map((n) => n.name))
            .flat());
        // console.debug('Dendrogram: Dendrogram.rendererOnSelectionChanged(), ' +
        //   `nodeNameSet = ${JSON.stringify([...nodeNameSet])}`);

        grid.dataFrame.selection.init(
          (rowI) => {
            const nodeName = !!leafCol ? leafCol.get(rowI) : `${rowI}`;
            return nodeNameSet.has(nodeName);
          },
          false);
      }

      const newSelection: DG.BitSet = grid.dataFrame.selection;
      let selectionChanged: boolean = oldSelection.length !== newSelection.length || // != -> true (changed)
        oldSelection.trueCount !== newSelection.trueCount;
      if (!selectionChanged) {
        for (let rowI: number = 0; rowI < oldSelection.length; rowI++) {
          if (oldSelection.get(rowI) !== newSelection.get(rowI)) {
            selectionChanged = true;
            break;
          }
        }
      }
      if (selectionChanged)
        grid.dataFrame.selection.fireChanged();
    }, 0 /* next event cycle*/);
  }

  function dataFrameOnCurrentRowChanged(_value: any) {
    const leafCol: DG.Column | null = !!leafColName ? grid.dataFrame.getCol(leafColName) : null;
    const idx: number = grid.dataFrame.currentRowIdx;
    const currentLeafName: string | null = idx == -1 ? null : !!leafCol ? leafCol.get(idx) : `${idx}`;

    const th: ITreeHelper = new TreeHelper();
    const currentLeaf: MarkupNodeType | null = !renderer.treeRoot ? null : th.getNodeList(renderer.treeRoot)
      .find((leaf) => currentLeafName == leaf.name) ?? null;
    const current: RectangleTreeHoverType<MarkupNodeType> | null = currentLeaf ? {
      node: currentLeaf,
      nodeHeight: placer!.getNodeHeight(renderer.treeRoot, currentLeaf)!,
    } : null;

    renderer.current = current;
  }

  function dataFrameOnMouseOverRowChanged(_value: any) {
    if (!renderer || !placer) return;

    const leafCol: DG.Column | null = !!leafColName ? grid.dataFrame.getCol(leafColName) : null;
    const idx: number = grid.dataFrame.mouseOverRowIdx;
    const mouseOverLeafName: string | null = idx == -1 ? null : !!leafCol ? leafCol.get(idx) : `${idx}`;

    const th: ITreeHelper = new TreeHelper();
    const mouseOverLeaf: MarkupNodeType | null = th.getLeafList(renderer.treeRoot)
      .find((leaf) => mouseOverLeafName == leaf.name) ?? null;
    const mouseOver: RectangleTreeHoverType<MarkupNodeType> | null = mouseOverLeaf ? {
      node: mouseOverLeaf,
      nodeHeight: placer!.getNodeHeight(renderer.treeRoot, mouseOverLeaf)!,
    } : null;

    renderer.mouseOver = mouseOver;
  }

  function dataFrameOnSelectionChanged(_value: any) {
    if (!renderer || !placer) return;

    const leafList: MarkupNodeType[] = th.getLeafList(renderer.treeRoot);
    const leafDict: { [name: string]: MarkupNodeType } = {};
    for (const node of leafList) {
      if (node.name in leafDict)
        throw new Error('Non unique key tree node name');
      leafDict[node.name] = node;
    }

    const rowDict: { [name: string]: number } = {};
    const leafCol: DG.Column | null = leafColName ? grid.dataFrame.getCol(leafColName) : null;
    const rowCount: number = grid.dataFrame.rowCount;
    for (let rowI = 0; rowI < rowCount; rowI++) {
      const leafName: string = leafCol ? leafCol.get(rowI) : `${rowI}`;
      rowDict[leafName] = rowI;
    }

    const selLeafNames: { [name: string]: number } = {};
    const selIndexes = grid.dataFrame.selection.getSelectedIndexes();
    for (const selRowI of selIndexes) {
      const leafName: string = leafCol ? leafCol.get(selRowI) : `${selRowI}`;
      selLeafNames[leafName] = selRowI;
    }
    const selNodeList: MarkupNodeType[] = th.getNodesByLeaves<MarkupNodeType>(renderer.treeRoot, selLeafNames);

    const selections: RectangleTreeHoverType<MarkupNodeType>[] = [];
    for (const selNode of selNodeList)
      selections.push({node: selNode, nodeHeight: placer.getNodeHeight(renderer.treeRoot, selNode)!});

    renderer.selections = selections;
  }

  // Variable to track if filter is changed and prevent the sorting change event
  let filterChangeCounter = 0;

  function dataFrameOnFilterChanged(_value: any) {
    // TODO: Filter newick tree
    console.debug('Dendrogram: injectTreeForGridUI2() grid.dataFrame.onFilterChanged()');
    filterChangeCounter += 1;
    // to prevent nested fire event in event handler
    window.setTimeout(() => {
      alignGridWithTree();
    }, 0);
  }

  function dfOnSortingChanged(_value?: any) {
    // If the reordering is caused by the filter change, return
    if (filterChangeCounter > 0) {
      filterChangeCounter -= 1;
      return;
    }
    const treeOverlay = ui.div();
    treeOverlay.style.width = treeNb.root!.style.width;
    treeOverlay.style.height = treeNb.root!.style.height;
    treeOverlay.classList.add('dendrogram-overlay');

    const sortInfoDiv = ui.div('Revert columns sort order to see Dendrogram Tree');
    const realignButton = ui.button('Revert sort', () => {
      alignGridWithTree();

      treeNb?.root?.removeChild(treeOverlay);
      sortingSub = grid.onRowsSorted.subscribe(dfOnSortingChanged);
    });

    const infoContainer = ui.divV(
      [sortInfoDiv, realignButton],
    );

    treeOverlay.appendChild(infoContainer);
    treeNb.root?.appendChild(treeOverlay);
    sortingSub.unsubscribe();
  }

  function treeNeighborOnClosed() {
    for (const sub of subs) { sub.unsubscribe(); }
  }

  const subs: Unsubscribable[] = [];
  subs.push(treeNb.onClosed.subscribe(treeNeighborOnClosed));
  subs.push(renderer.onCurrentChanged.subscribe(rendererOnCurrentChanged));
  subs.push(renderer.onMouseOverChanged.subscribe(rendererOnMouseOverChanged));
  subs.push(renderer.onSelectionChanged.subscribe(rendererOnSelectionChanged));

  let sortingSub = grid.onRowsSorted.subscribe(dfOnSortingChanged);
  subs.push(sortingSub);

  subs.push(grid.onRowsResized.subscribe(dataFrameOnFilterChanged));
  subs.push(grid.dataFrame.onCurrentRowChanged.subscribe(dataFrameOnCurrentRowChanged));
  subs.push(grid.dataFrame.onMouseOverRowChanged.subscribe(dataFrameOnMouseOverRowChanged));
  subs.push(grid.dataFrame.onSelectionChanged.subscribe(dataFrameOnSelectionChanged));
  subs.push(grid.dataFrame.onFilterChanged.subscribe(dataFrameOnFilterChanged));

  return treeNb;
}

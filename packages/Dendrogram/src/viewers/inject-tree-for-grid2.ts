import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import {TreeHelper} from '../../src/utils/tree-helper';
import {GridTreeRendererBase} from './tree-renderers/grid-tree-renderer-base';
import {ITreeHelper, NodeCuttedType, NodeType, TreeCutOptions} from '@datagrok-libraries/bio';
import {LeafRangeGridTreeRenderer} from '../../src/viewers/tree-renderers/grid-tree-renderer';
import {markupNode, MarkupNodeType} from './tree-renderers/markup';
import {attachDivToGrid} from '../utils';
import {
  PROPS as D_PROPS,
  PROPS_CATS as D_PROPS_CATS,
  TreeDefaultPalette,
  TreeColorNames,
  TRANS_ALPHA, IDendrogram
} from './dendrogram';
import {DendrogramTreeStyler} from './tree-renderers/dendrogram-tree-styler';
import {setAlpha, toRgba} from '@datagrok-libraries/utils/src/color';
import wu from 'wu';
import {RectangleTreeHoverType} from './tree-renderers/rectangle-tree-placer';
import {GridTreePlacer} from './tree-renderers/grid-tree-placer';
import {Unsubscribable} from 'rxjs';
import {render} from 'datagrok-api/ui';

export function injectTreeForGridUI2(
  grid: DG.Grid, newickRoot: NodeType, leafColName?: string, neighborWidth: number = 100, cut?: TreeCutOptions
): GridNeighbor {
  const th: ITreeHelper = new TreeHelper();

  const treeNb: GridNeighbor = attachDivToGrid(grid, neighborWidth);

  // const treeDiv = ui.div();
  // treeRoot.appendChild(treeDiv);
  // treeRoot.style.backgroundColor = '#FFF0F0';
  // treeRoot.style.setProperty('overflow-y', 'hidden', 'important');

  // TODO: adapt tree: bio.NodeType to MarkupNodeType
  markupNode(newickRoot);
  const totalLength: number = (newickRoot as MarkupNodeType).subtreeLength!;
  if (Number.isNaN(totalLength))
    throw new Error('Can not calculate totalLength for the tree.');

  const placer: GridTreePlacer<MarkupNodeType> = new GridTreePlacer<MarkupNodeType>(grid, totalLength);

  const renderer: GridTreeRendererBase<MarkupNodeType> =
    LeafRangeGridTreeRenderer.create(grid, newickRoot, placer);
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
    const treeHeight: number = (newickRoot as MarkupNodeType).subtreeLength;
    cutSlider = ui.sliderInput('', 0, 0, treeHeight);
    $(cutSlider.root).find('input').each((_, el) => {
      el.setAttribute('step', '0.01');
      el.style.width = '100%';
      el.style.height = `100%`;
    });
    cutSlider.root.style.position = 'absolute';
    cutSlider.root.style.top = `${0}px`;
    cutSlider.root.style.backgroundColor = '#FFF0F0';

    treeNb.root!.appendChild(cutSlider.root);

    cutSlider.onChanged(() => {
      // console.debug('Dendrogram: injectTreeToGrid() cutSlider.onChanged() ' + `${cutSlider!.value}`);
      // th.cutTreeToGrid(newickRoot, cutSlider!.value!, grid.dataFrame, leafColName, 'Cluster');

      const t1 = Date.now();
      th.treeCutAsTree(newickRoot, cutSlider!.value!, true);
      const t2 = Date.now();
      console.debug('Dendrogram: injectTreeForGrid() cutSlider.onChanged() treeCutAsTree() ' +
        `ET: ${((t2 - t1) / 1000).toString()}`);

      const newickRootCopy = JSON.parse(JSON.stringify(newickRoot));
      const newickRootCutted = th.treeCutAsTree(newickRootCopy, cutSlider!.value!);
      th.markClusters(newickRootCutted as NodeCuttedType, cut.dataDf, leafColName ?? null, cut.clusterColName);
      th.buildClusters(newickRootCutted as NodeCuttedType, cut.clusterDf, cut.clusterColName, leafColName);

      markupNode(newickRootCutted!);
      renderer.treeRoot = newickRootCutted as MarkupNodeType;
    });
  }

  function alignGridWithTree(): void {
    const [viewedRoot] = th.setGridOrder(newickRoot, grid, leafColName);
    markupNode(viewedRoot);
    const source = viewedRoot ? {type: 'biojs', data: viewedRoot} :
      {type: 'biojs', data: {name: 'NONE', branch_length: 1, children: []}};

    renderer.treeRoot = viewedRoot as MarkupNodeType;
  }

  // initial alignment tree with grid
  alignGridWithTree();

  // -- Inject properties --

  try {
    const lineWidthProperty = DG.Property.int(D_PROPS.lineWidth,
      (obj) => {
        let k = 11;
      },
      (obj, value) => {
        let k = 11;
      },
      1);
    lineWidthProperty.category = `Dendrogram ${D_PROPS_CATS.APPEARANCE}`;
    DG.Property.registerAttachedProperty('GridLook', lineWidthProperty);
  } catch (err: any) {
    console.warn(err);
  }

  // -- Handling events --

  function rendererOnCurrentChanged() {
    window.setTimeout(() => {
      if (!renderer || !placer) return;

      const oldCurrentRowIdx: number = grid.dataFrame.currentRowIdx;
      // only one current becomes first of sub leaves
      const currentLeaf: NodeType | null = renderer.current ? th.getLeafList(renderer.current.node)[0] : null;
      let newCurrentRowIdx: number = -1;
      if (currentLeaf) {
        if (leafColName) {
          const leafCol: DG.Column = grid.dataFrame.getCol(leafColName);
          const rowCount = grid.dataFrame.rowCount;
          for (let rowI: number = 0; rowI < rowCount; rowI++) {
            const rowLeafName: string = leafCol.get(rowI);
            if (rowLeafName == currentLeaf.name) {
              newCurrentRowIdx = rowI;
              break;
            }
          }
        } else {
          newCurrentRowIdx = parseInt(currentLeaf.name);
        }
      }
      if (newCurrentRowIdx != oldCurrentRowIdx) {
        grid.dataFrame.currentRowIdx = newCurrentRowIdx;
        grid.invalidate(); // fixing stall current on changed currentRowIdx to -1
      }
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
      if (newMouseOverRowIdx != oldMouseOverRowIdx) {
        grid.dataFrame.mouseOverRowIdx = newMouseOverRowIdx;
      }
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
        grid.dataFrame.selection.init((rowI) => { return false; }, false);
      } else {
        const leafCol: DG.Column | null = !!leafColName ? grid.dataFrame.getCol(leafColName) : null;
        const nodeNameSet = new Set(
          renderer.selectedNodes
            .map((sn) => th.getNodeList(sn).map((n) => n.name))
            .flat());
        console.debug('Dendrogram: Dendrogram.rendererOnSelectionChanged(), ' +
          `nodeNameSet = ${JSON.stringify([...nodeNameSet])}`);

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

  function dataFrameOnCurrentRowChanged(value: any) {
    const leafCol: DG.Column | null = !!leafColName ? grid.dataFrame.getCol(leafColName) : null;
    const idx: number = grid.dataFrame.currentRowIdx;
    const currentLeafName: string | null = idx == -1 ? null : !!leafCol ? leafCol.get(idx) : `${idx}`;

    const th: ITreeHelper = new TreeHelper();
    const currentLeaf: MarkupNodeType | null = th.getNodeList(renderer.treeRoot)
      .find((leaf) => currentLeafName == leaf.name) ?? null;
    const current: RectangleTreeHoverType<MarkupNodeType> | null = currentLeaf ? {
      node: currentLeaf,
      nodeHeight: placer!.getNodeHeight(renderer.treeRoot, currentLeaf)!
    } : null;

    renderer.current = current;
  }

  function dataFrameOnMouseOverRowChanged(value: any) {
    if (!renderer || !placer) return;

    const leafCol: DG.Column | null = !!leafColName ? grid.dataFrame.getCol(leafColName) : null;
    const idx: number = grid.dataFrame.mouseOverRowIdx;
    const mouseOverLeafName: string | null = idx == -1 ? null : !!leafCol ? leafCol.get(idx) : `${idx}`;

    const th: ITreeHelper = new TreeHelper();
    const mouseOverLeaf: MarkupNodeType | null = th.getLeafList(renderer.treeRoot)
      .find((leaf) => mouseOverLeafName == leaf.name) ?? null;
    const mouseOver: RectangleTreeHoverType<MarkupNodeType> | null = mouseOverLeaf ? {
      node: mouseOverLeaf,
      nodeHeight: placer!.getNodeHeight(renderer.treeRoot, mouseOverLeaf)!
    } : null;

    renderer.mouseOver = mouseOver;
  }

  function dataFrameOnSelectionChanged(value: any) {
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

  function dataFrameOnFilterChanged(value: any) {
    // TODO: Filter newick tree
    console.debug('Dendrogram: injectTreeForGridUI2() grid.dataFrame.onFilterChanged()');

    // to prevent nested fire event in event handler
    window.setTimeout(() => { alignGridWithTree(); }, 0);
  }

  const subs: Unsubscribable[] = [];
  subs.push(renderer.onCurrentChanged.subscribe(rendererOnCurrentChanged));
  subs.push(renderer.onMouseOverChanged.subscribe(rendererOnMouseOverChanged));
  subs.push(renderer.onSelectionChanged.subscribe(rendererOnSelectionChanged));

  subs.push(grid.dataFrame.onCurrentRowChanged.subscribe(dataFrameOnCurrentRowChanged));
  subs.push(grid.dataFrame.onMouseOverRowChanged.subscribe(dataFrameOnMouseOverRowChanged));
  subs.push(grid.dataFrame.onSelectionChanged.subscribe(dataFrameOnSelectionChanged));

  subs.push(grid.dataFrame.onFilterChanged.subscribe(dataFrameOnFilterChanged));

  return treeNb;
}


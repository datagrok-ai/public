import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {filter, map} from 'rxjs/operators';
import {interval, Unsubscribable} from 'rxjs';

import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import $ from 'cash-dom';
import {TreeHelper} from '../utils/tree-helper';

// const getBranchScaleOld = PhylocanvasGL.prototype.getBranchScale;
// PhylocanvasGL.prototype.getBranchScale = function(...args) {
//   return getBranchScaleOld(...args);
// };

/** Modifies
 * @param {string} leafColName Column name of grid.dataFrame to use as leaf name/key
 *                                 undefined - use row index as leaf name/key
 */
export function injectTreeToGridUI(
  grid: DG.Grid, newickText: string, leafColName?: string, neighborWidth: number = 100,
  cut?: { min: number, max: number, clusterColName: string }
): GridNeighbor {
  const th: bio.ITreeHelper = new TreeHelper();
  const subs: Unsubscribable[] = [];
  //const _grid = grid;
  const treeN = attachDivToGrid(grid, neighborWidth);
  const treeRoot = treeN.root!;
  const pcDiv = ui.div();

  const leftMargin: number = 0;

  treeRoot.appendChild(pcDiv);
  treeRoot.style.backgroundColor = '#FFF0F0';
  treeRoot.style.setProperty('overflow-y', 'hidden', 'important');

  pcDiv.style.backgroundColor = '#F0FFF0';
  pcDiv.style.position = 'absolute';
  pcDiv.style.left = `${leftMargin}px`;

  const newickRoot: bio.NodeType = bio.Newick.parse_newick(newickText);
  const [viewedRoot, warnings]: [bio.NodeType, string[]] = th.setGridOrder(newickRoot, grid, 'id');

  const pcViewer = new bio.PhylocanvasGL(pcDiv, {
    interactive: false,
    alignLabels: true,
    padding: 1, // required for top most joint
    nodeSize: 2,
    treeToCanvasRatio: 0.99, // 1 - hide leaf labels, 0.99 to fully viewed top most joint

    showLabels: true,
    showLeafLabels: true,
    size: {width: 400, height: 300},

    source: {type: 'biojs', data: viewedRoot},
  });
  pcViewer.view.style.backgroundImage = 'none';
  pcViewer.deck.setProps({useDevicePixels: true}); // fix blurred

  // const syncer = new TreeToGridSyncer(nDiv, newickRoot, pcViewer, grid, leafColName, true);

  let cutSlider: DG.InputBase<number | null> | null = null;
  if (cut) {
    // TODO: Get max from tree height
    //@ts-ignore
    const treeHeight: number = pcViewer.props.source.data.totalSubtreeLength;
    cutSlider = ui.sliderInput('', 0, 0, treeHeight);
    $(cutSlider.root).find('input').each((_, el) => {
      el.setAttribute('step', '0.01');
      el.style.width = '100%';
      el.style.height = `100%`;
    });
    // cutSlider.root.setAttribute('step', '0.01');
    cutSlider.root.style.position = 'absolute';
    cutSlider.root.style.backgroundColor = '#FFF0F0';

    treeRoot.appendChild(cutSlider.root);

    subs.push(cutSlider.onChanged(() => {
      console.debug('PhyloTreeViewer: injectTreeToGrid() cutSlider.onChanged() ' + `${cutSlider!.value}`);
      th.cutTreeToGrid(newickRoot, cutSlider!.value!, grid.dataFrame, 'id', 'Cluster');
    }));
  }

  // TODO: inject properties to grid Viewer

  pcCalcSize();

  function pcCalcSize() {
    const leafCount = grid.dataFrame.filter.trueCount;

    const width: number = treeRoot.clientWidth - leftMargin;
    let height: number = grid.props.rowHeight * leafCount;

    const firstVisibleRowIdx: number = Math.floor(grid.vertScroll.min);
    //const firstVisibleRowIdx: number = grid.vertScroll.min;

    console.debug('PhyloTreeViewer: injectTreeGridUI pcCalcSize() ' + `{ size: ${width} x ${height} }`);
    if (height == 0)
      height = 1;
    try {
      pcViewer.setProps({size: {width, height}});
    } catch (ex) {}
    pcDiv.style.top = `${grid.colHeaderHeight - firstVisibleRowIdx * grid.props.rowHeight}px`;

    if (cutSlider) {
      const cutDiv = cutSlider.root;
      cutDiv.style.top = `${0}px`;
      cutDiv.style.height = `${grid.colHeaderHeight}px`;
      cutDiv.style.width = `${width}px`;
    }
  }

  let colHeaderHeight = grid.colHeaderHeight;
  interval(200).pipe(
    map((_) => grid.colHeaderHeight),
    filter((hh) => hh != colHeaderHeight)
  ).subscribe(() => {
    colHeaderHeight = grid.colHeaderHeight;
    console.debug('PhyloTreeViewer: injectTreeGridUI() grid.onColHeaderHeightChanged() ' +
      `colHeaderHeight=${grid.colHeaderHeight}`);

    pcCalcSize();
  });

  let rowHeight = grid.props.rowHeight;
  interval(200).pipe(
    map((_) => grid.props.rowHeight),
    filter((rh) => rh != rowHeight)
  ).subscribe(() => {
    rowHeight = grid.props.rowHeight;
    console.debug('PhyloTreeViewer: injectTreeGridUI() grid.onRowHeightChanged() ' +
      `rowHeight=${grid.props.rowHeight}`);

    pcCalcSize();
  });


  grid.onRowsResized.subscribe(() => {
    pcCalcSize();
  });

  grid.vertScroll.onValuesChanged.subscribe(() => {
    pcCalcSize();
  });

  ui.onSizeChanged(treeRoot).subscribe(() => {
    pcCalcSize();
  });

  // grid.onBeforeDrawContent.subscribe(() => {
  //   pcCalcSize();
  // });

  grid.dataFrame.onFilterChanged.subscribe((args) => {
    // TODO: Filter newick tree
    console.debug('PhyloTreeViewer: injectTreeGridUI() grid.dataFrame.onFilterChanged()');

    // to prevent nested fire event in event handler
    window.setTimeout(() => {
      const [viewedRoot] = th.setGridOrder(newickRoot, grid, 'id');
      const source = viewedRoot ? {type: 'biojs', data: viewedRoot} :
        {type: 'biojs', data: {name: 'NONE', branch_length: 1, children: []}};

      try {
        pcViewer.setProps({source: source});
      } catch (ex) { } // ignore exception on empty source

      pcCalcSize();
    }, 0 /* next event cycle */);
  });

  return treeN;
}

/** By Dimitri Petrov */
export function attachDivToGrid(grid: DG.Grid, nWidth: number = 100): GridNeighbor {
  // const nRowCount = 100;
  // const nColCount = 5;
  // const dframe: DG.DataFrame = grok.data.demo.randomWalk(nRowCount, nColCount);
  // const colDate = dframe.columns.addNewDateTime('New Date');
  // for (let nR = 0; nR < nRowCount; ++nR) {
  //   colDate.set(nR, dayjs());
  // }

  let neighbor: GridNeighbor;
  const eDiv = ui.div();
  // let button = ui.button('Close', () => {
  //   neighbor.close();
  // });
  // eDiv.appendChild(button);
  neighbor = new GridNeighbor(eDiv, grid, nWidth);
  return neighbor;
}
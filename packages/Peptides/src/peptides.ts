import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {createPeptideSimilaritySpaceViewer} from './utils/peptide-similarity-space';
import {PeptidesModel} from './model';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {SARViewer, SARViewerVertical} from './viewers/sar-viewer';
import {SubstViewer} from './viewers/subst-viewer';

type viewerTypes = SARViewer | SARViewerVertical | SubstViewer;
export class Peptides {
  private helpUrl = '/help/domains/bio/peptides.md';
  /**
   * Class initializer
   *
   * @param {DG.Grid} tableGrid Working talbe grid.
   * @param {DG.TableView} view Working view.
   * @param {DG.DataFrame} currentDf Working table.
   * @param {StringDictionary} options SAR viewer options
   * @param {DG.Column} col Aligned sequences column.
   * @memberof Peptides
   */
  async init(
    tableGrid: DG.Grid, view: DG.TableView, currentDf: DG.DataFrame, options: StringDictionary, col: DG.Column,
    originalDfColumns: string[]) {
    function adjustCellSize(grid: DG.Grid) {
      const colNum = grid.columns.length;
      for (let i = 0; i < colNum; ++i) {
        const iCol = grid.columns.byIndex(i)!;
        iCol.width = isNaN(parseInt(iCol.name)) ? 50 : 40;
      }
      grid.props.rowHeight = 20;
    }

    for (let i = 0; i < tableGrid.columns.length; i++) {
      const aarCol = tableGrid.columns.byIndex(i);
      if (aarCol &&
          aarCol.name &&
          aarCol.column?.semType != 'aminoAcids'
      ) {
        //@ts-ignore
        tableGrid.columns.byIndex(i)?.visible = false;
      }
    }

    const originalDfName = currentDf.name;
    const dockManager = view.dockManager;

    PeptidesModel.getOrInit(currentDf);

    const sarViewer = await currentDf.plot.fromType('peptide-sar-viewer', options) as SARViewer;
    sarViewer.helpUrl = this.helpUrl;

    const sarViewerVertical = await currentDf.plot.fromType('peptide-sar-viewer-vertical') as SARViewerVertical;
    sarViewerVertical.helpUrl = this.helpUrl;

    const sarViewersGroup: viewerTypes[] = [sarViewer, sarViewerVertical];

    const peptideSpaceViewer = await createPeptideSimilaritySpaceViewer(
      currentDf, col, 't-SNE', 'Levenshtein', 100, view, `${options['activityColumnName']}Scaled`);
    dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.RIGHT, null, 'Peptide Space viewer');

    let nodeList = dockViewers(sarViewersGroup, DG.DOCK_TYPE.RIGHT, dockManager, DG.DOCK_TYPE.DOWN);

    const substViewer = await currentDf.plot.fromType(
      'substitution-analysis-viewer', {'activityColumnName': `${options['activityColumnName']}Scaled`}) as SubstViewer;
    const substViewersGroup = [substViewer];

    tableGrid.props.allowEdit = false;
    adjustCellSize(tableGrid);

    const hideIcon = ui.iconFA('window-close', () => {
      const viewers = [];
      for (const viewer of view.viewers) {
        if (viewer.type !== DG.VIEWER.GRID)
          viewers.push(viewer);
      }
      viewers.forEach((v) => v.close());

      const cols = (currentDf.columns as DG.ColumnList);
      for (const colName of cols.names()) {
        if (!originalDfColumns.includes(colName))
          cols.remove(colName);
      }

      currentDf.selection.setAll(false);
      currentDf.filter.setAll(true);

      tableGrid.setOptions({'colHeaderHeight': 20});
      tableGrid.columns.setVisible(originalDfColumns);
      tableGrid.props.allowEdit = true;
      currentDf.name = originalDfName;

      view.setRibbonPanels(ribbonPanels);
    }, 'Close viewers and restore dataframe');

    let isSA = false;
    const switchViewers = ui.iconFA('toggle-on', () => {
      $(switchViewers).toggleClass('fa-toggle-off').toggleClass('fa-toggle-on');
      nodeList?.forEach((node) => {
        view.dockManager.close(node);
        node.container.destroy();
      });
      const getCurrentViewerGroup = () => isSA ? substViewersGroup : sarViewersGroup;
      getCurrentViewerGroup().forEach((v) => v.removeFromView());
      isSA = !isSA;
      nodeList = dockViewers(getCurrentViewerGroup(), DG.DOCK_TYPE.LEFT, dockManager, DG.DOCK_TYPE.DOWN);
    }, 'Toggle viewer group');

    const ribbonPanels = view.getRibbonPanels();
    view.setRibbonPanels([[hideIcon, switchViewers]]);
  }
}

function dockViewers(
  viewerList: viewerTypes[], attachDirection: DG.DockType, dockManager: DG.DockManager,
  initialAttachDirection?: DG.DockType): DG.DockNode[] | null {
  const viewerListLength = viewerList.length;
  if (viewerListLength === 0)
    return null;

  let currentViewer = viewerList[0];
  const nodeList = [dockManager.dock(currentViewer, initialAttachDirection, null, currentViewer.name ?? '')];
  const ratio = 1 / viewerListLength;

  for (let i = 1; i < viewerListLength; i++) {
    currentViewer = viewerList[i];
    nodeList.push(dockManager.dock(currentViewer, attachDirection, nodeList[i - 1], currentViewer.name ?? '', ratio));
  }
  return nodeList;
}

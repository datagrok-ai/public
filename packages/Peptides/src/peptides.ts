import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {createPeptideSimilaritySpaceViewer} from './utils/peptide-similarity-space';
import {PeptidesModel} from './model';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import { SARViewer, SARViewerVertical } from './viewers/sar-viewer';
import { SubstViewer } from './viewers/subst-viewer';
// import $ from 'cash-dom';

/**
 * Peptides controller class.
 *
 * @export
 * @class Peptides
 */
export class Peptides {
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
    tableGrid: DG.Grid,
    view: DG.TableView,
    currentDf: DG.DataFrame,
    options: StringDictionary,
    col: DG.Column,
    originalDfColumns: string[],
  ) {

    function adjustCellSize(grid: DG.Grid) {
      for (let i = 0; i < grid.columns.length; ++i) {
        const col = grid.columns.byIndex(i)!;
        col.width = isNaN(parseInt(col.name)) ? 50 : 40;
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

    const initialFiter = currentDf.filter.clone();
    // const originalDfColumns = (currentDf.columns as DG.ColumnList).names();
    const originalDfName = currentDf.name;

    PeptidesModel.getOrInit(currentDf);

    // const substViewer = view.addViewer(
    //   'substitution-analysis-viewer', {'activityColumnName': `${options['activityColumnName']}Scaled`},
    // );
    // const substNode = view.dockManager.dock(substViewer, DG.DOCK_TYPE.RIGHT, null, 'Substitution Analysis');

    // const layout1 = view.saveLayout();
    // view.dockManager.close(substNode);

    const dockSAR = function (
        sarViewer: SARViewer, sarViewerVertical: SARViewerVertical, peptideSpaceViewer: DG.ScatterPlotViewer) {
      const sarNode = view.dockManager.dock(sarViewer, DG.DOCK_TYPE.DOWN, null, 'SAR Viewer');
      const sarVNode = view.dockManager.dock(sarViewerVertical, DG.DOCK_TYPE.RIGHT, sarNode, 'SAR Vertical Viewer');
      const psNode = view.dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.LEFT, sarNode, 'Peptide Space Viewer', 0.3);

      return [sarNode, sarVNode, psNode];
    };

    const dockSubst = function (substViewer: SubstViewer) {
      const substNode = view.dockManager.dock(substViewer, DG.DOCK_TYPE.RIGHT, null, 'Substitution Analysis');

      return [substNode];
    }

    const helpUrl = '/help/domains/bio/peptides.md';

    // currentDf.temp['viewerGroup'] = 'SAR';
    // const sarViewer = view.addViewer('peptide-sar-viewer', options);
    const sarViewer = await currentDf.plot.fromType('peptide-sar-viewer', options) as SARViewer;
    sarViewer.helpUrl = helpUrl;
    // const sarNode = view.dockManager.dock(sarViewer, DG.DOCK_TYPE.DOWN, null, 'SAR Viewer');

    // const sarViewerVertical = view.addViewer('peptide-sar-viewer-vertical');
    const sarViewerVertical = await currentDf.plot.fromType('peptide-sar-viewer-vertical') as SARViewerVertical;
    sarViewerVertical.helpUrl = helpUrl;
    // const sarVNode = view.dockManager.dock(sarViewerVertical, DG.DOCK_TYPE.RIGHT, sarNode, 'SAR Vertical Viewer');

    const peptideSpaceViewer = await createPeptideSimilaritySpaceViewer(
      currentDf, col, 't-SNE', 'Levenshtein', 100, view, `${options['activityColumnName']}Scaled`);
    // const psNode = view.dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.LEFT, sarNode, 'Peptide Space Viewer', 0.3);

    let viewerGroup1 = [sarViewer, sarViewerVertical, peptideSpaceViewer];

    // const layout1 = view.saveLayout();
    // const nodeList = [sarNode, sarVNode];
    let nodeList = dockSAR(sarViewer, sarViewerVertical, peptideSpaceViewer);

    const substViewer = await currentDf.plot.fromType(
      'substitution-analysis-viewer', {'activityColumnName': `${options['activityColumnName']}Scaled`}) as SubstViewer;
    // const substNode = view.dockManager.dock(substViewer, DG.DOCK_TYPE.RIGHT, null, 'Substitution Analysis');
    let viewerGroup2 = [substViewer];

    // const layout2 = view.saveLayout();
    // view.dockManager.close(substNode);

    // const viewerList = [sarViewer, sarViewerVertical, peptideSpaceViewer];
    tableGrid.props.allowEdit = false;
    adjustCellSize(tableGrid);

    const hideIcon = ui.iconFA('window-close', () => { //undo?, times?
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
    //TODO: fix layouts
    const switchViewers = ui.iconFA('toggle-on', () => {
      $(switchViewers).toggleClass('fa-toggle-off').toggleClass('fa-toggle-on');
      nodeList.forEach(node => {
        view.dockManager.close(node);
        node.container.destroy();
      });
      (isSA ? viewerGroup2 : viewerGroup1).forEach((v) => v.removeFromView());
      nodeList = isSA ? dockSAR(sarViewer, sarViewerVertical, peptideSpaceViewer) : dockSubst(substViewer);
      isSA = !isSA;
    }, 'Toggle viewer group');

    const ribbonPanels = view.getRibbonPanels();
    view.setRibbonPanels([[hideIcon, switchViewers]]);
  }
}

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {createPeptideSimilaritySpaceViewer} from './utils/peptide-similarity-space';
import {addViewerToHeader} from './viewers/stacked-barchart-viewer';
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
   * @param {{[key: string]: string}} options SAR viewer options
   * @param {DG.Column} col Aligned sequences column.
   * @param {string} activityColumnChoice Activity column name.
   * @memberof Peptides
   */
  async init(
    tableGrid: DG.Grid,
    view: DG.TableView,
    currentDf: DG.DataFrame,
    options: {[key: string]: string},
    col: DG.Column,
    activityColumnChoice: string,
  ) {
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

    const originalDfColumns = (currentDf.columns as DG.ColumnList).names();
    const originalDfName = currentDf.name;

    const substViewer = view.addViewer(
      'substitution-analysis-viewer', {'activityColumnName': options['activityColumnName']},
    );
    view.dockManager.dock(substViewer, DG.DOCK_TYPE.RIGHT, null, 'Substitution Analysis');

    const sarViewer = view.addViewer('peptide-sar-viewer', options);
    const sarNode = view.dockManager.dock(sarViewer, DG.DOCK_TYPE.DOWN, null, 'SAR Viewer');

    const sarViewerVertical = view.addViewer('peptide-sar-viewer-vertical');
    const sarVNode = view.dockManager.dock(sarViewerVertical, DG.DOCK_TYPE.RIGHT, sarNode, 'SAR Vertical Viewer');

    const peptideSpaceViewer = await createPeptideSimilaritySpaceViewer(
      currentDf,
      col,
      't-SNE',
      'Levenshtein',
      100,
      view,
      `${activityColumnChoice}Scaled`,
    );
    const psNode = view.dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.LEFT, sarNode, 'Peptide Space Viewer', 0.3);
    // const sarDockNodes = [sarNode, sarVNode, psNode];
    // const sarViewers = [sarViewer, sarViewerVertical, peptideSpaceViewer];

    const StackedBarchartProm = currentDf.plot.fromType('StackedBarChartAA');
    addViewerToHeader(tableGrid, StackedBarchartProm);
    tableGrid.props.allowEdit = false;

    const hideIcon = ui.iconFA('window-close', () => { //undo?, times?
      const viewers = [];
      for (const viewer of view.viewers) {
        if (viewer.type !== DG.VIEWER.GRID) {
          viewers.push(viewer);
        }
      }
      viewers.forEach((v) => v.close());

      const cols = (currentDf.columns as DG.ColumnList);
      for (const colName of cols.names()) {
        if (!originalDfColumns.includes(colName)) {
          cols.remove(colName);
        }
      }

      currentDf.selection.setAll(false);
      currentDf.filter.setAll(true);

      tableGrid.setOptions({'colHeaderHeight': 20});
      tableGrid.columns.setVisible(originalDfColumns);
      tableGrid.props.allowEdit = true;
      currentDf.name = originalDfName;

      view.setRibbonPanels(ribbonPanels);
    }, 'Close viewers and restore dataframe');

    // let substViewer: DG.Viewer | null = null;
    // let substNode: DG.DockNode | null = null;
    // let isSA = false;
    // let viewLayout1: DG.ViewLayout | null = null;
    // let viewLayout2: DG.ViewLayout | null = null;
    // const switchViewers = ui.iconFA('toggle-on', () => {
    //   if (isSA) {
    //     viewLayout2 = view.saveLayout();
    //     // view.dockManager.close(substNode!);
    //     substViewer?.close();

    //     view.loadLayout(viewLayout1!);

    //     $(switchViewers).removeClass('fa-toggle-off');
    //     $(switchViewers).addClass('fa-toggle-on');
    //   } else {
    //     viewLayout1 = view.saveLayout();
    //     // sarDockNodes.forEach((node) => view.dockManager.close(node));
    //     sarViewers.forEach((v) => v.close());

    //     if (viewLayout2 === null) {
    //       substViewer = view.addViewer(
    //         'substitution-analysis-viewer', {'activityColumnName': options['activityColumnName']},
    //       );
    //       substNode = view.dockManager.dock(substViewer, DG.DOCK_TYPE.DOWN, null, 'Substitution Analysis');
    //     } else {
    //       view.loadLayout(viewLayout2);
    //     }

    //     $(switchViewers).removeClass('fa-toggle-on');
    //     $(switchViewers).addClass('fa-toggle-off');
    //   }
    //   isSA = !isSA;
    // });

    const ribbonPanels = view.getRibbonPanels();
    view.setRibbonPanels([[hideIcon]]);
  }
}

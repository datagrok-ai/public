import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {createPeptideSimilaritySpaceViewer} from './utils/peptide-similarity-space';
import {addViewerToHeader} from './viewers/stacked-barchart-viewer';

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

    const sarViewer = view.addViewer('peptide-sar-viewer', options);
    const sarNode = view.dockManager.dock(sarViewer, DG.DOCK_TYPE.DOWN, null, 'SAR Viewer');

    const sarViewerVertical = view.addViewer('peptide-sar-viewer-vertical');
    view.dockManager.dock(sarViewerVertical, DG.DOCK_TYPE.RIGHT, sarNode, 'SAR Vertical Viewer');

    const peptideSpaceViewer = await createPeptideSimilaritySpaceViewer(
      currentDf,
      col,
      't-SNE',
      'Levenshtein',
      100,
      view,
      `${activityColumnChoice}Scaled`,
    );
    view.dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.LEFT, sarNode, 'Peptide Space Viewer', 0.3);

    const StackedBarchartProm = currentDf.plot.fromType('StackedBarChartAA');
    addViewerToHeader(tableGrid, StackedBarchartProm);

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

      view.setRibbonPanels(ribbonPanels);
    }, 'Close viewers and restore dataframe');

    const ribbonPanels = view.getRibbonPanels();
    view.setRibbonPanels([[hideIcon]]);
  }
}

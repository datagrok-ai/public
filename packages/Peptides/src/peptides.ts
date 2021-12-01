import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {createPeptideSimilaritySpaceViewer} from './utils/peptide-similarity-space';
import {addViewerToHeader} from './viewers/stacked-barchart-viewer';

export class Peptides {
  async init(
      tableGrid: DG.Grid,
      view: DG.TableView,
      currentDf: DG.DataFrame,
      options: {[key: string]: string},
      col: DG.Column,
      activityColumnChoice: string
    ) {
    for (let i = 0; i < tableGrid.columns.length; i++) {
      const col = tableGrid.columns.byIndex(i);
      if (col &&
          col.name &&
          col.column?.semType != 'aminoAcids'
      ) {
        //@ts-ignore
        tableGrid.columns.byIndex(i)?.visible = false;
      }
    }

    const originalDfColumns = (currentDf.columns as DG.ColumnList).names();

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

        const cols = (view.dataFrame.columns as DG.ColumnList);
        for (const colName of cols.names()) {
          if (!originalDfColumns.includes(colName)) {
            cols.remove(colName);
          }
        }

        tableGrid.setOptions({'colHeaderHeight': 20});
        tableGrid.columns.setVisible(originalDfColumns);

        view.setRibbonPanels(ribbonPanels);
    }, 'Close viewers and restore dataframe');
    
    const ribbonPanels = view.getRibbonPanels();
    view.setRibbonPanels([[hideIcon]]);
  }
}
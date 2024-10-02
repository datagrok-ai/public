import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {hierarchicalClusteringUI} from '../utils/hierarchical-clustering';
import {DistanceMetric} from '@datagrok-libraries/bio/src/trees';

export async function heatmapDemo() {
  try {
    const csv: string = await _package.files.readAsText('data/demog-short.csv');
    const df = DG.DataFrame.fromCsv(csv);
    const tv = grok.shell.addTableView(df, DG.DOCK_TYPE.FILL);
    await hierarchicalClusteringUI(df, ['AGE'], DistanceMetric.Euclidean, 'ward', 300, {tableView: tv});

    const adjustTreeHeight = () => {
      const rowCount = tv!.dataFrame.filter.trueCount;
      const rowsGridHeight: number = tv!.grid.root.clientHeight - tv!.grid.colHeaderHeight;
        tv!.grid.props.rowHeight = rowsGridHeight / (rowCount + 1);
    };


    const rootNode = tv.dockManager.rootNode;
    const extraViewerNode = tv.dockManager.dock(
      tv.addViewer(DG.VIEWER.SCATTER_PLOT), DG.DOCK_TYPE.DOWN, rootNode, DG.VIEWER.SCATTER_PLOT, 0.25,
    );
    tv.dockManager.dock(
      tv.addViewer(DG.VIEWER.HISTOGRAM), DG.DOCK_TYPE.RIGHT, extraViewerNode, DG.VIEWER.HISTOGRAM, 0.5,
    );
    // const viewerNode = tv.dockManager.dock(tv.grid, DG.DOCK_TYPE.TOP, null, DG.VIEWER.HEAT_MAP, 0.7);
    tv.dockManager.dock(tv.filters(), DG.DOCK_TYPE.LEFT, rootNode, DG.VIEWER.FILTERS, 0.2);
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp('/help/visualize/viewers/heat-map.md');
    tv!.grid.onBeforeDrawContent.subscribe(() => {
      adjustTreeHeight();
    });
    adjustTreeHeight();
    // tv.grid.props.showRowHeader = false;
    tv!.grid.props.isGrid = false;
    tv!.grid.props.isHeatmap = true;
    tv!.grid.props.showRowHeader = false;
    // this.tv!.grid.props.showAddNewRowIcon = false;
    setTimeout(() => {
      adjustTreeHeight();
        tv!.grid.invalidate();
    }, 100);
  } catch (err) {
    console.error(err);
  }
}

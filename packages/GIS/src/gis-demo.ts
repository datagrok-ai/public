import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


const MAP_PATH = 'geo/earthquakes.csv';
const MAP_NAME = 'Map';

export async function gisDemo(options?: object | null): Promise<void> {
  const df = await grok.data.getDemoTable(MAP_PATH);
  const tableView = grok.shell.addTableView(df);

  grok.shell.windows.showHelp = true;
  grok.shell.windows.help.syncCurrentObject = false;

  DG.debounce(df.onSemanticTypeDetected, 800).subscribe((_) => {
    const viewer = tableView.addViewer(MAP_NAME, options);
    grok.shell.windows.help.showHelp(viewer.helpUrl);

    const rootNode = tableView.dockManager.rootNode;
    const scatterplotNode = tableView.dockManager.dock(tableView.addViewer(DG.VIEWER.SCATTER_PLOT),
      DG.DOCK_TYPE.RIGHT, rootNode, DG.VIEWER.SCATTER_PLOT, 0.5);
    tableView.dockManager.dock(tableView.addViewer(DG.VIEWER.HISTOGRAM), DG.DOCK_TYPE.RIGHT,
      scatterplotNode, DG.VIEWER.HISTOGRAM, 0.3);
    const viewerNode = tableView.dockManager.dock(viewer, DG.DOCK_TYPE.TOP, null, MAP_NAME, 0.7);
    tableView.dockManager.dock(tableView.filters(), DG.DOCK_TYPE.LEFT, viewerNode, DG.VIEWER.FILTERS, 0.3);
  });
}

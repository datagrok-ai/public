import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';


const VIEWER_TABLES_PATH: {[key: string]: string} = {
  'Chord': 'energy_uk.csv',
  'Globe': 'geo/earthquakes.csv',
  'Radar': 'demog.csv',
  'Sankey': 'energy_uk.csv',
  'Sunburst': 'demog.csv',
  'SurfacePlot': 'files/surface-plot.csv',
  'Timelines': 'files/ae.csv',
  'Tree': 'demog.csv',
  'Word cloud': 'word_cloud.csv',
};


export async function viewerDemo(viewerName: string, options?: object | null) {
  const df = await (['SurfacePlot', 'Timelines'].includes(viewerName) ?
    grok.data.loadTable(`${_package.webRoot}${VIEWER_TABLES_PATH[viewerName]}`) :
    grok.data.getDemoTable(VIEWER_TABLES_PATH[viewerName]));

  const tableView = grok.shell.addTableView(df);

  grok.shell.windows.showHelp = true;
  grok.shell.windows.help.syncCurrentObject = false;

  if (viewerName === 'Globe') {
    DG.debounce(df.onSemanticTypeDetected, 1000).subscribe((_) => {
      const viewer = tableView.addViewer(viewerName, options);
      grok.shell.windows.help.showHelp(viewer.helpUrl);
      dockViewers(tableView, viewer, viewerName);
    });
    return;
  }

  const viewer = tableView.addViewer(viewerName, options);
  grok.shell.windows.help.showHelp(viewer.helpUrl);

  dockViewers(tableView, viewer, viewerName);
}

function dockViewers(tableView: DG.TableView, viewer: DG.Viewer, viewerName: string) {
  const rootNode = tableView.dockManager.rootNode;

  if (viewerName === DG.VIEWER.WORD_CLOUD) {
    tableView.dockManager.dock(tableView.filters(), DG.DOCK_TYPE.RIGHT, rootNode, 'Filters', 0.6);
    tableView.dockManager.dock(viewer, DG.DOCK_TYPE.TOP, null, viewerName, 0.7);
    return;
  }

  const scatterplotNode = tableView.dockManager.dock(tableView.addViewer('scatterplot'), DG.DOCK_TYPE.RIGHT,
    rootNode, 'Scatter plot', 0.5);
  tableView.dockManager.dock(tableView.addViewer('histogram'), DG.DOCK_TYPE.RIGHT, scatterplotNode, 'Histogram', 0.3);
  const viewerNode = tableView.dockManager.dock(viewer, DG.DOCK_TYPE.TOP, null, viewerName, 0.7);
  tableView.dockManager.dock(tableView.filters(), DG.DOCK_TYPE.LEFT, viewerNode, 'Filters', 0.3);
}

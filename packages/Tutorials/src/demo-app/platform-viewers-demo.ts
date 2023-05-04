import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';


// TODO: add customized datasets, fix the docking, etc.
const VIEWER_TABLES_PATH: {[key: string]: string} = {
  'Scatter plot': 'files/demog.csv',
  Histogram: 'files/demog.csv',
  'Line chart': 'sensors/eeg.csv',
  'Bar chart': 'files/demog.csv',
  'Pie chart': 'files/demog.csv',
  'Trellis plot': 'files/demog.csv',
  'Matrix plot': 'files/demog.csv',
  '3d scatter plot': 'files/demog.csv',
  'Density plot': 'files/demog.csv',
  'PC Plot': 'files/demog.csv',
  'Network diagram': 'got-s1-edges.csv',
  'Box plot': 'files/demog.csv',
  'Tree map': 'files/demog.csv',
  'Heat map': 'files/demog.csv',
  Statistics: 'files/demog.csv',
  'Correlation plot': 'sensors/eeg.csv',
  Calendar: 'files/demog.csv',
  Grid: 'files/demog.csv',
  Markup: 'files/demog.csv',
  'Tile Viewer': 'files/demog.csv',
  Form: 'files/demog.csv',
  'Shape Map': 'files/demog.csv',
  'Pivot table': 'files/demog.csv',
  Map: 'files/demog.csv',
};


export async function viewerDemo(viewerName: string, options?: object | null) {
  const df = ['Line chart', 'Network diagram', 'Correlation plot'].includes(viewerName) ?
    await grok.data.getDemoTable(VIEWER_TABLES_PATH[viewerName]) :
    await grok.data.loadTable(`${_package.webRoot}${VIEWER_TABLES_PATH[viewerName]}`);

  const tableView = grok.shell.addTableView(df);

  grok.shell.windows.showHelp = true;
  grok.shell.windows.help.syncCurrentObject = false;

  const viewer = tableView.addViewer(viewerName, options);
  grok.shell.windows.help.showHelp(viewer.helpUrl);

  dockViewers(tableView, viewer, viewerName);
}

function dockViewers(tableView: DG.TableView, viewer: DG.Viewer, viewerName: string) {
  const rootNode = tableView.dockManager.rootNode;

  // if (viewerName === 'WordCloud') {
  //   tableView.dockManager.dock(tableView.filters(), DG.DOCK_TYPE.RIGHT, rootNode, 'Filters', 0.6);
  //   tableView.dockManager.dock(viewer, DG.DOCK_TYPE.TOP, null, viewerName, 0.7);
  //   return;
  // }

  const scatterplotNode = tableView.dockManager.dock(tableView.addViewer('scatterplot'), DG.DOCK_TYPE.RIGHT,
    rootNode, 'Scatter plot', 0.5);
  tableView.dockManager.dock(tableView.addViewer('histogram'), DG.DOCK_TYPE.RIGHT, scatterplotNode, 'Histogram', 0.3);
  const viewerNode = tableView.dockManager.dock(viewer, DG.DOCK_TYPE.TOP, null, viewerName, 0.7);
  tableView.dockManager.dock(tableView.filters(), DG.DOCK_TYPE.LEFT, viewerNode, 'Filters', 0.3);
}

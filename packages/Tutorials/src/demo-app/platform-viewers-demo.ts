import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {delay} from '@datagrok-libraries/utils/src/test';


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
  Statistics: 'files/demog.csv',
  'Correlation plot': 'sensors/eeg.csv',
  Calendar: 'files/demog.csv',
  Grid: 'files/smiles.csv',
  Markup: 'files/demog.csv',
  'Tile Viewer': 'chem/sar_small.csv',
  Form: 'files/sar-small.csv',
  //'Shape Map': 'geo/us_2016_election_by_county.csv',
  'Pivot table': 'files/demog.csv',
  Filters: 'files/filters.csv',
};

const VIEWER_LAYOUTS_FILE_NAMES: {[key: string]: string} = {
  'Trellis plot': 'trellis-plot-viewer-layout.json',
  Form: 'form-viewer-layout.json',
  Grid: 'grid-layout.json',
  Filters: 'filters-layout.json',
};

const MARKUP_CONTENT = `# What's Markdown?

Markdown is a lightweight [markup language](https://en.wikipedia.org/wiki/Markdown)
used to create formatted text.

The **Markup viewer** allows you to display text using Markdown notation. 
You can edit the text using the properties panel of the viewer (press F4).

# Dynamic properties:

The text in the **Markup viewer** can include dynamic properties 
that automatically update when the data changes.

* Table name: **#{t.name}**
* Row count: **#{t.rowCount}**
* Selected: **#{t.selection.trueCount}**
* Filtered: **#{t.filter.trueCount}**
* Current row: **#{t.currentRow}**
* Current column: **#{t.currentCol}**
* Current cell value: **#{t.currentCell}**`;


async function getLayout(viewerName: string, df: DG.DataFrame): Promise<DG.ViewLayout> {
  const layout = await _package.files.readAsText(VIEWER_LAYOUTS_FILE_NAMES[viewerName]);
  return DG.ViewLayout.fromJson(layout.replaceAll("tableName", df.name));
}

async function loadViewerDemoLayout(tableView: DG.TableView, viewerName: string): Promise<void> {
  const layout = await getLayout(viewerName, tableView.dataFrame);
  tableView.loadLayout(layout);
  const viewer = Array.from(tableView.viewers)?.find((viewer) => viewer.type === viewerName);
  if (viewer)
    grok.shell.windows.help.showHelp(viewer.helpUrl);
}

export async function viewerDemo(viewerName: string, options?: object | null) {
  const df = ['Line chart', 'Network diagram', 'Correlation plot', 'Tile Viewer', 'Shape Map'].includes(viewerName) ?
    await grok.data.getDemoTable(VIEWER_TABLES_PATH[viewerName]) :
    await grok.data.loadTable(`${_package.webRoot}/${VIEWER_TABLES_PATH[viewerName]}`);
  await grok.data.detectSemanticTypes(df);

  const tableView = grok.shell.addTableView(df);

  grok.shell.windows.showContextPanel = false;
  grok.shell.windows.showHelp = true;
  grok.shell.windows.help.syncCurrentObject = false;

  if (['Form', 'Trellis plot', 'Grid', 'Filters'].includes(viewerName)) {
    await loadViewerDemoLayout(tableView, viewerName);
    return;
  }

  if (viewerName === DG.VIEWER.TILE_VIEWER) {
    DG.debounce(df.onSemanticTypeDetected, 800).subscribe((_) => {
      const viewer = tableView.addViewer(viewerName, options);
      grok.shell.windows.help.showHelp(viewer.helpUrl);
      dockViewers(tableView, viewer, viewerName);
    });
    return;
  }

  const viewer = tableView.addViewer(viewerName, options);

  if (viewerName === DG.VIEWER.NETWORK_DIAGRAM) {;
    viewer.props.node1ColumnName = 'Source';
    viewer.props.node2ColumnName = 'Target';
    (viewer.root.firstElementChild as HTMLElement).style.width = '100%';
  }

  if (viewerName === DG.VIEWER.MARKUP)
    viewer.props.content = MARKUP_CONTENT;

  grok.shell.windows.help.showHelp(viewer.helpUrl);

  dockViewers(tableView, viewer, viewerName);
}

function dockViewers(tableView: DG.TableView, viewer: DG.Viewer, viewerName: string) {
  const rootNode = tableView.dockManager.rootNode;

  if (viewerName === DG.VIEWER.MARKUP) {
    const viewerNode = tableView.dockManager.dock(viewer, DG.DOCK_TYPE.TOP, null, viewerName, 0.7);
    tableView.dockManager.dock(tableView.filters(), DG.DOCK_TYPE.LEFT, viewerNode, 'Filters', 0.3);
    return;
  }
  if (viewerName === DG.VIEWER.TILE_VIEWER) {
    tableView.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, viewerName, 0.75);
    return;
  }
  if (viewerName === DG.VIEWER.STATISTICS) {
    tableView.dockManager.dock(tableView.filters(), DG.DOCK_TYPE.RIGHT, rootNode, 'Filters', 0.6);
    tableView.dockManager.dock(viewer, DG.DOCK_TYPE.TOP, null, viewerName, 0.7);
    return;
  }

  const extraViewerName = viewerName === DG.VIEWER.SCATTER_PLOT ? DG.VIEWER.BAR_CHART : DG.VIEWER.SCATTER_PLOT;
  const extraViewer2Name = viewerName === DG.VIEWER.HISTOGRAM ? DG.VIEWER.BAR_CHART : DG.VIEWER.HISTOGRAM;

  const extraViewerNode = tableView.dockManager.dock(tableView.addViewer(extraViewerName), DG.DOCK_TYPE.RIGHT, rootNode, extraViewerName, 0.5);
  tableView.dockManager.dock(tableView.addViewer(extraViewer2Name), DG.DOCK_TYPE.RIGHT, extraViewerNode, 'histogram', 0.3);
  const viewerNode = tableView.dockManager.dock(viewer, DG.DOCK_TYPE.TOP, null, viewerName, 0.7);
  tableView.dockManager.dock(tableView.filters(), DG.DOCK_TYPE.LEFT, viewerNode, 'Filters', 0.3);
}

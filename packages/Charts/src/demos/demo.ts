import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';


const VIEWER_TABLES_PATH: {[key: string]: string} = {
  Chord: 'energy_uk.csv',
  Globe: 'geo/earthquakes.csv',
  GroupAnalysis: 'files/r-groups.csv',
  Radar: 'demog.csv',
  Sankey: 'energy_uk.csv',
  Sunburst: 'demog.csv',
  SurfacePlot: 'files/surface-plot.csv',
  Timelines: 'files/ae.csv',
  Tree: 'demog.csv',
  WordCloud: 'word_cloud.csv',
};


export async function viewerDemo(viewerName: string, options?: object | null) {
  const df = await (['GroupAnalysis', 'SurfacePlot', 'Timelines'].includes(viewerName) ?
    grok.data.loadTable(`${_package.webRoot}${VIEWER_TABLES_PATH[viewerName]}`) :
    grok.data.getDemoTable(VIEWER_TABLES_PATH[viewerName]));

  const tableView = grok.shell.addTableView(df);

  if (['Globe', 'GroupAnalysis'].includes(viewerName)) {
    DG.debounce(df.onSemanticTypeDetected, 800).subscribe((_) => {
      const viewer = tableView.addViewer(viewerName, options);
      if (viewerName === 'Globe')
        dockViewers(tableView, viewer, viewerName);
    });
    return;
  }

  const viewer = tableView.addViewer(viewerName, options);
  dockViewers(tableView, viewer, viewerName);
}

function dockViewers(tableView: DG.TableView, viewer: DG.Viewer, viewerName: string) {
  const rootNode = tableView.dockManager.rootNode;
  const node1 = tableView.dockManager.dock(tableView.addViewer('scatterplot'), 'right', rootNode,
    'Scatter plot', 0.5);
  tableView.dockManager.dock(tableView.addViewer('histogram'), 'right', node1, 'Histogram', 0.3);
  const node3 = tableView.dockManager.dock(viewer, 'up', null, viewerName, 0.7);
  tableView.dockManager.dock(tableView.filters(), 'left', node3, 'Filters', 0.3);
}

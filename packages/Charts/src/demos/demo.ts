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
        tableView.dockManager.dock(viewer, 'up', null, viewerName, 0.75);
      tableView.filters();
    });
    return;
  }

  const viewer = tableView.addViewer(viewerName, options);
  tableView.dockManager.dock(viewer, 'up', null, viewerName, 0.75);
  tableView.filters();
  //TODO: set grid instead of null, 'up' -> 'right' (for histogram and scatterplot)
}

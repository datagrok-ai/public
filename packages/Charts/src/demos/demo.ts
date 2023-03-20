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
  let df: DG.DataFrame;

  if (['GroupAnalysis', 'SurfacePlot', 'Timelines'].includes(viewerName))
    df = await grok.data.loadTable(`${_package.webRoot}${VIEWER_TABLES_PATH[viewerName]}`);
  else
    df = await grok.data.getDemoTable(VIEWER_TABLES_PATH[viewerName]);

  const tableView = grok.shell.addTableView(df);

  if (['Globe', 'GroupAnalysis'].includes(viewerName)) {
    DG.debounce(df.onSemanticTypeDetected, 300).subscribe((_) => tableView.addViewer(viewerName, options));
    return;
  }

  tableView.addViewer(viewerName, options);
}

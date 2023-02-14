import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';


const VIEWER_TABLES_PATH: {[key: string]: string} = {
  ChordViewer: 'energy_uk.csv',
  GlobeViewer: 'geo/earthquakes.csv',
  GroupAnalysisViewer: 'files/r-groups.csv',
  RadarViewer: 'demog.csv',
  SankeyViewer: 'energy_uk.csv',
  SunburstViewer: 'demog.csv',
  SurfacePlot: 'files/surface-plot.csv',
  TimelinesViewer: 'files/ae.csv',
  TreeViewer: 'demog.csv',
  WordCloudViewer: 'word_cloud.csv',
};


export async function viewerDemo(viewerName: string, options?: object | null) {
  let df: DG.DataFrame;

  if (['GroupAnalysisViewer', 'SurfacePlot', 'TimelinesViewer'].includes(viewerName))
    df = await grok.data.loadTable(`${_package.webRoot}${VIEWER_TABLES_PATH[viewerName]}`);
  else
    df = await grok.data.getDemoTable(VIEWER_TABLES_PATH[viewerName]);

  const tableView = grok.shell.addTableView(df);

  if (['GlobeViewer', 'GroupAnalysisViewer'].includes(viewerName)) {
    DG.debounce(df.onSemanticTypeDetected, 300).subscribe((_) => tableView.addViewer(viewerName, options));
    return;
  }

  tableView.addViewer(viewerName, options);
}

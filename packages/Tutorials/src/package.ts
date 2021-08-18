import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Track } from './tutorial';
import { ScatterPlotTutorial } from './tracks/eda/tutorials/scatter-plot';

export const _package = new DG.Package();

//name: Tutorials
//tags: app
export async function trackOverview() {
  const chem = new Track('Cheminformatics');
  const eda = new Track('Exploratory Data Analysis');
  const ml = new Track('Machine Learning');
  const tracks = DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'track', [chem, eda, ml].map((t) => t.name)),
    DG.Column.fromList('int', 'start', [10, 1, 15]),
    DG.Column.fromList('int', 'end', [29, 29, 29]),
  ]);
  const timelines = await tracks.plot.fromType('TimelinesViewer', {
    showZoomSliders: false,
    lineWidth: 20,
  });
  grok.shell.newView('Tracks', [timelines.root]);
  const sp = new ScatterPlotTutorial();
  await sp.run();
  eda.tutorials.push(sp);
  console.log(eda.name, eda.tutorials);
}

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Track, TutorialRunner } from './tutorial';
import { eda } from './tracks/eda';
import { ml } from './tracks/ml';


export const _package = new DG.Package();

//name: Tutorials
//tags: app
export async function trackOverview() {
  const chem = new Track('Cheminformatics');
  const tracks = DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'track', [chem, eda, ml].map((t) => t.name)),
    DG.Column.fromList('int', 'start', [10, 1, 15]),
    DG.Column.fromList('int', 'end', [29, 29, 29]),
  ]);
  const timelines = await tracks.plot.fromType('TimelinesViewer', {
    showZoomSliders: false,
    lineWidth: 20,
  });
  //grok.shell.newView('Tracks', [timelines.root]);
  const runEda = new TutorialRunner(eda);
  const runChem = new TutorialRunner(chem);
  const runMl = new TutorialRunner(ml);
  let root = ui.div([
    runEda.root,runChem.root,runMl.root,
    ui.panel([],{id:'tutorial-child-node'}),
  ], {style:{overflowY:'scroll'}});

  console.log(runEda);
  grok.shell.dockManager.dock(root, DG.DOCK_TYPE.RIGHT, null, 'Tutorials', 0.35);
}

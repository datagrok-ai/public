import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Track, TutorialRunner } from './tutorial';
import { eda } from './tracks/eda';
import { ml } from './tracks/ml';
import '../css/tutorial.css';


export const _package = new DG.Package();

//name: Tutorials
//tags: app
export async function trackOverview() {
  const chem = new Track('Cheminformatics');
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

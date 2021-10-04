import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Track, TutorialRunner } from './track';
import { chem } from './tracks/chem';
import { eda } from './tracks/eda';
import { da } from './tracks/data-access';
import { ml } from './tracks/ml';
import { TutorialWidget } from './widget';
import '../css/tutorial.css';


export const _package = new DG.Package();

//name: Tutorials
//tags: app
export async function trackOverview() {
  const tracks = [eda, chem, ml, da];
  let root = ui.div([
    ...tracks.map((track) => new TutorialRunner(track).root),
    ui.panel([],{id:'tutorial-child-node', style:{paddingTop:'10px'}}),
  ], 'tutorials-root');

  grok.shell.dockManager.dock(root, DG.DOCK_TYPE.RIGHT, null, 'Tutorials', 0.3);
  
}

//output: widget tutorial
//tags: dashboard
export function tutorialWidget(): DG.Widget {
  return new TutorialWidget();
}

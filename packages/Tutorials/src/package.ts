import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { TutorialRunner } from './tutorial-runner';
import { chem } from './tracks/chem';
import { eda } from './tracks/eda';
import { da } from './tracks/data-access';
import { ml } from './tracks/ml';
import { TutorialWidget } from './widget';
import '../css/tutorial.css';
import { Tutorial } from './tutorial';
import { Track } from './track';


export const _package = new DG.Package();

const externalTutorials: Tutorial[] = [];
const testTrack = new Track('Test Track', externalTutorials, '');
const tracks = [eda, chem, ml, da, testTrack];

//name: Tutorials
//tags: app
//top-menu: Help | Tutorials @Toolbox Help | Tutorials
export function trackOverview() {
  let root = ui.div([
    ...tracks.map((track) => new TutorialRunner(track).root),
    ui.panel([], { id: 'tutorial-child-node', style: { paddingTop: '10px' } }),
  ], 'tutorials-root');

  grok.shell.dockManager.dock(root, DG.DOCK_TYPE.RIGHT, null, 'Tutorials', 0.3);
}

//output: widget tutorial
export function tutorialWidget(): DG.Widget {
  return new TutorialWidget(...tracks.map((track) => new TutorialRunner(track)));
}

//tags: init
export async function tutorialsInit() {
  const tutorialFuncs = DG.Func.find({ tags: ['tutorial'] });
  for (const func of tutorialFuncs) {
    const tutorialList = await grok.functions.call(func.nqName);
    for (const t of tutorialList) {
      t.track = testTrack;
      externalTutorials.push(t);
    }
  }
}

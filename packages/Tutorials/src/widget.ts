/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Track, TutorialRunner } from './track';
import { Tutorial } from './tutorial';
import { chem } from './tracks/chem';
import { eda } from './tracks/eda';
import { ml } from './tracks/ml';
import '../css/tutorial.css';

export class TutorialWidget extends DG.Widget {
  caption: string;
  order: string;
  allTutorials: HTMLHeadingElement = ui.h2('');
  completeTutorials: HTMLHeadingElement = ui.h2('');
  totalProgress: HTMLHeadingElement = ui.h2('');

  constructor() {
    super(ui.panel([], 'tutorial-widget'));

    const runEda = new TutorialRunner(eda);
    const runChem = new TutorialRunner(chem);
    const runMl = new TutorialRunner(ml);

    let edaProgress = ui.element('progress');
    edaProgress.max = '0';
    edaProgress.value = '0';

    let chemProgress = ui.element('progress');
    chemProgress.max = '0';
    chemProgress.value = '0';

    let mlProgress = ui.element('progress');
    mlProgress.max = '0';
    mlProgress.value = '0';

    (async () => {

        console.clear();
        console.log(runEda.track)
        console.log(runEda.track.tutorials.length)
        console.log(await runEda.getCompleted(runEda.track.tutorials))

        let complete = await runEda.getCompleted(runEda.track.tutorials) + await runEda.getCompleted(runChem.track.tutorials) + await runMl.getCompleted(runMl.track.tutorials);
        let total = runEda.track.tutorials.length + runChem.track.tutorials.length + runMl.track.tutorials.length;
        let progress = 100/total*complete;

        this.completeTutorials.append(String(complete));
        this.allTutorials.append(String(total));
        this.totalProgress.append(String(Math.round(progress))+'%');

        edaProgress.max = runEda.track.tutorials.length;
        edaProgress.value = await runEda.getCompleted(runEda.track.tutorials);

        chemProgress.max = runChem.track.tutorials.length;
        chemProgress.value = await runChem.getCompleted(runChem.track.tutorials);
        
        mlProgress.max = runMl.track.tutorials.length;
        mlProgress.value = await runMl.getCompleted(runMl.track.tutorials);
    
    })();

    this.root.append(ui.divV([
        ui.divH([
            ui.divV([
                ui.label('Progress'),
                this.totalProgress
            ]),
            ui.divV([
                ui.label('Complete'),
                this.completeTutorials
            ]),
            ui.divV([
                ui.label('All Tutorials'),
                this.allTutorials
            ]),
        ]),
        ui.divV([
            ui.divH([
                ui.divText(runEda.track.name),
            ]),
            edaProgress
        ], 'tutorials-track'),
        ui.divV([
            ui.divH([
                ui.divText(runChem.track.name),
            ]),
            chemProgress
        ], 'tutorials-track'),
        ui.divV([
            ui.divH([
                ui.divText(runMl.track.name),
            ]),
            mlProgress
        ], 'tutorials-track'),
    ], 'tutorial'));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Tutorials');
    this.order = super.addProperty('order', DG.TYPE.STRING, '8');
  }
}

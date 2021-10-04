/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { Track, TutorialRunner } from './track';
import { Tutorial } from './tutorial';
import { chem } from './tracks/chem';
import { eda } from './tracks/eda';
import { ml } from './tracks/ml';
import '../css/tutorial.css';

export class TutorialWidget extends DG.Widget {
  caption: string;
  order: string;
  totalTracks: number = 0;
  totalTutorials:number = 0;
  totalCompleted:number = 0;
  totalProgress:number = 0;

  runners:TutorialRunner[];

  constructor(...runners:TutorialRunner[]) {
    super(ui.panel([], 'tutorial-widget'));
    
    this.runners = runners;
    let tracksRoot = ui.div([]);
    let tracks = runners.length;

    let tracksStas = ui.h1('');
    let tutorialsStats = ui.h1('');
    let completeStats = ui.h1('');
    let progressStats = ui.h1('');

    (async () => {
        this.totalTracks = runners.length;
        let i = 0;

        while(runners){
        
            let complete = await runners[i].getCompleted(runners[i].track.tutorials);
            let total = runners[i].track.tutorials.length;

            this.totalTutorials += total;
            this.totalCompleted += complete;
            this.totalProgress = 100/this.totalTutorials*this.totalCompleted;
            
            console.log('-------');

            let progressBar = ui.element('progress');
            progressBar.max = total
            progressBar.value = complete;
            
            tracksRoot.append(ui.divV([
                ui.divH([
                    ui.divText(runners[i].track.name),
                    ui.divText(String(complete)+' / '+String(total))
                ], 'widget-tutorials-track-details'),
                progressBar
            ], 'tutorials-track'));
            i++;

            if (i==runners.length)
                tracksStas.innerHTML = String(this.totalTracks);
            tutorialsStats.innerHTML = String(this.totalTutorials);
            completeStats.innerHTML = String(this.totalCompleted);
            progressStats.innerHTML = String(Math.round(this.totalProgress)+'%');
        }
        
    })();

    this.root.append(ui.divV([
        ui.divH([
            ui.divV([
                ui.label('Tracks'),
                tracksStas
            ]),
            ui.divV([
                ui.label('Tutorials'),
                tutorialsStats
            ]),
            ui.divV([
                ui.label('Complete'),
                completeStats
            ]),
            ui.divV([
                ui.label('Progress'),
                progressStats
            ]),
        ], 'widget-tutorials-summary'),
        tracksRoot,
    ], 'tutorial'));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Tutorials');
    this.order = super.addProperty('order', DG.TYPE.STRING, '8');
  }
}

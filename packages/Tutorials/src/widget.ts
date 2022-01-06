/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Track, TutorialRunner } from './track';
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

    (async () => {
        this.totalTracks = runners.length;
        let i = 0;

        while(i < this.totalTracks){
        
            let complete = await runners[i].getCompleted(runners[i].track.tutorials);
            let total = runners[i].track.tutorials.length;
            this.totalTutorials += total;
            this.totalCompleted += complete;
            this.totalProgress = 100/this.totalTutorials*this.totalCompleted;
            
            let root = runners[i].root;

            let trackroot = ui.divH([
                ui.divText(runners[i].track.name, 'widget-tutorial-track-name'),
                ui.divText(String(complete)+' / '+String(total), 'widget-tutorial-track-progress'),
            ], 'widget-tutorial-track-body');

            trackroot.addEventListener('click', function(){
                let dockRoot = ui.div([root,
                    ui.panel([],{id:'tutorial-child-node', style:{paddingTop:'10px'}}),
                  ], 'tutorials-root');
                grok.shell.dockManager.dock(dockRoot, DG.DOCK_TYPE.RIGHT, null, 'Tutorials', 0.3);
            })

            tracksRoot.append(ui.divV([
                trackroot,
                ui.div([],{
                    style:{
                        position: 'absolute',
                        width: String(Math.round(100/total*complete))+'%',
                        backgroundColor: 'rgba(32, 131, 213, 0.15)',
                        height: '100%',
                    }
                })
            ], 'widget-tutorial-track'));
            i++;
        }
        
    })();

    this.root.append(ui.divV([
        tracksRoot,
    ], 'tutorial'));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Tutorials');
    this.order = super.addProperty('order', DG.TYPE.STRING, '8');
  }
}
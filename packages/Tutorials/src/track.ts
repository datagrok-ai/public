import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { _package } from './package';
import { Tutorial } from './tutorial';
import { tutorials } from './tracks/chem';


/** A collection of tutorials */
export class Track {
  tutorials: Tutorial[];
  name: string;
  helpUrl:string
  completed: number = 0;

  constructor(name: string, tutorials: Tutorial[], helpUrl:string) {
    this.name = name;
    this.tutorials = tutorials;
    this.helpUrl = helpUrl;
    tutorials.forEach((t) => t.track = this);
  }
}

export class TutorialRunner {
  root: HTMLDivElement = ui.panel([], 'tutorial');
  track: Track;

  async run(t: Tutorial): Promise<void> {
    $('.tutorial').hide();
    $('#tutorial-child-node').html('');
    $('#tutorial-child-node').append(t.root);
    t.clearRoot();
    await t.run();
  }

  async getCompleted(tutorials: Tutorial[]) {
    let i = 0;
    let completed = 0;
    while (tutorials[i]) {
      await tutorials[i].updateStatus();
      if (tutorials[i].status == true) {
        completed++;
      }
      i++;
    }
    return completed;
  }

  constructor(track: Track, onStartTutorial?: (t: Tutorial) => Promise<void>) {
    let complete: number = 0;
    let total: number = track.tutorials.length;

    this.track = track;

    let progress = ui.element('progress');
    progress.max = '100';//track.tutorials.length;
    progress.value = '2';

    //Container for progress bar details (percents & completed tutorials)
    let progressDetails = ui.divH([], 'tutorials-track-details');

    (async () => {
      complete = await this.getCompleted(track.tutorials);
      track.completed = complete;
      if(total !=0){
        if(complete!=0){  
            progress.value = 100/total*complete 
            progressDetails.append(ui.divText(progress.value.toFixed() + '% complete'));
            progressDetails.append(ui.divText(complete+' / '+total));
        }else{   
        progressDetails.append(ui.divText('0% complete'));     
        progressDetails.append(ui.divText(complete+' / '+total));
        }
      }else{
        progressDetails.append(ui.divText('0% complete'));
        progressDetails.append(ui.divText(complete+' / '+total));
      }
    })();


    let trackRoot = ui.divV([
      ui.divH([ui.h1(track.name), ui.icons.info(()=>{})], 'tutorials-track-title'),
      ui.divH([progress]),
      progressDetails,
      ui.divV(track.tutorials.map((t) => {
        const el = new TutorialCard(t).root;
        el.addEventListener('click', () => {
          if (onStartTutorial == null) {
            this.run(t);
          } else {
            onStartTutorial(t);
          }
        });
        return el;
      })),
      ui.link('Clear storage', () => {
        let i = 0;
        while (track.tutorials[i]) {
          grok.dapi.userDataStorage.remove(Tutorial.DATA_STORAGE_KEY, track.tutorials[i].name);
          i++;
        }
        grok.shell.info('complete');
      }, '', '')
    ], 'tutorials-track');

    trackRoot.setAttribute('data-name', track.name);
    this.root.append(trackRoot)
  }

}

class TutorialCard {
  root: HTMLDivElement = ui.div();
  tutorial: Tutorial;

  constructor(tutorial: Tutorial) {
    this.tutorial = tutorial;

    let img = ui.image(`${_package.webRoot}images/${tutorial.name.toLowerCase().replace(/ /g, '-')}.png`, 90, 70);
    let icon = ui.div([], 'tutorials-card-status');
    icon.append(ui.iconFA('check', () => { }));
    $(icon).hide();
    $(icon).children().removeClass('fal');
    $(icon).children().addClass('far');

    (async () => {
      await tutorial.updateStatus();
      if (tutorial.status == true) {
        $(icon).show();
      }
    })();

    this.root = ui.divH([
      img,
      ui.tooltip.bind(
        ui.divV([
          ui.divText(tutorial.name, 'tutorials-card-title'),
          ui.divText(tutorial.description, 'tutorials-card-description')
        ]),
        `<b>${tutorial.name}</b><br>${tutorial.description}`),
      icon
    ], 'tutorials-card');
    this.root.setAttribute('data-name', tutorial.name);
  }
}

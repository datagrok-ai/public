import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { _package } from './package';
import { Track } from './track';
import { Tutorial } from './tutorial';


export class TutorialRunner {
  root: HTMLDivElement = ui.panel([], 'tutorial');
  track: Track;

  async run(t: Tutorial): Promise<void> {
    if (t instanceof TutorialSubstitute)
      t = await TutorialSubstitute.getTutorial(t, this.track);

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
      ui.divH([ui.h1(track.name), ui.button(ui.icons.help(()=>{}),()=>window.open(track.helpUrl,'_blank'), 'Read more about '+track.name+'\n'+track.helpUrl)], 'tutorials-track-title'),
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
      // ui.link('Clear storage', () => {
      //   let i = 0;
      //   while (track.tutorials[i]) {
      //     grok.dapi.userDataStorage.remove(Tutorial.DATA_STORAGE_KEY, track.tutorials[i].name);
      //     i++;
      //   }
      //   grok.shell.info('complete');
      // }, '', '')
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

    let img = ui.image(tutorial.imageUrl, 90, 70);
    let icon = ui.div([], 'tutorials-card-status');
    let title = ui.divText(tutorial.name, 'tutorials-card-title');
    let description = ui.divText(tutorial.description, 'tutorials-card-description');
    icon.append(ui.iconFA('check-circle', () => { }));
    $(icon).hide();
    $(icon).children().removeClass('fal');
    $(icon).children().addClass('fas');

    (async () => {
      await tutorial.updateStatus();
      if (tutorial.status == true) {
        $(icon).show();
        $(title).css('color','var(--grey-4)');
        $(description).css('color','var(--grey-4)');
        $(img).css('mix-blend-mode','luminosity');
        $(img).css('opacity','0.7');
      }
    })();

    this.root = ui.divH([
      img,
      ui.tooltip.bind(
        ui.divV([
          title,
          description
        ]),
        `<b>${tutorial.name}</b><br>${tutorial.description}`),
      icon
    ], 'tutorials-card');
    this.root.setAttribute('data-name', tutorial.name);
  }
}

export class TutorialSubstitute extends Tutorial {
  name: string;
  description: string;
  imageUrl: string;
  func: DG.Func;

  constructor(name: string, description: string, icon: string, func: DG.Func) {
    super();
    this.name = name;
    this.description = description;
    this.imageUrl = icon;
    this.func = func;
  }

  static async getTutorial(substitute: TutorialSubstitute, track?: Track): Promise<Tutorial> {
    const tutorial = await grok.functions.call(substitute.func.nqName);
    if (track)
      tutorial.track = track;
    return tutorial;
  }

  get steps() { return 1; }
  async _run() {}
}

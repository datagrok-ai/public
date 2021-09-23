import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { async, fromEvent, Observable } from 'rxjs';
import { filter, first, map } from 'rxjs/operators';
import { _package } from './package';
import { eda, tutorials } from './tracks/eda';


/** A base class for tutorials */
export abstract class Tutorial extends DG.Widget {
  abstract get name(): string;
  abstract get description(): string;

  track: Track | null = null;
  demoTable: string = 'demog.csv';
  get t(): DG.DataFrame {
    return grok.shell.t;
  }

  nextLink: HTMLAnchorElement = ui.link('next',
    '',
    'Go to the next tutorial', {
      classes: 'grok-tutorial-next',
      style: { display: 'none' },
    });
  header: HTMLHeadingElement = ui.h2('');
  subheader: HTMLHeadingElement = ui.h3('');
  activity: HTMLDivElement = ui.div([], 'grok-tutorial-description');
  status: boolean = false;

  static DATA_STORAGE_KEY: string = 'tutorials';

  async updateStatus(): Promise<void> {
    const info = await grok.dapi.userDataStorage.getValue(Tutorial.DATA_STORAGE_KEY, this.name);
    this.status = info ? true : false;
  }

  constructor() {
    super(ui.div([]));
    this.updateStatus();

    this.root.append(this.header);
    this.root.append(this.subheader);
    this.root.append(this.activity);
    this.root.append(this.nextLink);
  }

  protected abstract _run(): Promise<void>;

  async run(): Promise<void> {
    if (this.demoTable) {
      grok.shell.addTableView(await grok.data.getDemoTable(this.demoTable));
    }

    await this._run();

    this.title('Congratulations!');
    this.describe('You have successfully completed this tutorial.');
  
    await grok.dapi.userDataStorage.postValue(Tutorial.DATA_STORAGE_KEY, this.name, new Date().toUTCString());

    const tutorials = this.track?.tutorials;
    if (!tutorials) {
      console.error('The launched tutorial is not bound to any track.');
      return;
    }
    let id = this.track!.tutorials.indexOf(this);
    
    if (id < tutorials.length - 1){
      this.root.append(ui.div([
        ui.bigButton('Start', () => {
          $('#tutorial-child-node').html('');
          $('#tutorial-child-node').append(tutorials[++id].root);
          tutorials[id].run();
        }),
        ui.button('Cancel', () => {
          $('.tutorial').show();
          $('#tutorial-child-node').html('');
        })
      ]))
    } else if (id == tutorials.length - 1) {
      this.root.append(ui.div([ui.bigButton('Complete', () => {
        $('.tutorial').show();
        $('#tutorial-child-node').html('');
      })]))
    }
  }

  title(text: string): void {
    this.activity.append(ui.h2(text));
  }

  describe(text: string): void {
    const div = ui.div([], 'grok-tutorial-description-entry');
    div.innerHTML = text;
    this.activity.append(div);
    this._scroll();
  }

  async action(instructions: string, completed: Observable<any>,
    hint?: HTMLElement | null, hintSub?: DG.StreamSubscription | null): Promise<void> {
    hint?.classList.add('grok-tutorial-target-hint');
    const instructionDiv = ui.divText(instructions, 'grok-tutorial-entry-instruction');
    const entry = ui.div([instructionDiv], 'grok-tutorial-entry');
    this.activity.append(entry);
    this._scroll();
    await this.firstEvent(completed);
    instructionDiv.classList.add('grok-tutorial-entry-success');
    hintSub?.cancel();
    hint?.classList?.remove('grok-tutorial-target-hint');
  }

  protected _scroll(): void {
    this.root.scrollTop = this.root.scrollHeight;
  }

  clearRoot(): void {
    $(this.root).children().each((idx, el) => $(el).empty());
  }

  firstEvent(eventStream: Observable<any>): Promise<void> {
    return new Promise<void>((resolve, reject) => {
      eventStream.pipe(first()).subscribe((_: any) => resolve());
    });
  }

  protected async openPlot(name: string, check: (viewer: DG.Viewer) => boolean): Promise<DG.Viewer> {
    // TODO: Expand toolbox / accordion API coverage
    const getViewerIcon = (el: HTMLElement) => $(el).find(`i.svg-${name.replace(' ', '-')}`).get()[0];
    const view = grok.shell.v as DG.View;
    let viewer: DG.Viewer;

    await this.action(`Open ${name}`,
      grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
        const found = check(data.args.viewer);
        if (found) {
          viewer = data.args.viewer;
        }
        return found;
      })),
      view.type === 'TableView' ? getViewerIcon((<DG.TableView>view).toolboxPage.accordion.root) : null,
    );

    return viewer!;
  }

  protected async dlgInputAction(dlg: DG.Dialog, instructions: string, caption: string, value: string) {
    const inp = dlg.inputs.filter((input: DG.InputBase) => input.caption == caption)[0];
    if (inp == null) return;
    await this.action(instructions,
      new Observable((subscriber: any) => {
        if (inp.stringValue === value) subscriber.next(inp.stringValue);
        inp.onChanged(() => {
          if (inp.stringValue === value) subscriber.next(inp.stringValue);
        });
      }),
      inp.root,
    );
  }

  /** A helper method to access text inputs in a view. */
  protected async textInpAction(root: HTMLElement, instructions: string, caption: string, value: string) {
    const inputRoot = $(root)
      .find('div.ui-input-text.ui-input-root')
      .filter((idx, inp) => $(inp).find('label.ui-label.ui-input-label')[0]?.textContent === caption)[0];
    if (inputRoot == null) return;
    const input = $(inputRoot).find('input.ui-input-editor')[0] as HTMLInputElement;
    const source = fromEvent(input, 'input').pipe(map((_) => input.value), filter((val) => val === value));
    await this.action(instructions, source, inputRoot);
  }

  /** Prompts the user to open a view of the specified type, waits for it to open and returns it. */
  protected async openViewByType(instructions: string, type: string): Promise<DG.View> {
    let view: DG.View;

    await this.action(instructions, grok.events.onViewAdded.pipe(filter((v) => {
      if (v.type === type) {
        view = v;
        return true;
      }
      return false;
    })));

    return view!;
  }
}


/** A collection of tutorials */
export class Track {
  tutorials: Tutorial[];
  name: string;

  constructor(name: string, ...tutorials: Tutorial[]) {
    this.name = name;
    this.tutorials = tutorials;
    tutorials.forEach((t) => t.track = this);
  }
}

export class TutorialRunner {
  root: HTMLDivElement = ui.panel([],'tutorial');

  async run(t: Tutorial): Promise<void> {
    $('.tutorial').hide();
    $('#tutorial-child-node').append(t.root);

    await t.run();
  }

  async getCompleted(tutorials: Tutorial[]){
    let i = 0;
    let completed = 0;
    while (tutorials[i]){
      await tutorials[i].updateStatus();
      if (tutorials[i].status == true){
        completed++
      }
      i++
    }
    return completed
  }

  constructor(track: Track, onStartTutorial?: (t: Tutorial) => Promise<void>) {
    let iconList = ui.iconFA('list');
    iconList.style.color = 'var(--grey-4)';
    
    let listLengh = {
      style:{
        color:'var(--grey-6)',
        margin:'0px 0px 0px 5px'
      }
    }

    let trackTitle = ui.h1(track.name);
    trackTitle.style.margin = '0 0 10px 0';

    let progress = ui.element('progress');
    progress.max = track.tutorials.length;
    progress.value = '0'
    progress.style.margin = '5px 0';
    progress.style.width = '100%';

    (async () => {
        progress.value = await this.getCompleted(track.tutorials);   
    })();

    this.root.append(ui.divV([
      ui.divH([trackTitle,
        ui.button(ui.iconFA('sync',()=>{}),()=>{
          let i = 0;
          while(track.tutorials[i]){
           grok.dapi.userDataStorage.remove(Tutorial.DATA_STORAGE_KEY, track.tutorials[i].name);
           i++;
          }
          grok.shell.info('complete');
         }),
      ], {style:{alignItems:'flex-end'}}),
      ui.divH([iconList, ui.divText(' Tutorials: '+track.tutorials.length,listLengh),      ]),
      ui.divH([progress]),
      ui.divV(track.tutorials.map((t) => {
        const el = new TutorialCard(t).root;
        el.style.padding = '5px';
        el.addEventListener('mouseover', () => {
          el.style.cursor = 'pointer';
          el.style.background = 'var(--grey-1)';
        });
        el.addEventListener('mouseout', () => {
          el.style.cursor = 'normal';
          el.style.background = 'white';
        });
        el.addEventListener('click', () => {
          if (onStartTutorial == null) {
            this.run(t);
          } else {
            onStartTutorial(t);
          }
        });
        return el;
      }))
    ],{style:{marginBottom:'15px'}}));
  }

}

class TutorialCard {
  root: HTMLDivElement = ui.div();
  tutorial: Tutorial;

  constructor(tutorial: Tutorial) {
    this.tutorial = tutorial;
    
    let tutorialTitle = {
      style:{
        fontWeight:'bold',
        color:'var(--grey-6)',
        margin:'0px 0px 5px 0px'
      }
    };

    let tutorialDesciption = {
      style:{
        fontWeight:'normal',
        color:'var(--grey-4)',
        margin:'0px',
        maxHeight: '50px',
        textOverflow: 'ellipsis',
        display: '-webkit-box',
        '-webkit-line-clamp':'3',
        '-webkit-box-orient': 'vertical',
        overflow: 'hidden',
      }
    }

    let img = ui.image( `${_package.webRoot}images/${tutorial.name.toLowerCase().replace(/ /g, '-')}.png`,90, 70);
    img.style.border = '1px solid var(--grey-1)';
    img.style.backgroundSize = 'cover';

    let icon = ui.div([], {style:{maxWidth:'20px', position:'absolute', right:'15px'}});
    

    (async()=>{
      await tutorial.updateStatus();
      if(tutorial.status == true){
        icon.append(ui.iconFA('check', ()=>{console.log(tutorial.status)}));
        icon.style.color = 'var(--green-2)';
      }
    })();
    

    this.root = ui.divH([
      ui.div([img], {style:{minWidth:'90px', marginRight:'10px'}}),
      ui.tooltip.bind(
      ui.divV([
        ui.divText(tutorial.name, tutorialTitle),
        ui.divText(tutorial.description, tutorialDesciption)
      ]), 
      `<b>${tutorial.name}</b><br>${tutorial.description}`),
      icon
    ], {style:{position:'relative'}});
  }
}

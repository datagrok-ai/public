import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { fromEvent, Observable, Subject } from 'rxjs';
import { filter, first, map } from 'rxjs/operators';
import { _package } from './package';
import { Track, TutorialRunner } from './track';


/** A base class for tutorials */
export abstract class Tutorial extends DG.Widget {
  abstract get name(): string;
  abstract get description(): string;
  abstract get steps():number;

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
  header: HTMLHeadingElement = ui.h1('');
  headerDiv: HTMLDivElement = ui.divH([],'tutorials-root-header');
  subheader: HTMLHeadingElement = ui.h3('');
  activity: HTMLDivElement = ui.div([], 'tutorials-root-description');
  status: boolean = false;
  closed: boolean = false;
  activeHints: HTMLElement[] = [];
  progressDiv: HTMLDivElement = ui.divV([],'tutorials-root-progress');
  progress: HTMLProgressElement = ui.element('progress');
  progressSteps: HTMLDivElement = ui.divText('');

  static DATA_STORAGE_KEY: string = 'tutorials';

  async updateStatus(): Promise<void> {
    const info = await grok.dapi.userDataStorage.getValue(Tutorial.DATA_STORAGE_KEY, this.name);
    this.status = info ? true : false;
  }

  constructor() {

    super(ui.div([], 'tutorials-track'));
    this.updateStatus();
    this.root.append(this.headerDiv);
    this.root.append(this.progressDiv);
    this.progress.max = 0;
    this.progress.value = 1;
    this.root.append(this.subheader);
    this.root.append(this.activity);
    this.root.append(this.nextLink);
  }

  protected abstract _run(): Promise<void>;

  async run(): Promise<void> {
    this.progressDiv.append(this.progress);
    this.progress.max = this.steps;

    this.progressSteps = ui.divText('Step: '+String(this.progress.value)+' of '+this.steps);
    this.progressDiv.append(this.progressSteps);

    const tutorials = this.track?.tutorials;
    if (!tutorials) {
      console.error('The launched tutorial is not bound to any track.');
      return;
    }

    let id = tutorials.indexOf(this);
    
    
    let closeTutorial = ui.button(ui.iconFA('times-circle'),()=>{
      if (this.status != true){
        this.status = true;
        updateProgress(this.track);
      }  
      this.clearRoot();
      this.closed = true;
      this.onClose.next();
      // this._removeHints(this.activeHints);
      //grok.shell.tableView(this.t.name).close();
      $('.tutorial').show();
      $('#tutorial-child-node').html('');
    });

    closeTutorial.style.minWidth = '30px';
    this.headerDiv.append(this.header);
    this.headerDiv.append(closeTutorial);

    if (this.demoTable) {
      grok.shell.addTableView(await grok.data.getDemoTable(this.demoTable));
    }
    this.closed = false;
    await this._run();

    this.title('Congratulations!');
    this.describe('You have successfully completed this tutorial.');

    console.clear();
    console.log(id);

    await grok.dapi.userDataStorage.postValue(Tutorial.DATA_STORAGE_KEY, this.name, new Date().toUTCString());
    
    function updateProgress(track:any){
      track.completed++;
      $(`.tutorials-track[data-name ='${track?.name}']`)
      .find(`.tutorials-card[data-name='${track.tutorials[id].name}']`)
      .children('.tutorials-card-status').show();
      $(`.tutorials-track[data-name ='${track?.name}']`).find('progress').prop('value', String(100/track.tutorials.length*(track.completed)));
      $(`.tutorials-track[data-name ='${track?.name}'] > .tutorials-track-details`).children().first().text($(`.tutorials-track[data-name ='${track?.name}']`).find('progress').prop('value')+'% complete');
      $(`.tutorials-track[data-name ='${track?.name}'] > .tutorials-track-details`).children().last().text(String(track.completed+' / '+track.tutorials.length));
    }


    if (id < tutorials.length - 1){
      this.root.append(ui.divV([
        ui.divText('Next "'+tutorials[id+1].name+'"', {style:{margin:'5px 0'}}),
        ui.divH([
          ui.bigButton('Start', () => {
            console.log(id)
            if (this.status != true){
              this.status = true;
              updateProgress(this.track);
            }  
            this.clearRoot();
            $('#tutorial-child-node').html('');
            $('#tutorial-child-node').append(tutorials[id+1].root);
            tutorials[id+1].run();
          }),
          ui.button('Cancel', ()=>{
            console.log(id);
            if (this.status != true){
              this.status = true;
              updateProgress(this.track);
            }  
            //grok.shell.tableView(this.t.name).close();
            this.clearRoot();
            $('.tutorial').show();
            $('#tutorial-child-node').html('');
          })
        ], {style:{marginLeft:'-4px'}})
      ]))
    } else if (id == tutorials.length - 1) {
      this.root.append(ui.div([
        ui.divText(this.track?.name+' complete!'),
        ui.bigButton('Complete', ()=>{
          if (this.status != true){
            this.status = true;
            updateProgress(this.track);
          }  
          //grok.shell.tableView(this.t.name).close();
          this.clearRoot();
          $('.tutorial').show();
          $('#tutorial-child-node').html('');
        })
      ]))
    }
  }

  title(text: string): void {
    this.activity.append(ui.h3(text));
  }

  describe(text: string): void {
    const div = ui.div([]);
    div.innerHTML = text;
    this.activity.append(div);
    this._scroll();
  }

  _placeHint(hint: HTMLElement) {
    hint.classList.add('tutorials-target-hint');
    let hintIndicator = ui.element('div');
    hintIndicator.classList = 'blob';
    hint?.append(hintIndicator);

    let width = hint ? hint.clientWidth : 0;
    let height = hint ? hint.clientHeight : 0;

    let hintnode = hint?.getBoundingClientRect();
    let indicatornode = hintIndicator?.getBoundingClientRect();

    let hintPosition = $(hint).css('position');

    if(hintPosition == 'absolute'){
      $(hintIndicator).css('position','absolute');
      $(hintIndicator).css('left','0');
      $(hintIndicator).css('top','0');
    }
    if(hintPosition == "relative"){
      $(hintIndicator).css('position','absolute');
      $(hintIndicator).css('left','0');
      $(hintIndicator).css('top','0');
    }

    if(hintPosition == 'static'){
      $(hintIndicator).css('position','absolute');
      if(hintnode.left+1 == indicatornode.left)
        $(hintIndicator).css('margin-left',0)
      else
        $(hintIndicator).css('margin-left',-width)
      if(hintnode.top+1 == indicatornode.top)
        $(hintIndicator).css('margin-top',0)
      else
        $(hintIndicator).css('margin-top',-height) 
    }
  }

  _removeHints(hint: HTMLElement | HTMLElement[]) {
    const removeHint = (h: HTMLElement) => {
      $(h).find('div.blob')[0]?.remove();
      h.classList.remove('tutorials-target-hint');
    };
    if (hint instanceof HTMLElement) {
      removeHint(hint);
    } else if (Array.isArray(hint)) {
      hint.forEach((h) => removeHint(h));
    }
  }

  async action(instructions: string, description:string, completed: Observable<any> | Promise<void>,
    hint?: HTMLElement | HTMLElement[] | null): Promise<void> {
    if (this.closed) {
      return;
    }
    this.activeHints.length = 0;
    if (hint instanceof HTMLElement) {
      this.activeHints.push(hint);
      this._placeHint(hint);
    } else if (Array.isArray(hint)) {
      this.activeHints.push(...hint);
      hint.forEach((h) => this._placeHint(h));
    }

    const instructionDiv = ui.divText(instructions, 'grok-tutorial-entry-instruction');
    const descriptionDiv = ui.divText('', {style:{margin:'0px 0px 0px 15px'}});
    const instructionIndicator = ui.div([],'grok-tutorial-entry-indicator')
    const entry = ui.divH([
        instructionIndicator,
        instructionDiv
    ], 'grok-tutorial-entry');
    
    this.activity.append(entry);
    descriptionDiv.innerHTML = description;
    this.activity.append(descriptionDiv);
    this._scroll();

    const currentStep = completed instanceof Promise ? completed : this.firstEvent(completed);
    await currentStep;

    instructionDiv.classList.add('grok-tutorial-entry-success');
    instructionIndicator.classList.add('grok-tutorial-entry-indicator-success');
    $(descriptionDiv).hide();
    this.progress.value++;
    this.progressSteps.innerHTML = '';
    this.progressSteps.append('Step: '+String(this.progress.value)+' of '+this.steps);
    
    if (hint != null)
      this._removeHints(hint);
  }

  protected _scroll(): void {
    this.root.scrollTop = this.root.scrollHeight;
  }

  clearRoot(): void {
    this.progress.value = 0;
    $(this.root).children().each((idx, el) => $(el).empty());
  }

  firstEvent(eventStream: Observable<any>): Promise<void> {
    return new Promise<void>((resolve, reject) => {
      const eventSub = eventStream.pipe(first()).subscribe((_: any) => resolve());
      const closeSub = this.onClose.subscribe(() => {
        eventSub.unsubscribe();
        closeSub.unsubscribe();
        this._removeHints(this.activeHints);
        reject();
      });
    }).catch((_) => console.log('Closing tutorial', this.name));
  }

  _onClose: Subject<void> = new Subject();
  get onClose() { return this._onClose; }

  /** Prompts the user to open a viewer of the specified type and returns it. */
  protected async openPlot(name: string,description: string, check: (viewer: DG.Viewer) => boolean): Promise<DG.Viewer> {
    // TODO: Expand toolbox / accordion API coverage
    const getViewerIcon = (el: HTMLElement) => {
      const selector = name == 'filters' ? 'i.fa-filter' : `i.svg-${name.replace(' ', '-')}`;
      const icon = $(el).find(selector);
      return icon[0];
    }
    const view = grok.shell.v as DG.View;
    let viewer: DG.Viewer;

    await this.action(`Open ${name}`, description,
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

  /** Prompts the user to put the specified value into a dialog input. */
  protected async dlgInputAction(dlg: DG.Dialog, instructions: string, description:string, caption: string, value: string) {
    const inp = dlg.inputs.filter((input: DG.InputBase) => input.caption == caption)[0];
    if (inp == null) return;
    await this.action(instructions, description,
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
  protected async textInpAction(root: HTMLElement, instructions: string, description:string, caption: string, value: string) {
    const inputRoot = $(root)
      .find('div.ui-input-text.ui-input-root')
      .filter((idx, inp) => $(inp).find('label.ui-label.ui-input-label')[0]?.textContent === caption)[0];
    if (inputRoot == null) return;
    const input = $(inputRoot).find('input.ui-input-editor')[0] as HTMLInputElement;
    const source = fromEvent(input, 'input').pipe(map((_) => input.value), filter((val) => val === value));
    await this.action(instructions, description, source, inputRoot);
  }

  /** Prompts the user to open a view of the specified type, waits for it to open and returns it. */
  protected async openViewByType(instructions: string, description: string, type: string): Promise<DG.View> {
    let view: DG.View;

    await this.action(instructions, description, grok.events.onViewAdded.pipe(filter((v) => {
      if (v.type === type) {
        view = v;
        return true;
      }
      return false;
    })));

    return view!;
  }

  /** Prompts the user to open a dialog with the specified title, waits for it to open and returns it. */
  protected async openDialog(instructions: string, description: string, title: string, hint: HTMLElement | null = null): Promise<DG.Dialog> {
    let dialog: DG.Dialog;

    await this.action(instructions, description, grok.events.onDialogShown.pipe(filter((dlg) => {
      if (dlg.title === title) {
        dialog = dlg;
        return true;
      }
      return false;
    })), hint);

    return dialog!;
  }

  /** Prompts the user to select a menu item in the context menu. */
  protected async contextMenuAction(instructions: string, description:string, label: string, hint: HTMLElement | null = null): Promise<void> {
    const commandClick =  new Promise<void>((resolve, reject) => {
      const sub = grok.events.onContextMenu.subscribe((data) => {
        data.args.menu.onContextMenuItemClick.pipe(
          filter((mi) => (new DG.Menu(mi)).toString() === label),
          first()).subscribe((_: any) => {
            sub.unsubscribe();
            resolve();
          });
      });
    });

    await this.action(instructions, description, commandClick, hint);
  }
}

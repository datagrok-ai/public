import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { fromEvent, Observable } from 'rxjs';
import { filter, first, map } from 'rxjs/operators';
import { _package } from './package';
import { Track, TutorialRunner } from './track';
import { View } from 'datagrok-api/dg';
import { node } from 'webpack';


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
  subheader: HTMLHeadingElement = ui.h3('');
  activity: HTMLDivElement = ui.div([], 'tutorials-root-description');
  status: boolean = false;
  progress: HTMLProgressElement = ui.element('progress');

  static DATA_STORAGE_KEY: string = 'tutorials';

  async updateStatus(): Promise<void> {
    const info = await grok.dapi.userDataStorage.getValue(Tutorial.DATA_STORAGE_KEY, this.name);
    this.status = info ? true : false;
  }

  constructor() {

    super(ui.div([], 'tutorials-track'));
    this.updateStatus();
    this.root.append(this.header);
    this.root.append(this.subheader);
    this.root.append(this.progress);
    this.progress.max = 0;
    this.progress.value = 0;
    this.root.append(this.activity);
    this.root.append(this.nextLink);
  }

  protected abstract _run(): Promise<void>;

  async run(): Promise<void> {
    this.progress.max = this.steps;

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

    let id = tutorials.indexOf(this);

    function updateProgress(track:any){
      track.completed++;
      $(`.tutorials-track[data-name ='${track?.name}']`)
      .find(`.tutorials-card[data-name='${track.tutorials[id].name}']`)
      .children('.tutorials-card-status').show();
      $(`.tutorials-track[data-name ='${track?.name}']`).find('progress').prop('value', String(100/track.tutorials.length*(track.completed)));
      $(`.tutorials-track[data-name ='${track?.name}'] > .tutorials-track-details`).children().first().text($(`.tutorials-track[data-name ='${track?.name}']`).find('progress').prop('value')+'% complete');
      $(`.tutorials-track[data-name ='${track?.name}'] > .tutorials-track-details`).children().last().text(String(track.completed+' / '+track.tutorials.length));
    }

    const switchToMainView = () => {
      $('.tutorial').show();
      if (this.status != true){
        updateProgress(this.track);
      }  
      $('#tutorial-child-node').html('');
      grok.shell.tableView(this.t.name).close();
    };
    
    if (id < tutorials.length - 1){
      this.root.append(ui.div([
        ui.bigButton('Start', () => {
          if (this.status != true){
            console.log('update completed')
            updateProgress(this.track);
          }  
          $('#tutorial-child-node').html('');
          $('#tutorial-child-node').append(tutorials[++id].root);
          tutorials[id].run();
        }),
        ui.button('Cancel', switchToMainView)
      ]))
    } else if (id == tutorials.length - 1) {
      this.root.append(ui.div([ui.bigButton('Complete', switchToMainView)]))
    }
  }

  title(text: string): void {
    this.activity.append(ui.h1(text));
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

    console.clear();
    console.log('hint:'+$(hint).css('position'));
    console.log('hint top:'+hintnode.top+' left:'+hintnode.left);
    console.log('indicator top:'+indicatornode.top+' left:'+indicatornode.left);


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

    
   // $(hintIndicator).css('margin-left', width);
   // $(hintIndicator).css('margin-top', -height);
  }

  _removeHint(hint: HTMLElement) {
    $(hint).find('div.blob')[0]?.remove();
    hint.classList.remove('tutorials-target-hint');
  }

  async action(instructions: string, completed: Observable<any> | Promise<void>,
    hint?: HTMLElement | HTMLElement[] | null): Promise<void> {
    if (hint instanceof HTMLElement) {
      this._placeHint(hint);
    } else if (Array.isArray(hint)) {
      hint.forEach((h) => this._placeHint(h));
    }

    const instructionDiv = ui.divText(instructions, 'grok-tutorial-entry-instruction');
    const instructionIndicator = ui.div([],'grok-tutorial-entry-indicator')
    const entry = ui.divH([instructionIndicator,instructionDiv], 'grok-tutorial-entry');
    this.activity.append(entry);
    this._scroll();
    if (completed instanceof Promise) {
      await completed;
    } else {
      await this.firstEvent(completed);
    }
    instructionDiv.classList.add('grok-tutorial-entry-success');
    instructionIndicator.classList.add('grok-tutorial-entry-indicator-success')
    this.progress.value++;

    if (hint instanceof HTMLElement) {
      this._removeHint(hint);
    } else if (Array.isArray(hint)) {
      hint.forEach((h) => this._removeHint(h));
    }
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
      eventStream.pipe(first()).subscribe((_: any) => resolve());
    });
  }

  /** Prompts the user to open a viewer of the specified type and returns it. */
  protected async openPlot(name: string, check: (viewer: DG.Viewer) => boolean): Promise<DG.Viewer> {
    // TODO: Expand toolbox / accordion API coverage
    const getViewerIcon = (el: HTMLElement) => {
      const selector = name == 'filters' ? 'i.fa-filter' : `i.svg-${name.replace(' ', '-')}`;
      const icon = $(el).find(selector);
      return icon[0];
    }
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

  /** Prompts the user to put the specified value into a dialog input. */
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

  /** Prompts the user to open a dialog with the specified title, waits for it to open and returns it. */
  protected async openDialog(instructions: string, title: string, hint: HTMLElement | null = null): Promise<DG.Dialog> {
    let dialog: DG.Dialog;

    await this.action(instructions, grok.events.onDialogShown.pipe(filter((dlg) => {
      if (dlg.title === title) {
        dialog = dlg;
        return true;
      }
      return false;
    })), hint);

    return dialog!;
  }

  /** Prompts the user to select a menu item in the context menu. */
  protected async contextMenuAction(instructions: string, label: string, hint: HTMLElement | null = null): Promise<void> {
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

    await this.action(instructions, commandClick, hint);
  }
}

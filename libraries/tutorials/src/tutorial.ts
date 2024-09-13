import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {fromEvent, interval, Observable, Subject} from 'rxjs';
import {filter, first, map} from 'rxjs/operators';
import {Track} from './track';


/** A base class for tutorials */
export abstract class Tutorial extends DG.Widget {
  abstract get name(): string;
  abstract get description(): string;
  abstract get steps(): number;

  track: Track | null = null;
  prerequisites: TutorialPrerequisites = {};
  demoTable: string = 'demog.csv';
  currentSection: HTMLElement | undefined;
  // manualMode: boolean = false;
  private _t: DG.DataFrame | null = null;

  get t(): DG.DataFrame | null {
    return this._t;
  }

  set t(df: DG.DataFrame | null) {
    this._t = df;
  }

  get url(): string {
    const removeSpaces = (s: string) => s.split(' ').join('');
    const root = window.location.origin;
    return `${root}/apps/tutorials/Tutorials/${removeSpaces(this.track!.name)}/${removeSpaces(this.name)}`;
  }

  imageUrl: string = '';
  nextLink: HTMLAnchorElement = ui.link('next',
    '',
    'Go to the next tutorial', {
      classes: 'grok-tutorial-next',
      style: {display: 'none'},
    });
  mainHeader: HTMLDivElement = ui.div([], 'tutorials-main-header');
  header: HTMLHeadingElement = ui.h1('');
  headerDiv: HTMLDivElement = ui.divH([], 'tutorials-root-header');
  subheader: HTMLHeadingElement = ui.h3('');
  activity: HTMLDivElement = ui.div([], 'tutorials-root-description');
  status: boolean = false;
  closed: boolean = false;
  activeHints: HTMLElement[] = [];
  progressDiv: HTMLDivElement = ui.divV([], 'tutorials-root-progress');
  progress: HTMLProgressElement = ui.element('progress');
  progressSteps: HTMLDivElement = ui.divText('');

  static DATA_STORAGE_KEY: string = 'tutorials';
  static SERVICES: {[service: string]: string} = {
    'jupyter': 'Jupyter',
    'grokCompute': 'GrokCompute',
    'grokConnect': 'Grok Connect',
    'h2o': 'H2O',
  };

  async updateStatus(): Promise<void> {
    const info = await grok.userSettings.getValue(Tutorial.DATA_STORAGE_KEY, this.name);
    this.status = !!info;
  }

  constructor() {
    super(ui.div([], 'tutorials-track'));
    this.updateStatus();
    this.progress.max = 0;
    this.progress.value = 1;
    this.mainHeader.append(this.headerDiv, this.progressDiv);
    this.root.append(this.mainHeader);
    this.root.append(this.subheader);
    this.root.append(this.activity);
    this.root.append(this.nextLink);
  }

  protected abstract _run(): Promise<void>;

  async run(): Promise<void> {
    this._addHeader();

    const tutorials = this.track?.tutorials;
    if (!tutorials) {
      console.error('The launched tutorial is not bound to any track.');
      return;
    }

    if (this.prerequisites.packages && Array.isArray(this.prerequisites.packages) &&
      this.prerequisites.packages.length > 0) {
      const missingPackages = [];
      for (const p of this.prerequisites.packages) {
        const packages = await grok.dapi.packages.list({filter: `shortName = "${p}"`});
        if (!packages.length || !(packages[0] instanceof DG.Package))
          missingPackages.push(p);
      }
      if (missingPackages.length) {
        grok.shell.error(`Please install package${missingPackages.length === 1 ? '' : 's'} ${
          missingPackages.join(', ')} to start the tutorial`);
        this.close();
        return;
      }
    }

    const services = await grok.dapi.admin.getServiceInfos();

    for (const [service, flag] of Object.entries(this.prerequisites)) {
      if (service in Tutorial.SERVICES && flag === true) {
        const serviceAvailable = await this.checkService(Tutorial.SERVICES[service], services);
        if (!serviceAvailable)
          return;
      }
    }

    const id = tutorials.indexOf(this);

    if (this.demoTable) {
      this._t = await grok.data.getDemoTable(this.demoTable);
      grok.shell.addTableView(this._t);
    }
    this.closed = false;

    try {
      await this._run();
    } catch (error) {
      // If the tutorial was closed during execution, exit without error
      if (!this.closed) return Promise.reject(error);
    }    

    this.endSection();

    this.title('Congratulations!');
    this.describe('You have successfully completed this tutorial.');

    await grok.userSettings.add(Tutorial.DATA_STORAGE_KEY, this.name, new Date().toUTCString());
    const statusMap = await this.track?.updateStatus();

    function updateProgress(track:any) {
      const trackRoot = $(`.tutorials-track[data-name ='${track?.name}']`);
      trackRoot
        .find(`.tutorials-card[data-name='${track.tutorials[id].name}']`)
        .children('.tutorials-card-status').show();
      trackRoot
        .find('progress')
        .prop('value', (100/track.tutorials.length*(track.completed)).toFixed());
      const progressNodes = $(`.tutorials-track[data-name ='${track?.name}'] > .tutorials-track-details`).children();
      progressNodes.first().text(trackRoot.find('progress').prop('value')+'% complete');
      progressNodes.last().text(String(track.completed+' / '+track.tutorials.length));
    }

    if (statusMap && Object.values(statusMap).every((v) => v)) {
      this.root.append(ui.div([
        ui.divText(this.track?.name+'is complete!'),
        ui.bigButton('Complete', ()=>{
          updateProgress(this.track);
          this._closeAll();
          this.clearRoot();
          $('.tutorial').show();
          $('#tutorial-child-node').html('');
        }),
      ]));
    } else if (statusMap) {
      // Find the first uncompleted tutorial in the track. Give preference to the first tutorial
      // after the current one, if there are uncompleted tutorials both before and after it.
      let nextId: number | null = null;
      for (const [tutorialId, completed] of Object.entries(statusMap)) {
        const numId = +tutorialId;
        if (!completed) {
          if (nextId === null)
            nextId = numId;
          else {
            if (nextId < id && id < numId) {
              nextId = numId;
              break;
            }
          }
        }
      }
      if (nextId === null) {
        console.error('Corrupted status map.');
        nextId = id + 1;
      }
      const tutorialNode = $('#tutorial-child-node');
      const nextTutorial = tutorials[nextId];
      this.root.append(ui.divV([
        ui.divText(`Next "${nextTutorial.name}"`, {style: {margin: '5px 0'}}),
        ui.divH([
          ui.bigButton('Start', () => {
            updateProgress(this.track);
            this.clearRoot();
            tutorialNode.html('');
            tutorialNode.append(nextTutorial.root);
            nextTutorial.run();
          }),
          ui.button('Cancel', () => {
            updateProgress(this.track);
            this._closeAll();
            this.clearRoot();
            $('.tutorial').show();
            tutorialNode.html('');
          }),
        ], {style: {marginLeft: '-4px'}}),
      ]));
    }
  }

  async checkService(name: string, services?: DG.ServiceInfo[]): Promise<boolean> {
    if (!services)
      services = await grok.dapi.admin.getServiceInfos();
    const service = services.find((si) => si.name === name);
    const serviceAvailable = service == null ? false : service.enabled && service.status === 'Running';
    if (!serviceAvailable) {
      grok.shell.error(`Service "${name}" not available. Please try running this tutorial later`);
      this.close();
    }
    return serviceAvailable;
  }

  close(): void {
    this.clearRoot();
    this.closed = true;
    this.onClose.next();
    this._closeAll();
    $('.tutorial').show();
    $('#tutorial-child-node').html('');
  }

  _addHeader(): void {
    this.progressDiv.append(this.progress);
    this.progress.max = this.steps;

    this.progressSteps = ui.divText(`Step: ${this.progress.value} of ${this.steps}`);
    this.progressDiv.append(this.progressSteps);

    // const manualMode = ui.button(ui.iconFA('forward'), () => {
    //   this.manualMode = !this.manualMode;
    //   $(manualMode.firstChild).toggleClass('fal fas');
    // }, 'Self-paced mode');
    // if (this.manualMode)
    //   $(manualMode.firstChild).toggleClass('fal fas');
    const closeTutorial = ui.button(ui.iconFA('times-circle'), () => this.close());

    const linkIcon = ui.button(ui.iconFA('link'), () => {
      navigator.clipboard.writeText(this.url);
      grok.shell.info('Link copied to clipboard');
    }, `Copy the tutorial link`);

    // manualMode.style.minWidth = '30px';
    closeTutorial.style.minWidth = '30px';
    this.header.textContent = this.name;
    this.headerDiv.append(ui.divH([this.header, linkIcon], {style: {alignItems: 'center'}}));
    // this.headerDiv.append(ui.div([manualMode, closeTutorial]));
    this.headerDiv.append(closeTutorial);
  }

  title(text: string, startSection: boolean = false): void {
    const h3 = ui.h3(text);
    if (this.currentSection) {
      if (startSection) {
        this.endSection();
        this.currentSection = ui.div(ui.divH([h3], 'tutorials-section-header'));
        this.activity.append(this.currentSection);
      } else
        this.currentSection.append(h3);
    } else if (startSection) {
      this.currentSection = ui.div(ui.divH([h3], 'tutorials-section-header'));
      this.activity.append(this.currentSection);
    } else
      this.activity.append(h3);
  }

  describe(text: string): void {
    const div = ui.div();
    div.innerHTML = text;
    if (this.currentSection)
      this.currentSection.append(div);
    else
      this.activity.append(div);
    div.scrollIntoView();
  }

  endSection() {
    if (!this.currentSection) return;
    this.currentSection.classList.add('tutorials-done-section');
    const chevron = ui.iconFA('chevron-left');
    chevron.classList.add('tutorials-chevron');
    const s = this.currentSection;
    s.children[0].append(chevron);
    $(chevron).on('click', () => {
      $(chevron).toggleClass('tutorials-chevron-expanded');
      $(s).toggleClass('tutorials-done-section tutorials-done-section-expanded');
    });
    this.currentSection = undefined;
  }

  _placeHints(hint: HTMLElement | HTMLElement[]) {
    if (hint instanceof HTMLElement) {
      this.activeHints.push(hint);
      ui.hints.addHintIndicator(hint, false);
    } else if (Array.isArray(hint)) {
      this.activeHints.push(...hint);
      hint.forEach((h) => {
        if (h != null)
          ui.hints.addHintIndicator(h, false);
      });
    }
  }

  _setHintVisibility(hints: HTMLElement[], visibility: boolean) {
    hints.forEach((hint) => {
      if (hint != null)
        hint.style.visibility = visibility ? 'visible' : 'hidden';
    });
  }

  _removeHints(hint: HTMLElement | HTMLElement[]) {
    if (hint instanceof HTMLElement)
      ui.hints.remove(hint);
    else if (Array.isArray(hint)) {
      hint.forEach((h) => {
        if (h != null)
          ui.hints.remove(h);
      });
    }
  }

  async action(instructions: string, completed: Observable<any> | Promise<void>,
    hint: HTMLElement | HTMLElement[] | null = null, description: string = ''): Promise<void> {
    if (this.closed)
      return;

    this.activeHints.length = 0;
    if (hint != null)
      this._placeHints(hint);

    const view = grok.shell.v;
    const hints = Array.from(document.getElementsByClassName('ui-hint-blob')) as HTMLElement[];
    const sub = grok.events.onCurrentViewChanged.subscribe(() => {
      if (hint)
        this._setHintVisibility(hints, grok.shell.v === view);
    });

    const instructionDiv = ui.divText(instructions, 'grok-tutorial-entry-instruction');
    const descriptionDiv = ui.divText('', {classes: 'grok-tutorial-step-description', style: {
      margin: '0px 0px 0px 15px',
    }});
    const chevron = ui.iconFA('chevron-left');
    chevron.classList.add('tutorials-chevron');
    const instructionIndicator = ui.div([], 'grok-tutorial-entry-indicator');
    const entry = ui.divH([
      instructionIndicator,
      instructionDiv,
    ], 'grok-tutorial-entry');
    descriptionDiv.innerHTML = description;

    if (this.currentSection) {
      this.currentSection.append(entry);
      this.currentSection.append(descriptionDiv);
    } else {
      this.activity.append(entry);
      this.activity.append(descriptionDiv);
    }
    descriptionDiv.scrollIntoView();

    const currentStep = completed instanceof Promise ? completed : this.firstEvent(completed);
    await currentStep;

    instructionDiv.classList.add('grok-tutorial-entry-success');
    instructionIndicator.classList.add('grok-tutorial-entry-indicator-success');

    if (hint != null)
      this._removeHints(hint);
    sub.unsubscribe();

    // if (this.manualMode && manual !== false) {
    //   const nextStepIcon = ui.iconFA('forward', undefined, 'Next step');
    //   nextStepIcon.className = 'grok-icon fas fa-forward tutorials-next-step';
    //   entry.append(nextStepIcon);
    //   await this.firstEvent(fromEvent(nextStepIcon, 'click'));
    //   nextStepIcon.remove();
    // }

    $(descriptionDiv).hide();
    if (description.length != 0)
      entry.append(chevron);

    $(chevron).on('click', () => {
      $(chevron).toggleClass('tutorials-chevron-expanded');
      $(descriptionDiv).toggle();
    });
    ui.tooltip.bind(entry, description);

    this.progress.value++;
    this.progressSteps.innerHTML = '';
    this.progressSteps.append(`Step: ${this.progress.value} of ${this.steps}`);
  }

  clearRoot(): void {
    this.progress.value = 1;
    $(this.root).children().each((idx, el) => el.classList.contains('tutorials-main-header') ?
      ($(this.headerDiv).empty(), $(this.progressDiv).empty()) : $(el).empty());
  }

  firstEvent(eventStream: Observable<any>): Promise<void> {
    return new Promise<void>((resolve, reject) => {
      const eventSub = eventStream.pipe(first()).subscribe((_: any) => resolve());
      const closeSub = this.onClose.subscribe(() => {
        eventSub.unsubscribe();
        closeSub.unsubscribe();
        this._removeHints(this.activeHints);
        // eslint-disable-next-line
        reject();
      });
    }).catch((_) => console.log('Closing tutorial', this.name));
  }

  /** Closes all visual components that were added when working on tutorial, e.g., table views. */
  _closeAll(): void {
    // TODO: Take into account dialogs and other views
    if (this.t?.name) {
      grok.shell.tableView(this.t.name)?.close();
      grok.shell.closeTable(this.t);
    }
  }

  _onClose: Subject<void> = new Subject();
  get onClose() {return this._onClose;}

  private getElement(element: HTMLElement, selector: string,
    filter: ((idx: number, el: Element) => boolean) | null = null): EleLoose | null {
    const nodes = $(element).find(selector);
    return (filter ? nodes.filter(filter) : nodes)[0] ?? null;
  }

  protected get menuRoot(): HTMLElement {
    return grok.shell.windows.simpleMode ? grok.shell.v.ribbonMenu.root : grok.shell.topMenu.root;
  }

  protected getMenuItem(name: string, horizontalMenu?: boolean): HTMLElement | null {
    return this.getElement(this.menuRoot, `div.d4-menu-item.d4-menu-group${horizontalMenu ? '.d4-menu-item-horz' : ''}`,
      (idx, el) => Array.from(el.children).some((c) => c.textContent === name));
  }

  protected getSidebarHints(paneName: string, commandName: string): HTMLElement[] {
    const pane = grok.shell.sidebar.getPane(paneName);
    const command = this.getElement(pane.content, `div.d4-toggle-button[data-view=${commandName}]`) ??
      this.getElement(pane.content, 'div.d4-toggle-button', (idx, el) => el.textContent === commandName)!;
    return [pane.header, command];
  }

  /** Prompts the user to open a viewer of the specified type and returns it. */
  protected async openPlot(name: string, check: (viewer: DG.Viewer) => boolean,
    description: string = ''): Promise<DG.Viewer> {
    // TODO: Expand toolbox / accordion API coverage
    const getViewerIcon = (el: HTMLElement) => {
      const selector = name == 'filters' ? 'i.fa-filter' : `i.svg-${name.replace(' ', '-')}`;
      return this.getElement(el, selector);
    };
    const view = grok.shell.v as DG.View;
    let viewer: DG.Viewer;

    await this.action(`Open ${name}`,
      grok.events.onViewerAdded.pipe(filter((data: DG.EventData) => {
        const found = check(data.args.viewer);
        if (found)
          viewer = data.args.viewer;

        return found;
      })),
      view.type === 'TableView' ? getViewerIcon((<DG.TableView>view).toolboxPage.accordion.root) : null,
      description,
    );

    return viewer!;
  }

  /** Prompts the user to put the specified value into a dialog input. */
  protected async dlgInputAction(dlg: DG.Dialog, instructions: string, caption: string,
    value: string, description: string = '', historyHint: boolean = false, count: number = 0): Promise<void> {
    const inp = dlg.inputs.filter((input: DG.InputBase) => input.caption == caption)[count];
    if (inp == null) return;
    await this.action(instructions,
      new Observable((subscriber: any) => {
        if (inp.stringValue === value) subscriber.next(inp.stringValue);
        inp.onChanged.subscribe((inpValue) => {
          if (inpValue === value) subscriber.next(inpValue);
        });
      }),
      historyHint ? this.getElement(dlg.root, 'i.fa-history.d4-command-bar-icon') : inp.root,
      description,
    );
  }

  /** A helper method to access text inputs in a view. */
  protected async textInpAction(root: HTMLElement, instructions: string,
    caption: string, value: string, description: string = ''): Promise<void> {
    const inputRoot = this.getElement(root, 'div.ui-input-root', (idx, inp) =>
      $(inp).find('label.ui-label.ui-input-label')[0]?.textContent === caption);
    if (inputRoot == null) return;
    const input = this.getElement(inputRoot, 'input.ui-input-editor') as HTMLInputElement;
    const source = fromEvent(input, 'input').pipe(map((_) => input.value), filter((val) => val === value));
    await this.action(instructions, source, inputRoot, description);
  }

  /** A helper method to access choice inputs in a view. */
  protected async choiceInputAction(root: HTMLElement, instructions: string,
    caption: string, value: string, description: string = '') {
    let inputRoot = null;
    let select: HTMLSelectElement;
    $(root).find('.ui-input-root .ui-input-label span').each((idx, el) => {
      if (el.innerText === caption) {
        inputRoot = el.parentElement?.parentElement;
        if (inputRoot)
          select = this.getElement(inputRoot, 'select') as HTMLSelectElement;
      }
    });
    if (select! == null) return;
    const source = fromEvent(select, 'change');
    await this.action(instructions, select.value === value ?
      new Promise<void>((resolve) => resolve()) :
      source.pipe(map((_) => select.value), filter((v: string) => v === value)),
    inputRoot, description);
  };

  private async prepareColumnInpAction(root: HTMLElement, instructions: string, caption: string, columnName: string,
    description: string, inputSelector: string, valueSelector: string, intervalPeriod: number = 1000): Promise<void> {
    const columnInput = this.getElement(root, inputSelector, (idx, inp) =>
      this.getElement(inp as HTMLElement, 'label.ui-label.ui-input-label')?.textContent === caption);
    if (columnInput == null) return;
    const source = interval(intervalPeriod).pipe(
      map((_) => this.getElement(columnInput, valueSelector)?.textContent),
      filter((value) => value === columnName));
    return this.action(instructions, source, columnInput, description);
  }

  /** Prompts the user to choose a particular column in a column input with the specified caption. */
  protected async columnInpAction(root: HTMLElement, instructions: string,
    caption: string, columnName: string, description: string = '') {
    return this.prepareColumnInpAction(root, instructions, caption, columnName, description,
      'div.ui-input-root.ui-input-column', 'div.d4-column-selector-column');
  };

  /** Prompts the user to choose particular columns in a column input with the specified caption.
   * Column names should be given in the following format: `(3) AGE, HEIGHT, WEIGHT`. */
  protected async columnsInpAction(root: HTMLElement, instructions: string,
    caption: string, columnNames: string, description: string = '') {
    return this.prepareColumnInpAction(root, instructions, caption, columnNames, description,
      'div.ui-input-root.ui-input-columns', 'div.ui-input-editor > div.ui-input-column-names');
  };

  protected async buttonClickAction(root: HTMLElement, instructions: string,
    caption: string, description: string = '') {
    const btn = this.getElement(root, 'button.ui-btn', (idx, btn) => btn.textContent === caption);
    if (btn == null) return;
    const source = fromEvent(btn, 'click');
    await this.action(instructions, source, btn, description);
  };

  /** Prompts the user to open a view of the specified type, waits for it to open and returns it. */
  protected async openViewByType(instructions: string, type: string,
    hint: HTMLElement | HTMLElement[] | null = null, description: string = ''): Promise<DG.View> {
    let view: DG.View;

    // If the view was opened earlier, we find it and wait until it becomes current.
    for (const v of grok.shell.views) {
      if (v.type === type)
        view = v;
    }

    await this.action(instructions, view! == null ?
      grok.events.onViewAdded.pipe(filter((v) => {
        if (v.type === type) {
          view = v;
          return true;
        }
        return false;
      })) : grok.shell.v.type === view.type ?
        new Promise<void>((resolve, _) => resolve()) :
        grok.events.onCurrentViewChanged.pipe(filter((_) => grok.shell.v.type === view.type)),
    hint, description);

    return view!;
  }

  /** Prompts the user to open a dialog with the specified title, waits for it to open and returns it. */
  protected async openDialog(instructions: string, title: string,
    hint: HTMLElement | HTMLElement[] | null = null, description: string = ''): Promise<DG.Dialog> {
    let dialog: DG.Dialog;

    await this.action(instructions, grok.events.onDialogShown.pipe(filter((dlg) => {
      if (dlg.title === title) {
        dialog = dlg;
        return true;
      }
      return false;
    })), hint, description);

    return dialog!;
  }

  /** Prompts the user to open the "Add New Column" dialog, waits for it to open and returns it. */
  protected async openAddNCDialog(instructions: string = 'Open the "Add New Column" dialog',
    description: string = ''): Promise<DG.Dialog> {
    const addNCIcon = $('div.d4-ribbon-item').has('i.svg-add-new-column')[0];
    return await this.openDialog(instructions, 'Add New Column', addNCIcon, description);
  }

  /** Prompts the user to select a menu item in the context menu. */
  protected async contextMenuAction(instructions: string, label: string,
    hint: HTMLElement | HTMLElement[] | null = null, description: string = ''): Promise<void> {
    const commandClick = new Promise<void>((resolve, reject) => {
      const sub = grok.events.onContextMenu.subscribe((data) => {
        data.args.menu.onContextMenuItemClick.pipe(
          filter((mi) => (new DG.Menu(mi)).toString().toLowerCase() === label.toLowerCase()),
          first()).subscribe((_: any) => {
          sub.unsubscribe();
          resolve();
        });
      });
    });

    await this.action(instructions, commandClick, hint, description);
  }
}

type EleLoose = HTMLElement & Element & Node;

export interface TutorialPrerequisites {
  packages?: string[],
  jupyter?: boolean,
  grokCompute?: boolean,
  grokConnect?: boolean,
  h2o?: boolean,
}

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {delay} from '@datagrok-libraries/utils/src/test';


/** Type for {@link DemoScript} step */
export type Step = {
  name: string;
  func: () => Promise<void>;

  /** Step options: description and delay (in ms) after function ends */
  options?: {
    description?: string;
    delay?: number;
  }
};

/** Demo script class. Could be used for creating demo scripts to show the platform capabilities */
export class DemoScript {
  name: string = '';
  description: string = '';

  static currentObject: DemoScript | null = null;

  private readonly _isAutomatic: boolean = false;
  private readonly _autoStartFirstStep: boolean = false;
  private _currentStep: number = 0;
  private _isStopped: boolean = false;
  private _isCancelled: boolean = false;
  private _isStepProcessed: boolean = false;

  private _root: HTMLDivElement = ui.div([], {id: 'demo-script',
    classes: 'tutorials-root tutorials-track demo-app-script'});

  private _steps: Step[] = [];

  private _mainHeader: HTMLDivElement = ui.panel([], 'tutorials-main-header');
  private _header: HTMLHeadingElement = ui.h2('');
  private _headerDiv: HTMLDivElement = ui.divH([], 'tutorials-root-header');
  private _stopStartBtn: HTMLButtonElement = ui.button(ui.iconFA('pause'),
    async () => await this._changeStopState(), 'Play / pause');
  private _restartBtn: HTMLButtonElement = ui.button(ui.iconFA('redo'),
    async () => await this._restartScript(), 'Restart');
  private _nextStepBtn: HTMLButtonElement = ui.button(ui.iconFA('play'), async () => {
    if (!this._isStepProcessed)
      await this._nextStep();
  }, 'Next step');

  private _activity: HTMLDivElement = ui.panel([], 'tutorials-root-description');

  private _progressDiv: HTMLDivElement = ui.divV([], 'tutorials-root-progress');
  private _progress: HTMLProgressElement = ui.element('progress');
  private _progressSteps: HTMLDivElement = ui.divText('');

  private _node?: DG.DockNode;
  private _closeBtn: HTMLButtonElement = ui.button(ui.iconFA('chevron-left'), () => this._closeDock(), 'Back to demo');

  private _path?: string;
  DEMO_PATH: string = 'apps/Tutorials/Demo';

  get scriptDockNode(): DG.DockNode {
    const dockNode = Array.from(grok.shell.dockManager.rootNode.children)[0];
    return dockNode.container.containerElement === document.getElementsByClassName('panel-base splitter-container-horizontal')[0] ?
      dockNode : Array.from(dockNode.children)[0];
  }

  private _setBreadcrumbsInViewName(): void {
    const path = ['Home', 'Demo', ...((this._path ?? '').split('/'))];
    const breadcrumbs = ui.breadcrumbs(path);

    breadcrumbs.onPathClick.subscribe(async (value) => {
      const actualItem = value[value.length - 1];
      if (actualItem === breadcrumbs.path[breadcrumbs.path.length - 1])
        return;
      const tree = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Demo');
      tree.currentItem = actualItem === 'Demo' ? tree : tree.items.find((item) => item.text === actualItem)!;
    });

    if (grok.shell.v) {
      if (breadcrumbs.path.length !== 0 && breadcrumbs.path[0] === 'Home') { // integrate it to the actual breadcrumbs element
        const homeIcon = ui.iconFA('home', () => {
          grok.shell.v.close();
          grok.shell.v = DG.View.createByType(DG.VIEW_TYPE.HOME);
        });
        homeIcon.classList.add('demo-breadcrumbs-home-element');
        breadcrumbs.root.firstElementChild!.replaceWith(homeIcon);
      }
      const viewNameRoot = grok.shell.v.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
      if (viewNameRoot) {
        viewNameRoot.textContent = '';
        viewNameRoot.appendChild(breadcrumbs.root);
      }
    }
  }


  constructor(name: string, description: string, isAutomatic: boolean = false,
    options?: {autoStartFirstStep?: boolean, path?: string}) {
    this.name = name;
    this.description = description;
    this._isAutomatic = isAutomatic;
    this._autoStartFirstStep = options?.autoStartFirstStep ?? false;
    this._path = options?.path ?? '';

    this._progress.max = 0;
    this._progress.value = 1;

    DemoScript.currentObject = this;
  }

  /** Returns demo script steps */
  get steps(): Step[] {
    return this._steps;
  }

  /** Returns the amount of demo script steps */
  get stepNumber(): number {
    return this._steps.length;
  }


  /** Adds script header */
  private _addHeader(): void {
    this._createHeaderDiv();
    this._createProgressDiv();
    this._mainHeader.append(this._headerDiv, this._progressDiv);
  }

  /** Creates script header div */
  private _createHeaderDiv(): void {
    this._header.innerText = this.name;
    this._headerDiv.append(this._closeBtn);
    this._headerDiv.append(this._header);

    (this._nextStepBtn.firstChild as HTMLElement).className = 'grok-icon fas fa-play';

    this._headerDiv.append(this._isAutomatic ? this._stopStartBtn : this._nextStepBtn);
  }

  /** Creates script progress div */
  private _createProgressDiv(): void {
    this._progress.max = this.stepNumber;
    this._progressDiv.append(this._progress);
    this._progressSteps = ui.divText(`Step: ${this._progress.value} of ${this.stepNumber}`);

    this._progressDiv.append(this._progressSteps);
  }

  /** Adds description of the script */
  private _addDescription(): void {
    this._activity.append(ui.div(this.description, 'tutorials-root-description'));

    for (let i = 0; i < this.stepNumber; i++) {
      let instructionIndicator = ui.iconFA('clock');
      if (!this._isAutomatic) {
        if (i === 0) {
          instructionIndicator = ui.iconFA('play', () => this._nextStep(), 'Next step');
          instructionIndicator.className = 'grok-icon fas fa-play';
        }
      }
      const instructionDiv = ui.div(this._steps[i].name, 'grok-tutorial-entry-instruction');
      const currentStepDescription = ui.div(this._steps[i].options?.description,
        'grok-tutorial-step-description hidden');
      const entry = ui.divH([
        instructionIndicator,
        instructionDiv,
      ], 'grok-tutorial-entry');

      this._activity.append(entry, currentStepDescription);
    }
  }

  /** Initializes the root of the demo script */
  private _initRoot(): void {
    grok.shell.windows.showContextPanel = true;
    grok.shell.windows.showHelp = false;

    this._node = grok.shell.dockManager.dock(this._root, DG.DOCK_TYPE.FILL, this.scriptDockNode, this.name);
    this._node.container.containerElement.classList.add('tutorials-demo-script-container');

    this._addHeader();
    this._root.append(this._mainHeader);

    this._addDescription();
    this._root.append(this._activity);
  }

  private _setViewParams() {
    if (grok.shell.v) {
      grok.shell.v.name = this.name;
      grok.shell.v.path = `${this.DEMO_PATH}/${(this._path ?? '').replaceAll(' ', '-')}`;
      this._setBreadcrumbsInViewName();
    }
  }

  /** Processes next step */
  private async _nextStep(): Promise<void> {
    this._isStepProcessed = true;
    if (!this._isAutomatic) {
      this._nextStepBtn.classList.add('disabled');
      (this._nextStepBtn.firstChild as HTMLElement).classList.add('fa-disabled');
    }

    const entry = this._activity.getElementsByClassName('grok-tutorial-entry')[this._currentStep];
    const entryIndicator = this._activity.getElementsByClassName('grok-icon')[this._currentStep];
    const entryInstruction = this._activity.getElementsByClassName('grok-tutorial-step-description')[this._currentStep];

    entryIndicator.className = 'grok-icon far fa-spinner-third fa-spin';
    entryInstruction.classList.remove('hidden');
    entryInstruction.classList.add('visible');

    const currentStep = entry as HTMLDivElement;
    const stepDelay = this._steps[this._currentStep].options?.delay ?
      this._steps[this._currentStep].options?.delay! : 2000;

    try {
      this._setViewParams();
      await this._steps[this._currentStep].func();
      this._setViewParams();
    } catch (e) {
      console.error(e);
    }

    this._scrollTo(this._root, currentStep.offsetTop - this._mainHeader.offsetHeight);
    if (this._isAutomatic) {
      await this._countdown(entry as HTMLElement, entryIndicator as HTMLElement, stepDelay);
      await delay(stepDelay);
    }

    const newEntryIndicator = ui.iconFA('check');
    entryIndicator.replaceWith(newEntryIndicator);
    newEntryIndicator.className = 'grok-icon far fa-check';

    this._progress.value++;
    this._progressSteps.innerText = `Step: ${this._progress.value} of ${this.stepNumber}`;

    this._currentStep++;
    this._isStepProcessed = false;

    if (this._currentStep === this.stepNumber) {
      this._isAutomatic ? this._stopStartBtn.replaceWith(this._restartBtn) :
        this._nextStepBtn.replaceWith(this._restartBtn);
      return;
    }

    if (!this._isAutomatic) {
      const nextStepEntryIndicator = this._activity.getElementsByClassName('grok-icon')[this._currentStep];
      const startNextStepIcon = ui.iconFA('play', () => this._nextStep(), 'Next step');
      startNextStepIcon.className = 'grok-icon fas fa-play';
      nextStepEntryIndicator.replaceWith(startNextStepIcon);
      this._nextStepBtn.classList.remove('disabled');
      (this._nextStepBtn.firstChild as HTMLElement).classList.remove('fa-disabled');
    }
    if (grok.shell.v instanceof DG.TableView)
      await grok.data.detectSemanticTypes(grok.shell.tv.dataFrame);
  }

  /** Starts the demo script actions */
  private async _startScript(): Promise<void> {
    for (let i = this._currentStep; i < this.stepNumber; i++) {
      if (this._isStopped || this._isCancelled)
        break;

      await this._nextStep();
    }
  }

  /**
   * Scrolls to the current step
   * @param element - Current step element in root
   * @param y - y coordinate of the element
   */
  private _scrollTo(element: HTMLDivElement, y: number): void {
    element.focus();
    element.scrollTop = y;
  }

  /**
   * Adds an interactive delay indicator
   * @param element - Current step element
   * @param indicator - Current step indicator
   * @param time - Indicator animation time
   */
  private async _countdown(element: HTMLElement, indicator: HTMLElement, time: number): Promise<void> {
    const countdownDiv: HTMLDivElement = ui.div([], 'demo-script-countdown');

    indicator.classList.add('hidden');

    let countdown = time / 1000;
    const svg = this._createSVGIndicator(countdown);

    countdownDiv.append(svg);
    element.prepend(countdownDiv);

    const interval = setInterval(() => {
      countdown--;
      if (countdown === 0) {
        clearInterval(interval);
        countdownDiv.remove();

        indicator.classList.remove('hidden');
        indicator.classList.add('visible');
      }
    }, 1000);
  }

  /**
   * Creates SVG with countdown circle
   * @param countdown - countdown time
   * @returns SVG countdown indicator
   */
  private _createSVGIndicator(countdown: number): SVGSVGElement {
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    circle.setAttributeNS(null, 'cx', '7');
    circle.setAttributeNS(null, 'cy', '7');
    circle.setAttributeNS(null, 'r', '6');
    circle.setAttributeNS(null, 'style', `animation: countdown ${countdown}s linear infinite forwards`);
    svg.append(circle);

    return svg;
  }

  /** Changes the state of the demo script (stop/play) */
  private async _changeStopState(): Promise<void> {
    const icon = this._stopStartBtn.getElementsByClassName('grok-icon');
    icon[0].className = 'grok-icon fas fa-play';
    this._isStopped = !this._isStopped;

    if (!this._isStopped) {
      icon[0].className = 'grok-icon fal fa-pause';
      if (!this._isStepProcessed)
        await this._startScript();
    }
  }

  /** Restarts the script */
  private async _restartScript(): Promise<void> {
    grok.shell.dockManager.close(this._node!);
    this._clearRoot();
    this._setInitParams();
    await this.start();
  }

  /** Clears the root element */
  private _clearRoot(): void {
    this._root = ui.div([], {id: 'demo-script', classes: 'tutorials-root tutorials-track demo-app-script'});

    this._mainHeader = ui.panel([], 'tutorials-main-header');
    this._header = ui.h2('');
    this._headerDiv = ui.divH([], 'tutorials-root-header');

    this._activity = ui.panel([], 'tutorials-root-description');

    this._progressDiv = ui.divV([], 'tutorials-root-progress');
    this._progress = ui.element('progress');
    this._progressSteps = ui.divText('');

    this._progress.max = 0;
    this._progress.value = 1;
  }

  /** Sets initial parameters */
  private _setInitParams(): void {
    this._currentStep = 0;
    this._isStopped = false;
    this._isCancelled = false;

    const icon = this._stopStartBtn.getElementsByClassName('grok-icon');
    icon[0].className = 'grok-icon fal fa-pause';
    this._nextStepBtn.classList.remove('disabled');
  }

  /** Closes demo script dock */
  private _closeDock(): void {
    grok.shell.dockManager.close(this._node!);
    this.cancelScript();
  }

  /** Cancels the script */
  cancelScript(): void {
    this._isCancelled = true;
    DemoScript.currentObject = null;
  }

  /**
   * Adds a new step to script
   * @param name - Step name
   * @param func - Step function
   * @param options - Step options (description and delay after step ends)
   * @returns Returns the current demo script object
   */
  step(name: string, func: () => Promise<void>, options?: {description?: string, delay?: number}): this {
    this._steps[this.steps.length] = {
      name: name,
      func: func,
      options: options,
    };
    return this;
  }

  /** Starts the demo script */
  async start(): Promise<void> {
    this._initRoot();
    if (grok.shell.v.name === this.name) {
      grok.shell.v.close();
      this.scriptDockNode.container.setActiveChild(this._node!.container);
    }

    if (this._isAutomatic) {
      await this._startScript();
      return;
    }

    if (this._autoStartFirstStep)
      await this._nextStep();
  }
}

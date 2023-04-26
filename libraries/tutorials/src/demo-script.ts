import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {delay} from '@datagrok-libraries/utils/src/test';


/** Type for {@link DemoScript} step */
export type Step = {
  name: string;
  func: () => void;

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

  private _currentStep: number = 0;
  private _isStopped: boolean = false;
  private _isCancelled: boolean = false;

  private _root: HTMLDivElement = ui.div([], {id: 'demo-script',
    classes: 'tutorials-root tutorials-track demo-app-script'});

  private _steps: Step[] = [];

  private _mainHeader: HTMLDivElement = ui.panel([], 'tutorials-main-header');
  private _header: HTMLHeadingElement = ui.h1('');
  private _headerDiv: HTMLDivElement = ui.divH([], 'tutorials-root-header');
  private _stopStartBtn: HTMLButtonElement = ui.button(ui.iconFA('pause'),
    () => this._changeStopState(), 'Play / pause');
  private _restartBtn: HTMLButtonElement = ui.button(ui.iconFA('redo'), () => this._restartScript(), 'Restart');

  private _activity: HTMLDivElement = ui.panel([], 'tutorials-root-description');

  private _progressDiv: HTMLDivElement = ui.divV([], 'tutorials-root-progress');
  private _progress: HTMLProgressElement = ui.element('progress');
  private _progressSteps: HTMLDivElement = ui.divText('');


  constructor(name: string, description: string) {
    this.name = name;
    this.description = description;

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
    this._headerDiv.append(this._header);

    this._headerDiv.append(this._stopStartBtn);
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
      const instructionIndicator = ui.iconFA('clock');
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

    const node = grok.shell.dockManager.dock(this._root, DG.DOCK_TYPE.RIGHT, null, this.name, 0.3);
    node.container.containerElement.classList.add('tutorials-demo-script-container');

    this._addHeader();
    this._root.append(this._mainHeader);

    this._addDescription();
    this._root.append(this._activity);
  }

  /** Starts the demo script actions */
  private async _startScript(): Promise<void> {
    const entry = this._activity.getElementsByClassName('grok-tutorial-entry');
    const entryIndicators = this._activity.getElementsByClassName('grok-icon');
    const entryInstructions = this._activity.getElementsByClassName('grok-tutorial-step-description');

    for (let i = this._currentStep; i < this.stepNumber; i++) {
      if (this._isStopped || this._isCancelled)
        break;

      entryIndicators[i].className = 'grok-icon far fa-spinner-third fa-spin';
      entryInstructions[i].classList.remove('hidden');
      entryInstructions[i].classList.add('visible');

      const currentStep = entry[i] as HTMLDivElement;

      this._steps[i].func();
      this._scrollTo(this._root, currentStep.offsetTop - this._mainHeader.offsetHeight);
      await delay(this._steps[i].options?.delay ? this._steps[i].options?.delay! : 2000);

      entryIndicators[i].className = 'grok-icon far fa-check';
      this._progress.value++;
      this._progressSteps.innerText = `Step: ${this._progress.value} of ${this.stepNumber}`;

      this._currentStep++;
    }

    if (this._currentStep === this.stepNumber)
      this._stopStartBtn.replaceWith(this._restartBtn);
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

  /** Changes the state of the demo script (stop/play) */
  private _changeStopState(): void {
    const icon = this._stopStartBtn.getElementsByClassName('grok-icon');
    icon[0].className = 'grok-icon fas fa-play';
    this._isStopped = !this._isStopped;

    if (!this._isStopped) {
      icon[0].className = 'grok-icon fal fa-pause';
      this._startScript();
    }
  }

  /** Restarts the script */
  private _restartScript(): void {
    const scriptDockNode = Array.from(grok.shell.dockManager.rootNode.children)[1];
    if (scriptDockNode.container.containerElement.classList.contains('tutorials-demo-script-container')) {
      grok.shell.dockManager.close(scriptDockNode);
      grok.shell.closeAll();
    }
    this._clearRoot();
    this._setInitParams();
    this.start();
  }

  /** Clears the root element */
  private _clearRoot(): void {
    this._root = ui.div([], {id: 'demo-script', classes: 'tutorials-root tutorials-track demo-app-script'});

    this._mainHeader = ui.panel([], 'tutorials-main-header');
    this._header = ui.h1('');
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
  step(name: string, func: () => void, options?: {description?: string, delay?: number}): this {
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
    this._startScript();
  }
}

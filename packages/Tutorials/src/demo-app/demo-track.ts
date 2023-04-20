import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {delay} from '@datagrok-libraries/utils/src/test';


interface Step {
  name: string;
  func: () => void;

  options?: {
    description?: string;
    delay?: number;
  }
};

export class DemoScript {
  name: string = '';
  description: string = '';
  steps: Step[] = [];

  root: HTMLDivElement = ui.div([], 'tutorials-root tutorials-track');
  
  mainHeader: HTMLDivElement = ui.panel([], 'tutorials-main-header');
  header: HTMLHeadingElement = ui.h1('');
  headerDiv: HTMLDivElement = ui.divH([], 'tutorials-root-header');

  activity: HTMLDivElement = ui.panel([], 'tutorials-root-description');

  progressDiv: HTMLDivElement = ui.divV([], 'tutorials-root-progress');
  progress: HTMLProgressElement = ui.element('progress');
  progressSteps: HTMLDivElement = ui.divText('');


  constructor(name: string, description: string) {
    this.name = name;
    this.description = description;
    this.root.setAttribute('id', 'demo-script');

    this.progress.max = 0;
    this.progress.value = 1;
  }

  get stepNumber(): number {
    return this.steps.length;
  }

  _addHeader(): void {
    this.header.innerText = this.name;
    this.headerDiv.append(this.header);
    this.headerDiv.append(ui.button(ui.iconFA('times-circle'),()=>{grok.shell.info('stop demo')}));
    this.headerDiv.style.alignItems = 'center';

    this.progress.max = this.stepNumber;
    this.progressDiv.append(this.progress);
    this.progressSteps = ui.divText(`Step: ${this.progress.value} of ${this.stepNumber}`);

    this.progressDiv.append(this.progressSteps);

    this.mainHeader.append(this.headerDiv, this.progressDiv);
    this.mainHeader.style.paddingTop = '0px';
  }

  _addDescription(): void {
    this.activity.append(ui.div(this.description, 'tutorials-root-description'));
    for (let i = 0; i < this.stepNumber; i++) {
      const instructionIndicator = ui.iconFA('clock');
      instructionIndicator.style.color = 'var(--grey-2)';

      const instructionDiv = ui.div(this.steps[i].name, 'grok-tutorial-entry-instruction');
      instructionDiv.style.color = 'var(--grey-6)';
      const currentStepDescription = ui.div(this.steps[i].options?.description, 'grok-tutorial-step-description');
      currentStepDescription.style.paddingLeft = '22px';

      const entry = ui.divH([
        instructionIndicator,
        instructionDiv,
      ], 'grok-tutorial-entry');

      this.activity.append(entry, currentStepDescription);
    }
  }

  step(name: string, func: () => void, options?: {description?: string, delay?: number}): this {
    this.steps[this.steps.length] = {
      name: name,
      func: func,
      options: options
    };
    return this;
  }

  // TODO: add cancel button
  async start() {
    grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.RIGHT, null, this.name, 0.3);

    this._addHeader();
    this.root.append(this.mainHeader);

    this._addDescription();
    this.root.append(this.activity);

    const entryIndicators = this.activity.querySelectorAll('.tutorials-track .grok-icon');
    const entryInstructions = this.activity.getElementsByClassName('grok-tutorial-entry-instruction');

    for (let i = 0; i < this.stepNumber; i++) {
      entryIndicators[i].className = 'grok-icon far fa-spinner-third fa-spin';
      this.steps[i].func();
      entryIndicators[i].setAttribute('style','color: var(--blue-1)');

      await delay(this.steps[i].options?.delay ? this.steps[i].options?.delay! : 2000);;

      entryIndicators[i].className = 'grok-icon far fa-check';
      entryIndicators[i].setAttribute('style','color: var(--green-2)');
      //entryInstructions[i].classList.add('grok-tutorial-entry-success');

      this.progress.value++;
      this.progressSteps.innerText = `Step: ${this.progress.value} of ${this.stepNumber}`;
    }
  }
}

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

  root: HTMLDivElement = ui.div([ui.panel([], {id: 'demo-script'}),], 'tutorials-root');

  mainHeader: HTMLDivElement = ui.div([], 'tutorials-main-header');
  header: HTMLHeadingElement = ui.h1('');
  headerDiv: HTMLDivElement = ui.divH([], 'tutorials-root-header');

  activity: HTMLDivElement = ui.div([], 'tutorials-root-description');

  progressDiv: HTMLDivElement = ui.divV([], 'tutorials-root-progress');
  progress: HTMLProgressElement = ui.element('progress');
  progressSteps: HTMLDivElement = ui.divText('');


  constructor(name: string, description: string) {
    this.name = name;
    this.description = description;

    this.progress.max = 0;
    this.progress.value = 1;
  }

  get stepNumber(): number {
    return this.steps.length;
  }

  _addHeader(): void {
    this.header.innerText = this.name;
    this.headerDiv.append(ui.divH([this.header], {style: {alignItems: 'center'}}));

    this.progress.max = this.stepNumber;
    this.progressDiv.append(this.progress);
    this.progressSteps = ui.divText(`Step: ${this.progress.value} of ${this.stepNumber}`);

    this.progressDiv.append(this.progressSteps);

    this.mainHeader.append(this.headerDiv, this.progressDiv);
  }

  _addDescription(): void {
    this.activity.append(ui.h3(this.description));
    for (let i = 0; i < this.stepNumber; i++) {
      const instructionIndicator = ui.div([], 'grok-tutorial-entry-indicator');
      const instructionDiv = ui.div(this.steps[i].name, 'grok-tutorial-entry-instruction');
      const currentStepDescription = ui.div(this.steps[i].options?.description, 'grok-tutorial-step-description');

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

    const entryIndicators = this.activity.getElementsByClassName('grok-tutorial-entry-indicator');
    const entryInstructions = this.activity.getElementsByClassName('grok-tutorial-entry-instruction');

    for (let i = 0; i < this.stepNumber; i++) {
      this.steps[i].func();

      await delay(this.steps[i].options?.delay ? this.steps[i].options?.delay! : 2000);;

      entryIndicators[i].classList.add('grok-tutorial-entry-indicator-success');
      entryInstructions[i].classList.add('grok-tutorial-entry-success');

      this.progress.value++;
      this.progressSteps.innerText = `Step: ${this.progress.value} of ${this.stepNumber}`;
    }
  }
}

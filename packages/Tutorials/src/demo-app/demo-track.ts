// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';

import {delay} from '@datagrok-libraries/utils/src/test';
// import {filter} from 'rxjs/operators';


interface Step {
  func: Function;
  delay: number;
};

export class DemoScript {
  steps: Step[] = [];

  get stepNumber(): number {
    return this.steps.length;
  }

  addStep(step: Step) {
    this.steps[this.steps.length] = step;
  }

  async startScript() {
    for (let i = 0; i < this.stepNumber; i++) {
      await this.steps[i].func();
      await delay(this.steps[i].delay);
    }
  }
}

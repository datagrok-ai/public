// import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {delay} from '@datagrok-libraries/utils/src/test';
// import {filter} from 'rxjs/operators';


interface Step {
  // TODO: add step description in property panel
  description: string;
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
      const acc = ui.accordion();
      const accIcon = ui.element('i');
      accIcon.className = 'grok-icon svg-icon svg-view-layout';
      acc.addTitle(ui.span([accIcon, ui.label(`Demo test script`)]));
      acc.addPane('Test', () => ui.divText(this.steps[i].description), true);

      await this.steps[i].func();
      grok.shell.o = acc.root;
      await delay(this.steps[i].delay);
    }
  }
}

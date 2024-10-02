/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import {filter} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {interval, Observable, fromEvent} from 'rxjs';
import {DiffStudio} from '../app';
import {USE_CASES} from '../use-cases';
import {UI_TIME, TITLE} from '../ui-constants';
import {POPULATION_MODEL} from './constants';

export class DiffStudioTutorial extends Tutorial {
  get name() {
    return 'Diff Studio';
  }
  get description() {
    return 'In-browser solver of ordinary differential equations';
  }
  get steps() {return 11;}

  demoTable: string = '';
  helpUrl: string = 'https://datagrok.ai/help/compute/diff-studio';

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Diff Studio enables the modeling of processes defined by ordinary differential equations.');
    this.describe(ui.link('Learn more', this.helpUrl).outerHTML);
    this.title(`Earth's Population`);

    // 1. Run
    const browseView = grok.shell.view('Browse') as DG.BrowseView;
    const appsGroup = browseView.mainTree.getOrCreateGroup('Apps', null, false);
    const diffStudioTree = appsGroup.getOrCreateGroup('Diff Studio');
    await this.action(
      'Run Diff Studio',
      fromEvent(diffStudioTree.root, 'dblclick'),
      diffStudioTree.root,
      'Go to <b>Browse > Apps</b>, and double click <b>Diff Studio</b>',
    );

    await new Promise((resolve) => setTimeout(resolve, UI_TIME.APP_RUN_SOLVING * 2));

    const diffStudio = new DiffStudio();
    await diffStudio.runSolverApp(POPULATION_MODEL);

    await new Promise((resolve) => setTimeout(resolve, UI_TIME.APP_RUN_SOLVING));

    // 2. Play
    const finalInputAction = diffStudio.inputAction('Time', 'Final', 2200);
    await this.action(
      'Set "Final" to 2200',
      finalInputAction.promise,
      finalInputAction.root,
      'Simulate the population until 2200.',
    );

    const capacityInputAction = diffStudio.inputAction('Parameters', 'Carrying capacity', 30);
    await this.action(
      'Set "Carrying capacity" to 30',
      capacityInputAction.promise,
      capacityInputAction.root,
      'Move a slider to explore the impact of carrying capacity on Earth\'s population.',
    );

    // 3. Model
    const modelTabHeader = diffStudio.getTabHeaderRoot('Model');
    await this.action(
      'Click the Model tab',
      fromEvent(modelTabHeader, 'click'),
      modelTabHeader,
      'Go to the <b>Model</b> tab, and modify the underlying mathematical model.',
    );

    // 4. Add equation
    const equation = 'dR/dt = -P';
    const equationWithNoSpaces = equation.replaceAll(' ', '');
    await this.action(
      'Add equation',
      interval(1000).pipe(filter(() => diffStudio.getEquations().replaceAll(' ', '').includes(equationWithNoSpaces))),
      null,
      `Add the equation <b>${equation}</b> to the <b>#equations:</b>-block. It describes the <i>resource depletion (R)</i>.`,
    );

    // 5. Add initial value
    const initCondition = 'R = 3000';
    const initConditionWithNoSpaces = initCondition.replaceAll(' ', '');
    await this.action(
      'Set initial value',
      interval(1000).pipe(filter(() => diffStudio.getEquations().replaceAll(' ', '').includes(initConditionWithNoSpaces))),
      null,
      `Add <b>${initCondition}</b> to the <b>#inits:</b>-block. It defines the initial value of <b>R</b>.`,
    );

    // 6. Simulate
    const runTabHeader = diffStudio.getTabHeaderRoot('Run');
    await this.action(
      'Click the Run tab',
      fromEvent(runTabHeader, 'click'),
      runTabHeader,
      'Go to the <b>Run</b> tab, and explore the updated model.',
    );
  }
}

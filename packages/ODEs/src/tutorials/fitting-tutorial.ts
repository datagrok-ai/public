/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {interval, fromEvent} from 'rxjs';
import {DiffStudio} from '../app';
import {UI_TIME} from '../ui-constants';
import {POPULATION_MODEL_UPD} from './constants';
import {getElement, getView} from './utils';

/** Tutorial on solving differential equations */
export class FittingTutorial extends Tutorial {
  get name() {
    return 'Parameter optimization';
  }
  get description() {
    return 'Learn how to find the input conditions that lead to a specified output of the model';
  }
  get steps() {return 11;}

  demoTable: string = '';
  helpUrl: string = 'https://datagrok.ai/help/compute/function-analysis#parameter-optimization';

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Parameter optimization solves an inverse problem: finding the input conditions that lead to a specified output of the model.');
    this.describe(ui.link('Learn more', this.helpUrl).outerHTML);
    this.title('Ball flight');

    // 1. Run model catalog
    const browseView = grok.shell.view('Browse') as DG.BrowseView;
    grok.shell.v = browseView;
    const appsGroup = browseView.mainTree.getOrCreateGroup('Apps', null, false);
    appsGroup.expanded = true;
    await new Promise((resolve) => setTimeout(resolve, UI_TIME.APP_RUN_SOLVING));
    const modelCatalogNode = appsGroup.items.find((node) => node.text === 'Model Catalog');

    if (modelCatalogNode === undefined) {
      grok.shell.error('Cannot run this tutorial: the package Compute is not installed');
      return;
    }

    // 1. Model Catalog
    await this.action(
      'Open Model Catalog',
      fromEvent(modelCatalogNode.root, 'dblclick'),
      modelCatalogNode.root,
      'Go to <b>Browse > Apps</b>, and double click <b>Model Catalog</b>',
    );

    // 2. Run model
    const modelIconRoot = await getElement(grok.shell.v.root, 'span.d4-link-label[name="span-ballFlight"]');
    if (modelIconRoot === null) {
      grok.shell.error('Model Catalog run timeout exceeded');
      return;
    }

    await this.action(
      'Run the "Ball flight" model',
      fromEvent(modelIconRoot, 'dblclick'),
      modelIconRoot,
      'Double click <b>Ball flight</b>.',
    );

    // 3. Play
    const modelView = await getView('Ball flight');
    if (modelView === null) {
      grok.shell.error('Model run timeout exceeded');
      return;
    }

    await this.action(
      'Play',
      fromEvent(modelIconRoot, 'dblclick'),
      modelView.root,
      'Double click <b>Ball flight</b>.',
    );

    /*await new Promise((resolve) => setTimeout(resolve, UI_TIME.APP_RUN_SOLVING * 2));
    grok.shell.view('Template').close();

    const diffStudio = new DiffStudio();
    await diffStudio.runSolverApp(POPULATION_MODEL_UPD);

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
    const equationWithoutSpaces = equation.replaceAll(' ', '');
    await this.action(
      'Add equation',
      interval(1000).pipe(filter(() => diffStudio.getEquations().replaceAll(' ', '').includes(equationWithoutSpaces))),
      null,
      `Add the equation <b>${equation}</b> to the <b>#equations:</b>-block. It describes the <i>resource depletion (R)</i>.`,
    );

    // 5. Add initial value
    const initCondition = 'R = 3000';
    const initConditionWithoutSpaces = initCondition.replaceAll(' ', '');
    await this.action(
      'Set initial value',
      interval(1000).pipe(filter(() => diffStudio.getEquations().replaceAll(' ', '').includes(initConditionWithoutSpaces))),
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

    // 7. Back to model
    await this.action(
      'Click the Model tab',
      fromEvent(modelTabHeader, 'click'),
      modelTabHeader,
      'Diff Studio automatically generates the user interface. To enhance usability, navigate to the <b>Model</b> tab.',
    );

    // 8. Add annotation
    const annotation = '{caption: Resources; category: Initial values; min: 3000; max: 4000; units: mu}';
    const annotationWithouSpaces = annotation.replaceAll(' ', '');
    await this.action(
      'Annotate R',
      interval(1000).pipe(filter(() => diffStudio.getEquations().replaceAll(' ', '').includes(annotationWithouSpaces))),
      null,
      `Add <b>${annotation}</b> right after <b>${initCondition}</b> in the <b>#inits:</b>-block.`,
    );

    // 9. Update annotation
    const category = 'Initial values';
    const expected = `Population; category: ${category}`.replaceAll(' ', '');
    await this.action(
      'Update category of P',
      interval(1000).pipe(filter(() => diffStudio.getEquations().replaceAll(' ', '').includes(expected))),
      null,
      `Replace the current category of <b>P</b> with <b>${category}</b> in the <b>#inits:</b>-block. This will place <i>initial values</i> in the same inputs group.`,
    );

    // 10. Final checks
    await this.action(
      'Click the Run tab',
      fromEvent(runTabHeader, 'click'),
      runTabHeader,
      'Go to the <b>Run</b> tab, and check the updates.',
    );*/
  } // _run
} // DifferentialEquationsTutorial

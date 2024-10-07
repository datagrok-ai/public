/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter, map} from 'rxjs/operators';
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
    return 'Learn how to find the slider conditions that lead to a specified output of the model';
  }
  get steps() {return 11;}

  demoTable: string = '';
  helpUrl: string = 'https://datagrok.ai/help/compute/function-analysis#parameter-optimization';

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Parameter optimization solves an inverse problem: finding the slider conditions that lead to a specified output of the model.');
    this.describe(ui.link('Learn more', this.helpUrl).outerHTML);
    this.title('Ball flight');

    if (grok.shell.view('Browse') === undefined) {
      grok.shell.v = DG.View.createByType('browse');
      await new Promise((resolve) => setTimeout(resolve, UI_TIME.APP_RUN_SOLVING));
    }

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

    const formRoot = await getElement(modelView.root, 'div.d4-flex-row.ui-div.ui-form');
    const angleInputRoot = formRoot.children[5] as HTMLElement;

    await this.action(
      'Change "Angle"',
      fromEvent(angleInputRoot, 'click'),
      angleInputRoot,
      `This is a ball flight simulation. Move slider, and explore the impact of <b>Angle</b> on <b>Max distance</b>,
      <b>Max height</b> and ball\'trajectory.`,
    );

    // 4. Run fitting
    this.title('Fit max distance');

    const fitIcnRoot = document.querySelector('div.d4-ribbon-panel')
      .querySelector('i.grok-icon.fal.fa-chart-line') as HTMLElement;

    await this.action(
      'Click "Fit inputs"',
      fromEvent(fitIcnRoot, 'click'),
      fitIcnRoot,
      `How should the ball be thrown so that it flies exactly 10 meters? Run parameter optimization to answer this question.
      Click the "Fit inputs" icon on the top panel`,
    );

    // 5. Switch Velocity
    const fittingView = await getView('ballFlight - fitting');
    if (fittingView === null) {
      grok.shell.error('Fitting run timeout exceeded');
      return;
    }

    const fitFormRoot = await getElement(fittingView.root, 'div.ui-div.ui-form');
    const velocityInputRoot = fitFormRoot.children[9] as HTMLElement;
    const velocitySwitcher = velocityInputRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Switch on "Velocity"',
      fromEvent(velocitySwitcher, 'click'),
      velocitySwitcher,
      'Let\'s find the initial velocity and angle. Switch on <b>Velocity</b>.',
    );

    // 6. Switch Angle
    const angleFitInputRoot = fitFormRoot.children[12] as HTMLElement;
    const angleSwitcher = angleFitInputRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Switch on "Angle"',
      fromEvent(angleSwitcher, 'click'), //!!!
      angleSwitcher,
    );

    // 7. Target value
    const distFitInputRoot = fitFormRoot.children[16] as HTMLElement;
    const distSwitcher = distFitInputRoot.querySelector('input.ui-input-editor') as HTMLInputElement;
    const source = fromEvent(distSwitcher, 'input').pipe(map((_) => distSwitcher.value), filter((val) => val === '10'));

    await this.action(
      'Set "Max distance, m" to 10',
      source,
      distSwitcher,
      'Set target value for <b>Max distance, m</b>.',
    );

    // 8. Run
    const runIcnRoot = document.querySelector('div.d4-ribbon-panel')
      .querySelector('i.grok-icon.fal.fa-play.fas') as HTMLElement;

    await this.action(
      'Click "Run"',
      fromEvent(runIcnRoot, 'click'),
      runIcnRoot,
      `Click the <b>Run</b> icon on the top panel. This will launch fitting <b>Velocity</b> and <b>Angle</b>. 
      You can customize optimizer's settings in the <b>Using</b> block.`,
    );

    // 9. Results

    // get icon

    await this.action(
      'Click "Run"',
      fromEvent(modelIconRoot, 'dblclick'), //!!!
      modelView.root,
      `Click the <b>Run</b> icon on the top panel. This will launch fitting <b>Velocity</b> and <b>Angle</b>. 
      You can customize optimizer's settings in the <b>Using</b> block.`,
    );

    // 10. Explore

    // 11. Switch off
    this.describe('Fit ball\'s trajectory');

    // get icon

    await this.action(
      'Switch off "Angle"',
      fromEvent(modelIconRoot, 'dblclick'), //!!!
      modelView.root,
      `Click the <b>Run</b> icon on the top panel. This will launch fitting <b>Velocity</b> and <b>Angle</b>. 
      You can customize optimizer's settings in the <b>Using</b> block.`,
    );

    // 12. Check target trajectory

    // add tableview

    await this.action(
      'Click "Ball flight trajectory"',
      fromEvent(modelIconRoot, 'dblclick'), //!!!
      modelView.root,
      'How to throw a ball so that it follows a given trajectory? Click the <b>Ball flight target</b> view to explore it.',
    );

    // 13. Check target trajectory

    // add tableview

    await this.action(
      'Click "Ball flight trajectory"',
      fromEvent(modelIconRoot, 'dblclick'), //!!!
      modelView.root,
      'How to throw a ball so that it follows a given trajectory? Click the <b>Ball flight target</b> view to explore it.',
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

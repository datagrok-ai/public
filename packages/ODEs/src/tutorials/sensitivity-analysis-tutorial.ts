/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter, map} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {fromEvent} from 'rxjs';
import {UI_TIME} from '../ui-constants';
import {getElement, getView} from './utils';

/** Tutorial on sensitivity analysis */
export class SensitivityAnalysisTutorial extends Tutorial {
  get name() {
    return 'Sensitivity analysis';
  }
  get description() {
    return 'Learn how to analyze the relationship between inputs and outputs of the model';
  }
  get steps() {return 12;}

  helpUrl: string = 'https://datagrok.ai/help/compute/function-analysis#sensitivity-analysis';

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Sensitivity Analysis runs the computation multiple times with varying inputs, and analyzes the relationship between inputs and outputs.');
    this.describe(ui.link('Learn more', this.helpUrl).outerHTML);
    this.title('Model');
    this.describe('Consider ball flight simulation.');

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
      'Go to <b>Browse > Apps</b>, and double click <b>Model Catalog</b>.',
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
      'Move slider, and explore the impact of <b>Angle</b> on <b>Max distance</b>, <b>Max height</b> and ball\'s trajectory.',
    );

    // 4. Run sens.analysis
    this.title('Analysis');
    this.describe('How do Angle and other parameters affect the flight trajectory? Let\'s answer this question.');

    const senAnIcnRoot = document.querySelector('div.d4-ribbon-panel')
      .querySelector('i.grok-icon.fal.fa-analytics') as HTMLElement;

    await this.action(
      'Run sensitivity analysis',
      fromEvent(senAnIcnRoot, 'click'),
      senAnIcnRoot,
      'Click the "Run sensitivity analysis" icon on the top panel.',
    );

    // 5. Set samples
    const sensAnView = await getView('ballFlight - comparison');
    if (sensAnView === null) {
      grok.shell.error('Sensitivity analysis run timeout exceeded');
      return;
    }

    const sensAnFormRoot = await getElement(sensAnView.root, 'div.ui-div.ui-form');
    const children = sensAnFormRoot.children;

    // switch off trajectory
    const trajectoryDiv = children[23];
    const trajectorySwitcherWgt = trajectoryDiv.querySelector('div.ui-input-switch.ui-input-switch-on') as HTMLElement;
    trajectorySwitcherWgt.click();

    const samplesInputRoot = children[1] as HTMLElement;
    const samplesInputEditor = samplesInputRoot.querySelector('input.ui-input-editor') as HTMLInputElement;
    const samplesSource = fromEvent(samplesInputEditor, 'input').pipe(map((_) => samplesInputEditor.value), filter((val) => val === '100'));

    await this.action(
      'Set "Samples" to 100',
      samplesSource,
      samplesInputRoot,
      'Monte Carlo is the default method. Increase <b>Samples</b> to get more accurate results.',
    );

    // 6. Switch Angle
    const angleFitInputRoot = sensAnFormRoot.children[12] as HTMLElement;
    const angleSwitcher = angleFitInputRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Switch on "Angle"',
      fromEvent(angleSwitcher, 'click'),
      angleSwitcher,
    );

    // 7. Target value
    const distFitInputRoot = sensAnFormRoot.children[16] as HTMLElement;
    const distSwitcher = distFitInputRoot.querySelector('input.ui-input-editor') as HTMLInputElement;
    const numSource = fromEvent(distSwitcher, 'input').pipe(map((_) => distSwitcher.value), filter((val) => val === '10'));

    await this.action(
      'Set "Max distance" to 10',
      numSource,
      distSwitcher,
      'Set target value for <b>Max distance</b>.',
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

    this.describe(`Find fitting results in the grid rows. There are values of 
    <b>Velocity</b> and <b>Angle</b>, as well viewers visualizing the goodness of fit.`);

    const ballFlightTable = await grok.dapi.files.readCsv('System:AppData/DiffStudio/ball-flight-trajectory.csv');
    ballFlightTable.name = 'Ball trajectory';
    grok.shell.addTable(ballFlightTable);

    await new Promise((resolve) => setTimeout(resolve, 6000));

    // 9. Switch off Max distance
    this.title('Fit curve');
    this.describe('How to throw a ball so that it follows a given trajectory?\nYou may check the target in <b>Tables > Ball trajectory</b>.');

    const maxDistRoot = sensAnFormRoot.children[16] as HTMLElement;
    const maxDistSwitcher = maxDistRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Switch off "Max distance"',
      fromEvent(maxDistSwitcher, 'click'),
      maxDistSwitcher,
    );

    // 10. Switch on Trajectory
    const trajectoryRoot = sensAnFormRoot.children[20] as HTMLElement;
    const trajectorySwitcher = trajectoryRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Switch on "Trajectory"',
      fromEvent(trajectorySwitcher, 'click'),
      trajectorySwitcher,
    );

    // 11. Select table
    const tableInputRoot = sensAnFormRoot.querySelector('div.ui-input-choice.ui-input-table.ui-input-root');
    const tableChoiceRoot = tableInputRoot.querySelector('select.ui-input-editor.d4-invalid') as HTMLSelectElement;
    const dfSource = fromEvent(tableChoiceRoot, 'input').pipe(map((_) => tableChoiceRoot.value), filter((val) => val === 'Ball trajectory'));

    await this.action(
      'Set "Trajectory" to "Ball trajectory"',
      dfSource,
      tableChoiceRoot,
    );

    // 12. Run
    await this.action(
      'Click "Run"',
      fromEvent(runIcnRoot, 'click'),
      runIcnRoot,
    );

    await new Promise((resolve) => setTimeout(resolve, 3000));

    this.describe('Compare the simulated and target trajectories. The first row in the grid presents the best values for <b>Velocity</b> and <b>Angle</b>.');
  } // _run
} // SensitivityAnalysisTutorial

/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter, map} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {fromEvent} from 'rxjs';
import {UI_TIME} from '../ui-constants';
import {getElement, getView} from './utils';

/** Viewers description */
const viewerInfo = [
  `# Model runs
  
  Input values of **Angle** and output values of **Max Distance** and **Max Height** for each of the 100 model evaluations.`,
  `# Correlations
  
  A positive value indicates a direct correlation, while a negative value indicates an inverse correlation:
  
  * the larger **Angle**, the shorter **Max distance**
  * the larger **Angle**, the greater **Max height**
  * the larger **Max distance**, the shorter **Max height**`,
  `# Variations
  
   Multidimensional data visualization illustrates the relationship between **Angle** and the simulated values of **Max Distance** and **Max Height**.`,
  `# Graphs
  
  Line charts showing the dependence of **Max distance** and **Max height** on **Angle**.`,
];

/** Tutorial on sensitivity analysis */
export class SensitivityAnalysisTutorial extends Tutorial {
  get name() {
    return 'Sensitivity analysis';
  }
  get description() {
    return 'Learn how to analyze the relationship between inputs and outputs of the model';
  }
  get steps() {return 12;}

  demoTable: string = '';
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
      grok.shell.warning('Cannot run this tutorial: the package Compute is not installed');
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
      grok.shell.warning('Model Catalog run timeout exceeded');
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
      grok.shell.warning('Model run timeout exceeded');
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
    this.describe('How does Angle affect the flight trajectory? Let\'s answer this question.');

    const senAnIcnRoot = document.querySelector('div.d4-ribbon-panel')
      .querySelector('i.grok-icon.fal.fa-analytics') as HTMLElement;

    await this.action(
      'Run sensitivity analysis',
      fromEvent(senAnIcnRoot, 'click'),
      senAnIcnRoot,
      'Click the "Run sensitivity analysis" icon on the top panel.',
    );

    // 5. Set samples
    const sensAnView = await getView('ballFlight - comparison') as DG.TableView;
    if (sensAnView === null) {
      grok.shell.warning('Sensitivity analysis run timeout exceeded');
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
    const angleFitInputRoot = children[16] as HTMLElement;
    const angleSwitcher = angleFitInputRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Switch on "Angle"',
      fromEvent(angleSwitcher, 'click'),
      angleSwitcher,
    );

    // 7. Run sens.analysis
    const runIcnRoot = document.querySelector('div.d4-ribbon-panel').querySelector('i.grok-icon.fal.fa-play.fas') as HTMLElement;
    const runBtnRoot = children[children.length - 1].querySelector('button.ui-btn') as HTMLButtonElement;

    let resolve: (value: void | PromiseLike<void>) => void;
    const runSensAnPromise = new Promise<void>((res, rej) => resolve = res);

    runIcnRoot.addEventListener('click', (e) => resolve());
    runBtnRoot.addEventListener('click', (e) => resolve());

    await this.action(
      'Run sensitivity analysis',
      runSensAnPromise,
      [runIcnRoot, runBtnRoot],
      `Click the <b>Run</b> button or the <b>Run</b> icon on the top panel.`,
    );

    // 8. Explore viewers

    // Wait for computations complete
    const pcPlotRoot = await getElement(sensAnView.root, 'div.d4-pc-plot');
    if (pcPlotRoot === null) {
      grok.shell.warning('Sensitivity analysis run timeout exceeded');
      return;
    }

    const viewerRoots = [...sensAnView.viewers].map((v) => v.root);

    let idx = 0;
    let hint: HTMLElement;
    let msg: HTMLDivElement;
    let popup: HTMLDivElement;
    const nextBtn = ui.button('next', () => hint.click(), 'Go to the next viewer');
    const prevBtn = ui.button('prev', () => {
      idx -= 2;
      hint.click();
    }, 'Go to the previous viewer');
    const doneBtn = ui.button('done', () => hint.click(), 'Go to the next step');
    const btnsDiv = ui.divH([prevBtn, nextBtn, doneBtn]);
    btnsDiv.style.marginLeft = 'auto';
    btnsDiv.style.marginRight = '0';

    const step = () => {
      if (idx < viewerRoots.length) {
        msg = ui.divV([ui.markdown(viewerInfo[idx]), btnsDiv]);
        popup = ui.hints.addHint(viewerRoots[idx], msg, 'left');
        doneBtn.hidden = (idx < viewerRoots.length - 1);
        nextBtn.hidden = (idx === viewerRoots.length - 1);
        prevBtn.hidden = (idx < 1);
        hint = ui.hints.addHintIndicator(popup, undefined, 4000);
        hint.onclick = () => {
          popup.remove();
          ++idx;
          step();
        };
      }
    };

    step();

    await this.action(
      'Explore each viewer',
      fromEvent(doneBtn, 'click'),
      undefined,
      'Press "Next" to switch to the next viewer',
    );

    // 9. Optimization
    this.title('Optimization');
    this.describe('Solve extremum problems');

    const sliders = pcPlotRoot.querySelectorAll('div.d4-range-selector');
    const maxDistSlider = sliders[sliders.length - 2] as HTMLElement;

    await this.action(
      'Move slider',
      fromEvent(maxDistSlider, 'mousedown'),
      maxDistSlider,
      'Move the bottom slider to the top. Find the maximum of <b>Max distance</b> and the corresponding <b>Angle</b> in the filtered grid.',
    );

    await new Promise((resolve) => setTimeout(resolve, 1000));

    const clearBtn = ui.button('clear', () => hint.click(), 'Go to the next step');
    const btnDiv = ui.divH([clearBtn]);
    msg = ui.divV([
      ui.markdown(`# Maximum\n\n**Angle** providing the highest **Max distance**`),
      btnDiv,
    ]);
    popup = ui.hints.addHint(sensAnView.grid.root, msg, 'left');
    btnDiv.style.marginLeft = 'auto';
    btnDiv.style.marginRight = '0';
    hint = ui.hints.addHintIndicator(popup, undefined, 4000);
    hint.onclick = () => {
      popup.remove();
      ++idx;
      step();
    };

    await this.action(
      'Explore solution',
      fromEvent(clearBtn, 'click'),
      undefined,
      'Press "Clear" to go to the next step',
    );

    this.describe(`The Monte Carlo and Sobol method studies a model at randomly taken points.
      Set <b>Method</b> to "Grid" to get exploration at non-random points.`);
  } // _run
} // SensitivityAnalysisTutorial

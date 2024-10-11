/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter, map} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {fromEvent} from 'rxjs';
import {getElement, getView} from './utils';

/** Monte Carlo viewers description */
const monteCarloViewersInfo = [
  `# Model runs
  
  Input and output values for each of the 100 model evaluations.`,
  `# Correlations
  
  A positive value shows direct correlation; a negative, inverse:
  
  * the larger **Angle**, the greater **Max height**
  * the larger **Max distance**, the shorter **Max height**`,
  `# Variations
  
   Multidimensional visualization of the relationship between **Angle** and the simulated values of **Max Distance** and **Max Height**.`,
  `# Graphs
  
  The dependence of **Max distance** and **Max height** on **Angle**.`,
];

/** Sobol viewers description */
const sobolViewersInfo = [
  `# Scatter
  
  Here, **Max distance** specifies markers size and color. Change them to visualize **Max height**.`,
  `# Max distance
  
  The contribution of varying inputs alone: **Angle** has the highest impact.

  (to explore **Max height**, change value in the value selector)`,
  `# Max distance

  Overall impact of parameters, including their interactions with each other: **Angle** produces greater contribution than **Velocity**.`,
];

/** Describe viewers and return the Done button */
function describeViewers(viewerRoots: HTMLElement[], description: string[]): HTMLButtonElement {
  if (viewerRoots.length !== description.length)
    throw new Error('Non-equal size of viewer roots and descritions');

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
      msg = ui.divV([ui.markdown(description[idx]), btnsDiv]);
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

  return doneBtn;
}

/** Description of a single element */
function singleDescription(root: HTMLElement, description: string, tooltip: string): HTMLButtonElement {
  const clearBtn = ui.button('clear', () => hint.click(), tooltip);
  const btnDiv = ui.divH([clearBtn]);
  const msg = ui.divV([
    ui.markdown(description),
    btnDiv,
  ]);
  const popup = ui.hints.addHint(root, msg, 'left');
  btnDiv.style.marginLeft = 'auto';
  btnDiv.style.marginRight = '0';
  const hint = ui.hints.addHintIndicator(popup, undefined, 4000);
  hint.onclick = () => popup.remove();

  return clearBtn;
}

/** Tutorial on sensitivity analysis */
export class SensitivityAnalysisTutorial extends Tutorial {
  get name() {
    return 'Sensitivity analysis';
  }
  get description() {
    return 'Learn how to analyze the relationship between inputs and outputs of the model';
  }
  get steps() {return 15;}

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
      await new Promise((resolve) => setTimeout(resolve, 100));
    }

    // 1. Run model catalog
    const browseView = grok.shell.view('Browse') as DG.BrowseView;
    grok.shell.v = browseView;
    const appsGroup = browseView.mainTree.getOrCreateGroup('Apps', null, false);
    appsGroup.expanded = true;
    await new Promise((resolve) => setTimeout(resolve, 100));
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
    const angleInputRoot = formRoot!.children[5] as HTMLElement;

    await this.action(
      'Change "Angle"',
      fromEvent(angleInputRoot, 'click'),
      angleInputRoot,
      'Move slider, and explore the impact of <b>Angle</b> on <b>Max distance</b>, <b>Max height</b> and ball\'s trajectory.',
    );

    // 4. Run sens.analysis
    this.title('Analysis');
    this.describe('How does Angle affect the flight trajectory? Let\'s answer this question.');

    const senAnIcnRoot = document.querySelector('i.grok-icon.fal.fa-analytics') as HTMLElement;

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
    const children = sensAnFormRoot!.children;

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
      'Monte Carlo is the default method, and <b>Samples</b> defines the number of model evaluations. Increase <b>Samples</b> to get more accurate results.',
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
    const runIcnRoot = document.querySelector('i.grok-icon.fal.fa-play.fas') as HTMLElement;
    const runBtnRoot = children[children.length - 1].querySelector('button.ui-btn') as HTMLButtonElement;

    let resolve: (value: void | PromiseLike<void>) => void;
    let runSensAnPromise = new Promise<void>((res, rej) => resolve = res);

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

    let viewerRoots = [...sensAnView.viewers].map((v) => v.root);
    let doneBtn = describeViewers(viewerRoots, monteCarloViewersInfo);

    await this.action(
      'Explore each viewer',
      fromEvent(doneBtn, 'click'),
      undefined,
      'Click "Next" to switch to the next viewer.',
    );

    // 9. Optimization
    this.title('Optimization');
    this.describe(`Use features of ${ui.link('PC plot', 'https://datagrok.ai/help/visualize/viewers/pc-plot').outerHTML} to solve extreme problems.`);

    const sliders = pcPlotRoot.querySelectorAll('div.d4-range-selector');
    const maxDistSlider = sliders[sliders.length - 2] as HTMLElement;

    await this.action(
      'Move slider',
      fromEvent(maxDistSlider, 'mousedown'),
      maxDistSlider,
      'Move the bottom slider to the top. Find the maximum of <b>Max distance</b> and the corresponding <b>Angle</b> in the filtered grid.',
    );

    await new Promise((resolve) => setTimeout(resolve, 1000));

    // 10. Check grid
    let clearBtn = singleDescription(
      sensAnView.grid.root,
      '# Maximum\n\n**Angle** providing the highest **Max distance**',
      'Go to the next step',
    );

    await this.action(
      'Explore solution',
      fromEvent(clearBtn, 'click'),
      undefined,
      'Click "Clear" to go to the next step.',
    );

    // 11. Parameters' impact
    this.title('Parameters\' impact');
    this.describe('Explore which of the throw parameters has the most significant impact on <b>Max distance</b> and <b>Max height</b>.');

    const methodInputRoot = children[0] as HTMLElement;
    const methodChoiceRoot = methodInputRoot.querySelector('select.ui-input-editor') as HTMLSelectElement;
    const methodSource = fromEvent(methodChoiceRoot, 'input').pipe(map((_) => methodChoiceRoot.value), filter((val) => val === 'Sobol'));

    await this.action(
      'Set "Method" to "Sobol"',
      methodSource,
      methodInputRoot,
      'The Sobol method provides a quantitative assessment of the parameters\' impact.',
    );

    // 12. Switch Velocity
    const velocityFitInputRoot = children[12] as HTMLElement;
    const velocitySwitcher = velocityFitInputRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Switch on "Velocity"',
      fromEvent(velocitySwitcher, 'click'),
      velocitySwitcher,
    );

    // 13. Run sens.analysis
    runSensAnPromise = new Promise<void>((res, rej) => resolve = res);
    await this.action(
      'Run sensitivity analysis',
      runSensAnPromise,
      runIcnRoot,
    );

    // 14. Explore viewers
    // Wait for computations complete
    const barChartRoot = await getElement(sensAnView.root, 'div.d4-bar-chart');
    if (pcPlotRoot === null) {
      grok.shell.warning('Sensitivity analysis run timeout exceeded');
      return;
    }

    viewerRoots = [...sensAnView.viewers].map((v) => v.root).slice(3);
    doneBtn = describeViewers(viewerRoots, sobolViewersInfo);

    await this.action(
      'Explore each viewer',
      fromEvent(doneBtn, 'click'),
      undefined,
      'Click "Next" to switch to the next viewer',
    );

    // 15. The Grid  method note
    const methodLabelRoot = children[0].querySelector('label.ui-label') as HTMLElement;
    clearBtn = singleDescription(
      methodLabelRoot,
      '# Method\n\nThe Monte Carlo and Sobol methods study a model at randomly taken points. Select "Grid" to get exploration at non-random points.',
      'Complete this tutorial',
    );

    await this.action(
      'Click "Clear"',
      fromEvent(clearBtn, 'click'),
    );
  } // _run
} // SensitivityAnalysisTutorial

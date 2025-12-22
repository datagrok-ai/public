/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter, map} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {fromEvent} from 'rxjs';
import {getElement, getView, describeElements, singleDescription, closeWindows, PAUSE, getLegendDiv, getBallFlightModelLegend} from './utils';
import {runDescriber, Tour, DescriptionPage} from './ui-describer';
import '../../../../css/ui-describer.css';

/** Monte Carlo viewers description */
const monteCarloViewersInfo = [
  `# Model runs 🧮\n\nInput and output values for each of the 100 model evaluations.`,
  getLegendDiv(
    '# Correlations\n\nA positive value shows direct correlation; a negative, inverse:',
    [
      '🔼 the larger **Angle**, the greater **Max height**',
      '🔽 the larger **Max distance**, the shorter **Max height**',
    ],
  ),  
  `# Variations 🔀
  
   Multidimensional visualization of the relationship between **Angle** and the simulated values of **Max Distance** and **Max Height**.`,
  `# Graphs 📈
  
  The dependence of **Max distance** and **Max height** on **Angle**.`,
];

/** Sobol viewers description */
const sobolViewersInfo = [
  `# Scatter\n\nHere, **Max distance** specifies markers size and color. Change them to visualize **Max height**.`,
  getLegendDiv(
    '# Max distance\n\nThe contribution of varying inputs alone:',
    ['↗️ **Angle** has the highest impact.'],
  ),
  getLegendDiv(
    '# Max distance\n\nOverall impact of parameters, including their interactions with each other:',
    ['⚡ **Angle** produces greater contribution than **Velocity**.'],
  ),
];

/** Help links */
enum LINK {
  SENS_AN = 'https://datagrok.ai/help/compute/function-analysis#sensitivity-analysis',
  MONTE_CARLO = 'https://en.wikipedia.org/wiki/Monte_Carlo_method',
  SOBOL = 'https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis',
  DIF_STUDIO = 'https://datagrok.ai/help/compute/diff-studio',
};

/** Tutorial on sensitivity analysis */
export class SensitivityAnalysisTutorial extends Tutorial {
  get name() {
    return 'Sensitivity analysis';
  }
  get description() {
    return 'Analyze the relationship between inputs and outputs of the model';
  }
  get steps() {return 16;}

  get icon() {
    return '📊🔍';
  }

  demoTable: string = '';
  helpUrl: string = LINK.SENS_AN;

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Sensitivity Analysis runs the computation multiple times with varying inputs, and analyzes the relationship between inputs and outputs.');
    this.describe(ui.link('Learn more', this.helpUrl).outerHTML);
    this.title('Model');
    this.describe('Consider ball flight simulation.');
    closeWindows();

    // 1. Open Apps
    let browseHeader = document.querySelector('div[class="panel-titlebar disable-selection grok-browse-header"]');
    let browseIcon = document.querySelector('div[name="Browse"]') as HTMLElement;
    if (browseHeader === null)      
      browseIcon.click();

    const browsePanel = grok.shell.browsePanel;    
    const appsGroupRoot = await getElement(browsePanel.root, 'div[name="tree-Apps"]');
    if (appsGroupRoot === null) {
      grok.shell.warning('Failed to open Apps');
      return;
    }

    const appView = grok.shell.view('Apps');
    if ((appView !== null) && (appView !== undefined))
      appView.close();
    
    await this.action(
      'Open Apps',
      fromEvent(appsGroupRoot, 'click'),
      appsGroupRoot,
      'Go to <b>Browse</b> and click <b>Apps</b>',
    );

    await new Promise((resolve) => setTimeout(resolve, PAUSE));
    
    appsGroupRoot.dispatchEvent(new Event("dblclick", { bubbles: true, cancelable: true }));

    // 2. Run Model catalog
    const galleryGrid = await getElement(document,'div[class="grok-gallery-grid"]');
    if (galleryGrid === null) {
      grok.shell.warning('Failed to open apps');
      return;
    }

    browseIcon.click();

    let name = 'Model-Hub';;
    const modelCatalogIcn = await getElement(galleryGrid,`div[name="div-${name}"]`);
    name = name.replace('-',' ');

    if (modelCatalogIcn === null) {
      grok.shell.warning(`${name} not found: install the Compute package`);
      return;
    }

    await this.action(
      `Run ${name}`,
      fromEvent(modelCatalogIcn, 'dblclick'),
      modelCatalogIcn,
      `Double-click the ${name} icon`,
    );

    // 3. Run model
    const rootDiv = document.querySelector('div[id="rootDiv"]') as HTMLElement;
    const modelIconRoot = await getElement(rootDiv, 'span.d4-link-label[name="span-ballFlight"]');
      if (modelIconRoot === null) {
        grok.shell.warning(`${name} run timeout exceeded`);
        return;
    }

    await this.action(
      'Run the "Ball flight" model',
      fromEvent(modelIconRoot, 'dblclick'),
      modelIconRoot,
      'Double click <b>Ball flight</b>.',
    );

    // 4. Play
    const modelView = await getView('Ball flight');
    if (modelView === null) {
      grok.shell.warning('Model run timeout exceeded');
      return;
    }

    await new Promise((resolve) => setTimeout(resolve, 500));

    const modelRoot = await getElement(modelView.root, 'div.d4-line-chart');
    if (modelRoot === null) {
      grok.shell.warning('Sensitivity analysis run timeout exceeded');
      return;
    }

    let btnToClick = runDescriber({
      pages: [{
        root: modelRoot,
        description: getBallFlightModelLegend(),
        position: 'left',
        elements: {major: modelView.root},
      }],
      btnsText: {done: 'OK', next: '', prev: ''},
    });

    await this.action('Click "OK"', fromEvent(btnToClick, 'click'));

    // 5. Run sens.analysis
    this.title('Analysis');
    this.describe('How does Angle affect the flight trajectory? Let\'s answer this question.');

    const ribbonPannels = modelView.getRibbonPanels();
    if (ribbonPannels.length < 1) {
      grok.shell.warning('Failed to run model analysis features');
      return;      
    }
    
    const rightPanel = ribbonPannels[ribbonPannels.length - 1];
    if (rightPanel.length < 2) {
      grok.shell.warning('Failed to load model analysis features');
      return;      
    }
    
    const senAnIcnRoot = rightPanel[rightPanel.length - 1];

    await this.action(
      'Run sensitivity analysis',
      fromEvent(senAnIcnRoot, 'click'),
      senAnIcnRoot,
      'Click the "Run sensitivity analysis" icon on the top panel.',
    );

    // 6. Set samples
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

    this.describe(`${ui.link('Monte Carlo', LINK.MONTE_CARLO).outerHTML} is the default method, and the <b>Samples</b> value defines the number of model evaluations.`);

    await this.action(
      'Set "Samples" to 100',
      samplesSource,
      samplesInputRoot,
      'Increase <b>Samples</b> to get more accurate results.',
    );

    // 7. Switch Angle
    const angleFitInputRoot = children[16] as HTMLElement;
    const angleSwitcher = angleFitInputRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Toggle the "Angle" parameter',
      fromEvent(angleSwitcher, 'click'),
      angleSwitcher,
    );

    // 8. Run sens.analysis
    const runIcnRoot = document.querySelector('i.grok-icon.fal.fa-play.fas') as HTMLElement;
    const runBtnRoot = children[children.length - 1].querySelector('button.ui-btn') as HTMLButtonElement;
    if (runBtnRoot !== null)
      runBtnRoot.hidden = true;

    await this.action(
      'Run sensitivity analysis',
      fromEvent(runIcnRoot, 'click'),//runSensAnPromise,
      runIcnRoot,
      `Click the <b>Run</b> button or the <b>Run</b> icon on the top panel.`,
    );

    // 9. Explore viewers

    // Wait for computations complete
    const pcPlotRoot = await getElement(sensAnView.root, 'div.d4-pc-plot');
    if (pcPlotRoot === null) {
      grok.shell.warning('Sensitivity analysis run timeout exceeded');
      return;
    }

    let viewerRoots = [...sensAnView.viewers].map((v) => v.root);
    btnToClick = runDescriber({
      pages: viewerRoots.map((root, idx) => {
        return {
          root: root,
          description: monteCarloViewersInfo[idx],
          position: 'left',
          elements: {major: root},
        };
      }),
    });

    await this.action('Explore each viewer', fromEvent(btnToClick, 'click'));

    // 10. Optimization
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

    await new Promise((resolve) => setTimeout(resolve, 1500));

    // 11. Check grid
    btnToClick = runDescriber({
      pages: [{
        root: sensAnView.grid.root,
        description: '# Maximum ↔️\n\n**Angle** providing the highest **Max distance**',
        position: 'left',
        elements: {major: sensAnView.grid.root},
      }],
      btnsText: {done: 'ok', next: '', prev: ''},
    });

    await this.action('Explore the solution', fromEvent(btnToClick, 'click'));

    // 12. Parameters' impact

    // Align elements
    const panelRoot = sensAnView.root.querySelector('div.panel-base') as HTMLElement;

    this.title('Parameters\' impact');
    this.describe(`Explore which of the throw parameters has the most significant impact on
    <b>Max distance</b> and <b>Max height</b>. The ${ui.link('Sobol', LINK.SOBOL).outerHTML} method provides a quantitative assessment of the parameters' impact.`);

    const methodInputRoot = children[0] as HTMLElement;
    const methodChoiceRoot = methodInputRoot.querySelector('select.ui-input-editor') as HTMLSelectElement;
    const methodSource = fromEvent(methodChoiceRoot, 'input').pipe(map((_) => methodChoiceRoot.value), filter((val) => val === 'Sobol'));

    await this.action(
      'Set "Method" to "Sobol"',
      methodSource,
      methodInputRoot,
    );

    // 13. Switch Velocity
    const velocityFitInputRoot = children[12] as HTMLElement;
    const velocitySwitcher = velocityFitInputRoot.querySelector('div.ui-input-editor') as HTMLElement;

    await this.action(
      'Toggle the "Velocity" parameter',
      fromEvent(velocitySwitcher, 'click'),
      velocitySwitcher,
    );

    // 14. Run sens.analysis
    await this.action(
      'Run sensitivity analysis',
      fromEvent(runIcnRoot, 'click'),//runSensAnPromise,
      runIcnRoot,
    );

    // 15. Explore viewers
    // Wait for computations complete
    const barChartRoot = await getElement(sensAnView.root, 'div.d4-bar-chart');
    if (barChartRoot === null) {
      grok.shell.warning('Sensitivity analysis run timeout exceeded');
      return;
    }

    viewerRoots = [...sensAnView.viewers].map((v) => v.root).slice(3);
    btnToClick = runDescriber({
      pages: [...sensAnView.viewers].map((v) => v.root).slice(3).map((root, idx) => {
        return {
          root: root,
          position: 'left',
          description: sobolViewersInfo[idx],
          elements: {major: root},
        };
      }),
    });

    await this.action('Explore each viewer', fromEvent(btnToClick, 'click'));

    await new Promise((resolve) => setTimeout(resolve, 2000));

    // 16. The Grid  method note
    panelRoot.style.width = '300px';
    const methodLabelRoot = children[0] as HTMLElement;//.querySelector('label.ui-label') as HTMLElement;
    btnToClick = runDescriber({
      pages: [{
        root: methodLabelRoot,
        position: 'left',
        description: '# Method 📖\n\n * The Monte Carlo and Sobol methods study a model at randomly taken points.\n\n* Use "Grid" to get exploration at non-random points.',
        elements: {major: methodLabelRoot},
      }],
      btnsText: {done: 'clear', next: '', prev: ''},
    });

    await this.action('Click "Clear"', fromEvent(btnToClick, 'click'));

    this.describe(`Apply ${ui.link('Sensitivity Analysis', LINK.SENS_AN).outerHTML} to both ${name} and 
    ${ui.link('Diff Studio', LINK.DIF_STUDIO).outerHTML} models.`);
  } // _run   
} // SensitivityAnalysisTutorial

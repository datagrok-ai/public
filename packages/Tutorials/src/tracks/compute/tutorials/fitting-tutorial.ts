/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter, map} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {fromEvent} from 'rxjs';
import {getElement, getView, closeWindows, describeElements, PAUSE, getBallFlightModelLegend, buildToggleOverlay} from './utils';
import { runDescriber } from './ui-describer';

/** Fitting results info */
const fittingInfo = [
  `# Grid 🧮

  The first row presents the best fit:

  * the computed <b>Velocity</b> and <b>Angle</b>
  * a visualization illustrating the goodness of fit`,
  `# Max distance ↔️

  The bar chart illustrates the difference between the target and computed values of <b>Max distance</b>.`,
  `# Loss 📉

  The line chart shows the reduction of the loss function over iterations.`,
];

/** Help links */
enum LINK {
  DIF_STUDIO = 'https://datagrok.ai/help/compute/diff-studio',
  FITTING = 'https://datagrok.ai/help/compute/function-analysis#parameter-optimization',
};

/** Tutorial on parameters optimization */
export class FittingTutorial extends Tutorial {
  get name() {
    return 'Parameter optimization';
  }
  get description() {
    return 'Find the initial conditions that produce a desired model output';
  }
  get steps() {return 15;}

  get icon() {
    return '⚙️📈';
  }

  demoTable: string = '';
  helpUrl: string = LINK.FITTING;

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Parameter optimization solves an inverse problem: finding the input conditions that lead to a specified output of the model.');
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

    let name = 'Model-Hub';
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

    await this.action('Run the "Ball flight" model',
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
      grok.shell.warning('Fitting run timeout exceeded');
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

    // 5. Run fitting
    this.title('Fit scalar output');
    this.describe('How should the ball be thrown to make it fly exactly 10 meters? Let\'s find the answer.');

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

    const fitIcnRoot = rightPanel[rightPanel.length - 2];

    await this.action(
      'Click "Fit inputs"',
      fromEvent(fitIcnRoot, 'click'),
      fitIcnRoot,
      'Click the "Fit inputs" icon on the top panel.',
    );

    // 6. Switch Velocity
    const fittingView = await getView('ballFlight - fitting') as DG.TableView;
    if (fittingView === null) {
      grok.shell.warning('Fitting run timeout exceeded');
      return;
    }

    const fitFormRoot = await getElement(fittingView.root, 'div.ui-form');
    await new Promise((resolve) => setTimeout(resolve, 1000));

    const switchers = fitFormRoot!.querySelectorAll('div.ui-input-bool-switch');
    const velocitySwitcher = switchers[3] as HTMLElement;
    const velocityToggle = velocitySwitcher.querySelector('div.ui-input-editor') as HTMLElement;

    this.describe('Let\'s find the initial velocity and angle.');

    // Build 'fake' overlay to prevent wide hint for the case of narrow input form
    const velocityOverlay = buildToggleOverlay(velocityToggle);

    await this.action(
      'Toggle the "Velocity" parameter',
      fromEvent(velocityToggle, 'click'),
      velocityOverlay,
    );

    velocityOverlay.remove();

    // 7. Switch Angle    
    const angleSwitcher = switchers[4] as HTMLElement;
    const angleToggle = angleSwitcher.querySelector('div.ui-input-editor') as HTMLElement;

    // Build 'fake' overlay to prevent wide hint for the case of narrow input form
    const angleOverlay = buildToggleOverlay(angleToggle);

    await this.action(
      'Toggle the "Angle" parameter',
      fromEvent(angleToggle, 'click'),
      angleOverlay,
    );

    angleOverlay.remove();

    // 8. Target value
    const distFitInputRoot = fitFormRoot!.children[38] as HTMLElement;
    const distSwitcher = distFitInputRoot.querySelector('input.ui-input-editor') as HTMLInputElement;
    const numSource = fromEvent(distSwitcher, 'input').pipe(map((_) => distSwitcher.value), filter((val) => val === '10'));

    await this.action(
      'Set "Max distance" to 10',
      numSource,
      distSwitcher,
      'Set the target value for <b>Max distance</b>.',
    );

    // 9. Run
    const runIcnRoot = document.querySelector('i.grok-icon.fal.fa-play.fas') as HTMLElement;

    await this.action(
      'Click "Run"',
      fromEvent(runIcnRoot, 'click'),
      runIcnRoot,
      `Click the <b>Run</b> icon on the top panel to launch fitting <b>Velocity</b> and <b>Angle</b>.
      Customize optimizer's settings in the <b>Using</b> block.`,
    );

    // 10. Explore
    const barChartRoot = await getElement(fittingView.root, 'div.d4-bar-chart');
    if (barChartRoot === null) {
      grok.shell.warning('Parameter optimization run timeout exceeded');
      return;
    }

    this.describe(`To fit <b>Velocity</b> and <b>Angle</b>, Datagrok iteratively minimizes the ${ui.link('loss function', 'https://en.wikipedia.org/wiki/Loss_function').outerHTML}.`);

    await new Promise((resolve) => setTimeout(resolve, 2000));

    const grid = fittingView.grid;
    const gridRoot = grid.root;
    const gridCellRoots = gridRoot.querySelectorAll('div.d4-g-cell') as unknown as HTMLElement[];
    const rootsToDescribe = [
      gridRoot,
      gridCellRoots[0] ?? gridRoot,
      gridCellRoots[grid.table.rowCount] ?? gridRoot,
    ];

    btnToClick = runDescriber({
      pages: rootsToDescribe.map((root, idx) => {
        return {
          root: root,
          description: fittingInfo[idx],
          position: 'left',
          elements: {major: root},
        };
      }),
    });

    await this.action('Explore results', fromEvent(btnToClick, 'click'));

    const ballFlightTable = await grok.dapi.files.readCsv('System:AppData/Tutorials/ball-flight-trajectory.csv');
    ballFlightTable.name = 'Ball trajectory';
    grok.shell.addTable(ballFlightTable);

    // 11. Switch off Max distance
    this.title('Fit curve');
    this.describe('How to throw a ball so that it follows a given trajectory?\nYou may check the target in <b>Tables > Ball trajectory</b>.');

    const maxDistRoot = fitFormRoot!.children[37] as HTMLElement;
    const maxDistSwitcher = maxDistRoot.querySelector('div.ui-input-editor') as HTMLElement;

    // Build 'fake' overlay to prevent wide hint for the case of narrow input form
    const maxDistOverlay = buildToggleOverlay(maxDistSwitcher);
    
    await this.action(
      'Disable "Max distance"',
      fromEvent(maxDistSwitcher, 'click'),
      maxDistOverlay,
    );

    maxDistOverlay.remove();

    // 12. Switch on Trajectory
    const trajectoryRoot = fitFormRoot!.children[41] as HTMLElement;
    const trajectorySwitcher = trajectoryRoot.querySelector('div.ui-input-editor') as HTMLElement;

    // Build 'fake' overlay to prevent wide hint for the case of narrow input form
    const trajectoryOverlay = buildToggleOverlay(trajectorySwitcher);

    await this.action(
      'Toggle "Trajectory"',
      fromEvent(trajectorySwitcher, 'click'),
      trajectoryOverlay,
    );

    trajectoryOverlay.remove();

    // 13. Select table
    const tableInputRoot = fitFormRoot!.querySelector('div.ui-input-choice.ui-input-table.ui-input-root') as HTMLElement;
    const tableChoiceRoot = tableInputRoot.querySelector('select.ui-input-editor') as HTMLSelectElement;
    tableChoiceRoot.value = '';
    const dfSource = fromEvent(tableChoiceRoot, 'input').pipe(map((_) => tableChoiceRoot.value), filter((val) => val === 'Ball trajectory'));

    await this.action(
      'Set "Trajectory" to "Ball trajectory"',
      dfSource,
      tableInputRoot,
    );

    // 14. Run
    await this.action(
      'Click "Run"',
      fromEvent(runIcnRoot, 'click'),
      runIcnRoot,
    );

    await new Promise<void>((resolve) => {grid.onAfterDrawContent.subscribe((e) => resolve());});

    // 15. Check trajectory
    const bestVelocity = grid.table.get('Velocity', 0) as number;
    const bestAngle = grid.table.get('Angle', 0) as number;

    btnToClick = runDescriber({
      pages: [{
        root: grid.cell('[Height]', 0).element,
        description: '# The best fit 🎯\n\nTo make the ball follow the specified trajectory, use the following throw parameters:\n\n' +
          `* <b>Velocity</b> = ${bestVelocity.toFixed(2)} m/sec\n\n` + `* <b>Angle</b> = ${bestAngle.toFixed(2)} deg`,
        position: 'left',
        elements: {major: grid.cell('[Height]', 0).element}
      }],
      btnsText: {done: 'clear', next: '', prev: ''},
    });

    await this.action('Explore the fitted trajectory', fromEvent(btnToClick, 'click'));

    this.describe(`Apply ${ui.link('Parameter Optimization', LINK.FITTING).outerHTML} to both ${name} and
    ${ui.link('Diff Studio', LINK.DIF_STUDIO).outerHTML} models.`);
  } // _run
} // FittingTutorial

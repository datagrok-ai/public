/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {interval, fromEvent} from 'rxjs';
import {describeElements, singleDescription, closeWindows, getElement, getViewWithElement} from './utils';

import '../../../../css/tutorial.css';

enum DS_CONSTS {
  EDIT_RIBBON_IDX = 2,
  EDIT_TOGGLE_IDX = 2,
  REFRESH_RIBBON_IDX = 1,
  REFRESH_ICN_IDX = 0,
};

enum LINKS {
  DIF_STUDIO = 'https://datagrok.ai/help/compute/diff-studio',
  COMPS_SYNTAX = 'https://datagrok.ai/help/compute/diff-studio#model-components-and-syntax',
  ADVANCES = 'https://datagrok.ai/help/compute/diff-studio#advanced-features',
  MODELS = 'https://datagrok.ai/help/compute/models',
  LOTKA_VOLT = 'https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations',
};

/** Earth's population modeling */
export const POPULATION_MODEL = `#name: Lotka-Volterra
#equations: 
  dx/dt = alpha * x
  dy/dt = -gamma * y

#argument: t
  initial = 0 {min: 0; max: 2; caption: Start; category: Time}
  final = 15  {min: 2; max: 150; caption: Finish; category: Time}
  step = 0.1  {min: 0.1; max: 1; caption: Step; category: Time}

#inits:  
  x = 20 {min: 2; max: 40; category: Seed; caption: Prey}
  y = 2  {min: 2; max: 10; category: Seed; caption: Predator}

#parameters:
  alpha = 1.1 {min: 0.1; max: 1.5; category: Parameters} [The maximum prey per capita growth rate]
  beta = 0.4  {min: 0.1; max: 1; category: Parameters} [The effect of the presence of predators on the prey death rate]
  gamma = 1.1 {min: 0.1; max: 1.5; category: Parameters} [The predator's per capita death rate]
  delta = 0.4 {min: 0.1; max: 1; category: Parameters} [The effect of the presence of prey on the predator's growth rate]

#output:
  t {caption: Time}
  x {caption: Prey}
  y {caption: Predator}`;

/** Diff Studio UI info */
const uiInfo = [
  `# Model
  
  You can enter the model inputs here.`,
  `# Graphs
  
  The model is incomplete, leading to the following effects:

  * The prey population grows infinitely
  * The predator population steadily declines`,
  `# Dataframe

  The outcomes of the numerical simulation.`,
];

/** Diff Studio editor info */
const editorInfo = [
  `# Equations
   
  Place differential equations in this block.`,
  `# Argument
  
  Define the argument in this block: 
  * the initial value
  * the final value
  * the grid step`,
  `# Annotation
  
  Diff Studio automatically generates the user interface. To enhance usability, specify these annotations inside braces \`{}\` after each argument or parameter:

  * \`min\` and \`max\`: Creates a slider control with defined range for quick model adjustments
  * \`caption\`: Adds a descriptive label
  * \`category\`: Groups related arguments together
  
  Additionally, you can add tooltips inside brackets \`[]\`.`,
  `# Inits
  
  Define initial values here.`,
  `# Parameters
  
  If your model has parameters, specify them here.`,
  `# Output

  The computation output is a dataframe. Customize its columns here.`,
];

/** Tutorial on solving differential equations */
export class DifferentialEquationsTutorial extends Tutorial {
  get name() {
    return 'Differential equations';
  }
  get description() {
    return 'Learn how to model processes defined by ordinary differential equations with Diff Studio';
  }
  get steps() {return 14;}

  demoTable: string = '';
  helpUrl: string = LINKS.DIF_STUDIO;

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Diff Studio lets you create models declaratively using a simple, intuitive syntax.');
    this.describe(ui.link('Learn more', this.helpUrl).outerHTML);
    this.title('Models definition review');
    this.describe(`Let\'s implement the Lotka-Volterra predator-prey ${ui.link('model', LINKS.LOTKA_VOLT).outerHTML}.`);
    closeWindows();

    if (grok.shell.view('Browse') === undefined) {
      grok.shell.v = DG.View.createByType('browse');
      await new Promise((resolve) => setTimeout(resolve, 100));
    }

    // 1. Open Apps
    const browseView = grok.shell.view('Browse') as DG.BrowseView;
    grok.shell.v = browseView;
    browseView.showTree = true;

    const appsGroupRoot = await getElement(browseView.root, 'div[name="tree-Apps"]');
    if (appsGroupRoot === null) {
      grok.shell.warning('Failed to open Apps');
      return;
    }

    await this.action(
      'Open Apps',
      fromEvent(appsGroupRoot, 'click'),
      appsGroupRoot,
      'Go to <b>Browse</b> and click <b>Apps</b>',
    );

    // 2. Run Diff Studio
    const diffStudIcon = await getElement(browseView.root,'div[name="div-Diff-Studio"]');

    if (diffStudIcon === null) {
      grok.shell.warning('Diff Studio not found: install the Diff Studio package');
      return;
    }

    await this.action(
      'Run Diff Studio',
      fromEvent(diffStudIcon, 'dblclick'),
      diffStudIcon,
      'Double-click the Diff Studio icon',
    );

    const viewToClose = await getViewWithElement('div.d4-line-chart');
    if (viewToClose === null) {
      grok.shell.warning('Failed to run Diff Studio');
      return;      
    }
    viewToClose.close();  

    const openModelFunc: DG.Func = await grok.functions.eval('DiffStudio:ivpFileHandler');
    const openModelFuncCall = openModelFunc.prepare({'content': POPULATION_MODEL});
    await openModelFuncCall.call();

    await new Promise((resolve) => setTimeout(resolve, 300));
    
    // 3. Explore elements
    this.describe('We start with an incomplete model.');

    const dsView = grok.shell.v as DG.TableView;
    const dsViewRoot = dsView.root;
    let uiFormRoot = dsViewRoot.querySelector('div.ui-form') as HTMLElement;
    let inputRoots = uiFormRoot.querySelectorAll('div.ui-input.ui-input-root.ui-input-float');
    const gridRoot = dsView.grid.root;
    const lineChartRoot = dsViewRoot.querySelector('div.d4-line-chart') as HTMLElement;

    let doneBtn = describeElements([inputRoots[0] as HTMLElement, lineChartRoot, gridRoot], uiInfo);
    await this.action(
      'Explore the interface',
      fromEvent(doneBtn, 'click'),
      undefined,
      'Click "Next" to go to the next item.',
    );

    // 4. Play with inputs    
    let finishEditor = inputRoots[1].querySelector('input[class="ui-input-editor"]') as HTMLInputElement;
    await this.action(
      'Set "Finish" to 150',
      interval(100).pipe(filter(() => finishEditor.value == '150')),
      finishEditor,
      'You can also enter "150" (instead of dragging a slider) to achieve the goal.',
    );

    // 5. Go to ODEs
    this.describe('Explore the underlying mathematical model.');    
    const ribbonPanels = dsView.getRibbonPanels();    
    const editToggle = ribbonPanels[DS_CONSTS.EDIT_RIBBON_IDX][DS_CONSTS.EDIT_TOGGLE_IDX];
    await this.action(
      'Open equations editor',
      interval(100).pipe(filter(() => dsViewRoot.querySelector('div.cm-line') != null)),
      editToggle,
      'Turn on the <b>Edit</b> toggle',
    );

    // 6. Explore equations
    const editorRoot = dsViewRoot.querySelector('div.panel-base.splitter-container-horizontal') as HTMLElement;
    let ignore = await getElement(editorRoot, 'div.cm-line');
    editorRoot.style.width = '515px';
    const lineRoots = editorRoot.querySelectorAll('div[class="cm-line"]') as unknown as HTMLElement[];

    doneBtn = describeElements([
      lineRoots[0],
      lineRoots[4],
      lineRoots[5],
      lineRoots[9],
      lineRoots[13],
      lineRoots[20],
    ], editorInfo);

    await this.action(
      'Explore editor',
      fromEvent(doneBtn, 'click'),
      undefined,
      'Click "Next" to go to the next item.',
    );

    this.describe(`Diff Studio enables creating models declaratively using a simple ${ui.link('syntax', LINKS.COMPS_SYNTAX).outerHTML}.`);

    // 7. Complete 1st equation  
    this.title('Improving models');
    this.describe('Let\'s modify the model so that it takes into account the interaction between predator and prey.');

    const tutorialPanelRoot = document.querySelector('div.tutorials-root-description.ui-div');

    let equation = 'dx/dt = alpha * x - beta * x * y';
    let rawEquation = equation.replaceAll(' ', '');
    let codeDiv = ui.divV([
      ui.label('Complete the first equation to'),
      ui.divH([
        ui.div(equation, 'tutorials-code-section'),
        ui.div(ui.iconFA('copy', () => {
          grok.shell.info('Copied to clipboard');
          navigator.clipboard.writeText(equation);
        }, 'Copy to clipboard'), 'tutorials-copy-icon'),
      ]),
    ]);

    setTimeout(() => tutorialPanelRoot?.append(codeDiv), 100);
    
    await this.action(
      'Complete the first equation',
      interval(100).pipe(filter(() => lineRoots[1].textContent?.replaceAll(' ', '') == rawEquation)),
    );

    codeDiv.hidden = true;

    // 8. Complete 2nd equation    
    equation = 'dy/dt = -gamma * y + delta * x * y';
    rawEquation = equation.replaceAll(' ', '');
    codeDiv = ui.divV([
      ui.label('Complete the second equation to'),
      ui.divH([
        ui.div(equation, 'tutorials-code-section'),
        ui.div(ui.iconFA('copy', () => {
          grok.shell.info('Copied to clipboard');
          navigator.clipboard.writeText(equation);
        }, 'Copy to clipboard'), 'tutorials-copy-icon'),
      ]),
    ]);

    setTimeout(() => tutorialPanelRoot?.append(codeDiv), 100);

    await this.action(
      'Complete the second equation',
      interval(100).pipe(filter(() => lineRoots[2].textContent?.replaceAll(' ', '') == rawEquation)),
    );

    codeDiv.hidden = true;

    this.describe(`Find more capabilities at the ${ui.link('link', LINKS.ADVANCES).outerHTML}.`);

    // 9. Refresh
    const refreshIcn = ribbonPanels[DS_CONSTS.REFRESH_RIBBON_IDX][DS_CONSTS.REFRESH_ICN_IDX];
    await this.action(
      'Apply changes',
      fromEvent(refreshIcn, 'click'),
      refreshIcn,
      'Click the Refresh icon',
    );

    // 10. Check meaning
    const description = '# Updates\n\nNow, the model takes into account:' + 
      '\n* the effect of the presence of predators on the prey death rate' +
      '\n* the effect of the presence of prey on the predator\'s growth rate';
    
    let okBtn = singleDescription(lineChartRoot, description, 'Go to the next step');

    await this.action(
      'Click "OK"',
      fromEvent(okBtn, 'click'),
      undefined,
      'Check the updates. Click "OK" to go to the next step.',
    );

    // 11. Close editor
    this.title('Exploration');

    await this.action(
      'Close equations editor',
      interval(100).pipe(filter(() => dsViewRoot.querySelector('div.cm-line') == null)),
      //fromEvent(editToggle.querySelector('div.ui-input-editor')!, 'click'),
      editToggle,
      'Turn off the <b>Edit</b> toggle',
    );

    editorRoot.style.width = '241px';

    // 12. Play with inputs
    uiFormRoot = await getElement(dsViewRoot, 'div.ui-form') as HTMLElement;
    inputRoots = uiFormRoot.querySelectorAll('div.ui-input.ui-input-root.ui-input-float');

    const preyEditor = inputRoots[3].querySelector('input[class="ui-input-editor"]') as HTMLInputElement;
    await this.action(
      'Set "Prey" to 2',
      interval(100).pipe(filter(() => preyEditor.value == '2')),
      preyEditor,
      'Reduce the initial value of the prey population.',
    );    

    const deltaEditor = inputRoots[8].querySelector('input[class="ui-input-editor"]') as HTMLInputElement;
    await this.action(
      'Set "Delta" to 0.1',
      interval(100).pipe(filter(() => deltaEditor.value == '0.1')),
      deltaEditor,
      'Reduce the effect of preys on the predator\'s growth rate.',
    );

    // 13. Play with inputs
    finishEditor = inputRoots[1].querySelector('input[class="ui-input-editor"]') as HTMLInputElement;
    await this.action(
      'Set "Finish" to 150',
      interval(100).pipe(filter(() => finishEditor.value == '150')),
      finishEditor,
      'Get a simulation over a longer time period.',
    );

    this.describe(`Find useful Diff Studio ${ui.link('models', LINKS.MODELS).outerHTML}.`);
  } // _run
} // DifferentialEquationsTutorial

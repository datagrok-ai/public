/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {filter} from 'rxjs/operators';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {interval, fromEvent} from 'rxjs';
import {closeWindows, getElement, getViewWithElement, PAUSE, getTextWithSlider, simulateMouseEventsWithMove,
  DELAY,
  getLegendDiv} from './utils';

import '../../../../css/tutorial.css';
import '../../../../css/ui-describer.css';
import { DescriptionPage, getRect, runDescriber } from './ui-describer';

enum DS_CONSTS {
  EDIT_RIBBON_IDX = 2,
  EDIT_TOGGLE_IDX = 3,
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

enum SIZE_LIMITS {
  WINDOW_MIN_HEIGHT = 580,
  WINDOW_MIN_WIDTH = 930,
  TUTORIAL_PANEL_MIN_WIDTH = 190,  
};

/** Earth's population modeling */
export const POPULATION_MODEL = `#name: Lotka-Volterra
#equations: 
  dx/dt = alpha * x
  dy/dt = -gamma * y

#argument: t
  initial = 0 {min: 0; max: 2; caption: Start; category: Time}    [Simulation start.]
  final = 15  {min: 2; max: 150; caption: Finish; category: Time} [Simulation end.]
  step = 0.1  {min: 0.1; max: 1; caption: Step; category: Time}   [Grid step.]

#inits:  
  x = 20 {min: 2; max: 40; category: Seed; caption: Prey}
  y = 2  {min: 2; max: 10; category: Seed; caption: Predator}

#parameters:
  alpha = 1.1 {min: 0.1; max: 1.5; category: Parameters} [Maximum prey growth rate.]
  beta = 0.4  {min: 0.1; max: 1; category: Parameters}   [Predator impact on prey death rate.]
  gamma = 1.1 {min: 0.1; max: 1.5; category: Parameters} [Predator mortality rate.]
  delta = 0.4 {min: 0.1; max: 1; category: Parameters}   [Effect of prey on predator growth.]

#output:
  t {caption: Time}
  x {caption: Prey}
  y {caption: Predator}`;

/** Diff Studio editor info */
const editorInfo = new Map([  
  ['eqs', `# Equations
   
  Place differential equations in this block.`],
  ['arg', `# Argument
  
  Define the argument in this block: 
  * the initial value
  * the final value
  * the grid step`],
  ['annot', `# Annotation
  
  Diff Studio automatically generates the user interface. To enhance usability, specify these annotations inside braces \`{}\` after each argument or parameter:

  * \`min\` and \`max\`: Creates a slider control with defined range for quick model adjustments
  * \`caption\`: Adds a descriptive label
  * \`category\`: Groups related arguments together.`],
  ['hint', `# Hints
   
  Additionally, you can add tooltips inside brackets \`[]\`.`],
  ['inits', `# Inits
  
  Define initial values here.`],
  ['params', `# Parameters
  
  If your model has parameters, specify them here.`],
  ['out', `# Output

  The computation output is a dataframe. Customize its columns here.`],
]);

/** Tutorial on solving differential equations */
export class DifferentialEquationsTutorial extends Tutorial {
  get name() {
    return 'Differential equations';
  }
  get description() {
    return 'Learn how to model processes defined by ordinary differential equations with Diff Studio';
  }
  get steps() {return 13;}

  get icon() {
    return '📈🧮';
  }

  demoTable: string = '';
  helpUrl: string = LINKS.DIF_STUDIO;

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Diff Studio lets you create models declaratively using a simple, intuitive syntax.');
    this.describe(ui.link('Learn more', this.helpUrl).outerHTML);
    this.title('Models definition review');
    this.describe(`Let\'s implement the Lotka-Volterra predator-prey ${ui.link('model', LINKS.LOTKA_VOLT).outerHTML}.`);
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

    // 2. Run Diff Studio
    const galleryGrid = await getElement(document,'div[class="grok-gallery-grid"]');
    if (galleryGrid === null) {
      grok.shell.warning('Failed to open apps');
      return;
    }

    browseIcon.click();

    const diffStudIcon = await getElement(galleryGrid, 'div[name="div-Diff-Studio"]');
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

    await new Promise((resolve) => setTimeout(resolve, DELAY));

    // 3. Explore elements
    const dsView = grok.shell.v as DG.TableView;
    const dsViewRoot = dsView.root;
    const inputsPanel = dsViewRoot.querySelector('div[class="panel-base splitter-container-horizontal"]') as HTMLElement;
    let uiFormRoot = dsViewRoot.querySelector('div.ui-form') as HTMLElement;
    let inputRoots = uiFormRoot.querySelectorAll('div.ui-input.ui-input-root.ui-input-float');
    const gridRoot = dsView.grid.root;
    const lineChartRoot = dsViewRoot.querySelector('div.d4-line-chart') as HTMLElement;
    
    const isTutorialPanelWide = (dsViewRoot.getBoundingClientRect().x >= SIZE_LIMITS.TUTORIAL_PANEL_MIN_WIDTH);

    // Set column format for better layout
    const col = dsView.grid.col('Prey');
    if (col != null)
      col.format = '0.00E0';

    let doneBtn = runDescriber({pages: [
      { // Diff Studio view
        root: dsViewRoot.parentElement ?? dsViewRoot,
        description: '# App\n\nThis is the main view of Diff Studio.', 
        position: 'left',
        elements: {
          major: dsViewRoot.parentElement ?? dsViewRoot,
        },
      },
      { // Line charts
        root: lineChartRoot,
        description: this.getLegend(),
        position: 'left',
        elements: {
          major: lineChartRoot,
        },
      },
      { // Grid
        root: gridRoot,
        description: '# Dataframe 🧮\n\nThe outcomes of the numerical simulation.',
        position: 'top',
        elements: {
          major: gridRoot,
        },
      },
      { // Inputs panel
        root: inputsPanel,
        description: '# Inputs 3️⃣7️⃣\n\nYou can enter the model inputs here.',
        position: isTutorialPanelWide ? 'left' : 'right',
        elements: {
          major: inputsPanel,
        },
      },
      { // Single input
        root: isTutorialPanelWide ? (inputRoots[1] as HTMLElement) : (inputRoots[1].querySelector('label') ?? inputRoots[1] as HTMLElement),
        description: this.getInputDescription(), 
        position: isTutorialPanelWide ? 'left' : 'top', 
        elements: {
          major: getRect(inputRoots[1] as HTMLElement, {paddingBottom: 3}),
          extra: lineChartRoot,
        },
      },
    ]});
    await this.action(
      'Explore the interface',
      fromEvent(doneBtn, 'click'),
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
    const splitBar = dsViewRoot.querySelector('div.splitbar-vertical') as HTMLElement;    
    simulateMouseEventsWithMove(splitBar, 515, 0);    
    const lineRoots = editorRoot.querySelectorAll('div[class="cm-line"]') as unknown as HTMLElement[];

    await new Promise((resolve) => setTimeout(resolve, DELAY));

    doneBtn = runDescriber({pages: this.getCodeDescriptionPages(lineRoots, dsViewRoot)});

    await this.action(
      'Explore editor',
      fromEvent(doneBtn, 'click'),
    );

    this.describe(`Diff Studio enables creating models declaratively using a simple ${ui.link('syntax', LINKS.COMPS_SYNTAX).outerHTML}.`);

    // 7. Complete 1st equation
    this.title('Improving models');
    this.describe('Let\'s modify the model so that it takes into account the interaction between predator and prey.');

    const tutorialPanelRoot = document.querySelector('div.tutorials-root-description.ui-div');

    let equation = 'dx/dt = alpha * x - beta * x * y';
    let rawEquation = equation.replaceAll(' ', '');
    let codeDiv = ui.divV([
      ui.label('Update the first equation (🐰) to obtain'),
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
      ui.label('Update the second equation (🦊) to obtain'),
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
    const clearBtn = runDescriber({
      pages: [{
        root: lineChartRoot,
        description: this.getUpdatesDescription(),
        position: 'left',
        elements: {
          major: lineChartRoot,
        },
      }],
      btnsText: {done: 'clear', next: '', prev: ''},
    });  
    await this.action(
      'Check the updates',
      fromEvent(clearBtn, 'click'),
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

    simulateMouseEventsWithMove(splitBar, -450, 0);

    // 12. Play with input
    uiFormRoot = await getElement(dsViewRoot, 'div.ui-form') as HTMLElement;
    if (uiFormRoot == null)
      return;

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
    let finishEditor = inputRoots[1].querySelector('input[class="ui-input-editor"]') as HTMLInputElement;
    await this.action(
      'Set "Finish" to 150',
      interval(100).pipe(filter(() => finishEditor.value == '150')),
      finishEditor,
      'Get a simulation over a longer time period.',
    );

    this.describe(`Find useful Diff Studio ${ui.link('models', LINKS.MODELS).outerHTML}.`);
  } // _run

  private getLegend(): HTMLElement {
    return getLegendDiv('# Graphs\n\nThe model is incomplete, leading to the following effects:', [
      '🐰 The prey population grows infinitely.',
      '🦊 The predator population steadily declines.'      
    ]);
  }

  private getUpdatesDescription(): HTMLElement {
    return getLegendDiv('# Updates\n\nNow, the model takes into account:', [
      '🦊➠🐰 the effect of the presence of predators on the prey death rate',
      '🐰➠🦊 the effect of the presence of prey on the predator\'s growth rate'      
    ]);
  }

  private getInputDescription(): HTMLElement {
    return ui.divV([
      ui.markdown('# Interact 🖱️\n\n Set "Finish" to 150 and explore the results.'),
      getTextWithSlider('Move the slider', 'and check the updates.')
    ]);
  }

  private getCodeDescriptionPages(lineRoots: HTMLElement[], viewRoot: HTMLElement): DescriptionPage[] {
    const isWide = window.innerWidth >= SIZE_LIMITS.WINDOW_MIN_WIDTH;

    const pages: DescriptionPage[] = [
      { // Equations block
        root: lineRoots[1],
        description: editorInfo.get('eqs')!,
        position: 'left',
        elements: {
          major: isWide ? getRect(lineRoots[0], {width: 162, height: 57}) : viewRoot
        },
      },
      { // Argument block
        root: lineRoots[4],
        description: editorInfo.get('arg')!,
        position: 'left',
        elements: {
          major: isWide ? getRect(lineRoots[4], {width: 105, height: 76}) : viewRoot,
        },
      },      
    ];
    
    if (isWide) {
      pages.push(...[
        { // Annotation block
          root: lineRoots[4],
          description: editorInfo.get('annot')!,
          position: 'left',
          elements: {
            major: isWide ? getRect(lineRoots[4], {width: 477, height: 76}) : viewRoot
          },
        },
        { // Annotation with hint block
          root: lineRoots[4],
          description: editorInfo.get('hint')!,
          position: 'left',
          elements: {
            major: isWide ? getRect(lineRoots[4], {width: 627, height: 76}) : viewRoot,
          },
        },
      ]);
    }
    
    pages.push(...[
      { // Initial conditions block
        root: lineRoots[9],
        description: editorInfo.get('inits')!,
        position: 'left',
        elements: {
          major: isWide ? getRect(lineRoots[9], {width: 455, height: 57}) : viewRoot,
        },
      },
      { // Parameters block
        root: lineRoots[13],
        description: editorInfo.get('params')!,
        position: 'left',
        elements: {
          major: isWide ? getRect(lineRoots[13], {width: 412, height: 95}) : viewRoot,
        },
      }
    ]);  
    
    if (window.innerHeight >= SIZE_LIMITS.WINDOW_MIN_HEIGHT) {
      pages.push({ // Outputs block
        root: lineRoots[19],
        description: editorInfo.get('out')!,
        position: 'left',
        elements: {
          major: isWide ? getRect(lineRoots[19], {width: 185, height: 76}) : viewRoot,
        },
      });
    }

    return pages;
  } // getFullDescriptionPages
} // DifferentialEquationsTutorial

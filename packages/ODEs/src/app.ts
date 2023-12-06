// Application for solving initial value problems (IVP)

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import {basicSetup, EditorView} from "codemirror";
import {EditorState} from "@codemirror/state";
import {python} from "@codemirror/lang-python";
import {autocompletion} from "@codemirror/autocomplete";

import {DF_NAME, CONTROL_EXPR, TEMPLATES, MAX_LINE_CHART, DEMO_TEMPLATE} from './constants';
import {getIVP, getScriptLines, getScriptParams, IVP} from './scripting-tools';

/** State of IVP code editor */
enum EDITOR_STATE {
  CLEAR = 0,
  BASIC_TEMPLATE = 1,
  ADVANCED_TEMPLATE = 2,  
  FROM_FILE = 3,
  EXTENDED_TEMPLATE = 4,
  CHEM_REACT = 5,
  ROBERT = 6,
  FERM = 7,
  PKPD = 8,
  ACID_PROD = 9,
  NIMOTUZUMAB = 10,
};

/** Context help links */
enum HELP_LINKS {
  SOLVER = '/help/compute/solver.md',
  QUICK_START = '/help/explore/dim-reduction.md', // link to "Quick start" page
  FROM_SCRATCH = '/help/explore/anova.md', // link to "From scratch" page
  EXTENSIONS = '/help/visualize/viewers/3d-scatter-plot.md', // link to "Extensions" page
};

/** Get problem with respect to IVP editor state. */
function getProblem(state: EDITOR_STATE): string {
  switch (state) {
    case EDITOR_STATE.BASIC_TEMPLATE:
      return TEMPLATES.BASIC;

    case EDITOR_STATE.ADVANCED_TEMPLATE:
      return TEMPLATES.ADVANCED;
      
    case EDITOR_STATE.EXTENDED_TEMPLATE:
      return TEMPLATES.EXTENDED;
    
    case EDITOR_STATE.CHEM_REACT:
      return TEMPLATES.CHEM_REACT;

    case EDITOR_STATE.ROBERT:
      return TEMPLATES.ROBERTSON;    

    case EDITOR_STATE.FERM:
      return TEMPLATES.FERMENTATION;

    case EDITOR_STATE.PKPD:
      return TEMPLATES.PK_PD;

    case EDITOR_STATE.ACID_PROD:
      return TEMPLATES.ACID_PROD;

    case EDITOR_STATE.NIMOTUZUMAB:
      return TEMPLATES.NIMOTUZUMAB;

    default:
      return TEMPLATES.EMPTY;
  }
}

/** Completions of control expressions */
const completions = [
  {label: `${CONTROL_EXPR.NAME}: `, type: "keyword", info: "name of the problem"},  
  {label: `${CONTROL_EXPR.TAGS}: `, type: "keyword", info: "scripting tags"},
  {label: `${CONTROL_EXPR.DESCR}: `, type: "keyword", info: "descritpion of the problem"},
  {label: `${CONTROL_EXPR.DIF_EQ}:\n  `, type: "keyword", info: "block of differential equation(s) specification"},
  {label: `${CONTROL_EXPR.EXPR}:\n  `, type: "keyword", info: "block of auxiliary expressions & computations"},
  {label: `${CONTROL_EXPR.ARG}: `, type: "keyword", info: "independent variable specification"},
  {label: `${CONTROL_EXPR.INITS}:\n  `, type: "keyword", info: "initial values of the problem"},
  {label: `${CONTROL_EXPR.PARAMS}:\n  `, type: "keyword", info: "parameters of the problem"},
  {label: `${CONTROL_EXPR.CONSTS}:\n  `, type: "keyword", info: "constants definition"},
  {label: `${CONTROL_EXPR.TOL}: `, type: "keyword", info: "tolerance of numerical solution"},
  {label: `${CONTROL_EXPR.LOOP}:\n  `, type: "keyword", info: "loop feature"},
  {label: `${CONTROL_EXPR.UPDATE}:\n  `, type: "keyword", info: "update model feature"},
  {label: `${CONTROL_EXPR.OUTPUT}:\n  `, type: "keyword", info: "output specification"},
];

/** Control expressions completion utilite */
function contrCompletions(context: any) {
  let before = context.matchBefore(/[#]\w*/)
  if (!context.explicit && !before) return null
  return {
    from: before ? before.from : context.pos,
    options: completions,
    validFor: /^\w*$/
  }
}

/** Return options of line chart */
function lineChartOptions(ivp: IVP): Object {
  const outputColsCount = (ivp.outputs) ? ivp.outputs.size : ivp.inits.size;

  return {    
    showTitle: true,
    sharex: true, 
    multiAxis: outputColsCount > MAX_LINE_CHART, //true,
    multiAxisLegendPosition: "RightTop",
  };
}

/** Run solver application */
export async function runSolverApp() {

  /** Get JS-script for solving the current IVP */
  const exportToJS = () => {
    try {
      const scriptText = getScriptLines(getIVP(editorView.state.doc.toString())).join('\n');      
      const script = DG.Script.create(scriptText);
      const sView = DG.ScriptView.create(script);
      grok.shell.addView(sView);
    }
    catch (err) {
      if (err instanceof Error) {
        const b = new DG.Balloon();
        b.error(err.message);
      }
    }
  };

  /** Solve the current IVP */
  const solve = async () => {  
    try {  
      const ivp = getIVP(editorView.state.doc.toString());
      const scriptText = getScriptLines(ivp).join('\n');    
      const script = DG.Script.create(scriptText);    
      const params = getScriptParams(ivp);    
      const call = script.prepare(params);

      await call.call();
      solutionTable = call.outputs[DF_NAME];
      solverView.dataFrame = call.outputs[DF_NAME];
      solverView.name = solutionTable.name;

      if (!solutionViewer)
        solutionViewer = DG.Viewer.lineChart(solutionTable, lineChartOptions(ivp));
      
      solverView.dockManager.dock(solutionViewer, 'right');

    } catch (err) {
        if (err instanceof Error) {
          const b = new DG.Balloon();
          b.error(err.message);
        }
  }}; 
   
  let solutionTable = DG.DataFrame.create();
  let solverView = grok.shell.addTableView(solutionTable);
  let solutionViewer: DG.Viewer | null = null;
  solverView.name = 'Template';
  let div = ui.divV([]);   

  /** Code editor for IVP specifying */
  let editorView = new EditorView({
    doc: TEMPLATES.BASIC,
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    parent: div
  });

  let editorState: EDITOR_STATE = EDITOR_STATE.BASIC_TEMPLATE;
  let toShowWarning = true;

  /** Load IVP from file */
  const loadFn = async () => {
    let text = '';
    const dlg = ui.dialog('Open a file');
    const fileInp = document.createElement('input');
    fileInp.type = 'file';
    fileInp.onchange = () => {
      //@ts-ignore
      const [file] = document.querySelector("input[type=file]").files;
      const reader = new FileReader();
      reader.addEventListener("load", () => { 
        text = reader.result as string; 
        setState(EDITOR_STATE.FROM_FILE, text);
        dlg.close();
      }, false);
          
      if (file) 
        reader.readAsText(file);
    }     

    dlg.add(fileInp);        
    fileInp.click();
  };

  /** Save the current IVP to file */
  const saveFn = async () => {
    const link = document.createElement("a");
    const file = new Blob([editorView.state.doc.toString()], {type: 'text/plain'});
    link.href = URL.createObjectURL(file);
    link.download = "equations.txt";
    link.click();
    URL.revokeObjectURL(link.href);
  };

  /** Set IVP code editor state */
  const setState = (state: EDITOR_STATE, text?: string | undefined) => {
    editorState = state;
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;
    
    if (solutionViewer) {
      solutionViewer.close();
      solutionViewer = null;
    }

    const newState = EditorState.create({
      doc: text ?? getProblem(state), 
      extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    });

    editorView.setState(newState);
  };

  /** Overwrite the editor content */
  const overwrite = async (state?: EDITOR_STATE, fn?: () => Promise<void>) => {
    if (toShowWarning) {      
      const boolInput = ui.boolInput('Show this warning', true, () => toShowWarning = !toShowWarning);      
      const dlg = ui.dialog({title: 'Overwrite?', helpUrl: HELP_LINKS.SOLVER});
      solverView.append(dlg);

      dlg
        .add(ui.label('This will overwrite the current project.'))
        .add(ui.label('Do you want to go on?'))
        .add(ui.divH([boolInput.root]))
        .onCancel(() => dlg.close())
        .onOK(async () => {
          if (fn)
            await fn();
          else
            setState(state ?? EDITOR_STATE.CLEAR);          
        })
        .show();
    }
    else if (fn)
      await fn();
    else
      setState(state ?? EDITOR_STATE.CLEAR);
  };

  editorView.dom.addEventListener<"contextmenu">("contextmenu", (event) => {
    event.preventDefault();
    DG.Menu.popup()
      .item('Load...', async () => await overwrite(undefined, loadFn), undefined, {description: 'Load problem from local file'})
      .item('Save...', saveFn, undefined, {description: 'Save problem to local file'})
      .item('Clear...', async () => await overwrite(EDITOR_STATE.CLEAR), undefined, {description: 'Clear problem'})
      .separator()
      .group('Templates')
      .item('Basic...', async () => await overwrite(EDITOR_STATE.BASIC_TEMPLATE), undefined, {description: 'Open basic template'})
      .item('Advanced...', async () => await overwrite(EDITOR_STATE.ADVANCED_TEMPLATE), undefined, {description: 'Open advanced template'})
      .item('Extended...', async () => await overwrite(EDITOR_STATE.EXTENDED_TEMPLATE), undefined, {description: 'Open extended template'})
      .endGroup()
      .group('Use cases')
      .item('Chem reactions...', async () => await overwrite(EDITOR_STATE.CHEM_REACT), undefined, {description: 'Mass-action kinetics illustration'})
      .item("Robertson's model...", async () => await overwrite(EDITOR_STATE.ROBERT), undefined, {description: "Robertson's chemical reaction model"})
      .item('Fermentation...', async () => await overwrite(EDITOR_STATE.FERM), undefined, {description: 'Fermentation process simulation'})
      //.separator()
      .item('PK-PD...', async () => await overwrite(EDITOR_STATE.PKPD), undefined, {description: 'Pharmacokinetic-pharmacodynamic model'})
      .item('Acid production...', async () => await overwrite(EDITOR_STATE.ACID_PROD), undefined, {description: 'Gluconic acid production model'})
      .item('Nimotuzumab...', async () => await overwrite(EDITOR_STATE.NIMOTUZUMAB), undefined, {description: 'Nimotuzumab disposition model'})
      .endGroup()
      .show();    
  });
  
  editorView.dom.style.overflow = 'auto';
  editorView.dom.style.height = '100%';

  solverView.dockManager.dock(div, 'left');
  solverView.helpUrl = HELP_LINKS.SOLVER;

  const exportIcon = ui.iconFA('file-import', exportToJS, 'Export to JavaScript script');
  exportIcon.classList.add("fal");

  const playIcon = ui.iconFA('play', solve, 'Solve (F5)');  
  playIcon.style.color = "var(--green-2)";
  playIcon.classList.add("fas"); 
  
  solverView.root.addEventListener('keydown', async (e) => {
    if (e.key === "F5") {
      e.stopImmediatePropagation();
      e.preventDefault();
      await solve();
    }
  });

  solverView.setRibbonPanels([[exportIcon], [playIcon]]);
} // runSolverApp

/** Run solver demo application */
export async function runSolverDemoApp() {
  let solutionTable: DG.DataFrame;
  let solverView: DG.TableView;
  let solutionViewer: DG.Viewer | null = null;  
  let div = ui.div([]);   

  /** Code editor for IVP specifying */
  let editorView = new EditorView({
    doc: DEMO_TEMPLATE[0],
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    parent: div
  });

  /** Return new state of editor */
  const newState = (lines: string[]) => EditorState.create({
    doc: lines.join('\n'),
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
  });

  /** Solve the current IVP */
  const solve = async () => {  
    try {  
      const ivp = getIVP(editorView.state.doc.toString());
      const scriptText = getScriptLines(ivp).join('\n');    
      const script = DG.Script.create(scriptText);    
      const params = getScriptParams(ivp);    
      const call = script.prepare(params);

      await call.call();
      solutionTable = call.outputs[DF_NAME];
      solverView.dataFrame = call.outputs[DF_NAME];

      if (!solutionViewer)
        solutionViewer = DG.Viewer.lineChart(solutionTable, lineChartOptions(ivp));
      
      solverView.dockManager.dock(solutionViewer, DG.DOCK_TYPE.RIGHT);
    } catch (err) {
        if (err instanceof Error) {
          const b = new DG.Balloon();
          b.error(err.message);
        }
  }};

  const demoScript = new DemoScript('Solving ODEs',
    'Interactive solution of complex problems given in a declarative form');

  await demoScript
    .step('Project', async () => {      
      solutionTable = DG.DataFrame.create();
      solverView = grok.shell.addTableView(solutionTable);      
      solverView.name = 'Demo';
      editorView.dom.style.overflow = 'auto';
      editorView.dom.style.height = '100%';
      solverView.dockManager.dock(div, DG.DOCK_TYPE.LEFT);
      solverView.helpUrl = HELP_LINKS.SOLVER;
      solverView.setRibbonPanels([]);
    }, {description: 'Each project starts with the name.', delay: 0})
    .step('Equations', async () => {  
      editorView.setState(newState(DEMO_TEMPLATE.slice(0, 4)));
    }, {description: 'Declarative formulas define differential equations.', delay: 0})
    .step('Argument', async () => {  
      editorView.setState(newState(DEMO_TEMPLATE.slice(0, 9)));      
    }, {description: 'The variable, its initial & final values and modeling step.', delay: 0})
    .step('Initials', async () => {  
      editorView.setState(newState(DEMO_TEMPLATE));      
    }, {description: 'Initial values of the functions to be found.', delay: 0})
    .step('Solution', async () => { 
      grok.shell.windows.context.visible = false;
      await solve();      
    }, {description: 'Numerical results & their visualization.', delay: 0})
    .step('Try it', async () => {

      editorView.dom.addEventListener<"contextmenu">("contextmenu", (event) => {
        event.preventDefault();
        DG.Menu.popup()
          .item('Load...', loadFn, undefined, {description: 'Load problem from local file'})
          .item('Save...', saveFn, undefined, {description: 'Save problem to local file'})
          .item('Clear...', async () => setState(EDITOR_STATE.CLEAR), undefined, {description: 'Clear problem'})
          .separator()
          .group('Templates')
          .item('Basic...', async () => setState(EDITOR_STATE.BASIC_TEMPLATE), undefined, {description: 'Open basic template'})
          .item('Advanced...', async () => setState(EDITOR_STATE.ADVANCED_TEMPLATE), undefined, {description: 'Open advanced template'})
          .item('Extended...', async () => setState(EDITOR_STATE.EXTENDED_TEMPLATE), undefined, {description: 'Open extended template'})
          .endGroup()
          .group('Use cases')
          .item('Chem reactions...', async () => setState(EDITOR_STATE.CHEM_REACT), undefined, {description: 'Mass-action kinetics illustration'})
          .item("Robertson's model...", async () => setState(EDITOR_STATE.ROBERT), undefined, {description: "Robertson's chemical reaction model"})
          .item('Fermentation...', async () => setState(EDITOR_STATE.FERM), undefined, {description: 'Fermentation process simulation'})
          //.separator()
          .item('PK-PD...', async () => setState(EDITOR_STATE.PKPD), undefined, {description: 'Pharmacokinetic-pharmacodynamic model'})
          .item('Acid production...', async () => setState(EDITOR_STATE.ACID_PROD), undefined, {description: 'Gluconic acid production model'})
          .item('Nimotuzumab...', async () => setState(EDITOR_STATE.NIMOTUZUMAB), undefined, {description: 'Nimotuzumab disposition model'})
          .endGroup()
          .show();    
      });
    
      const exportIcon = ui.iconFA('file-import', exportToJS, 'Export to JavaScript script');
      exportIcon.classList.add("fal");
    
      const playIcon = ui.iconFA('play', solve, 'Solve (F5)');  
      playIcon.style.color = "var(--green-2)";
      playIcon.classList.add("fas"); 
      
      solverView.root.addEventListener('keydown', async (e) => {
        if (e.key === "F5") {
          e.stopImmediatePropagation();
          e.preventDefault();
          await solve();
        }
      });
    
      solverView.setRibbonPanels([[exportIcon], [playIcon]]);
    }, {description: 'Modify formulas and press Run button. Check use cases via context menu.', delay: 0})
    .start();
  
  /** Get JS-script for solving the current IVP */
  const exportToJS = () => {
    try {
      const scriptText = getScriptLines(getIVP(editorView.state.doc.toString())).join('\n');      
      const script = DG.Script.create(scriptText);
      const sView = DG.ScriptView.create(script);
      grok.shell.addView(sView);
    }
    catch (err) {
      if (err instanceof Error) {
        const b = new DG.Balloon();
        b.error(err.message);
      }
    }
  };  

  /** Load IVP from file */
  const loadFn = async () => {
    let text = '';
    const dlg = ui.dialog('Open a file');
    const fileInp = document.createElement('input');
    fileInp.type = 'file';
    fileInp.onchange = () => {
      //@ts-ignore
      const [file] = document.querySelector("input[type=file]").files;
      const reader = new FileReader();
      reader.addEventListener("load", () => { 
        text = reader.result as string; 
        setState(EDITOR_STATE.FROM_FILE, text);
        dlg.close();
      }, false);
          
      if (file) 
        reader.readAsText(file);
    }     

    dlg.add(fileInp);        
    fileInp.click();
  };

  /** Save the current IVP to file */
  const saveFn = async () => {
    const link = document.createElement("a");
    const file = new Blob([editorView.state.doc.toString()], {type: 'text/plain'});
    link.href = URL.createObjectURL(file);
    link.download = "equations.txt";
    link.click();
    URL.revokeObjectURL(link.href);
  };

  /** Set IVP code editor state */
  const setState = async (state: EDITOR_STATE, text?: string | undefined) => {
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;
    
    if (solutionViewer) {
      solutionViewer.close();
      solutionViewer = null;
    }

    const newState = EditorState.create({
      doc: text ?? getProblem(state), 
      extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    });

    editorView.setState(newState);

    if (state != EDITOR_STATE.CLEAR)
      await solve();
  };
} // runSolverDemoApp

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {basicSetup, EditorView} from "codemirror";
import {EditorState} from "@codemirror/state";
import {python} from "@codemirror/lang-python";
import {autocompletion} from "@codemirror/autocomplete";

import {getIVP, getScriptLines, getScriptParams, DF_NAME, CONTROL_EXPR} from './scripting-tools';

/** */
const TEMPLATE_BASIC = `#name: Basic 
#differential equations:
  dy/dt = -y + sin(t) / t

#argument: t
  initial = 0.01
  final = 15.0
  step = 0.001

#initial values:  
  y = 0`;

/** */
const TEMPLATE_ADVANCED = `This is an advanced template. Modify it. 
Use multi-line formulas if needed.
Add new equations, expressions, constants & parameters.
Edit these description lines if required.

#name: Advanced
#differential equations:
  dx/dt = E1 * y + sin(t)

  dy/dt = E2 * x - pow(t, 5)

#expressions:
  E1 = C1 * exp(-t) + P1
  E2 = C2 * cos(2 * t) + P2

#argument: t
  start = 0.0
  finish = 2.0
  step = 0.01

#initial values:
  x = 2.0
  y = 0.0

#constants:
  C1 = 1.0
  C2 = 3.0

#parameters:
  P1 = 1.0
  P2 = -1.0

#tolerance: 0.00005`;

/** */
enum EDITOR_STATE {
  CLEAR = 0,
  BASIC_TEMPLATE = 1,
  ADVANCED_TEMPLATE = 2,
};

/** */
function getSampleIVP(state: EDITOR_STATE): string {
  switch (state) {
    case EDITOR_STATE.BASIC_TEMPLATE:
      return TEMPLATE_BASIC;

    case EDITOR_STATE.ADVANCED_TEMPLATE:
      return TEMPLATE_ADVANCED;
    
    default:
      return '';
  }
};

/** Completions with of control */
//const completions = Object.values(CONTROL_EXPR).map((val) => {return {label: `${val}: `, type: "keyword"}});
const completions = [
  {label: `${CONTROL_EXPR.NAME}: `, type: "keyword", info: "name of the problem"},
  {label: `${CONTROL_EXPR.DIF_EQ}:\n  `, type: "keyword", info: "block of differential equation(s) specification"},
  {label: `${CONTROL_EXPR.EXPR}:\n  `, type: "keyword", info: "block of auxiliary expressions & computations"},
  {label: `${CONTROL_EXPR.ARG}: `, type: "keyword", info: "independent variable specification"},
  {label: `${CONTROL_EXPR.INITS}:\n  `, type: "keyword", info: "initial values of the problem"},
  {label: `${CONTROL_EXPR.PARAMS}:\n  `, type: "keyword", info: "parameters of the problem"},
  {label: `${CONTROL_EXPR.CONSTS}:\n  `, type: "keyword", info: "constants definition"},
  {label: `${CONTROL_EXPR.TOL}: `, type: "keyword", info: "tolerance of numerical solution"},
];

/** */
function contrCompletions(context: any) {
  let before = context.matchBefore(/[#]/)
  if (!context.explicit && !before) return null
  return {
    from: before ? before.from : context.pos,
    options: completions,
    validFor: /^\w*$/
  }
}

/** */
export async function runSolverApp() {
  const exportToJS = () => {    
    const scriptText = getScriptLines(getIVP(editorView.state.doc.toString())).join('\n');      
    const script = DG.Script.create(scriptText);
    const sView = DG.ScriptView.create(script);
    grok.shell.addView(sView);
  };

  const exportBtn = ui.button('export', exportToJS, 'Export to JavaScript script');

  const solve = async () => {    
    const ivp = getIVP(editorView.state.doc.toString());
    const scriptText = getScriptLines(ivp).join('\n');    
    const script = DG.Script.create(scriptText);    
    const params = getScriptParams(ivp);    
    const call = script.prepare(params);

    try {
      await call.call();
      solutionTable = call.outputs[DF_NAME];
      solverView.dataFrame = call.outputs[DF_NAME];
      solverView.name = solutionTable.name;

      if (solutionViewer)
        solutionViewer.close();  
      
      solutionViewer = DG.Viewer.lineChart(solutionTable, {
        showTitle: true,
        autoLayout: false,
        sharex: true, 
        multiAxis: true,
        multiAxisLegendPosition: "RightTop",
      });
      
      solverView.dockManager.dock(solutionViewer, 'right');
      
    } catch (err) {
        if (err instanceof Error) {
          const b = new DG.Balloon();
          b.error(err.message);
        }
  }};
  
  const solveBtn = ui.bigButton('solve', solve, 'Solve the problem');
   
  let solutionTable = DG.DataFrame.create();
  let solverView = grok.shell.addTableView(solutionTable);
  let solutionViewer: DG.Viewer | null = null;
  solverView.name = 'Template';
  let div = ui.divV([]);   

  let editorView = new EditorView({
    doc: TEMPLATE_BASIC,
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    parent: div
  });

  let editorState: EDITOR_STATE = EDITOR_STATE.BASIC_TEMPLATE;
  const editorTooltip = ui.tooltip.bind(editorView.dom, () => {
    switch (editorState) {
      case EDITOR_STATE.BASIC_TEMPLATE:      
        return 'Modify the problem. Right-click and open complex template';

      case EDITOR_STATE.CLEAR:        
        return 'Define initial value problem here. Right-click and select a template';
      
      default:        
        return 'Modify the problem. Add new equations, expressions, constants & parameters. Change name, argument & tolerance if needed.';
  }});

  const openFn = () => {};

  const saveFn = () => {};

  const setState = (state: EDITOR_STATE) => {
    editorState = state;
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;
    
    if (solutionViewer) {
      solutionViewer.close();
      solutionViewer = null;
    }

    const newState = EditorState.create({
      doc: getSampleIVP(state), 
      extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    });

    editorView.setState(newState);
  };

  editorView.dom.addEventListener<"contextmenu">("contextmenu", (event) => {
    event.preventDefault();
    DG.Menu.popup()
      .item('Load...', openFn, undefined, {description: 'Load problem from local file'})
      .item('Save...', saveFn, undefined, {description: 'Save problem to local file'})
      .item('Clear...', () => setState(EDITOR_STATE.CLEAR), undefined, {description: 'Clear problem'})
      .separator()
      .item('Basic...', () => setState(EDITOR_STATE.BASIC_TEMPLATE), undefined, {description: 'Open basic template'})
      .item('Advanced...', () => setState(EDITOR_STATE.ADVANCED_TEMPLATE), undefined, {description: 'Open advanced template'})
      .show();    
  });
  
  editorView.dom.style.overflow = 'auto';

  solverView.dockManager.dock(div, 'left');
  
  div.appendChild(ui.h3(' '));

  const buttons = ui.buttonsInput([exportBtn, solveBtn]);
  //buttons.style.alignSelf = 'end';  

  div.appendChild(buttons);

  //solve();
} // runSolverApp
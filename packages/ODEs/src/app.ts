import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {basicSetup, EditorView} from "codemirror";
import {python} from "@codemirror/lang-python";
import {autocompletion} from "@codemirror/autocomplete";

import {getIVP, getScriptLines, getScriptParams, DF_NAME, TEMPLATE, CONTROL_EXPR} from './scripting-tools';

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
  const exportFn = () => {    
    const scriptText = getScriptLines(getIVP(newView.state.doc.toString())).join('\n');      
    const script = DG.Script.create(scriptText);
    const sView = DG.ScriptView.create(script);
    grok.shell.addView(sView);
  };

  const exportBtn = ui.button('export', exportFn, 'Export to JavaScript script');

  const solveFn = async () => {    
    const ivp = getIVP(newView.state.doc.toString());
    const scriptText = getScriptLines(ivp).join('\n');    
    const script = DG.Script.create(scriptText);    
    const params = getScriptParams(ivp);    
    const call = script.prepare(params);

    try {
      await call.call();
      df = call.outputs[DF_NAME];
      view.dataFrame = call.outputs[DF_NAME];
      view.name = df.name;
  
      if (!viewer) {
        viewer = DG.Viewer.lineChart(df, {
          showTitle: true,
          autoLayout: false,
          sharex: true, 
          multiAxis: true,
          multiAxisLegendPosition: "RightTop",
        });
        view.dockManager.dock(viewer, 'right');
      }
    } catch (err) {
        if (err instanceof Error) {
          const b = new DG.Balloon();
          b.error(err.message);
        }
  }};
  
  const solveBtn = ui.bigButton('solve', solveFn, 'Solve the problem');
   
  let df = DG.DataFrame.create();
  let view = grok.shell.addTableView(df);
  let viewer: DG.Viewer | null = null;
  view.name = 'Template';
  let div = ui.divV([]);

  let newView = new EditorView({
    doc: TEMPLATE,
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    parent: div
  });

  let toolTipIdx = 0;
  const editorTooltip = ui.tooltip.bind(newView.dom, () => {
    switch (toolTipIdx) {
      case 0:
        ++toolTipIdx;
        return 'Modify name & argument if needed';
      
      default:
        ++toolTipIdx;
        return 'Start typing with "#" to add the desired required blocks';
  }});
  
  newView.dom.style.overflow = 'auto';

  view.dockManager.dock(div, 'left');
  
  div.appendChild(ui.h3(' '));

  const buttons = ui.buttonsInput([exportBtn, solveBtn]);
  //buttons.style.alignSelf = 'end';

  div.appendChild(buttons);  
} // runSolverApp
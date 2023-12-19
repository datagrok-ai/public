// Application for solving initial value problems (IVP) & demo app

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import {basicSetup, EditorView} from "codemirror";
import {EditorState} from "@codemirror/state";
import {python} from "@codemirror/lang-python";
import {autocompletion} from "@codemirror/autocomplete";

import {DF_NAME, CONTROL_EXPR, MAX_LINE_CHART} from './constants';
import {TEMPLATES, DEMO_TEMPLATE} from './templates';
import { USE_CASES } from './use-cases';
import {HINT, TITLE, LINK, HOT_KEY, ERROR_MSG, INFO, WARNING, MISC} from './ui-constants';
import {getIVP, getScriptLines, getScriptParams} from './scripting-tools';

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
      return USE_CASES.CHEM_REACT;

    case EDITOR_STATE.ROBERT:
      return USE_CASES.ROBERTSON;    

    case EDITOR_STATE.FERM:
      return USE_CASES.FERMENTATION;

    case EDITOR_STATE.PKPD:
      return USE_CASES.PK_PD;

    case EDITOR_STATE.ACID_PROD:
      return USE_CASES.ACID_PROD;

    case EDITOR_STATE.NIMOTUZUMAB:
      return USE_CASES.NIMOTUZUMAB;

    default:
      return TEMPLATES.EMPTY;
  }
} // getProblem

/** Completions of control expressions */
const completions = [
  {label: `${CONTROL_EXPR.NAME}: `, type: "keyword", info: INFO.NAME},  
  {label: `${CONTROL_EXPR.TAGS}: `, type: "keyword", info: INFO.TAGS},
  {label: `${CONTROL_EXPR.DESCR}: `, type: "keyword", info: INFO.DESCR},
  {label: `${CONTROL_EXPR.DIF_EQ}:\n  `, type: "keyword", info: INFO.DIF_EQ},
  {label: `${CONTROL_EXPR.EXPR}:\n  `, type: "keyword", info: INFO.EXPR},
  {label: `${CONTROL_EXPR.ARG}: `, type: "keyword", info: INFO.ARG},
  {label: `${CONTROL_EXPR.INITS}:\n  `, type: "keyword", info: INFO.INITS},
  {label: `${CONTROL_EXPR.PARAMS}:\n  `, type: "keyword", info: INFO.PARAMS},
  {label: `${CONTROL_EXPR.CONSTS}:\n  `, type: "keyword", info: INFO.CONSTS},
  {label: `${CONTROL_EXPR.TOL}: `, type: "keyword", info: INFO.TOL},
  {label: `${CONTROL_EXPR.LOOP}:\n  `, type: "keyword", info: INFO.LOOP},
  {label: `${CONTROL_EXPR.UPDATE}:\n  `, type: "keyword", info: INFO.UPDATE},
  {label: `${CONTROL_EXPR.OUTPUT}:\n  `, type: "keyword", info: INFO.OUTPUT},
  {label: `${CONTROL_EXPR.COMMENT}: `, type: "keyword", info: INFO.COMMENT},
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
function lineChartOptions(colNames: string[]): Object {
  const count = colNames.length;

  return {
    xColumnName: colNames[0],
    yColumnNames: (count > 1) ? colNames.slice(1) : colNames[0],     
    showTitle: true,
    sharex: true, 
    multiAxis: count > MAX_LINE_CHART,
    multiAxisLegendPosition: "RightTop",
  };
}

/** Run solver application */
export async function runSolverApp() {

  /** Get JS-script for solving the current IVP */
  const exportToJS = async () => {
    try {
      const ivp = getIVP(editorView.state.doc.toString());
      const scriptText = getScriptLines(ivp).join('\n');      
      const script = DG.Script.create(scriptText);

      // try to call computations - correctness check
      const params = getScriptParams(ivp);    
      const call = script.prepare(params);
      await call.call();

      const sView = DG.ScriptView.create(script);
      grok.shell.addView(sView);
    }
    catch (err) {
      if (err instanceof Error)
        grok.shell.error(`${ERROR_MSG.EXPORT_TO_SCRIPT_FAILS}: ${err.message}`);
      else
        grok.shell.error(`${ERROR_MSG.EXPORT_TO_SCRIPT_FAILS}: ${ERROR_MSG.CORE_ISSUE}`);
  }};

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

      if (!solutionViewer) {
        solutionViewer = DG.Viewer.lineChart(solutionTable, lineChartOptions(solutionTable.columns.names()));
        viewerDockNode = grok.shell.dockManager.dock(
          solutionViewer, 
          DG.DOCK_TYPE.TOP, 
          solverView.dockManager.
          findNode(solverView.grid.root
        ));
      }
      else {
        solutionViewer.dataFrame = solutionTable;

        if (toChangeSolutionViewerProps) {
          solutionViewer.setOptions(lineChartOptions(solutionTable.columns.names()));
          toChangeSolutionViewerProps = false;
        }
      }

    } catch (err) {
        if (err instanceof Error) 
          grok.shell.error(`${ERROR_MSG.SOLVING_FAILS}: ${err.message}`);
        else
          grok.shell.error(`${ERROR_MSG.SOLVING_FAILS}: ${ERROR_MSG.CORE_ISSUE}`);
  }};

  /** Run model application */
  const runModelApp = async () => {
    try {
      const ivp = getIVP(editorView.state.doc.toString());
      const scriptText = getScriptLines(ivp).join('\n');    
      const script = DG.Script.create(scriptText);
      const savedScript = await grok.dapi.scripts.save(script);
      const params = getScriptParams(ivp);

      // try to call computations - correctness check
      const testCall = script.prepare(params);
      await testCall.call();

      const call = savedScript.prepare(params);
      call.edit();
    } 
    catch (err) {
      if (err instanceof Error) 
        grok.shell.error(`${ERROR_MSG.APP_CREATING_FAILS}: ${err.message}`);        
      else
        grok.shell.error(`${ERROR_MSG.APP_CREATING_FAILS}: ${ERROR_MSG.CORE_ISSUE}`);
  }};
   
  let solutionTable = DG.DataFrame.create();
  let solverView = grok.shell.addTableView(solutionTable);
  let solutionViewer: DG.Viewer | null = null;
  let viewerDockNode: DG.DockNode | null = null;
  let toChangeSolutionViewerProps = false;
  solverView.name = MISC.VIEW_DEFAULT_NAME;
  let div = ui.divV([]);   

  /** Code editor for IVP specifying */
  let editorView = new EditorView({
    doc: TEMPLATES.BASIC,
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    parent: div
  });
  
  let editorState: EDITOR_STATE = EDITOR_STATE.BASIC_TEMPLATE;
  let toShowWarning = true;
  let isChanged = false;

  editorView.dom.addEventListener('keydown', async (e) => {
    if (e.key !== HOT_KEY.RUN)
      isChanged = true;    
  });

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
    link.download = MISC.FILE_DEFAULT_NAME;
    link.click();
    URL.revokeObjectURL(link.href);
  };

  /** Set IVP code editor state */
  const setState = async (state: EDITOR_STATE, text?: string | undefined) => {
    toChangeSolutionViewerProps = true;
    isChanged = false;
    editorState = state;
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;

    const newState = EditorState.create({
      doc: text ?? getProblem(state), 
      extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    });

    editorView.setState(newState);

    if (state != EDITOR_STATE.CLEAR)
      await solve();
    else
      if (solutionViewer && viewerDockNode) {
        grok.shell.dockManager.close(viewerDockNode);
        solutionViewer = null;
      }
  };

  /** Overwrite the editor content */
  const overwrite = async (state?: EDITOR_STATE, fn?: () => Promise<void>) => {
    if (toShowWarning && isChanged) {      
      const boolInput = ui.boolInput(WARNING.CHECK, true, () => toShowWarning = !toShowWarning);      
      const dlg = ui.dialog({title: WARNING.TITLE, helpUrl: LINK.DIF_STUDIO_MD});
      solverView.append(dlg);

      dlg
        .add(ui.label(WARNING.MES))
        .add(ui.label(WARNING.QUE))
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
      .item(TITLE.LOAD, async () => await overwrite(undefined, loadFn), undefined, {description: HINT.LOAD})
      .item(TITLE.SAVE, saveFn, undefined, {description: HINT.SAVE})
      .separator()         
      .group(TITLE.TEMPL)
      .item(TITLE.BASIC, async () => await overwrite(EDITOR_STATE.BASIC_TEMPLATE), undefined, {description: HINT.BASIC})
      .item(TITLE.ADV, async () => await overwrite(EDITOR_STATE.ADVANCED_TEMPLATE), undefined, {description: HINT.ADV})
      .item(TITLE.EXT, async () => await overwrite(EDITOR_STATE.EXTENDED_TEMPLATE), undefined, {description: HINT.EXT})
      .endGroup()
      .group(TITLE.CASES)
      .item(TITLE.CHEM, async () => await overwrite(EDITOR_STATE.CHEM_REACT), undefined, {description: HINT.CHEM})
      .item(TITLE.ROB, async () => await overwrite(EDITOR_STATE.ROBERT), undefined, {description: HINT.ROB})
      .item(TITLE.FERM, async () => await overwrite(EDITOR_STATE.FERM), undefined, {description: HINT.FERM})      
      .item(TITLE.PKPD, async () => await overwrite(EDITOR_STATE.PKPD), undefined, {description: HINT.PKPD})
      .item(TITLE.ACID, async () => await overwrite(EDITOR_STATE.ACID_PROD), undefined, {description: HINT.ACID})
      .item(TITLE.NIM, async () => await overwrite(EDITOR_STATE.NIMOTUZUMAB), undefined, {description: HINT.NIM})
      .endGroup()
      .separator()
      .item(TITLE.CLEAR, async () => await overwrite(EDITOR_STATE.CLEAR), undefined, {description: HINT.CLEAR})
      .show();    
  });
  
  editorView.dom.style.overflow = 'auto';
  editorView.dom.style.height = '100%';

  solverView.dockManager.dock(div, 'left');
  solverView.helpUrl = LINK.DIF_STUDIO_MD;

  const helpIcon = ui.iconFA('question', () => {window.open(LINK.DIF_STUDIO, '_blank')}, HINT.HELP);

  const exportButton = ui.button(TITLE.TO_JS, exportToJS, HINT.TO_JS);

  const playIcon = ui.iconFA('play', solve, HINT.RUN);  
  playIcon.style.color = "var(--green-2)";
  playIcon.classList.add("fas");

  const appButton = ui.bigButton(TITLE.APP, runModelApp, HINT.APP);
  
  solverView.root.addEventListener('keydown', async (e) => {
    if (e.key === HOT_KEY.RUN) {
      e.stopImmediatePropagation();
      e.preventDefault();
      await solve();
    }
  });

  const openMenu = DG.Menu.popup()
    .item(TITLE.FROM_FILE, async () => await overwrite(undefined, loadFn), undefined, {description: HINT.LOAD})
    .group(TITLE.TEMPL)
    .item(TITLE.BASIC, async () => await overwrite(EDITOR_STATE.BASIC_TEMPLATE), undefined, {description: HINT.BASIC})
    .item(TITLE.ADV, async () => await overwrite(EDITOR_STATE.ADVANCED_TEMPLATE), undefined, {description: HINT.ADV})
    .item(TITLE.EXT, async () => await overwrite(EDITOR_STATE.EXTENDED_TEMPLATE), undefined, {description: HINT.EXT})
    .endGroup()
    .group(TITLE.CASES)
    .item(TITLE.CHEM, async () => await overwrite(EDITOR_STATE.CHEM_REACT), undefined, {description: HINT.CHEM})
    .item(TITLE.ROB, async () => await overwrite(EDITOR_STATE.ROBERT), undefined, {description: HINT.ROB})
    .item(TITLE.FERM, async () => await overwrite(EDITOR_STATE.FERM), undefined, {description: HINT.FERM})   
    .item(TITLE.PKPD, async () => await overwrite(EDITOR_STATE.PKPD), undefined, {description: HINT.PKPD})
    .item(TITLE.ACID, async () => await overwrite(EDITOR_STATE.ACID_PROD), undefined, {description: HINT.ACID})
    .item(TITLE.NIM, async () => await overwrite(EDITOR_STATE.NIMOTUZUMAB), undefined, {description: HINT.NIM})
    .endGroup();

  const openIcon = ui.iconFA('folder-open', () => openMenu.show(), HINT.OPEN);
  const saveIcon = ui.iconFA('save', async () => {await saveFn()}, HINT.SAVE);

  solverView.setRibbonPanels([[openIcon, saveIcon], [helpIcon], [exportButton, appButton], [playIcon]]);
} // runSolverApp

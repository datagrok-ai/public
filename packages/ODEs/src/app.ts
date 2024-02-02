// Application for solving initial value problems (IVP) & demo app

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {basicSetup, EditorView} from "codemirror";
import {EditorState} from "@codemirror/state";
import {python} from "@codemirror/lang-python";
import {autocompletion} from "@codemirror/autocomplete";

import {DF_NAME, CONTROL_EXPR, MAX_LINE_CHART} from './constants';
import {TEMPLATES, DEMO_TEMPLATE} from './templates';
import { USE_CASES } from './use-cases';
import {HINT, TITLE, LINK, HOT_KEY, ERROR_MSG, INFO, WARNING, MISC, demoInfo, INPUT_TYPE, PATH} from './ui-constants';
import {getIVP, getScriptLines, getScriptParams, IVP, Input, SCRIPTING,
  BRACE_OPEN, BRACE_CLOSE, BRACKET_OPEN, BRACKET_CLOSE, ANNOT_SEPAR, CONTROL_SEP} from './scripting-tools';

/** State of IVP code editor */
enum EDITOR_STATE {
  EMPTY = 'empty',
  BASIC_TEMPLATE = 'template',
  ADVANCED_TEMPLATE = 'advanced',
  FROM_FILE = 'from-file',
  EXTENDED_TEMPLATE = 'extended',
  CHEM_REACT = 'chem-react',
  ROBERT = 'robertson',
  FERM = 'fermentation',
  PKPD = 'pk-pd',
  ACID_PROD = 'ga-production',
  NIMOTUZUMAB = 'nimotuzumab',
};

/** Models & templates */
const MODELS: string[] = [ EDITOR_STATE.BASIC_TEMPLATE,
  EDITOR_STATE.ADVANCED_TEMPLATE,
  EDITOR_STATE.EXTENDED_TEMPLATE,
  EDITOR_STATE.CHEM_REACT,
  EDITOR_STATE.ROBERT,
  EDITOR_STATE.FERM,
  EDITOR_STATE.PKPD,
  EDITOR_STATE.ACID_PROD,
  EDITOR_STATE.NIMOTUZUMAB,
];

/** Get problem with respect to IVP editor state */
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

/** Return help link with respect to IVP editor state */
function getLink(state: EDITOR_STATE): string {
  switch (state) {
    case EDITOR_STATE.CHEM_REACT:
      return LINK.CHEM_REACT;

    case EDITOR_STATE.ROBERT:
      return LINK.ROBERTSON;    

    case EDITOR_STATE.FERM:
      return LINK.FERMENTATION;

    case EDITOR_STATE.PKPD:
      return LINK.PKPD;

    case EDITOR_STATE.ACID_PROD:
      return LINK.GA_PRODUCTION;

    case EDITOR_STATE.NIMOTUZUMAB:
      return LINK.NIMOTUZUMAB;

    default:
      return LINK.DIF_STUDIO_REL;
  }
} // getLink

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
function getLineChartOptions(colNames: string[]): Object {
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
export async function runSolverApp(content?: string) {

  /** Get JS-script for solving the current IVP */
  const exportToJS = async () => {
    try {
      const ivp = getIVP(editorView.state.doc.toString());
      const scriptText = getScriptLines(ivp, true, true).join('\n');      
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
        grok.shell.error(`${ERROR_MSG.EXPORT_TO_SCRIPT_FAILS}: ${ERROR_MSG.SCRIPTING_ISSUE}`);
  }};

  /** Solve IVP */
  const solve = async (ivp: IVP, inputsPath: string) => {
    try {
      solverView.path = `${solverMainPath}${PATH.PARAM}${inputsPath}`;

      const start = ivp.arg.initial.value;
      const finish = ivp.arg.final.value;
      const step = ivp.arg.step.value;
  
      if (start >= finish)
        return;
  
      if ((step <= 0) || (step > finish - start))
        return;
  
      const scriptText = getScriptLines(ivp).join('\n');    
      const script = DG.Script.create(scriptText);
      const params = getScriptParams(ivp);    
      const call = script.prepare(params);
  
      await call.call();
        
      solutionTable = call.outputs[DF_NAME];
      solverView.dataFrame = call.outputs[DF_NAME];
      solverView.name = solutionTable.name;
  
      if (!solutionViewer) {
        solutionViewer = DG.Viewer.lineChart(solutionTable, getLineChartOptions(solutionTable.columns.names()));
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
          solutionViewer.setOptions(getLineChartOptions(solutionTable.columns.names()));
          toChangeSolutionViewerProps = false;
        }
      }

      isSolvingSuccess = true;
    } catch (error) {
      clearSolution();

      if (error instanceof Error) 
          grok.shell.error(error.message);
      else
          grok.shell.error(ERROR_MSG.SCRIPTING_ISSUE);      
    }
  }; // solve

  /** Run solving the current IVP */
  const runSolving = async (showApp: boolean) => {
    if (prevInputsNode !== null)
        inputsDiv.removeChild(prevInputsNode);

    try {
      const ivp = getIVP(editorView.state.doc.toString());      

      prevInputsNode = inputsDiv.appendChild(await getInputsUI(ivp, solve, startingInputs));

      runPane.header.hidden = !isSolvingSuccess;

      if (isSolvingSuccess) {
        toChangeInputs = false;      
        tabControl.currentPane = (showApp && isSolvingSuccess)? runPane : modelPane;
      }
      else
        tabControl.currentPane = modelPane;

      } catch (error) {
          prevInputsNode = null;
          runPane.header.hidden = true;
          tabControl.currentPane = modelPane;
          clearSolution();

          if (error instanceof Error) 
            grok.shell.error(error.message);
          else
            grok.shell.error(ERROR_MSG.UI_ISSUE);
      }          
  }; // runSolving

  /** Clear solution table & viewer */
  const clearSolution = () => {
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;

    if (solutionViewer && viewerDockNode) {
      grok.shell.dockManager.close(viewerDockNode);
      solutionViewer = null;
      solverView.path = PATH.EMPTY;
    }
  } // clearSolution
   
  let solutionTable = DG.DataFrame.create();
  const startingPath = window.location.href;
  let startingInputs: Map<string, number> | null = null;
  let solverView = grok.shell.addTableView(solutionTable);
  let solverMainPath: string = PATH.CUSTOM;
  let solutionViewer: DG.Viewer | null = null;
  let viewerDockNode: DG.DockNode | null = null;
  let toChangeSolutionViewerProps = false;
  solverView.name = MISC.VIEW_DEFAULT_NAME;
  let modelDiv = ui.divV([]);
  let inputsDiv = ui.divV([]);
  let prevInputsNode: Node | null = null;
  const tabControl = ui.tabControl();

  let editorState: EDITOR_STATE = EDITOR_STATE.BASIC_TEMPLATE;
  let toShowWarning = true;
  let isModelChanged = false;
  let toChangeInputs = false;
  let isSolvingSuccess = false;

  const modelPane = tabControl.addPane(TITLE.MODEL, () => modelDiv);
  const runPane = tabControl.addPane(TITLE.IPUTS, () => inputsDiv);

  tabControl.onTabChanged.subscribe(async (_) => { 
    if ((tabControl.currentPane === runPane) && toChangeInputs) 
      await runSolving(true);
  });

  /** Code editor for IVP specifying */
  let editorView = new EditorView({
    doc: content ?? TEMPLATES.BASIC,
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    parent: modelDiv
  }); 


  editorView.dom.addEventListener('keydown', async (e) => {
    if (e.key !== HOT_KEY.RUN) {
      isModelChanged = true;
      toChangeInputs = true;
      solverView.path = PATH.CUSTOM;
      solverMainPath = PATH.CUSTOM;
      startingInputs = null;
      solverView.helpUrl = LINK.DIF_STUDIO_REL;
      isSolvingSuccess = false;
      runPane.header.hidden = false;
    }
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
        setState(EDITOR_STATE.FROM_FILE, true, text);
        dlg.close();
      }, false);
          
      if (file) 
        reader.readAsText(file);
    }     

    dlg.add(fileInp);        
    fileInp.click();
  }; // loadFn

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
  const setState = async (state: EDITOR_STATE, toClearStartingInputs: boolean = true, text?: string | undefined) => {
    toChangeSolutionViewerProps = true;
    isModelChanged = false;
    editorState = state;
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;
    solverView.helpUrl = getLink(state);

    if (toClearStartingInputs)
      startingInputs = null;

    const newState = EditorState.create({
      doc: text ?? getProblem(state), 
      extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    });

    editorView.setState(newState);
    
    // set path
    if (state === EDITOR_STATE.FROM_FILE)
      solverMainPath = PATH.CUSTOM;
    else
      solverMainPath = `${PATH.MODEL}${state}`;

    switch(state) {
      case EDITOR_STATE.EMPTY:
        clearSolution();
        break;

      case EDITOR_STATE.BASIC_TEMPLATE:
      case EDITOR_STATE.ADVANCED_TEMPLATE:
      case EDITOR_STATE.EXTENDED_TEMPLATE:
        await runSolving(false);
        break;
      
      default:
        await runSolving(true);
        break;
    }
  }; // setState

  /** Overwrite the editor content */
  const overwrite = async (state?: EDITOR_STATE, fn?: () => Promise<void>) => {
    if (toShowWarning && isModelChanged) {      
      const boolInput = ui.boolInput(WARNING.CHECK, true, () => toShowWarning = !toShowWarning);      
      const dlg = ui.dialog({title: WARNING.TITLE, helpUrl: LINK.DIF_STUDIO_REL});
      solverView.append(dlg);

      dlg
        .add(ui.label(WARNING.MES))
        .add(boolInput.root)
        .onCancel(() => dlg.close())
        .onOK(async () => {
          if (fn)
            await fn();
          else
            setState(state ?? EDITOR_STATE.EMPTY);          
        })
        .show();
    }
    else if (fn)
      await fn();
    else
      setState(state ?? EDITOR_STATE.EMPTY);
  }; // overwrite

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
      .item(TITLE.CLEAR, async () => await overwrite(EDITOR_STATE.EMPTY), undefined, {description: HINT.CLEAR})
      .show();    
  });
  
  editorView.dom.style.overflow = 'auto';
  editorView.dom.style.height = '100%';
  const node = solverView.dockManager.dock(tabControl.root, 'left');
  if (node.container.dart.elementTitle)
    node.container.dart.elementTitle.hidden = true;
  solverView.helpUrl = LINK.DIF_STUDIO_REL;
  
  // routing
  if (content) {    
    await runSolving(false);
  }
  else {
    const modelIdx = startingPath.indexOf(PATH.MODEL);
    const paramsIdx = startingPath.indexOf(PATH.PARAM);

    if (modelIdx > -1) {
      const model = startingPath.slice(modelIdx + PATH.MODEL.length, (paramsIdx > -1) ? paramsIdx : undefined);
      
      if (MODELS.includes(model)) {
        startingInputs = new Map<string, number>();

        if (modelIdx < paramsIdx)
          try {
            startingPath.slice(paramsIdx + PATH.PARAM.length).split(PATH.AND).forEach((equality) => {
              const eqIdx = equality.indexOf(PATH.EQ);
              startingInputs?.set(equality.slice(0, eqIdx).toLowerCase(), Number(equality.slice(eqIdx + 1)));
            });
          } catch (error) {
            startingInputs = null;      
          }
        
        await setState(model as EDITOR_STATE, false);
      }
      else 
        await setState(EDITOR_STATE.BASIC_TEMPLATE);      
    } 
    else 
      await setState(EDITOR_STATE.BASIC_TEMPLATE);
  }

  const helpIcon = ui.iconFA('question', () => {window.open(LINK.DIF_STUDIO, '_blank')}, HINT.HELP);

  const exportButton = ui.bigButton(TITLE.TO_JS, exportToJS, HINT.TO_JS); 
  
  solverView.root.addEventListener('keydown', async (e) => {
    if (e.key === HOT_KEY.RUN) {
      e.stopImmediatePropagation();
      e.preventDefault();

      if (tabControl.currentPane === modelPane) 
        if ( toChangeInputs )
          await runSolving(true);
        else
          tabControl.currentPane = runPane;
  }});

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

  solverView.setRibbonPanels([[openIcon, saveIcon, exportButton, helpIcon]]);
} // runSolverApp

export async function runSolverAppNew(content?: string)  {

  /** Get JS-script for solving the current IVP */
  const exportToJS = async () => {
    try {
      const ivp = getIVP(editorView.state.doc.toString());
      const scriptText = getScriptLines(ivp, true, true).join('\n');      
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
        grok.shell.error(`${ERROR_MSG.EXPORT_TO_SCRIPT_FAILS}: ${ERROR_MSG.SCRIPTING_ISSUE}`);
  }};

  /** Solve IVP */
  const solve = async (ivp: IVP, inputsPath: string) => {
    try {
      solverView.path = `${solverMainPath}${PATH.PARAM}${inputsPath}`;

      const start = ivp.arg.initial.value;
      const finish = ivp.arg.final.value;
      const step = ivp.arg.step.value;
  
      if (start >= finish)
        return;
  
      if ((step <= 0) || (step > finish - start))
        return;
  
      const scriptText = getScriptLines(ivp).join('\n');    
      const script = DG.Script.create(scriptText);
      const params = getScriptParams(ivp);    
      const call = script.prepare(params);
  
      await call.call();
        
      solutionTable = call.outputs[DF_NAME];
      solverView.dataFrame = call.outputs[DF_NAME];
      solverView.name = solutionTable.name;
  
      if (!solutionViewer) {
        solutionViewer = DG.Viewer.lineChart(solutionTable, getLineChartOptions(solutionTable.columns.names()));
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
          solutionViewer.setOptions(getLineChartOptions(solutionTable.columns.names()));
          toChangeSolutionViewerProps = false;
        }
      }

      isSolvingSuccess = true;
    } catch (error) {
      clearSolution();

      if (error instanceof Error) 
          grok.shell.error(error.message);
      else
          grok.shell.error(ERROR_MSG.SCRIPTING_ISSUE);      
    }
  }; // solve

  /** Run solving the current IVP */
  const runSolving = async (showApp: boolean) => {
    if (prevInputsNode !== null)
        inputsDiv.removeChild(prevInputsNode);

    try {
      const ivp = getIVP(editorView.state.doc.toString());      

      prevInputsNode = inputsDiv.appendChild(await getInputsUI(ivp, solve, startingInputs));

      runPane.header.hidden = !isSolvingSuccess;

      if (isSolvingSuccess) {
        toChangeInputs = false;      
        tabControl.currentPane = (showApp && isSolvingSuccess)? runPane : modelPane;
      }
      else
        tabControl.currentPane = modelPane;

      } catch (error) {
          prevInputsNode = null;
          runPane.header.hidden = true;
          tabControl.currentPane = modelPane;
          clearSolution();

          if (error instanceof Error) 
            grok.shell.error(error.message);
          else
            grok.shell.error(ERROR_MSG.UI_ISSUE);
      }          
  }; // runSolving

  /** Clear solution table & viewer */
  const clearSolution = () => {
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;

    if (solutionViewer && viewerDockNode) {
      grok.shell.dockManager.close(viewerDockNode);
      solutionViewer = null;
      solverView.path = PATH.EMPTY;
    }
  } // clearSolution
   
  let solutionTable = DG.DataFrame.create();
  const startingPath = window.location.href;
  let startingInputs: Map<string, number> | null = null;
  let solverView = grok.shell.addTableView(solutionTable);
  let solverMainPath: string = PATH.CUSTOM;
  let solutionViewer: DG.Viewer | null = null;
  let viewerDockNode: DG.DockNode | null = null;
  let toChangeSolutionViewerProps = false;
  solverView.name = MISC.VIEW_DEFAULT_NAME;
  let modelDiv = ui.divV([]);
  let inputsDiv = ui.divV([]);
  let prevInputsNode: Node | null = null;
  const tabControl = ui.tabControl();

  let editorState: EDITOR_STATE = EDITOR_STATE.BASIC_TEMPLATE;
  let toShowWarning = true;
  let isModelChanged = false;
  let toChangeInputs = false;
  let isSolvingSuccess = false;

  const modelPane = tabControl.addPane(TITLE.MODEL, () => modelDiv);
  const runPane = tabControl.addPane(TITLE.IPUTS, () => inputsDiv);

  tabControl.onTabChanged.subscribe(async (_) => { 
    if ((tabControl.currentPane === runPane) && toChangeInputs) 
      await runSolving(true);
  });

  /** Code editor for IVP specifying */
  let editorView = new EditorView({
    doc: content ?? TEMPLATES.BASIC,
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    parent: modelDiv
  }); 


  editorView.dom.addEventListener('keydown', async (e) => {
    if (e.key !== HOT_KEY.RUN) {
      isModelChanged = true;
      toChangeInputs = true;
      solverView.path = PATH.CUSTOM;
      solverMainPath = PATH.CUSTOM;
      startingInputs = null;
      solverView.helpUrl = LINK.DIF_STUDIO_REL;
      isSolvingSuccess = false;
      runPane.header.hidden = false;
    }
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
        setState(EDITOR_STATE.FROM_FILE, true, text);
        dlg.close();
      }, false);
          
      if (file) 
        reader.readAsText(file);
    }     

    dlg.add(fileInp);        
    fileInp.click();
  }; // loadFn

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
  const setState = async (state: EDITOR_STATE, toClearStartingInputs: boolean = true, text?: string | undefined) => {
    toChangeSolutionViewerProps = true;
    isModelChanged = false;
    editorState = state;
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;
    solverView.helpUrl = getLink(state);

    if (toClearStartingInputs)
      startingInputs = null;

    const newState = EditorState.create({
      doc: text ?? getProblem(state), 
      extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    });

    editorView.setState(newState);
    
    // set path
    if (state === EDITOR_STATE.FROM_FILE) {
      solverView.basePath = PATH.APPS_DS;
      solverMainPath = PATH.CUSTOM;
    }
    else if (state === EDITOR_STATE.EMPTY) {
      solverView.basePath = PATH.APPS_DS;
      solverMainPath = PATH.EMPTY;
    }
    else {
      solverView.basePath = PATH.APP_DATA_DS;
      solverMainPath = `/${state}.${MISC.FILE_EXT}`;
    }

    switch(state) {
      case EDITOR_STATE.EMPTY:
        clearSolution();
        break;

      case EDITOR_STATE.BASIC_TEMPLATE:
      case EDITOR_STATE.ADVANCED_TEMPLATE:
      case EDITOR_STATE.EXTENDED_TEMPLATE:
        await runSolving(false);
        break;
      
      default:
        await runSolving(true);
        break;
    }
  }; // setState

  /** Overwrite the editor content */
  const overwrite = async (state?: EDITOR_STATE, fn?: () => Promise<void>) => {
    if (toShowWarning && isModelChanged) {      
      const boolInput = ui.boolInput(WARNING.CHECK, true, () => toShowWarning = !toShowWarning);      
      const dlg = ui.dialog({title: WARNING.TITLE, helpUrl: LINK.DIF_STUDIO_REL});
      solverView.append(dlg);

      dlg
        .add(ui.label(WARNING.MES))
        .add(boolInput.root)
        .onCancel(() => dlg.close())
        .onOK(async () => {
          if (fn)
            await fn();
          else
            setState(state ?? EDITOR_STATE.EMPTY);          
        })
        .show();
    }
    else if (fn)
      await fn();
    else
      setState(state ?? EDITOR_STATE.EMPTY);
  }; // overwrite

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
      .item(TITLE.CLEAR, async () => await overwrite(EDITOR_STATE.EMPTY), undefined, {description: HINT.CLEAR})
      .show();    
  });
  
  editorView.dom.style.overflow = 'auto';
  editorView.dom.style.height = '100%';
  const node = solverView.dockManager.dock(tabControl.root, 'left');
  if (node.container.dart.elementTitle)
    node.container.dart.elementTitle.hidden = true;
  solverView.helpUrl = LINK.DIF_STUDIO_REL;
  
  // routing
  if (content) {    
    await runSolving(false);
  }
  else {
    let eqIdx = startingPath.indexOf(PATH.EQ);
    const andIdx = startingPath.indexOf(PATH.AND);

    if ((eqIdx > -1) && (eqIdx < andIdx - 1)) {
      const model = startingPath.slice(eqIdx + 1, andIdx);

      if (MODELS.includes(model)) {
        try {
          startingInputs = new Map<string, number>();

          startingPath.slice(andIdx + 1).split(PATH.AND).forEach((equality) => {
            eqIdx = equality.indexOf(PATH.EQ);
            startingInputs?.set(equality.slice(0, eqIdx).toLowerCase(), Number(equality.slice(eqIdx + 1)));
          });
        } catch (error) {
          startingInputs = null;      
        }        

        await setState(model as EDITOR_STATE, false);        
      }
    }
    else 
      await setState(EDITOR_STATE.BASIC_TEMPLATE);
  }

  const helpIcon = ui.iconFA('question', () => {window.open(LINK.DIF_STUDIO, '_blank')}, HINT.HELP);

  const exportButton = ui.bigButton(TITLE.TO_JS, exportToJS, HINT.TO_JS); 
  
  solverView.root.addEventListener('keydown', async (e) => {
    if (e.key === HOT_KEY.RUN) {
      e.stopImmediatePropagation();
      e.preventDefault();

      if (tabControl.currentPane === modelPane) 
        if ( toChangeInputs )
          await runSolving(true);
        else
          tabControl.currentPane = runPane;
  }});

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

  solverView.setRibbonPanels([[openIcon, saveIcon, exportButton, helpIcon]]);
} // runSolverApp

/** Run solver demo application */
export async function runSolverDemoApp() { 

  /** Get JS-script for solving the current IVP */
  const exportToJS = async () => {
    try {
      const ivp = getIVP(editorView.state.doc.toString());
      const scriptText = getScriptLines(ivp, true, true).join('\n');      
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
        grok.shell.error(`${ERROR_MSG.EXPORT_TO_SCRIPT_FAILS}: ${ERROR_MSG.SCRIPTING_ISSUE}`);
  }};

  /** Solve IVP */
  const solve = async (ivp: IVP) => {
    try {      
      const start = ivp.arg.initial.value;
      const finish = ivp.arg.final.value;
      const step = ivp.arg.step.value;
  
      if (start >= finish)
        return;
  
      if ((step <= 0) || (step > finish - start))
        return;
  
      const scriptText = getScriptLines(ivp).join('\n');    
      const script = DG.Script.create(scriptText);
      const params = getScriptParams(ivp);    
      const call = script.prepare(params);
  
      await call.call();
        
      solutionTable = call.outputs[DF_NAME];
      solverView.dataFrame = call.outputs[DF_NAME];
      solverView.name = solutionTable.name;
  
      if (!solutionViewer) {
        solutionViewer = DG.Viewer.lineChart(solutionTable, getLineChartOptions(solutionTable.columns.names()));
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
          solutionViewer.setOptions(getLineChartOptions(solutionTable.columns.names()));
          toChangeSolutionViewerProps = false;
        }
      }

      isSolvingSuccess = true;
    } catch (error) {
      clearSolution();

      if (error instanceof Error) 
          grok.shell.error(error.message);
      else
          grok.shell.error(ERROR_MSG.SCRIPTING_ISSUE);      
    }
  }; // solve

  /** Run solving the current IVP */
  const runSolving = async (showApp: boolean) => {
    if (prevInputsNode !== null)
        inputsDiv.removeChild(prevInputsNode);

    try {
      const ivp = getIVP(editorView.state.doc.toString());      

      prevInputsNode = inputsDiv.appendChild(await getInputsUI(ivp, solve, null));

      runPane.header.hidden = !isSolvingSuccess;

      if (isSolvingSuccess) {
        toChangeInputs = false;      
        tabControl.currentPane = (showApp && isSolvingSuccess)? runPane : modelPane;
      }
      else
        tabControl.currentPane = modelPane;

      } catch (error) {
          prevInputsNode = null;
          runPane.header.hidden = true;
          tabControl.currentPane = modelPane;
          clearSolution();

          if (error instanceof Error) 
            grok.shell.error(error.message);
          else
            grok.shell.error(ERROR_MSG.UI_ISSUE);
      }          
  }; // runSolving

  /** Clear solution table & viewer */
  const clearSolution = () => {
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;

    if (solutionViewer && viewerDockNode) {
      grok.shell.dockManager.close(viewerDockNode);
      solutionViewer = null;
      solverView.path = PATH.EMPTY;
    }
  } // clearSolution
   
  let solutionTable = DG.DataFrame.create();  
  let solverView = grok.shell.addTableView(solutionTable);
  let solutionViewer: DG.Viewer | null = null;
  let viewerDockNode: DG.DockNode | null = null;
  let toChangeSolutionViewerProps = false;
  solverView.name = MISC.VIEW_DEFAULT_NAME;
  let modelDiv = ui.divV([]);
  let inputsDiv = ui.divV([]);
  let prevInputsNode: Node | null = null;
  const tabControl = ui.tabControl();

  let editorState: EDITOR_STATE = EDITOR_STATE.BASIC_TEMPLATE;
  let toShowWarning = false;
  let isModelChanged = false;
  let toChangeInputs = false;
  let isSolvingSuccess = false;

  const modelPane = tabControl.addPane(TITLE.MODEL, () => modelDiv);
  const runPane = tabControl.addPane(TITLE.IPUTS, () => inputsDiv);

  tabControl.onTabChanged.subscribe(async (_) => { 
    if ((tabControl.currentPane === runPane) && toChangeInputs) 
      await runSolving(true);
  });

  /** Code editor for IVP specifying */
  let editorView = new EditorView({
    doc: DEMO_TEMPLATE,
    extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    parent: modelDiv
  }); 

  editorView.dom.addEventListener('keydown', async (e) => {
    if (e.key !== HOT_KEY.RUN) {
      isModelChanged = true;
      toChangeInputs = true;
      solverView.helpUrl = LINK.DIF_STUDIO_REL;
      isSolvingSuccess = false;
      runPane.header.hidden = false;
    }
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
        setState(EDITOR_STATE.FROM_FILE, true, text);
        dlg.close();
      }, false);
          
      if (file) 
        reader.readAsText(file);
    }     

    dlg.add(fileInp);        
    fileInp.click();
  }; // loadFn

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
  const setState = async (state: EDITOR_STATE, toClearStartingInputs: boolean = true, text?: string | undefined) => {
    toChangeSolutionViewerProps = true;
    isModelChanged = false;
    editorState = state;
    solutionTable = DG.DataFrame.create();
    solverView.dataFrame = solutionTable;
    solverView.helpUrl = getLink(state);    

    const newState = EditorState.create({
      doc: text ?? getProblem(state), 
      extensions: [basicSetup, python(), autocompletion({override: [contrCompletions]})],
    });

    editorView.setState(newState);

    switch(state) {
      case EDITOR_STATE.EMPTY:
        clearSolution();
        break;

      case EDITOR_STATE.BASIC_TEMPLATE:
      case EDITOR_STATE.ADVANCED_TEMPLATE:
      case EDITOR_STATE.EXTENDED_TEMPLATE:
        await runSolving(false);
        break;
      
      default:
        await runSolving(true);
        break;
    }
  }; // setState

  /** Overwrite the editor content */
  const overwrite = async (state?: EDITOR_STATE, fn?: () => Promise<void>) => {
    if (toShowWarning && isModelChanged) {      
      const boolInput = ui.boolInput(WARNING.CHECK, true, () => toShowWarning = !toShowWarning);      
      const dlg = ui.dialog({title: WARNING.TITLE, helpUrl: LINK.DIF_STUDIO_REL});
      solverView.append(dlg);

      dlg
        .add(ui.label(WARNING.MES))
        .add(boolInput.root)
        .onCancel(() => dlg.close())
        .onOK(async () => {
          if (fn)
            await fn();
          else
            setState(state ?? EDITOR_STATE.EMPTY);          
        })
        .show();
    }
    else if (fn)
      await fn();
    else
      setState(state ?? EDITOR_STATE.EMPTY);
  }; // overwrite

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
      .item(TITLE.CLEAR, async () => await overwrite(EDITOR_STATE.EMPTY), undefined, {description: HINT.CLEAR})
      .show();    
  });
  
  editorView.dom.style.overflow = 'auto';
  editorView.dom.style.height = '100%';
  const node = solverView.dockManager.dock(tabControl.root, 'left');
  node.container.dart.elementTitle.hidden = true;
  solverView.helpUrl = LINK.DIF_STUDIO_REL;  

  const helpIcon = ui.iconFA('question', () => {window.open(LINK.DIF_STUDIO, '_blank')}, HINT.HELP);

  const exportButton = ui.bigButton(TITLE.TO_JS, exportToJS, HINT.TO_JS); 
  
  solverView.root.addEventListener('keydown', async (e) => {
    if (e.key === HOT_KEY.RUN) {
      e.stopImmediatePropagation();
      e.preventDefault();

      if (tabControl.currentPane === modelPane) 
        if (toChangeInputs)
          await runSolving(true);
        else
          tabControl.currentPane = runPane;
  }});

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

  solverView.setRibbonPanels([[openIcon, saveIcon, exportButton, helpIcon]]);

  grok.shell.windows.context.visible = false;
  grok.shell.windows.help.visible = false;
  grok.shell.windows.help.syncCurrentObject = true;

  const helpMD = ui.markdown(demoInfo);
  const divHelp = ui.div([helpMD]);  
  divHelp.style.padding = '10px';
  divHelp.style.overflow = 'auto';
  helpMD.style.fontWeight = 'lighter';
  solverView.dockManager.dock(divHelp, 'right', undefined, undefined, 0.3);

  await runSolving(false);
} // runSolverDemoApp

/** Return model inputs UI */
async function getInputsUI(ivp: IVP, solveFn: (ivp: IVP, inputsPath: string) => Promise<void>, startingInputs: Map<string, number> | null): Promise<HTMLDivElement> {
  /**  String to value */
  const strToVal = (s: string) => {
    let num = Number(s);
  
    if (!isNaN(num))
      return num;
    
    if (s === 'true')
      return true;
  
    if (s === 'false')
      return false;
  
    return s;
  };

  /** Return options with respect to the model input specification */
  const getOptions = (name: string, modelInput: Input, modelBlock: string) => {
    let options: DG.PropertyOptions = { 
      name: name,
      defaultValue: modelInput.value,
      type: DG.TYPE.FLOAT,
      inputType: INPUT_TYPE.FLOAT,
    };

    if (modelInput.annot !== null) {
      let annot = modelInput.annot;
      let descr: string | undefined = undefined;

      let posOpen = annot.indexOf(BRACKET_OPEN);
      let posClose = annot.indexOf(BRACKET_CLOSE);

      if (posOpen !== -1) {    
        if (posClose === -1)
          throw new Error(`${ERROR_MSG.MISSING_CLOSING_BRACKET}, see '${name}' in ${modelBlock}-block`);
    
        descr = annot.slice(posOpen + 1, posClose);
    
        annot = annot.slice(0, posOpen);
      }

      posOpen = annot.indexOf(BRACE_OPEN);
      posClose = annot.indexOf(BRACE_CLOSE);

      if (posOpen >= posClose)
        throw new Error(`${ERROR_MSG.INCORRECT_BRACES_USE}, see '${name}' in ${modelBlock}-block`);
    
      let pos: number;
      let key: string;
      let val;

      annot.slice(posOpen + 1, posClose).split(ANNOT_SEPAR).forEach((str) => {
        pos = str.indexOf(CONTROL_SEP);
      
        if (pos === -1)
          throw new Error(`${ERROR_MSG.MISSING_COLON}, see '${name}' in ${modelBlock}-block`);
      
        key = str.slice(0, pos).trim();
        val = str.slice(pos + 1).trim();

        // @ts-ignore
        options[key] = strToVal(val);
      });

      options.description = descr ?? '';
      options.name = options.caption ?? options.name;
      options.caption = options.name;
    }

    if (startingInputs) {
      options.defaultValue = startingInputs.get(options.name!.replace(' ', '').toLowerCase()) ?? options.defaultValue;
      modelInput.value = options.defaultValue;
    }

    return options;
  }; // getOptions
  
  const inputsByCategories = new Map<string, DG.InputBase[]>();
  inputsByCategories.set(TITLE.MISC, []);
  let options: DG.PropertyOptions;

  /** Pull input to appropriate category & add tooltip */
  const categorizeInput = (options: DG.PropertyOptions, input: DG.InputBase) => {
    let category = options.category;

    if (category === undefined)
      inputsByCategories.get(TITLE.MISC)?.push(input);
    else if (inputsByCategories.has(category))
      inputsByCategories.get(category)!.push(input)
    else
      inputsByCategories.set(category, [input]);

    input.setTooltip(options.description!);
  };

  /** Return line with inputs names & values */
  const getInputsPath = () => {
    let line = '';
    
    inputsByCategories.forEach((inputs, cat) => {
      if (cat !== TITLE.MISC)
        inputs.forEach((input) => {line += `${PATH.AND}${input.caption.replace(' ', '')}${PATH.EQ}${input.value}`});
    });

    inputsByCategories.get(TITLE.MISC)!.forEach((input) => {line += `${PATH.AND}${input.caption.replace(' ', '')}${PATH.EQ}${input.value}`});

    return line.slice(1); // we ignore 1-st '&'
  };

  // Inputs for argument
  for (const key in ivp.arg)
    if (key !== SCRIPTING.ARG_NAME) {
      //@ts-ignore
      options = getOptions(key, ivp.arg[key], CONTROL_EXPR.ARG);
      const input = ui.input.forProperty(DG.Property.fromOptions(options));
      input.onChanged(async () => {
        if (input.value !== null) {
          //@ts-ignore
          ivp.arg[key].value = input.value;
          await solveFn(ivp, getInputsPath());
        }
      });

      categorizeInput(options, input);
    }

  // Inputs for initial values
  ivp.inits.forEach((val, key, map) => {
    options = getOptions(key, val, CONTROL_EXPR.INITS);
    const input = ui.input.forProperty(DG.Property.fromOptions(options));
    input.onChanged(async () => {
      if (input.value !== null) {
        ivp.inits.get(key)!.value = input.value;
        await solveFn(ivp, getInputsPath());
      }
    });

    categorizeInput(options, input);
  });

  // Inputs for parameters
  if (ivp.params !== null)
    ivp.params.forEach((val, key, map) => {
      options = getOptions(key, val, CONTROL_EXPR.PARAMS);
      const input = ui.input.forProperty(DG.Property.fromOptions(options));
      input.onChanged(async () => {
        if (input.value !== null) {
          ivp.params!.get(key)!.value = input.value;
          await solveFn(ivp, getInputsPath());
        }
      });

      categorizeInput(options, input);
    });

  // Inputs for loop
  if (ivp.loop !== null) {
    options = getOptions(SCRIPTING.COUNT, ivp.loop.count, CONTROL_EXPR.LOOP);
    options.inputType = INPUT_TYPE.INT; // since it's an integer
    options.type = DG.TYPE.INT; // since it's an integer
    const input = ui.input.forProperty(DG.Property.fromOptions(options));
    input.onChanged(async () => {
      if (input.value !== null) {
        ivp.loop!.count.value = input.value;
        await solveFn(ivp, getInputsPath());
      }
    });

    categorizeInput(options, input);    
  }

  // Inputs form
  const form = ui.form([]);

  if (inputsByCategories.size === 1)
    inputsByCategories.get(TITLE.MISC)!.forEach((input) => form.append(input.root));
  else {
    inputsByCategories.forEach((inputs, category) => {
      if (category !== TITLE.MISC) {
        form.append(ui.h2(category));
        inputs.forEach((inp) => form.append(inp.root));
      }
    });

    if (inputsByCategories.get(TITLE.MISC)!.length > 0) {
      form.append(ui.h2(TITLE.MISC));
      inputsByCategories.get(TITLE.MISC)!.forEach((input) => form.append(input.root));     
    }
  } 

  form.style.overflowY = 'auto';
  form.style.padding = '5px';

  await solveFn(ivp, getInputsPath());
  
  return form;
} // getInputsUI

/** Return file preview view */
export async function getFilePreview(file: DG.FileInfo): Promise<DG.View> {
  const browseView = grok.shell.view(TITLE.BROWSE);
  const mainPath = browseView!.path;  
  const equations = await file.readAsString();
  const startingPath = window.location.href;
  let startingInputs: Map<string, number> | null = null;

  const paramsIdx = startingPath.indexOf(PATH.PARAM);

  if (paramsIdx > -1) {
    try {
      startingInputs = new Map<string, number>();
      
      console.log(startingPath.slice(paramsIdx + PATH.PARAM.length));

      startingPath.slice(paramsIdx + PATH.PARAM.length).split(PATH.AND).forEach((equality) => {
        let eqIdx = equality.indexOf(PATH.EQ);
        startingInputs!.set(equality.slice(0, eqIdx).toLowerCase(), Number(equality.slice(eqIdx + 1)));
      });

      console.log(startingInputs);
    } catch (error) {
      startingInputs = null;      
    }    
  }

  const modelDiv = ui.divV([]);
  const editorView = new EditorView({
    doc: equations,
    extensions: [basicSetup, python()],
    parent: modelDiv
  });

  editorView.dom.style.overflow = 'auto';
  editorView.dom.style.height = '100%';

  let toShowInfo = true;
  const openBtn = ui.bigButton(TITLE.OPEN, async () => { await runSolverApp(equations); }, HINT.OPEN_DS);

  try {
    const ivp = getIVP(equations);
    const scriptText = getScriptLines(getIVP(equations)).join('\n');    
    const script = DG.Script.create(scriptText);
    const params = getScriptParams(ivp);    
    const call = script.prepare(params);
  
    await call.call();
    
    let solutionTable = DG.DataFrame.create();
    const solverView = DG.TableView.create(solutionTable);
    solverView.setRibbonPanels([]);
    let graph: DG.Viewer;

    /** Solve IVP */
    const solve = async (ivp: IVP, inputsPath: string) => {
      try {      
        const start = ivp.arg.initial.value;
        const finish = ivp.arg.final.value;
        const step = ivp.arg.step.value;
  
        if (start >= finish)
          return;
   
        if ((step <= 0) || (step > finish - start))
          return;
    
        const scriptText = getScriptLines(ivp).join('\n');    
        const script = DG.Script.create(scriptText);
        const params = getScriptParams(ivp);    
        const call = script.prepare(params);
    
        await call.call();
          
        solutionTable = call.outputs[DF_NAME];
        solverView.dataFrame = solutionTable;
        graph.dataFrame = solutionTable;
        browseView!.path = `${mainPath}${PATH.PARAM}${inputsPath}`;
      } catch (error) {  
        if (error instanceof Error) 
          grok.shell.error(error.message);
        else
          grok.shell.error(ERROR_MSG.SCRIPTING_ISSUE);      
      }
    }; // solve

    solutionTable = call.outputs[DF_NAME];
    graph = DG.Viewer.lineChart(solutionTable, getLineChartOptions(solutionTable.columns.names()));
    solverView.dockManager.dock(graph, DG.DOCK_TYPE.TOP);    

    const inputsDiv = await getInputsUI(ivp, solve, startingInputs);    
    
    const tabCtrl = ui.tabControl();
    const modelPane = tabCtrl.addPane(TITLE.MODEL, () => modelDiv);
    const runPane = tabCtrl.addPane(TITLE.SOLUTION, () => inputsDiv);
    tabCtrl.currentPane = runPane;
    const node = solverView.dockManager.dock(tabCtrl.root, DG.DOCK_TYPE.LEFT);
    node.container.dart.elementTitle.hidden = true;    

    editorView.dom.addEventListener('keydown', async (e) => {
      if (e.key === HOT_KEY.RUN) {
        e.stopImmediatePropagation();
        e.preventDefault();
  
        if (tabCtrl.currentPane === modelPane)         
          tabCtrl.currentPane = runPane;
      }
      else if (toShowInfo) {
        toShowInfo = false;        
        grok.shell.info(ui.divV([ERROR_MSG.OPEN_IN_DIF_STUD, ui.divH([openBtn])]));
      }
    });

    return solverView;
  }
  catch (error) {
    const view = DG.View.create();
    view.helpUrl = LINK.DIF_STUDIO_REL;
    view.append(modelDiv).style.padding = '0px';

    editorView.dom.addEventListener('keydown', async (e) => {
      if (toShowInfo) {
        toShowInfo = false;        
        grok.shell.info(ui.divV([ERROR_MSG.OPEN_IN_DIF_STUD, ui.divH([openBtn])]));
      }
    });

    return view;
  }  
} // getFilePreview

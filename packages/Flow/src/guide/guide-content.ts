/** The actual guides: 4 multi-step tutorials and a set of how-to answers.
 *
 *  Every step highlights a CONCRETE element — a specific browser item, canvas
 *  node, ribbon icon, or context-panel field (addressed by the `data-testid` /
 *  `data-*` / `data-param` attributes across the UI) — and waits for a real
 *  action before advancing. Where a value must be entered, the step copies it
 *  to the clipboard so the user can simply paste it (and the literal value is in
 *  the instruction text as a fallback).
 *
 *  The concrete targets here were verified empirically against a live server
 *  (OpenFile/AddNewColumn signatures, the `fullPath`/`name`/`expression` param
 *  rows, and that `System:DemoFiles/demog.csv` exists with an `AGE` column). */

import {
  Guide, GuideStep, GuideContext,
  byTid, bySel, byNodeFunc, byNodeType, byBrowserFunc, byParam, paramFieldSelector, socketOf,
  byFileTreeConn, byFileTreeFile,
  untilClick, untilNodeType, untilMoreNodes, untilMoreConnections, untilFuncNode,
  untilFewerNodes, untilValueContains, untilValueMatches, untilValueNonEmpty, untilNodeRightOf,
  untilNodeSelected, untilNodeSelectedOfFunc, untilMoreCollapsed,
  untilFileTreeConnExpanded, untilScrolledIntoView, untilFuncNodeWithInput, isScrolledIntoView,
  untilExists, copyToClipboard, prefillSearch, hasFuncNode, hasNodeType,
} from './guide-model';
import {createNode} from '../rete/node-factory';

const TOUR_NODE_TYPE = 'Inputs/Table Input';

const DEMO_FORMULA = '${AGE} * 12';
const NEW_COL_NAME = 'My New Column';
const SEARCH_SEL = '[data-testid="ff-browser-search"]';
const TABLE_OUTPUT_TYPE = 'Outputs/Table Output';

// The KNIME-style Files browser: scientists bring data in by grabbing a file
// from a connection, not by typing a path. We steer every "load data" step to
// the Demo connection's demog.csv (top-level, below the folders → needs a scroll).
const DEMO_CONN = 'Demo';
const DEMOG_FILE = 'demog.csv';
const DEMOG_FILE_SEL = '[data-testid="ff-files-file-demog-csv"]';

const delay = (ms: number): Promise<void> => new Promise((r) => setTimeout(r, ms));

/** Reveal the toolbox and make sure the Files pane is expanded (click its header
 *  if a previous session left it collapsed). */
const openFiles = async (ctx: GuideContext): Promise<void> => {
  ctx.host.showFunctionBrowser();
  await delay(80);
  const header = document.querySelector(
    '[data-testid="ff-browser-files"] .funcflow-section-header') as HTMLElement | null;
  if (header?.classList.contains('collapsed')) header.click();
};

/** True once demog.csv exists in the Files tree AND is scrolled into view. */
const demogVisible = (ctx: GuideContext): boolean => {
  const fileEl = byFileTreeFile(DEMOG_FILE)(ctx);
  return !!fileEl && isScrolledIntoView(fileEl);
};

/** The reusable "bring data in from the Files browser" sequence: open Files →
 *  expand the Demo connection → scroll to demog.csv → double-click / drag it.
 *  Ends with an OpenFile node already pointing at the file (no path typing).
 *  Pass `skipIf` to skip the whole sequence when data is already loaded. */
const loadDemogViaFiles = (skipIf?: (ctx: GuideContext) => boolean): GuideStep[] => {
  const steps: GuideStep[] = [
    {
      title: 'Open the Files browser',
      text: 'In the toolbox on the left, the Files pane lists your data connections — this is the ' +
        'easiest way to bring data in. Click Next.',
      setup: openFiles,
      target: byTid('browser-files'),
    },
    {
      title: 'Open the Demo connection',
      text: `Double-click the “${DEMO_CONN}” connection (highlighted) — or click its ▸ triangle — to ` +
        'expand it and list its files.',
      setup: openFiles,
      target: byFileTreeConn(DEMO_CONN),
      until: untilFileTreeConnExpanded(DEMO_CONN),
    },
    {
      title: 'Scroll to demog.csv',
      text: `The files are listed below the folders. Scroll the Files list down until “${DEMOG_FILE}” ` +
        'comes into view.',
      // Highlight the whole Files pane while scrolling (the file itself is off-screen
      // and would otherwise pulse over unrelated rows); switch to the file only once
      // it's actually in view.
      target: (ctx) => demogVisible(ctx) ? byFileTreeFile(DEMOG_FILE)(ctx) : byTid('browser-files')(ctx),
      highlights: (ctx) => [demogVisible(ctx) ? byFileTreeFile(DEMOG_FILE)(ctx) : byTid('browser-files')(ctx)],
      until: untilScrolledIntoView(DEMOG_FILE_SEL),
    },
    {
      title: 'Add it to the canvas',
      text: `Double-click “${DEMOG_FILE}” (highlighted) — or drag it onto the canvas — to drop an ` +
        'Open File node already pointing at that file. No path typing needed.',
      target: byFileTreeFile(DEMOG_FILE),
      until: untilFuncNodeWithInput('OpenFile', 'fullPath', 'demog'),
    },
  ];
  return skipIf ? steps.map((s) => ({...s, skipIf})) : steps;
};

/** Reveal the function list with a cleared search box, so the user types the
 *  query themselves (the tour gates on what they type). */
const openSearch = async (ctx: GuideContext): Promise<void> => {
  ctx.host.showFunctionBrowser();
  await delay(40);
  prefillSearch('');
};

/** Reveal the list and pre-filter it to a function (used by quick how-to's). */
const findInBrowser = (term: string) => async (ctx: GuideContext): Promise<void> => {
  ctx.host.showFunctionBrowser();
  await delay(60);
  prefillSearch(term);
};

/** Copy a value so the next step can be "paste it (Ctrl+V)". */
const putOnClipboard = (text: string) => async (): Promise<void> => {
  await copyToClipboard(text);
};

// ============================ TUTORIALS ============================

/** Flagship, hands-on: open a real dataset and compute a new column — touches
 *  the browser, the canvas, connections, the context panel (two fields), the
 *  clipboard, Run, and inspect. This is what the Start-panel "tour" launches. */
const loadDataAddColumn: Guide = {
  id: 'load-data-add-column',
  kind: 'tutorial',
  title: 'Load data and add a column',
  summary: 'Open a demo dataset, compute a new column, and run it.',
  steps: [
    {
      title: 'Let\'s build a real flow',
      text: 'In a few minutes you\'ll load a demographics dataset, compute a new column, send it to ' +
        'an output, and run the flow. Click Next to start.',
    },
    ...loadDemogViaFiles(),
    {
      title: 'Search for “Add New Column”',
      text: 'Click the search box and type add new column (or addnewcolumn) — the function that ' +
        'computes a column.',
      setup: openSearch,
      target: byTid('browser-search'),
      until: untilValueMatches(SEARCH_SEL, 'addnewcolumn'),
    },
    {
      title: 'Add the Add New Column node',
      text: 'Double-click “Add New Column” (highlighted) to add it.',
      target: byBrowserFunc('AddNewColumn'),
      until: untilFuncNode('AddNewColumn'),
    },
    {
      title: 'Move it clear of Open File',
      text: 'New nodes land in the center, so this one overlaps Open File. Drag the “Add New Column” ' +
        'node to the right until it sits clear, with room to wire them together.',
      target: byNodeFunc('AddNewColumn'),
      position: 'top',
      until: untilNodeRightOf(byNodeFunc('AddNewColumn'), byNodeFunc('OpenFile'), 220),
    },
    {
      title: 'Connect the data',
      text: 'Drag from Open File\'s result output dot (right) to Add New Column\'s table input dot ' +
        '(left). Both dots are highlighted — matching colors mean compatible types.',
      target: byNodeFunc('AddNewColumn'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx),
        socketOf(byNodeFunc('AddNewColumn'), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Select Add New Column',
      text: 'Click the “Add New Column” node so its settings open on the right.',
      target: byNodeFunc('AddNewColumn'),
      until: untilNodeSelectedOfFunc('AddNewColumn'),
    },
    {
      title: 'Name the new column',
      text: `In the panel on the right, type ${NEW_COL_NAME} into the Name field (exactly that, any case).`,
      target: byParam('name'),
      position: 'left',
      until: untilValueContains(paramFieldSelector('name'), NEW_COL_NAME.toLowerCase()),
    },
    {
      title: 'Enter the formula',
      text: `We copied a formula to your clipboard. Paste it (Ctrl+V) into the Expression field: ` +
        `${DEMO_FORMULA} — it turns the AGE column into months.`,
      setup: putOnClipboard(DEMO_FORMULA),
      target: byParam('expression'),
      position: 'left',
      until: untilValueContains(paramFieldSelector('expression'), 'age'),
    },
    {
      title: 'Search for “Table Output”',
      text: 'Now mark the result. Click the search box and type table output (or tableoutput).',
      setup: openSearch,
      target: byTid('browser-search'),
      until: untilValueMatches(SEARCH_SEL, 'tableoutput'),
    },
    {
      title: 'Add the Table Output node',
      text: 'Double-click “Table Output” (highlighted) — it marks a table as your flow\'s result.',
      target: byTid('browser-item', TABLE_OUTPUT_TYPE),
      until: untilNodeType(TABLE_OUTPUT_TYPE),
    },
    {
      title: 'Move it clear of Add New Column',
      text: 'Drag the “Table Output” node to the right of Add New Column so you can wire them up.',
      target: byNodeType(TABLE_OUTPUT_TYPE),
      position: 'top',
      until: untilNodeRightOf(byNodeType(TABLE_OUTPUT_TYPE), byNodeFunc('AddNewColumn'), 220),
    },
    {
      title: 'Connect the result',
      text: 'Drag from Add New Column\'s “table →” pass-through output dot (highlighted) to the Table ' +
        'Output\'s input dot (highlighted).',
      target: byNodeType(TABLE_OUTPUT_TYPE),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFunc('AddNewColumn'), 'output', 'table__pt')(ctx),
        socketOf(byNodeType(TABLE_OUTPUT_TYPE), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Run the flow',
      text: 'Click Run (the ▶ icon) in the ribbon. Each node lights up as it executes and reports ' +
        'its row × column counts.',
      target: byTid('ribbon', 'run'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'run')),
    },
    {
      title: 'You built a working flow! 🎉',
      text: 'Inputs → transform → output, wired and run. To peek at any node\'s data, right-click its ' +
        'output dot and choose “Run up to here & preview”.',
    },
  ],
};

const findFunctions: Guide = {
  id: 'find-functions',
  kind: 'tutorial',
  title: 'Find the right function',
  summary: 'Browse by category, search by name, and add functions.',
  steps: [
    {
      title: 'The function list',
      text: 'Every Datagrok function is here, grouped by what it does. Let\'s learn to find things ' +
        'fast. Click Next.',
      setup: (ctx) => ctx.host.showFunctionBrowser(),
      target: byTid('browser'),
    },
    {
      title: 'Group by what it does',
      text: 'This dropdown buckets functions into Data Sources, Combine Tables, Transform Tables, ' +
        'Column Operations, and more — Data Sources first. Click Next when you\'ve seen it.',
      target: byTid('browser-groupby'),
    },
    {
      title: 'Search by name',
      text: 'Type join into the search box to filter the whole catalog down to matching functions.',
      target: byTid('browser-search'),
      until: untilValueMatches(SEARCH_SEL, 'join'),
    },
    {
      title: 'Add “Join Tables”',
      text: 'Double-click “Join Tables” (highlighted) to drop it on the canvas.',
      target: byBrowserFunc('JoinTables'),
      until: untilFuncNode('JoinTables'),
    },
    {
      title: 'Let Flow suggest the next step',
      text: 'Drag from Join Tables\' result output dot (highlighted) into empty canvas. Flow pops up ' +
        'compatible next functions, common ones first — pick any one to add it, already connected.',
      target: byNodeFunc('JoinTables'),
      position: 'top',
      highlights: (ctx) => [socketOf(byNodeFunc('JoinTables'), 'output', 'result')(ctx)],
      until: untilMoreNodes(),
    },
  ],
};

const organizeCanvas: Guide = {
  id: 'organize-canvas',
  kind: 'tutorial',
  title: 'Organize your canvas',
  summary: 'Add nodes, move them apart, collapse, tidy, undo/redo, and navigate.',
  steps: [
    {
      title: 'Add a Table Input',
      text: 'Let\'s get something to arrange. Double-click “Table Input” (highlighted) to add it.',
      setup: findInBrowser('Table Input'),
      target: byTid('browser-item', 'Inputs/Table Input'),
      until: untilNodeType('Inputs/Table Input'),
    },
    {
      title: 'Add a Table Output',
      text: 'Now double-click “Table Output” (highlighted). It lands in the center, on top of the ' +
        'first node — we\'ll fix that next.',
      setup: findInBrowser('Table Output'),
      target: byTid('browser-item', TABLE_OUTPUT_TYPE),
      until: untilNodeType(TABLE_OUTPUT_TYPE),
    },
    {
      title: 'Move it apart',
      text: 'Drag the “Table Output” node to the right until it no longer overlaps the Table Input.',
      target: byNodeType(TABLE_OUTPUT_TYPE),
      position: 'top',
      until: untilNodeRightOf(byNodeType(TABLE_OUTPUT_TYPE), byNodeType('Inputs/Table Input'), 220),
    },
    {
      title: 'Collapse a node',
      text: 'Click the ▾ caret in a node\'s title bar (highlighted) to fold it to just the title. ' +
        'Click ▸ to expand it again.',
      target: bySel('.ff-node:not(.ff-node-collapsed) .ff-node-caret'),
      until: untilMoreCollapsed(),
    },
    {
      title: 'Tidy the layout',
      text: 'Click “Tidy up layout” (highlighted) to auto-arrange everything left-to-right.',
      target: byTid('ribbon', 'layout'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'layout')),
    },
    {
      title: 'Undo it',
      text: 'Changed your mind? Click Undo (highlighted) — or press Ctrl+Z — to revert the layout.',
      target: byTid('ribbon', 'undo'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'undo')),
    },
    {
      title: 'Redo it',
      text: 'Click Redo (highlighted) to put it back. Undo/redo cover every edit.',
      target: byTid('ribbon', 'redo'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'redo')),
    },
    {
      title: 'Use the overview',
      text: 'The overview (bottom-right) appears once you have nodes and shows the whole graph. ' +
        'Click or drag inside it to jump around.',
      target: byTid('minimap'),
      position: 'left',
      until: untilClick(byTid('minimap')),
    },
    {
      title: 'Zoom to fit',
      text: 'Click “Zoom to fit” (highlighted) to frame the entire flow in the viewport.',
      target: byTid('ribbon', 'zoom-fit'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'zoom-fit')),
    },
  ],
};

const reuseScript: Guide = {
  id: 'reuse-script',
  kind: 'tutorial',
  title: 'See & reuse the generated script',
  summary: 'Flow is a glass box — view, run, and save the real script.',
  steps: [
    {
      title: 'Flow is a glass box',
      text: 'Your visual flow compiles to a real, readable Datagrok script. Let\'s look at it. ' +
        'Click Next.',
    },
    {
      title: 'See the steps',
      text: 'Click “See the steps” — the 👁 (eye) icon — to view the generated script for your flow.',
      target: byTid('ribbon', 'view-script'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'view-script')),
    },
    {
      title: 'Open or export it',
      text: 'From that dialog you can copy the script, export a .js file, or open it in the Script ' +
        'editor. Close the dialog, then click Next.',
    },
    {
      title: 'Save your flow',
      text: 'Click Save (highlighted) to keep your flow as a .ffjson file you can reopen and share.',
      target: byTid('ribbon', 'save'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'save')),
    },
    {
      title: 'No black box',
      text: 'You can always see exactly what a flow does and hand the script to a colleague. Done!',
    },
  ],
};

/** Ensure a sample node exists on the canvas so the canvas / context-panel parts
 *  of the interface tour have something concrete to point at. Added silently and
 *  programmatically (no detour to hunt for it in the toolbox). */
const ensureTourNode = async (ctx: GuideContext): Promise<void> => {
  const flow = ctx.host.getFlow();
  if (!flow) return;
  if (flow.getNodes().some((n) => n.dgTypeName === TOUR_NODE_TYPE)) return;
  const node = createNode(TOUR_NODE_TYPE);
  if (node) await flow.addNodeAtCenter(node);
};

/** A guided tour of the whole interface: every toolbox pane, every ribbon group,
 *  the canvas + a node's anatomy, the overview/status bar, and the context panel.
 *  Mostly "read & click Next" steps; a couple are interactive (selecting a node).
 *  A sample node is auto-added as a prerequisite for the canvas/panel sections. */
const interfaceTour: Guide = {
  id: 'interface-tour',
  kind: 'tutorial',
  title: 'Tour the interface',
  summary: 'Meet every control: toolbox, ribbon, canvas, and context panel.',
  steps: [
    {
      title: 'A tour of Flow\'s interface',
      text: 'We\'ll walk through every part of the screen — the toolbox, the ribbon, the canvas, and ' +
        'the context panel — and what each control does. Click Next to begin.',
      setup: (ctx) => ctx.host.showFunctionBrowser(),
    },

    // ---------------- Toolbox (left) ----------------
    {
      title: 'The toolbox',
      text: 'On the left is the toolbox — your palette of building blocks. Files, queries, built-in ' +
        'nodes, and every Datagrok function live here. Drag or double-click anything to add it.',
      setup: (ctx) => ctx.host.showFunctionBrowser(),
      target: byTid('browser'),
    },
    {
      title: 'Search',
      text: 'Type here to filter the whole catalog by name — it matches both the display name ' +
        '(“Open File”) and the raw name (“OpenFile”).',
      target: byTid('browser-search'),
    },
    {
      title: 'Group by',
      text: 'Choose how functions are organized: by what they do (default), by role, by tags, or by ' +
        'package.',
      target: byTid('browser-groupby'),
    },
    {
      title: 'Files',
      text: 'A file browser of your data connections. Expand a connection, then double-click — or drag ' +
        '— a file onto the canvas to load it as an Open File node.',
      target: byTid('browser-files'),
    },
    {
      title: 'Queries',
      text: 'Database queries, grouped into sub-sections by their data connection. Double-click or drag ' +
        'one to add it as a node.',
      target: byTid('browser-queries'),
    },
    {
      title: 'Built-in building blocks',
      text: 'Below the data, collapsible sections hold the wiring blocks: Inputs, Outputs, Constants, ' +
        'Comparisons, Utilities, and Debug. Click a header to expand it.',
      target: byTid('browser-section', 'Inputs'),
    },
    {
      title: 'Functions by category',
      text: 'Then come Datagrok\'s functions, grouped by what they do — Data Sources, Combine Tables, ' +
        'Transform Tables, Column Operations, Compute Values, Visualize — Data Sources first.',
      target: byTid('browser-section', 'Data Sources'),
    },

    // ---------------- Ribbon (top) ----------------
    {
      title: 'Run',
      text: 'The ▶ Run button executes your flow with live visualization — each node lights up as it ' +
        'runs and reports its row × column counts.',
      target: byTid('ribbon', 'run'),
      position: 'bottom',
    },
    {
      title: 'Debug, Continue, Stop',
      text: 'Next to Run: Debug (🐞) runs but pauses at Breakpoint nodes, Continue (▶▶) resumes, and ' +
        'Stop (■) halts a run.',
      target: byTid('ribbon', 'debug'),
      position: 'bottom',
    },
    {
      title: 'See the script',
      text: 'The 👁 eye opens the real, readable Datagrok script your flow compiles to — copy it, export ' +
        'a .js file, or open it in the Script editor. Flow is a glass box.',
      target: byTid('ribbon', 'view-script'),
      position: 'bottom',
    },
    {
      title: 'Save & Open',
      text: 'Save (💾) downloads your flow as a .ffjson file; Open (📂, just right of it) loads one back. ' +
        'Share the file with a colleague to hand off the whole flow.',
      target: byTid('ribbon', 'save'),
      position: 'bottom',
    },
    {
      title: 'Undo & Redo',
      text: 'Undo (↶) and Redo (↷) cover every edit — adding, moving, connecting, deleting. Ctrl+Z / ' +
        'Ctrl+Shift+Z work too.',
      target: byTid('ribbon', 'undo'),
      position: 'bottom',
    },
    {
      title: 'Tidy up layout',
      text: 'This auto-arranges the whole graph left-to-right along the flow of data — a one-click clean ' +
        'layout.',
      target: byTid('ribbon', 'layout'),
      position: 'bottom',
    },
    {
      title: 'Zoom',
      text: 'Zoom in (🔍+), zoom out (🔍−), and Zoom to fit (⤢) — the last frames the entire flow in the ' +
        'viewport. You can also scroll to zoom and drag to pan.',
      target: byTid('ribbon', 'zoom-fit'),
      position: 'bottom',
    },
    {
      title: 'Show / hide the toolbox',
      text: 'Toggles the toolbox panel, giving the canvas more room when you don\'t need the catalog.',
      target: byTid('ribbon', 'toggle-browser'),
      position: 'bottom',
    },
    {
      title: 'Help & tutorials',
      text: 'The 🎓 cap (and the floating help button) opens this menu of tutorials and how-to answers ' +
        'anytime.',
      target: byTid('ribbon', 'help'),
      position: 'bottom',
    },

    // ---------------- Canvas + node anatomy (needs a node) ----------------
    {
      title: 'The canvas',
      text: 'The center is the canvas, where your flow lives. Drag empty space to pan, scroll to zoom, ' +
        'and drag a node to move it. We added a sample node so you can see its parts.',
      setup: ensureTourNode,
      target: byTid('canvas'),
      position: 'top',
    },
    {
      title: 'A node\'s title bar',
      text: 'Every node shows a colored title bar — the color encodes its role. Below it, a plain-language ' +
        'summary tells you what the node does.',
      target: bySel('.ff-node [data-testid="ff-node-title"]'),
      position: 'right',
    },
    {
      title: 'Collapse / expand',
      text: 'The ▾ caret folds a node down to just its title (and back). Handy for tidying a busy canvas.',
      target: bySel('.ff-node [data-testid="ff-node-caret"]'),
      position: 'right',
    },
    {
      title: 'Status indicator',
      text: 'This dot shows the node\'s run state — idle, running, done, or error — and turns amber when a ' +
        'required input is still missing.',
      target: bySel('.ff-node [data-testid="ff-node-status"]'),
      position: 'right',
    },
    {
      title: 'Sockets',
      text: 'The colored dots are sockets — inputs on the left, outputs on the right. Drag between two ' +
        'compatible (same-colored) dots to connect nodes. The small gray squares at the corners are ' +
        'execution-order ports.',
      target: bySel('.ff-node .ff-socket'),
      position: 'right',
    },
    {
      title: 'The overview',
      text: 'Bottom-right is the overview (minimap) — it appears once you have nodes and shows the whole ' +
        'graph. Click or drag inside it to jump around; click its header to minimize it.',
      target: byTid('minimap'),
      position: 'left',
    },
    {
      title: 'The status bar',
      text: 'Along the bottom, the status bar reports your node and connection counts and any validation ' +
        'problems with the flow.',
      target: byTid('statusbar'),
      position: 'top',
    },

    // ---------------- Context panel (right) ----------------
    {
      title: 'Open a node\'s settings',
      text: 'Click the sample node on the canvas (highlighted). Its settings open in the context panel ' +
        'on the right.',
      setup: ensureTourNode,
      target: byNodeType(TOUR_NODE_TYPE),
      until: untilNodeSelected(),
    },
    {
      title: 'The context panel',
      text: 'On the right, the context panel edits whatever node is selected — its name, type, and every ' +
        'parameter.',
      target: byTid('property-panel'),
      position: 'left',
    },
    {
      title: 'Rename a node',
      text: 'The title row lets you rename the node — the new name flows through to the generated script.',
      target: byTid('property-title-row'),
      position: 'left',
    },
    {
      title: 'Node type',
      text: 'This badge shows the node\'s kind (input, output, utility, or function).',
      target: byTid('property-type-badge'),
      position: 'left',
    },
    {
      title: 'Parameters',
      text: 'Below, each of the node\'s parameters is an editable field you can type or paste into. ' +
        'Function nodes also show a Connections pane listing how each input/output is wired.',
      target: byTid('property-content'),
      position: 'left',
    },
    {
      title: 'That\'s the whole interface! 🎉',
      text: 'Toolbox to add, ribbon to act, canvas to compose, context panel to configure. Try the ' +
        '“Load data and add a column” tutorial next to build a real flow.',
    },
  ],
};

export const TUTORIALS: Guide[] = [loadDataAddColumn, findFunctions, organizeCanvas, reuseScript, interfaceTour];

// ============================ HOW-TO QUESTIONS ============================

function q(id: string, title: string, steps: GuideStep | GuideStep[]): Guide {
  return {id, kind: 'question', title, summary: title, steps: Array.isArray(steps) ? steps : [steps]};
}

// ---- prerequisite step builders (skipped when already satisfied) ----

/** Ensure a DG-function node exists (adds it via the browser if missing). */
const ensureFuncNode = (funcName: string, friendly: string): GuideStep => ({
  title: `First, add “${friendly}”`,
  text: `This needs a “${friendly}” node. Double-click it (highlighted) to add one.`,
  skipIf: hasFuncNode(funcName),
  setup: findInBrowser(funcName),
  target: byBrowserFunc(funcName),
  until: untilFuncNode(funcName),
});

/** Ensure a built-in node (e.g. Table Input/Output) exists. */
const ensureBuiltin = (typeName: string, friendly: string): GuideStep => ({
  title: `First, add “${friendly}”`,
  text: `This needs a “${friendly}” node. Double-click it (highlighted) to add one.`,
  skipIf: hasNodeType(typeName),
  setup: findInBrowser(friendly),
  target: byTid('browser-item', typeName),
  until: untilNodeType(typeName),
});

/** True once some Open File node has a non-empty path set. */
const openFileHasPath = (ctx: GuideContext): boolean =>
  (ctx.host.getFlow()?.getNodes() ?? []).some((n) =>
    (n.dgFuncName ?? '').toLowerCase().includes('openfile') &&
    !!String((n.inputValues ?? {})['fullPath'] ?? '').trim());

export const QUESTIONS: Guide[] = [
  q('how-add-function', 'How do I add a function?', {
    title: 'Add a function',
    text: 'Open the function list on the left, type in the search box (highlighted) or browse a ' +
      'category, then double-click a result — or drag it onto the canvas.',
    setup: (ctx) => ctx.host.showFunctionBrowser(),
    target: byTid('browser-search'),
    until: untilMoreNodes(),
  }),
  q('how-add-data', 'How do I bring data in?', loadDemogViaFiles()),
  q('how-add-column', 'How do I add a calculated column?', [
    ensureFuncNode('AddNewColumn', 'Add New Column'),
    {
      title: 'Open its settings',
      text: 'Click the “Add New Column” node (highlighted) so its fields open on the right. ' +
        '(Wire a table into its table input so the formula has data to work on.)',
      target: byNodeFunc('AddNewColumn'),
      until: untilNodeSelectedOfFunc('AddNewColumn'),
    },
    {
      title: 'Name the column',
      text: 'Type a name into the Name field (highlighted).',
      target: byParam('name'),
      position: 'left',
      until: untilValueNonEmpty(paramFieldSelector('name')),
    },
    {
      title: 'Write the formula',
      text: 'Type an expression into the Expression field — for example ${AGE} * 12, referencing ' +
        'columns of the incoming table with ${ColumnName}.',
      target: byParam('expression'),
      position: 'left',
      until: untilValueNonEmpty(paramFieldSelector('expression')),
    },
  ]),
  q('how-set-param', 'How do I edit a node\'s settings?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'Select the node',
      text: 'Click a node (the highlighted Table Input). Its parameters open in the panel on the ' +
        'right, each an editable field you can type or paste into.',
      target: byNodeType('Inputs/Table Input'),
      until: untilNodeSelected(),
    },
  ]),
  q('how-connect', 'How do I connect two nodes?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    ensureBuiltin('Outputs/Table Output', 'Table Output'),
    {
      title: 'Drag between the dots',
      text: 'Drag from the Table Input\'s output dot (right, highlighted) to the Table Output\'s ' +
        'input dot (left, highlighted). Matching colors mean compatible types.',
      target: byNodeType('Outputs/Table Output'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeType('Inputs/Table Input'), 'output', 'table')(ctx),
        socketOf(byNodeType('Outputs/Table Output'), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
  ]),
  q('how-run', 'How do I run my flow?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'Press Run',
      text: 'Click Run (the ▶ icon, highlighted) in the ribbon. Nodes light up as they execute; each ' +
        'reports its row × column counts underneath.',
      target: byTid('ribbon', 'run'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'run')),
    },
  ]),
  q('how-preview', 'How do I preview a node\'s data?', [
    ...loadDemogViaFiles(openFileHasPath),
    {
      title: 'Run up to here & preview',
      text: 'Right-click Open File\'s result output dot (highlighted) and choose “Run up to here & ' +
        'preview”. Flow runs just that slice and shows the data — no full run or output node needed.',
      target: byNodeFunc('OpenFile'),
      position: 'top',
      highlights: (ctx) => [socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx)],
      until: untilExists('[data-testid="ff-port-preview"], [data-testid="ff-output-panel"]'),
    },
  ]),
  q('how-delete', 'How do I delete a node?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'Delete it',
      text: 'Click the node (highlighted) to select it, then press Delete or Backspace — or ' +
        'right-click it and choose Delete.',
      target: byNodeType('Inputs/Table Input'),
      until: untilFewerNodes(),
    },
  ]),
  q('how-layout', 'How do I tidy up the layout?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    ensureBuiltin('Outputs/Table Output', 'Table Output'),
    {
      title: 'Tidy the layout',
      text: 'Click “Tidy up layout” (highlighted) in the ribbon to auto-arrange every node ' +
        'left-to-right along the flow of data.',
      target: byTid('ribbon', 'layout'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'layout')),
    },
  ]),
  q('how-categories', 'How do I find functions by category?', {
    title: 'Browse by category',
    text: 'Open the function list and use the “Group by” dropdown (highlighted). Functions bucket by ' +
      'what they do — Data Sources, Combine Tables, Transform Tables, Column Operations, Compute ' +
      'Values, Visualize — with Data Sources first. Click Finish when you\'ve seen it.',
    setup: (ctx) => ctx.host.showFunctionBrowser(),
    target: byTid('browser-groupby'),
  }),
  q('how-collapse', 'How do I collapse a node?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'Click the caret',
      text: 'Click the ▾ caret in a node\'s title bar (highlighted) to fold it down to just the ' +
        'title; click ▸ to expand it again. (The status dot is only an indicator — it no longer ' +
        'collapses.)',
      target: bySel('.ff-node:not(.ff-node-collapsed) .ff-node-caret'),
      until: untilMoreCollapsed(),
    },
  ]),
  q('how-undo', 'How do I undo a change?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'Undo',
      text: 'Click Undo (highlighted) in the ribbon, or press Ctrl+Z — this reverts the last edit ' +
        '(e.g. adding that node). Redo sits right next to it.',
      target: byTid('ribbon', 'undo'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'undo')),
    },
  ]),
  q('how-save', 'How do I save or share my flow?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'Save / share',
      text: 'Click Save (highlighted) to download a .ffjson file. Open it later from anywhere in ' +
        'Datagrok — or hand the file to a colleague — to restore the flow exactly.',
      target: byTid('ribbon', 'save'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'save')),
    },
  ]),
  q('how-view-script', 'How do I see the generated script?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'See the steps',
      text: 'Click “See the steps” — the 👁 (eye) icon, highlighted. Your visual flow compiles to a ' +
        'real, editable Datagrok script you can copy, export, or open in the Script editor.',
      target: byTid('ribbon', 'view-script'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'view-script')),
    },
  ]),
  q('how-open', 'How do I open a saved flow?', {
    title: 'Open a flow',
    text: 'Click the Open (folder) icon (highlighted) in the ribbon and pick a .ffjson file.',
    target: byTid('ribbon', 'open'),
    position: 'bottom',
    until: untilClick(byTid('ribbon', 'open')),
  }),
  q('how-navigate', 'How do I navigate a large flow?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'Use the overview',
      text: 'The overview (bottom-right, highlighted — it appears once you have nodes) shows the ' +
        'whole graph. Click or drag inside it to move; scroll the canvas to zoom; or use Zoom to ' +
        'fit in the ribbon.',
      target: byTid('minimap'),
      position: 'left',
      until: untilClick(byTid('minimap')),
    },
  ]),
];

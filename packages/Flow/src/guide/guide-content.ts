/** The actual guides: 6 multi-step tutorials and a set of how-to answers.
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
  byTid, bySel, byNodeFunc, byNodeFuncNth, byNodeType, byBrowserFunc, byParam, paramFieldSelector, socketOf,
  byFileTreeConn, byFileTreeFile, preferDialog, openDialogEl,
  untilClick, untilNodeType, untilMoreNodes, untilMoreConnections, untilFuncNode,
  untilFewerNodes, untilValueContains, untilValueMatches, untilValueNonEmpty, untilNodeRightOf,
  nodeIsRightOf, nodeIsApart, untilNodeApart, untilNodeMovedBy, untilVisible, untilNodeOfTypeSelected,
  untilNodeSelected, untilNodeSelectedOfFunc, untilMoreCollapsed, untilColumnCountAtLeast,
  untilFileTreeConnExpanded, untilScrolledIntoView, untilFuncNodeWithInput, isScrolledIntoView,
  untilExists, copyToClipboard, prefillSearch, hasFuncNode, hasNodeType, poll, el,
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
// System:DemoFiles/demog.csv has 11 columns (USUBJID, AGE, SEX, RACE, DIS_POP,
// HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY) and USUBJID is its key —
// the join how-to gates on selecting USUBJID and on picking every column.
const DEMOG_KEY = 'USUBJID';
const DEMOG_COLUMN_COUNT = 11;

const delay = (ms: number): Promise<void> => new Promise((r) => setTimeout(r, ms));

/** Frame the whole graph before a wiring step — with several nodes added one
 *  after another, the drag SOURCE (e.g. Open File) may have drifted off-screen.
 *  The small delay lets the previous step's own add-pan finish first — the
 *  step gate fires the moment the node exists, before its pan settles, and a
 *  late pan would override the fit. */
const fitGraph = async (ctx: GuideContext): Promise<void> => {
  await delay(350);
  try {
    await ctx.host.getFlow()?.zoomToFit();
  } catch {/* editor not ready */}
};

/** True while a node running `funcName` is selected (its settings are already
 *  open) — lets "click the node" steps skip silently instead of flashing. */
const funcNodeSelected = (funcName: string) => (): boolean => {
  const lc = funcName.toLowerCase();
  return (Array.from(document.querySelectorAll('.ff-node[data-selected="true"]')) as HTMLElement[])
    .some((n) => (n.dataset.func ?? '').toLowerCase().includes(lc));
};

/** Reveal the toolbox and make sure the Files tab is active (a previous
 *  session may have left another top tab selected). */
const openFiles = async (ctx: GuideContext): Promise<void> => {
  ctx.host.showFunctionBrowser();
  await delay(80);
  ctx.host.showToolboxTab('Files');
};

/** True once demog.csv exists in the Files tree AND is scrolled into view. */
const demogVisible = (ctx: GuideContext): boolean => {
  const fileEl = byFileTreeFile(DEMOG_FILE)(ctx);
  return !!fileEl && isScrolledIntoView(fileEl);
};

/** True once the Demo connection row is expanded (its files are listed). */
const demoConnExpanded = (ctx: GuideContext): boolean =>
  !!byFileTreeConn(DEMO_CONN)(ctx)?.querySelector('.d4-tree-view-tri-expanded');

/** The reusable "bring data in from the Files browser" sequence: open Files →
 *  expand the Demo connection → scroll to demog.csv → double-click / drag it.
 *  Ends with an OpenFile node already pointing at the file (no path typing).
 *  Pass `skipIf` to skip the whole sequence when data is already loaded.
 *  Steps whose action is already done (connection expanded, file in view) skip
 *  silently instead of flashing a card that instantly self-advances. */
const loadDemogViaFiles = (skipIf?: (ctx: GuideContext) => boolean): GuideStep[] => {
  const steps: GuideStep[] = [
    {
      title: 'Open the Files browser',
      text: 'In the toolbox on the left, the Files tab lists your data connections — this is the ' +
        'easiest way to bring data in. Click Next.',
      setup: openFiles,
      target: byTid('browser-files'),
    },
    {
      title: 'Open the Demo connection',
      text: `Double-click the “${DEMO_CONN}” connection (highlighted) — or click its ▸ triangle — to ` +
        'expand it and list its files.',
      skipIf: demoConnExpanded,
      setup: openFiles,
      target: byFileTreeConn(DEMO_CONN),
      until: untilFileTreeConnExpanded(DEMO_CONN),
    },
    {
      title: 'Scroll to demog.csv',
      text: `The files are listed below the folders. Scroll the Files list down until “${DEMOG_FILE}” ` +
        'comes into view.',
      skipIf: demogVisible,
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
  // An outer skipIf ("data already loaded") composes with each step's own.
  return skipIf ?
    steps.map((s) => ({...s, skipIf: (ctx: GuideContext) => skipIf(ctx) || (s.skipIf?.(ctx) ?? false)})) :
    steps;
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
      // Nodes land in free space now, so this repair step almost always skips.
      // The skipIf treats "node not there yet" as skipped too — otherwise the
      // step-count estimate includes it at guide start and the visible total
      // shrinks mid-run when it self-skips.
      title: 'Move it clear of Open File',
      text: 'Give it room: drag the “Add New Column” node by its title bar to the right of Open ' +
        'File, so there\'s space to wire them together.',
      skipIf: (ctx) => !byNodeFunc('AddNewColumn')(ctx) ||
        nodeIsApart(byNodeFunc('AddNewColumn'))(ctx),
      target: byNodeFunc('AddNewColumn'),
      position: 'top',
      until: untilNodeApart(byNodeFunc('AddNewColumn')),
    },
    {
      title: 'Connect the data',
      text: 'Drag from Open File\'s result output dot (right) to Add New Column\'s table input dot ' +
        '(left). Both dots are highlighted — matching colors mean compatible types. A line between ' +
        'the nodes means it worked.',
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
      skipIf: funcNodeSelected('AddNewColumn'),
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
      text: 'Double-click “Table Output” (highlighted) — it marks a table as your flow\'s result. ' +
        'Outputs don\'t float on the canvas: it docks into the OUTPUTS strip on the right edge.',
      target: byTid('browser-item', TABLE_OUTPUT_TYPE),
      until: untilNodeType(TABLE_OUTPUT_TYPE),
    },
    {
      title: 'Connect the result',
      text: 'Drag from Add New Column\'s “table →” pass-through output dot (highlighted) to the Table ' +
        'Output\'s dot at the very right edge, in the OUTPUTS strip (also highlighted).',
      target: byNodeType(TABLE_OUTPUT_TYPE),
      position: 'left',
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
      title: 'Open the result tab',
      text: 'Look at the bottom status bar — your result got its own tab next to “Canvas” ' +
        '(highlighted). Its dot turns green when the result is ready. Click it to open the result ' +
        'as a full table view.',
      // A prior step may have been skipped, leaving the tab empty forever — target
      // and gate on any output tab so the guide never dead-ends (skip-tolerance).
      target: (ctx) => bySel('.ff-view-tab[data-state="ready"]')(ctx) ?? bySel('.ff-view-tab[data-node-id]')(ctx),
      position: 'top',
      until: untilExists('.ff-view-tab[data-node-id][data-active="true"]'),
    },
    {
      title: 'A real table view',
      text: 'This is a full Datagrok table view. Scroll right to find “My New Column” — the column ' +
        'your formula computed. Add viewers, filter, reorder — whatever you arrange here is saved ' +
        'with the flow. Click “Canvas” (highlighted) to go back to the graph.',
      target: byTid('view-tab', 'canvas'),
      position: 'top',
      until: untilExists('.ff-view-tab[data-param="canvas"][data-active="true"]'),
    },
    {
      title: 'You built a working flow! 🎉',
      text: 'Data in → transform → result out, wired and run. Click any completed node to see its ' +
        'data in the bottom panel. Tip: the ⚡ bolt makes the flow rerun by itself after every ' +
        'change. When you\'re ready to share, click Save.',
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
      text: 'The FUNCTIONS list in the lower half of the left panel holds every Datagrok function, ' +
        'grouped by what it does. Let\'s learn to find things fast. Click Next.',
      setup: (ctx) => ctx.host.showFunctionBrowser(),
      target: byTid('browser-tree'),
    },
    {
      title: 'Group by what it does',
      text: 'The “by: …” button switches how the list is organized — by what functions do (Data ' +
        'Sources, Combine Tables, Transform Tables, …), by role, tags, or package. Click Next.',
      target: byTid('browser-groupby'),
    },
    {
      title: 'Search by name',
      text: 'Type "join" into the search box — the list below narrows to matching functions as ' +
        'you type.',
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
      text: 'Drag from the small circle next to “Result” on the Join Tables node (highlighted) into ' +
        'empty canvas. Flow pops up compatible next functions, common ones first — pick any one to ' +
        'add it, already connected.',
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
      title: 'Add a String Input',
      text: 'Now double-click “String Input” (highlighted). It lands next to the first node, in ' +
        'free space.',
      setup: findInBrowser('String Input'),
      target: byTid('browser-item', 'Inputs/String Input'),
      until: untilNodeType('Inputs/String Input'),
    },
    {
      title: 'Move a node',
      text: 'Drag the “String Input” node (highlighted) anywhere — grab it by its title bar. Nodes ' +
        'stay wherever you put them.',
      target: byNodeType('Inputs/String Input'),
      position: 'top',
      until: untilNodeMovedBy(byNodeType('Inputs/String Input'), 60),
    },
    {
      title: 'Collapse a node',
      text: 'Click the ▾ caret in a node\'s title bar (highlighted) to fold it to just the title — ' +
        'it shrinks to a single row. Click ▸ to expand it again.',
      target: bySel('.ff-node:not(.ff-node-collapsed) .ff-node-caret'),
      until: untilMoreCollapsed(),
    },
    {
      title: 'Tidy the layout',
      text: 'Click “Tidy up layout” (highlighted) — the nodes snap into a neat left-to-right ' +
        'arrangement.',
      target: byTid('ribbon', 'layout'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'layout')),
    },
    {
      title: 'Undo it',
      text: 'Changed your mind? Click Undo (highlighted) — or press Ctrl+Z — and the nodes jump ' +
        'back to where you had them.',
      target: byTid('ribbon', 'undo'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'undo')),
    },
    {
      title: 'Redo it',
      text: 'Click Redo (highlighted) to apply the tidy layout again. Undo/redo cover every edit.',
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
      // An empty canvas has no script to show and Save stays disabled — make
      // sure there's at least one node before pointing at those buttons.
      title: 'First, put something on the canvas',
      text: 'An empty flow has nothing to compile. Double-click “Table Input” (highlighted) to add ' +
        'a node — any node makes a script.',
      skipIf: (ctx) => (ctx.host.getFlow()?.getNodeCount() ?? 0) > 0,
      setup: findInBrowser('Table Input'),
      target: byTid('browser-item', 'Inputs/Table Input'),
      until: untilNodeType('Inputs/Table Input'),
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
      text: 'Click Save (highlighted) and give your flow a name — it\'s stored on the platform so ' +
        'you can reopen it from anywhere. (To hand someone a file instead: Flow → Export .ffjson.)',
      target: byTid('ribbon', 'save'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'save')),
    },
    {
      title: 'No black box',
      text: 'You can always see exactly what a flow does and hand the script to a colleague. Saved ' +
        'flows also appear in the toolbox\'s Workflows tab (highlighted) — drop one into another ' +
        'flow like any function. Done!',
      setup: (ctx) => ctx.host.showToolboxTab('Workflows'),
      target: byTid('browser-tab', 'Workflows'),
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
      text: 'Type here to filter everything below — functions, queries, workflows, and favorites. ' +
        'Each tab shows how many of its items match.',
      target: byTid('browser-search'),
    },
    {
      title: 'Files',
      text: 'The tabs on top hold your data collections. Files is a browser of your data connections — ' +
        'expand a connection, then double-click or drag a file onto the canvas to load it.',
      setup: (ctx) => ctx.host.showToolboxTab('Files'),
      target: byTid('browser-files'),
    },
    {
      title: 'Queries',
      text: 'The Queries tab lists database queries, grouped by their data connection. Double-click or ' +
        'drag one to add it as a node.',
      setup: (ctx) => ctx.host.showToolboxTab('Queries'),
      target: byTid('browser-queries'),
    },
    {
      title: 'Workflows',
      text: 'The Workflows tab holds flows you saved — each one can be reused as a single node inside ' +
        'another flow.',
      setup: (ctx) => ctx.host.showToolboxTab('Workflows'),
      target: byTid('browser-workflows'),
    },
    {
      title: 'Favorites',
      text: 'Hover any node in the toolbox and click its ★ to pin it to the Favorites tab — your ' +
        'personal shortlist of go-to steps.',
      setup: (ctx) => ctx.host.showToolboxTab('Favorites'),
      target: byTid('browser-favorites'),
    },
    {
      title: 'Group by',
      text: 'The Functions list below is grouped by what each function does (default). Click ' +
        '“by: …” to organize it by role, tags, or package instead.',
      target: byTid('browser-groupby'),
    },
    {
      title: 'Functions by category',
      text: 'The FUNCTIONS list groups everything by what it does — categories like Data Sources, ' +
        'Combine Tables, Transform Tables, and Visualize. Click a header to expand it.',
      target: byTid('browser-section', 'Data Sources'),
    },
    {
      title: 'Built-in building blocks',
      text: 'Further down, collapsible sections hold the wiring blocks: Inputs, Outputs, Constants, ' +
        'Utilities, and Debug.',
      target: byTid('browser-section', 'Inputs'),
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
      title: 'Autorun',
      text: 'The ⚡ bolt toggles autorun: after any change — a new wire, an edited parameter — the ' +
        'flow reruns by itself, and only the nodes your change affected. It\'s grey when off and ' +
        'lights up when on.',
      target: byTid('ribbon', 'autorun'),
      position: 'bottom',
    },
    {
      title: 'See the script',
      text: 'The 👁 eye shows the recipe behind your flow — the exact steps, as a script you can ' +
        'copy or share. Nothing is hidden.',
      target: byTid('ribbon', 'view-script'),
      position: 'bottom',
    },
    {
      title: 'Save & Open',
      text: 'Save stores your flow on the platform with a name, so you can reopen it from anywhere ' +
        '— saved flows also appear in the toolbox\'s Workflows tab for reuse. Open (📂) loads a ' +
        'saved flow back.',
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
      text: 'The 🎓 cap opens this menu of tutorials and how-to answers anytime.',
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
      text: 'This dot shows whether the node is idle, running, done, or failed — amber means it\'s ' +
        'still waiting for an input. When you change something, only the affected nodes rerun.',
      target: bySel('.ff-node [data-testid="ff-node-status"]'),
      position: 'bottom',
      avoid: (ctx) => [bySel('.ff-node')(ctx)],
    },
    {
      title: 'Sockets',
      text: 'The colored dots are sockets — inputs on the left, outputs on the right. Drag between two ' +
        'compatible (same-colored) dots to connect nodes.',
      // The first .ff-socket in DOM order is a hover-only order port hidden
      // under the title bar — highlight a REAL data socket in a socket row.
      target: bySel('.ff-node .ff-socket-row .ff-socket'),
      position: 'bottom',
      avoid: (ctx) => [bySel('.ff-node')(ctx)],
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
    {
      title: 'Result tabs',
      text: 'At the left of the status bar live the view tabs: “Canvas” is the graph, and after a ' +
        'run every table output gets its own tab — green dot means ready, amber means out of date.',
      target: byTid('view-tab', 'canvas'),
      position: 'top',
    },
    {
      title: 'The output panel',
      text: 'Outputs dock along this right-edge strip. After you run the flow, an output panel also ' +
        'opens just above the status bar — click any completed node to see its data there: tables ' +
        'as real grids, plots as live charts.',
      target: byTid('output-strip'),
      position: 'left',
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
      text: 'The title row renames the node, and the badge next to it shows the node\'s kind ' +
        '(input, output, utility, or function).',
      target: byTid('property-title-row'),
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

// ============================ HOW-TO QUESTIONS ============================

function q(id: string, title: string, steps: GuideStep | GuideStep[]): Guide {
  return {id, kind: 'question', title, summary: title, steps: Array.isArray(steps) ? steps : [steps]};
}

// ---- prerequisite step builders (skipped when already satisfied) ----

/** Ensure a DG-function node exists (adds it via the browser if missing).
 *  Searches by the friendly name — that's what a user would type. */
const ensureFuncNode = (funcName: string, friendly: string, opts: {title?: string; text?: string} = {},
): GuideStep => ({
  title: opts.title ?? `First, add “${friendly}”`,
  text: opts.text ?? `This needs a “${friendly}” node to demonstrate on. Double-click it ` +
    '(highlighted) to add one.',
  skipIf: hasFuncNode(funcName),
  setup: findInBrowser(friendly),
  target: byBrowserFunc(funcName),
  until: untilFuncNode(funcName),
});

/** Ensure a built-in node (e.g. Table Input/Output) exists. */
const ensureBuiltin = (typeName: string, friendly: string, opts: {title?: string; text?: string} = {},
): GuideStep => ({
  title: opts.title ?? `First, add “${friendly}”`,
  text: opts.text ?? `This needs a “${friendly}” node to demonstrate on. Double-click it ` +
    '(highlighted) to add one.',
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

/** True once some Add New Column node has a table wired into it. */
const ancTableConnected = (ctx: GuideContext): boolean => {
  const flow = ctx.host.getFlow();
  if (!flow) return false;
  return flow.getNodes().some((n) =>
    (n.dgFuncName ?? '').toLowerCase().includes('addnewcolumn') && flow.isInputConnected(n.id, 'table'));
};

/** True once some Table Output node has a table wired into it. */
const tableOutputConnected = (ctx: GuideContext): boolean => {
  const flow = ctx.host.getFlow();
  if (!flow) return false;
  return flow.getNodes().some((n) =>
    n.dgTypeName === TABLE_OUTPUT_TYPE && flow.isInputConnected(n.id, 'table'));
};

/** Dashboard publishing end-to-end: run a flow, explore the result tab, then
 *  Save → Create dashboard → the platform's standard Save-project dialog.
 *  Defined here (not with the other tutorials) because it reuses the
 *  prerequisite builders above. */
const publishDashboard: Guide = {
  id: 'publish-dashboard',
  kind: 'tutorial',
  title: 'Publish your results as a dashboard',
  summary: 'Run a flow, arrange the result view, and publish it as a project.',
  steps: [
    {
      title: 'From flow to dashboard',
      text: 'Any flow with table outputs can publish them as a Datagrok dashboard — a project ' +
        'colleagues open without ever seeing the graph. Let\'s build one. Click Next.',
    },
    ...loadDemogViaFiles(openFileHasPath),
    ensureBuiltin(TABLE_OUTPUT_TYPE, 'Table Output'),
    {
      title: 'Connect the table',
      text: 'Drag from Open File\'s result output dot (highlighted) to the Table Output\'s dot at ' +
        'the very right edge — outputs dock in the OUTPUTS strip there (also highlighted).',
      skipIf: tableOutputConnected,
      target: byNodeType(TABLE_OUTPUT_TYPE),
      position: 'left',
      highlights: (ctx) => [
        socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx),
        socketOf(byNodeType(TABLE_OUTPUT_TYPE), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Run the flow',
      text: 'Click Run (the ▶ icon) so the output gets a value to publish.',
      target: byTid('ribbon', 'run'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'run')),
    },
    {
      title: 'Open the result tab',
      text: 'In the bottom status bar, a tab named after your table gets a green dot when the run ' +
        'finishes. Click it (highlighted) to open the result as a full table view. (Empty dot? ' +
        'The output isn\'t wired in or the run didn\'t reach it — the tab tells you what to do.)',
      // Skip-tolerant: highlight and gate on any output tab, so a skipped connect
      // step upstream can't dead-end the guide here.
      target: (ctx) => bySel('.ff-view-tab[data-state="ready"]')(ctx) ?? bySel('.ff-view-tab[data-node-id]')(ctx),
      position: 'top',
      until: untilExists('.ff-view-tab[data-node-id][data-active="true"]'),
    },
    {
      title: 'Make it look like a dashboard',
      text: 'This is a real table view — the panel on the left lists viewers you can add, and you ' +
        'can reorder columns or set up filters right here. Whatever you arrange ships with the ' +
        'dashboard, even if you never open this tab again. Click Next when it looks the way ' +
        'colleagues should see it.',
    },
    {
      title: 'Save & publish',
      text: 'Click Save (highlighted) in the ribbon.',
      target: byTid('ribbon', 'save'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'save')),
    },
    {
      title: 'The Save dialog',
      text: 'Name the flow. Below, the Dashboard section lists your computed tables — keep ' +
        '“Create dashboard” checked and click Save. (No tables listed? The flow hasn\'t run — the ' +
        'dialog offers a Run button right there. A name clash? Pick another name.)',
      target: preferDialog(bySel('.ff-save-dash')),
      position: 'left',
      // Advance only when the Save dialog is actually done (saved or closed) —
      // never narrate the NEXT dialog while this one is still on screen.
      until: (ctx) => poll(() => el('.ff-save-dash') == null, ctx.signal),
    },
    {
      title: 'The Save-project dialog',
      text: 'Next, the platform\'s Save-project dialog opens with your output tables and layouts. ' +
        'Name the project and click OK — the dashboard will re-run your flow when opened, so it ' +
        'always shows fresh data.',
      target: preferDialog(bySel('.ff-save-dash')),
      position: 'left',
    },
    {
      title: 'Published! 🎉',
      text: 'Once the project is saved, your dashboard is a regular Datagrok project — find it in ' +
        'Browse > Projects and share it like any other. The flow remembers it: the next Save ' +
        'updates the same project instead of creating a new one.',
    },
  ],
};

export const TUTORIALS: Guide[] =
  [loadDataAddColumn, findFunctions, organizeCanvas, reuseScript, publishDashboard, interfaceTour];

export const QUESTIONS: Guide[] = [
  q('how-add-function', 'How do I add a function?', {
    title: 'Add a function',
    text: 'Open the function list on the left, type in the search box (highlighted) or browse a ' +
      'category, then double-click a result — or drag it onto the canvas.',
    setup: (ctx) => ctx.host.showFunctionBrowser(),
    target: byTid('browser-search'),
    until: untilMoreNodes(),
  }),
  q('how-add-data', 'How do I bring data in?', [
    ...loadDemogViaFiles(),
    {
      title: 'Your own files work too',
      text: 'To load a file from your computer, click the folder icon at the top of the Files tab ' +
        '(highlighted) — or simply drag the file from your computer onto the canvas. It becomes a ' +
        'node just the same, and is stored with the flow when you save.',
      setup: openFiles,
      target: byTid('browser-upload'),
    },
  ]),
  q('how-add-column', 'How do I add a calculated column?', [
    // The formula references a real column, so the answer builds on real data
    // (all skipped when the flow already has it).
    ...loadDemogViaFiles(openFileHasPath),
    ensureFuncNode('AddNewColumn', 'Add New Column'),
    {
      title: 'Wire the table in',
      text: 'Drag from Open File\'s result output dot (highlighted) to Add New Column\'s table ' +
        'input dot (highlighted) — the formula needs data to work on.',
      skipIf: ancTableConnected,
      target: byNodeFunc('AddNewColumn'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx),
        socketOf(byNodeFunc('AddNewColumn'), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Open its settings',
      text: 'Click the “Add New Column” node (highlighted) so its fields open on the right.',
      skipIf: funcNodeSelected('AddNewColumn'),
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
      text: 'Type an expression into the Expression field — for example ${AGE} * 12. ' +
        '${ColumnName} refers to a column of the incoming table. Run the flow to see the new ' +
        'column appear.',
      target: byParam('expression'),
      position: 'left',
      until: untilValueNonEmpty(paramFieldSelector('expression')),
    },
  ]),
  q('how-upload-file', 'How do I use a file from my computer?', {
    title: 'Open a local file',
    text: 'Click the folder icon at the top of the Files tab (highlighted) and pick a file — or ' +
      'simply drag one from your computer onto the canvas. Either way it becomes a node, and the ' +
      'file is stored with the flow when you save, so colleagues can rerun it.',
    setup: openFiles,
    target: byTid('browser-upload'),
  }),
  q('how-group-nodes', 'How do I group nodes?', {
    title: 'Group nodes into a frame',
    text: 'Select several nodes (drag a box around them, or Ctrl-click one by one), then press ' +
      'Ctrl+G — they get a titled frame that moves as one. The frame\'s ▾ caret folds the whole ' +
      'group into a single card; Ctrl+Shift+G ungroups. Grouping is purely visual — it never ' +
      'changes what the flow computes.',
  }),
  q('how-set-param', 'How do I edit a node\'s settings?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    {
      title: 'Select the node',
      text: 'Click the highlighted node. Its settings open in the panel on the right — each ' +
        'parameter is an editable field you can type or paste into.',
      target: byNodeType('Inputs/Table Input'),
      until: untilNodeOfTypeSelected('Inputs/Table Input'),
    },
  ]),
  q('how-connect', 'How do I connect two nodes?', [
    ensureBuiltin('Inputs/Table Input', 'Table Input'),
    ensureBuiltin('Outputs/Table Output', 'Table Output', {
      title: 'Also add “Table Output”',
      text: 'Now double-click “Table Output” (highlighted). Outputs don\'t float on the canvas — ' +
        'it docks into the OUTPUTS strip on the right edge.',
    }),
    {
      title: 'Drag between the dots',
      text: 'Drag from the Table Input\'s output dot (highlighted) to the Table Output\'s dot at ' +
        'the very right edge — outputs dock in the OUTPUTS strip there (also highlighted). ' +
        'Matching colors mean compatible types.',
      target: byNodeType('Outputs/Table Output'),
      position: 'left',
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
      text: 'Click Run (the ▶ icon, highlighted) at the top. Nodes light up as they execute; each ' +
        'reports its row × column counts underneath.',
      target: byTid('ribbon', 'run'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'run')),
    },
    {
      title: 'If Datagrok asks for values',
      text: 'A flow with unset inputs (like a Table Input with no table wired in) asks for them in ' +
        'a dialog when you press Run — pick or type the values and click OK, and the run proceeds.',
      target: preferDialog(byTid('ribbon', 'run')),
      position: 'bottom',
    },
  ]),
  q('how-autorun', 'How do I rerun automatically after every change?', [
    {
      title: 'Toggle Autorun',
      text: 'Click the ⚡ bolt (highlighted) in the ribbon. While it\'s on (colored), any change — a ' +
        'new wire, an edited parameter — reruns the flow by itself after a short pause. Click it ' +
        'again anytime to turn it off.',
      target: byTid('ribbon', 'autorun'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'autorun')),
    },
    {
      title: 'Only what changed',
      text: 'Autorun is smart about it: only the steps affected by your change are re-run — ' +
        'earlier results are reused. Flows that would ask for input values are left alone — run ' +
        'those with ▶.',
    },
  ]),
  q('how-out-of-date', 'Why do nodes say “Out of date”?', [
    {
      title: '“Out of date” = a change affects this node',
      text: 'When you edit a setting or rewire a connection, that node and every step after it ' +
        'lose their last result — they show “Out of date”. Steps before it keep theirs. ' +
        'Click Next.',
    },
    {
      title: 'Rerun to refresh',
      text: 'Rerun with ▶ Run (highlighted) — or toggle the ⚡ bolt next to it so the flow reruns ' +
        'by itself after every change. Right-clicking one node offers “Rerun this node only”.',
      target: byTid('ribbon', 'run'),
      position: 'bottom',
    },
    {
      title: 'The result tabs show it too',
      text: 'The tabs at the bottom (highlighted) mirror the state — an amber dot means the table ' +
        'you see is from the previous run; green means fresh.',
      target: (ctx) => bySel('.ff-view-tab[data-node-id]')(ctx) ?? byTid('view-tab', 'canvas')(ctx),
      position: 'top',
    },
  ]),
  q('how-func-editor', 'How do I edit parameters in the function\'s own dialog?', [
    ...loadDemogViaFiles(openFileHasPath),
    ensureFuncNode('AddNewColumn', 'Add New Column'),
    {
      // A freshly added node half-covers the previous one, hiding the very
      // dots the next step highlights — get it clear first.
      title: 'Move it clear of Open File',
      text: 'The new node landed overlapping Open File. Drag the “Add New Column” node aside (by ' +
        'its title bar) until both nodes are fully visible.',
      skipIf: (ctx) => !byNodeFunc('AddNewColumn')(ctx) || ancTableConnected(ctx) ||
        nodeIsApart(byNodeFunc('AddNewColumn'))(ctx),
      target: byNodeFunc('AddNewColumn'),
      position: 'top',
      until: untilNodeApart(byNodeFunc('AddNewColumn')),
    },
    {
      title: 'Wire the table in',
      text: 'The editor needs real data. Drag from Open File\'s result output dot (highlighted) to ' +
        'Add New Column\'s table input dot (highlighted).',
      skipIf: ancTableConnected,
      target: byNodeFunc('AddNewColumn'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx),
        socketOf(byNodeFunc('AddNewColumn'), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Select the node',
      text: 'Click the “Add New Column” node so its settings open on the right.',
      target: byNodeFunc('AddNewColumn'),
      until: untilNodeSelectedOfFunc('AddNewColumn'),
    },
    {
      title: 'Open the function\'s editor',
      text: 'In the parameters pane header (titled with the function name), click “Open editor” (highlighted). Flow opens the ' +
        'function\'s own dialog seeded with the real upstream table — running the flow up to that ' +
        'point first if it hasn\'t run yet.',
      target: byTid('prop-func-editor'),
      position: 'left',
      until: untilExists('.d4-flow-function-funccall-editor'),
    },
    {
      title: 'Edit and confirm',
      text: 'Configure the column in the dialog — real column pickers, live preview — then click OK. ' +
        'The values are written back into the node; with Autorun on, the affected nodes rerun ' +
        'right after.',
    },
  ]),
  q('how-rerun-node', 'How do I rerun just one node?', {
    title: 'Rerun this node only',
    text: 'After a run, right-click a node and choose “Rerun this node only”. It reuses the data ' +
      'that flowed into it during the last run, so the earlier steps don\'t run again — handy for ' +
      'tweaking one step of a slow flow. (Available once the node has run; rewiring clears the ' +
      'stored data until the next run.)',
  }),
  q('how-preview', 'How do I preview a node\'s data?', [
    ...loadDemogViaFiles(openFileHasPath),
    {
      title: 'Run up to here & preview',
      text: 'Right-click Open File\'s result output dot (highlighted) and choose “Run up to here & ' +
        'preview”. Flow runs just that slice and shows the data — no full run or output node needed.',
      target: byNodeFunc('OpenFile'),
      position: 'top',
      highlights: (ctx) => [socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx)],
      // untilExists would be fooled by a permanently-mounted-but-hidden panel.
      until: untilVisible('[data-testid="ff-port-preview"], [data-testid="ff-output-panel"]'),
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
    text: 'The FUNCTIONS list is grouped by what each function does — Data Sources, Combine ' +
      'Tables, Transform Tables, and more. The “by: …” button (highlighted) switches the ' +
      'grouping. Click Finish when you\'ve seen it.',
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
      text: 'Click Save (highlighted) and give your flow a name — it\'s stored on the platform so ' +
        'you can reopen it from anywhere, and it appears in the toolbox\'s Workflows tab for reuse. ' +
        '(To hand a colleague a file instead: Flow → Export .ffjson.)',
      target: byTid('ribbon', 'save'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'save')),
    },
  ]),
  q('how-reuse-flow', 'How do I reuse a saved flow inside another flow?', {
    title: 'Saved flows are building blocks',
    text: 'Save a flow (the Save button) and it appears in the toolbox\'s Workflows tab ' +
      '(highlighted) — every saved flow is listed there. Double-click or drag one onto the ' +
      'canvas to use it as a single step inside the flow you\'re building.',
    setup: async (ctx) => {
      ctx.host.showFunctionBrowser();
      await delay(60);
      ctx.host.showToolboxTab('Workflows');
    },
    target: byTid('browser-workflows'),
  }),
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
    text: 'Click the Open (folder) icon (highlighted) in the ribbon and pick one of your saved flows ' +
      'from the platform. (A .ffjson file from a colleague? Use Flow → Import .ffjson… instead.)',
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
  q('how-visualize', 'How do I add visualization nodes?', [
    ...loadDemogViaFiles(openFileHasPath),
    {
      title: 'Add a Scatter Plot',
      text: 'In the toolbox search box type scatter, then double-click “Scatter Plot” (highlighted) ' +
        'under Viewers to add a scatter-plot node.',
      skipIf: hasNodeType('Viewers/Scatter Plot'),
      setup: findInBrowser('Scatter Plot'),
      target: byTid('browser-item', 'Viewers/Scatter Plot'),
      until: untilNodeType('Viewers/Scatter Plot'),
    },
    {
      title: 'Move it clear of Open File',
      text: 'Drag the “Scatter Plot” node to the right until it no longer overlaps Open File.',
      skipIf: nodeIsRightOf(byNodeType('Viewers/Scatter Plot'), byNodeFunc('OpenFile'), 200),
      target: byNodeType('Viewers/Scatter Plot'),
      position: 'top',
      until: untilNodeRightOf(byNodeType('Viewers/Scatter Plot'), byNodeFunc('OpenFile'), 200),
    },
    {
      title: 'Feed data into the chart',
      text: 'Drag from Open File\'s result output dot (highlighted) to the Scatter Plot\'s table ' +
        'input dot (highlighted).',
      target: byNodeType('Viewers/Scatter Plot'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx),
        socketOf(byNodeType('Viewers/Scatter Plot'), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Open the chart\'s options',
      text: 'Click the “Scatter Plot” node so its options open in the panel on the right.',
      target: byNodeType('Viewers/Scatter Plot'),
      until: untilNodeOfTypeSelected('Viewers/Scatter Plot'),
    },
    {
      title: 'Choose the X column',
      text: 'Type age into the X field — that becomes the horizontal axis. (You can also wire a ' +
        'column straight into the X socket on the node.)',
      target: byParam('X'),
      position: 'left',
      until: untilValueContains(paramFieldSelector('X'), 'age'),
    },
    {
      title: 'Choose the Y column',
      text: 'Type height into the Y field for the vertical axis.',
      target: byParam('Y'),
      position: 'left',
      until: untilValueContains(paramFieldSelector('Y'), 'height'),
    },
    {
      title: 'Add a Bar Chart',
      text: 'Now search bar and double-click “Bar Chart” (highlighted) — every viewer type exposes ' +
        'its own options.',
      skipIf: hasNodeType('Viewers/Bar Chart'),
      setup: findInBrowser('Bar Chart'),
      target: byTid('browser-item', 'Viewers/Bar Chart'),
      until: untilNodeType('Viewers/Bar Chart'),
    },
    {
      title: 'Move the Bar Chart clear',
      text: 'If the new chart overlaps another node, drag the “Bar Chart” aside until it overlaps ' +
        'nothing.',
      skipIf: (ctx) => !byNodeType('Viewers/Bar Chart')(ctx) ||
        nodeIsApart(byNodeType('Viewers/Bar Chart'))(ctx),
      target: byNodeType('Viewers/Bar Chart'),
      position: 'top',
      until: untilNodeApart(byNodeType('Viewers/Bar Chart')),
    },
    {
      title: 'Wire the Bar Chart up',
      text: 'Drag from Open File\'s result output dot (highlighted) to the Bar Chart\'s table input ' +
        'dot (highlighted).',
      setup: fitGraph,
      target: byNodeType('Viewers/Bar Chart'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx),
        socketOf(byNodeType('Viewers/Bar Chart'), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Add a Pie Chart',
      text: 'Finally search pie and double-click “Pie Chart” (highlighted).',
      skipIf: hasNodeType('Viewers/Pie Chart'),
      setup: findInBrowser('Pie Chart'),
      target: byTid('browser-item', 'Viewers/Pie Chart'),
      until: untilNodeType('Viewers/Pie Chart'),
    },
    {
      title: 'Move the Pie Chart clear',
      text: 'Same again — if it overlaps, drag the “Pie Chart” aside until it overlaps nothing.',
      skipIf: (ctx) => !byNodeType('Viewers/Pie Chart')(ctx) ||
        nodeIsApart(byNodeType('Viewers/Pie Chart'))(ctx),
      target: byNodeType('Viewers/Pie Chart'),
      position: 'top',
      until: untilNodeApart(byNodeType('Viewers/Pie Chart')),
    },
    {
      title: 'Wire the Pie Chart up',
      text: 'Drag from Open File\'s result output dot (highlighted) to the Pie Chart\'s table input ' +
        'dot (highlighted).',
      setup: fitGraph,
      target: byNodeType('Viewers/Pie Chart'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFunc('OpenFile'), 'output', 'result')(ctx),
        socketOf(byNodeType('Viewers/Pie Chart'), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Run to see the charts',
      text: 'Click Run (the ▶ icon). When every node reports Done, click any chart node to render ' +
        'it in the bottom panel — use the gear in a chart\'s corner to fine-tune every setting.',
      target: byTid('ribbon', 'run'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'run')),
    },
  ]),
  q('how-join', 'How do I join two tables?', [
    ...loadDemogViaFiles(openFileHasPath),
    {
      title: 'Open the Demo connection',
      text: `Make sure the “${DEMO_CONN}” connection is expanded — double-click it (highlighted) if ` +
        'it isn\'t.',
      setup: openFiles,
      target: byFileTreeConn(DEMO_CONN),
      until: untilFileTreeConnExpanded(DEMO_CONN),
    },
    {
      title: 'Scroll to demog.csv',
      text: `Scroll the Files list down until “${DEMOG_FILE}” comes into view again (it sits below the ` +
        'folders).',
      setup: openFiles,
      // Highlight the whole Files pane while the file is off-screen, then snap to
      // the file once it's visible (same pattern as the data-loading steps).
      target: (ctx) => demogVisible(ctx) ? byFileTreeFile(DEMOG_FILE)(ctx) : byTid('browser-files')(ctx),
      highlights: (ctx) => [demogVisible(ctx) ? byFileTreeFile(DEMOG_FILE)(ctx) : byTid('browser-files')(ctx)],
      until: untilScrolledIntoView(DEMOG_FILE_SEL),
    },
    {
      title: 'Add the second table',
      text: `Double-click “${DEMOG_FILE}” (highlighted) to drop a second Open File node. (We use ` +
        'the same file twice so the columns are guaranteed to match — with your data these would ' +
        'be two different tables.)',
      target: byFileTreeFile(DEMOG_FILE),
      until: untilMoreNodes(),
    },
    {
      title: 'Move the new file aside',
      text: 'If the second Open File overlaps anything, drag it (by its title bar) below the ' +
        'first one so both are fully visible.',
      skipIf: (ctx) => !byNodeFuncNth('OpenFile', 1)(ctx) ||
        nodeIsApart(byNodeFuncNth('OpenFile', 1))(ctx),
      target: byNodeFuncNth('OpenFile', 1),
      position: 'top',
      until: untilNodeApart(byNodeFuncNth('OpenFile', 1)),
    },
    {
      title: 'Add Join Tables',
      text: 'Search join and double-click “Join Tables” (highlighted) to add it.',
      skipIf: hasFuncNode('JoinTables'),
      setup: findInBrowser('join'),
      target: byBrowserFunc('JoinTables'),
      until: untilFuncNode('JoinTables'),
    },
    {
      title: 'Move it into open space',
      text: 'Drag the “Join Tables” node to the right of both Open File nodes, clear of everything, ' +
        'so you can wire both tables into it.',
      skipIf: (ctx) => !byNodeFunc('JoinTables')(ctx) || nodeIsApart(byNodeFunc('JoinTables'))(ctx),
      target: byNodeFunc('JoinTables'),
      position: 'top',
      until: untilNodeApart(byNodeFunc('JoinTables')),
    },
    {
      title: 'Connect the first table',
      text: 'Drag from the first Open File node\'s result output dot (highlighted) to Join Tables\' ' +
        'table1 input dot (highlighted).',
      setup: fitGraph,
      target: byNodeFunc('JoinTables'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFuncNth('OpenFile', 0), 'output', 'result')(ctx),
        socketOf(byNodeFunc('JoinTables'), 'input', 'table1')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Connect the second table',
      text: 'Now drag from the second Open File node\'s result output dot (highlighted) to Join ' +
        'Tables\' table2 input dot (highlighted).',
      setup: fitGraph,
      target: byNodeFunc('JoinTables'),
      position: 'top',
      highlights: (ctx) => [
        socketOf(byNodeFuncNth('OpenFile', 1), 'output', 'result')(ctx),
        socketOf(byNodeFunc('JoinTables'), 'input', 'table2')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Open Join Tables\' settings',
      text: 'Click the “Join Tables” node. On the right you\'ll see keys1, keys2, values1 and values2 ' +
        '— the columns to match on, and the columns to carry over from each table.',
      target: byNodeFunc('JoinTables'),
      until: untilNodeSelectedOfFunc('JoinTables'),
    },
    {
      title: `Pick the first key — ${DEMOG_KEY}`,
      text: `Next to keys1, click the list icon (highlighted). Flow loads the real columns from the ` +
        `first table (running the flow up to that point if needed). In the dialog, select the ` +
        `${DEMOG_KEY} column and click OK.`,
      // While the dialog is open it sits behind the card — re-anchor the card to
      // the dialog (and stop pulsing the now-hidden icon) so it's not obscured.
      target: preferDialog(byTid('prop-pick-columns', 'keys1')),
      highlights: (ctx) => openDialogEl() ? [] : [byTid('prop-pick-columns', 'keys1')(ctx)],
      position: 'left',
      until: untilValueContains(paramFieldSelector('keys1'), DEMOG_KEY),
    },
    {
      title: `Pick the matching key — ${DEMOG_KEY}`,
      text: `Do the same for keys2 — click its list icon (highlighted) and select ${DEMOG_KEY} from ` +
        'the second table, then OK. Matching keys are what the join lines up rows on.',
      target: preferDialog(byTid('prop-pick-columns', 'keys2')),
      highlights: (ctx) => openDialogEl() ? [] : [byTid('prop-pick-columns', 'keys2')(ctx)],
      position: 'left',
      until: untilValueContains(paramFieldSelector('keys2'), DEMOG_KEY),
    },
    {
      title: 'Choose what to carry over — values1',
      text: 'values1 and values2 are required too — the columns to bring from each table. Click the ' +
        'list icon next to values1 (highlighted) and select every column (use the header checkbox to ' +
        'select all), then OK.',
      target: preferDialog(byTid('prop-pick-columns', 'values1')),
      highlights: (ctx) => openDialogEl() ? [] : [byTid('prop-pick-columns', 'values1')(ctx)],
      position: 'left',
      until: untilColumnCountAtLeast(paramFieldSelector('values1'), DEMOG_COLUMN_COUNT),
    },
    {
      title: 'And from the second table — values2',
      text: 'Finally, click the list icon next to values2 (highlighted) and select all of the second ' +
        'table\'s columns, then OK.',
      target: preferDialog(byTid('prop-pick-columns', 'values2')),
      highlights: (ctx) => openDialogEl() ? [] : [byTid('prop-pick-columns', 'values2')(ctx)],
      position: 'left',
      until: untilColumnCountAtLeast(paramFieldSelector('values2'), DEMOG_COLUMN_COUNT),
    },
    {
      title: 'Mark the result',
      text: 'Search table output and double-click “Table Output” (highlighted) — it marks the ' +
        'joined table as the flow\'s result and docks into the OUTPUTS strip on the right edge.',
      skipIf: hasNodeType(TABLE_OUTPUT_TYPE),
      setup: findInBrowser('Table Output'),
      target: byTid('browser-item', TABLE_OUTPUT_TYPE),
      until: untilNodeType(TABLE_OUTPUT_TYPE),
    },
    {
      title: 'Wire up the result',
      text: 'Drag from Join Tables\' result output dot (highlighted) to the Table Output\'s dot in ' +
        'the OUTPUTS strip at the very right edge (also highlighted).',
      setup: fitGraph,
      target: byNodeType(TABLE_OUTPUT_TYPE),
      position: 'left',
      highlights: (ctx) => [
        socketOf(byNodeFunc('JoinTables'), 'output', 'result')(ctx),
        socketOf(byNodeType(TABLE_OUTPUT_TYPE), 'input', 'table')(ctx),
      ],
      until: untilMoreConnections(),
    },
    {
      title: 'Run the join',
      text: 'Click Run (the ▶ icon). When Join Tables reports Done with its row × column counts, ' +
        'click it to preview the joined table in the bottom panel.',
      target: byTid('ribbon', 'run'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'run')),
    },
  ]),
  q('how-open-result', 'How do I see a result as a full table?', {
    title: 'Result tabs',
    text: 'Every table output of your flow gets its own tab next to “Canvas” in the bottom status ' +
      'bar. Run the flow, then click the tab once its dot turns green — the result opens as a ' +
      'full Datagrok table view where you can add viewers, filter, and rearrange. An amber dot ' +
      'means the result is out of date — run again to refresh it.',
    // Point at a real output tab when one exists; a highlighted "Canvas" tab
    // would be a decoy on a flow with no outputs yet.
    target: (ctx) => bySel('.ff-view-tab[data-node-id]')(ctx),
    position: 'top',
  }),
  q('how-table-layouts', 'Are my result views saved with the flow?', {
    title: 'Layouts persist',
    text: 'Yes. Whatever you arrange in a result tab — viewers, filters, column order — is saved ' +
      'with the flow and restored when you reopen it. Published dashboards use the same layouts, ' +
      'even for tabs you never opened.',
    target: (ctx) => bySel('.ff-view-tab[data-node-id]')(ctx),
    position: 'top',
  }),
  q('how-publish-dashboard', 'How do I publish my results as a dashboard?', [
    {
      title: 'Run, then Save',
      text: 'Run the flow, then click Save (highlighted — it stays greyed out until the canvas ' +
        'has something on it). In the Save dialog, the Dashboard section lists your computed ' +
        'tables — keep “Create dashboard” checked and click Save.',
      target: byTid('ribbon', 'save'),
      position: 'bottom',
    },
    {
      title: 'Then name the project',
      text: 'A second dialog opens with your output tables and their layouts — name the project ' +
        'and click OK. Colleagues open the dashboard without ever seeing the flow graph. The ' +
        '“Publish your results as a dashboard” tutorial walks this end to end.',
    },
  ]),
  q('how-update-dashboard', 'How do I update a published dashboard?', {
    title: 'Re-publish — same project',
    text: 'Just Save again with “Create dashboard” checked: the dashboard remembers its flow, so ' +
      'every publish updates the same project in place — tables, views, and layouts. To start a ' +
      'fresh project instead, click “publish as new” in the Save dialog\'s Dashboard section.',
  }),
];

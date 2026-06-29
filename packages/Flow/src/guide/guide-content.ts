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
  untilClick, untilNodeType, untilMoreNodes, untilMoreConnections, untilFuncNode,
  untilFewerNodes, untilValueContains, untilValueMatches, untilNodeRightOf,
  untilNodeSelected, untilNodeSelectedOfFunc, untilNodeCollapsed,
  untilNodeCountAtLeast, untilExists, copyToClipboard, prefillSearch,
} from './guide-model';

const DEMO_FILE = 'System:DemoFiles/demog.csv';
const DEMO_FORMULA = '${AGE} * 12';
const NEW_COL_NAME = 'My New Column';
const SEARCH_SEL = '[data-testid="ff-browser-search"]';
const TABLE_OUTPUT_TYPE = 'Outputs/Table Output';

const delay = (ms: number): Promise<void> => new Promise((r) => setTimeout(r, ms));

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
    {
      title: 'Search for “Open File”',
      text: 'Click the search box and type open file (or openfile) to find the data source that ' +
        'loads a file.',
      setup: openSearch,
      target: byTid('browser-search'),
      until: untilValueMatches(SEARCH_SEL, 'openfile'),
    },
    {
      title: 'Add the Open File node',
      text: 'Double-click “Open File” (highlighted) to drop it on the canvas.',
      target: byBrowserFunc('OpenFile'),
      until: untilFuncNode('OpenFile'),
    },
    {
      title: 'Select the Open File node',
      text: 'Click the “Open File” node on the canvas. Its settings open in the panel on the right.',
      target: byNodeFunc('OpenFile'),
      until: untilNodeSelectedOfFunc('OpenFile'),
    },
    {
      title: 'Paste the file path',
      text: `We copied a demo path to your clipboard. Click the File path field and paste it ` +
        `(Ctrl+V): ${DEMO_FILE}`,
      setup: putOnClipboard(DEMO_FILE),
      target: byParam('fullPath'),
      position: 'left',
      until: untilValueContains(paramFieldSelector('fullPath'), 'demog.csv'),
    },
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
      text: 'Type “join” into the search box to filter the whole catalog to matching functions.',
      target: byTid('browser-search'),
      until: untilValueContains('[data-testid="ff-browser-search"]', 'join'),
    },
    {
      title: 'Add a function',
      text: 'Double-click any result to drop it on the canvas (or drag it where you want it).',
      target: byTid('browser-tree'),
      until: untilMoreNodes(),
    },
    {
      title: 'Tip: let Flow suggest the next step',
      text: 'You can also drag from a node\'s output dot into empty canvas — Flow offers compatible ' +
        'next functions, the common ones first. Click Finish.',
      target: byTid('canvas'),
      position: 'top',
    },
  ],
};

const organizeCanvas: Guide = {
  id: 'organize-canvas',
  kind: 'tutorial',
  title: 'Organize your canvas',
  summary: 'Collapse nodes, tidy the layout, and navigate with the overview.',
  steps: [
    {
      title: 'Add a couple of nodes',
      text: 'Let\'s get something to arrange. Double-click two items in the function list so the ' +
        'canvas has at least two nodes.',
      setup: (ctx) => ctx.host.showFunctionBrowser(),
      target: byTid('browser-tree'),
      until: untilNodeCountAtLeast(2),
    },
    {
      title: 'Collapse a node',
      text: 'Click the ▾ caret in a node\'s title bar (highlighted) to fold it down to just the ' +
        'title. Click ▸ to expand it again.',
      target: bySel('.ff-node .ff-node-caret'),
      until: untilNodeCollapsed(),
    },
    {
      title: 'Tidy the layout',
      text: 'Click “Tidy up layout” (highlighted) to auto-arrange everything left-to-right.',
      target: byTid('ribbon', 'layout'),
      position: 'bottom',
      until: untilClick(byTid('ribbon', 'layout')),
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

export const TUTORIALS: Guide[] = [loadDataAddColumn, findFunctions, organizeCanvas, reuseScript];

// ============================ HOW-TO QUESTIONS ============================

function q(id: string, title: string, step: GuideStep): Guide {
  return {id, kind: 'question', title, summary: title, steps: [step]};
}

export const QUESTIONS: Guide[] = [
  q('how-add-data', 'How do I bring data in?', {
    title: 'Add a data source',
    text: `The “Open File” data source loads a file. We filtered the list to it and copied a demo ` +
      `path (${DEMO_FILE}) to your clipboard — double-click Open File (highlighted) to add it, ` +
      `then select it and paste the path into File path.`,
    setup: async (ctx) => {
      await findInBrowser('OpenFile')(ctx);
      await copyToClipboard(DEMO_FILE);
    },
    target: byBrowserFunc('OpenFile'),
    until: untilFuncNode('OpenFile'),
  }),
  q('how-add-column', 'How do I add a calculated column?', {
    title: 'Add a calculated column',
    text: 'Add the “Add New Column” function (highlighted), wire a table into it, then select it ' +
      'and fill in the Name and Expression fields in the panel on the right (e.g. ${AGE} * 12).',
    setup: findInBrowser('AddNewColumn'),
    target: byBrowserFunc('AddNewColumn'),
    until: untilFuncNode('AddNewColumn'),
  }),
  q('how-set-param', 'How do I edit a node\'s settings?', {
    title: 'Edit a node in the context panel',
    text: 'Click a node to select it — its parameters appear in the panel on the right, each as an ' +
      'editable field. Type or paste values there.',
    target: byTid('canvas'),
    until: untilNodeSelected(),
  }),
  q('how-add-function', 'How do I add a function?', {
    title: 'Add a function',
    text: 'Open the list on the left, search or browse, then double-click a function (or drag it ' +
      'onto the canvas).',
    setup: (ctx) => ctx.host.showFunctionBrowser(),
    target: byTid('browser-search'),
    until: untilMoreNodes(),
  }),
  q('how-connect', 'How do I connect two nodes?', {
    title: 'Connect nodes',
    text: 'Drag from a node\'s output dot (right side) to another node\'s input dot (left side). ' +
      'Matching colors = compatible types.',
    target: byTid('canvas'),
    position: 'top',
    until: untilMoreConnections(),
  }),
  q('how-run', 'How do I run my flow?', {
    title: 'Run the flow',
    text: 'Click Run (the ▶ icon, highlighted) in the ribbon. Nodes light up as they execute; each ' +
      'reports its row/column counts.',
    target: byTid('ribbon', 'run'),
    position: 'bottom',
    until: untilClick(byTid('ribbon', 'run')),
  }),
  q('how-preview', 'How do I preview a node\'s data?', {
    title: 'Inspect any node',
    text: 'Right-click a node\'s output dot and choose “Run up to here & preview”. Flow runs just ' +
      'that slice and shows the data — no full run or output node needed.',
    target: bySel('.ff-node .ff-socket-row-output'),
    until: untilExists('[data-testid="ff-port-preview"], [data-testid="ff-output-panel"]'),
  }),
  q('how-delete', 'How do I delete a node?', {
    title: 'Delete a node',
    text: 'Select a node and press Delete (or Backspace), or right-click it and choose Delete.',
    target: byTid('canvas'),
    until: untilFewerNodes(),
  }),
  q('how-layout', 'How do I tidy up the layout?', {
    title: 'Tidy the layout',
    text: 'Click “Tidy up layout” (highlighted) in the ribbon to auto-arrange nodes left-to-right.',
    target: byTid('ribbon', 'layout'),
    position: 'bottom',
    until: untilClick(byTid('ribbon', 'layout')),
  }),
  q('how-categories', 'How do I find functions by category?', {
    title: 'Browse by category',
    text: 'Use the “Group by” dropdown (highlighted). Functions bucket by what they do — Data ' +
      'Sources, Combine Tables, Transform Tables, Column Operations, Compute Values, Visualize.',
    setup: (ctx) => ctx.host.showFunctionBrowser(),
    target: byTid('browser-groupby'),
  }),
  q('how-collapse', 'How do I collapse a node?', {
    title: 'Collapse a node',
    text: 'Click the ▾ caret in a node\'s title bar (highlighted) to fold it; click ▸ to expand. ' +
      '(The status dot is just an indicator — it no longer collapses.)',
    target: bySel('.ff-node .ff-node-caret'),
    until: untilNodeCollapsed(),
  }),
  q('how-undo', 'How do I undo a change?', {
    title: 'Undo',
    text: 'Click Undo (highlighted) in the ribbon, or press Ctrl+Z. Redo is right next to it.',
    target: byTid('ribbon', 'undo'),
    position: 'bottom',
    until: untilClick(byTid('ribbon', 'undo')),
  }),
  q('how-save', 'How do I save or share my flow?', {
    title: 'Save / share',
    text: 'Click Save (highlighted) to download a .ffjson file. Open it later from anywhere in ' +
      'Datagrok to restore the flow.',
    target: byTid('ribbon', 'save'),
    position: 'bottom',
    until: untilClick(byTid('ribbon', 'save')),
  }),
  q('how-view-script', 'How do I see the generated script?', {
    title: 'See the script',
    text: 'Click “See the steps” — the 👁 (eye) icon, highlighted. Your visual flow compiles to a ' +
      'real, editable Datagrok script.',
    target: byTid('ribbon', 'view-script'),
    position: 'bottom',
    until: untilClick(byTid('ribbon', 'view-script')),
  }),
  q('how-open', 'How do I open a saved flow?', {
    title: 'Open a flow',
    text: 'Click the Open (folder) icon, highlighted, and pick a .ffjson file.',
    target: byTid('ribbon', 'open'),
    position: 'bottom',
    until: untilClick(byTid('ribbon', 'open')),
  }),
  q('how-navigate', 'How do I navigate a large flow?', {
    title: 'Navigate',
    text: 'Use the overview (bottom-right, appears once you have nodes) — click or drag to move; ' +
      'scroll to zoom; or use Zoom to fit in the ribbon.',
    target: byTid('minimap'),
    position: 'left',
    until: untilClick(byTid('minimap')),
  }),
];

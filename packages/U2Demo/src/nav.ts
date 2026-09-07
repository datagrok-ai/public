/* The demo's navigation vocabulary: one 2-level registry (groups → leaves) feeding the in-app
   tree, the URL path (`/<groupId>/<leafId>`) and the Browse sub-tree. */
import type {Signal} from '@datagrok-libraries/u2';
import {convergencePage} from './convergence';
import {funcConvergencePage} from './func-convergence';
import {funcsPage} from './funcs';
import {funcHistoryPage} from './func-history';
import {basicInputsPage, rangeSliderPage, multiSelectPage, asyncPage} from './pages/inputs';
import {containersPage, popupsPage} from './pages/containers';
import {listsPage, treesPage} from './pages/collections';
import {cardsPage, feedbackPage, tablesPage, sectionsPage, messagingPage} from './pages/display';
import {overviewPage} from './pages/overview';
import {formPage, propertyGridPage, objectFormPage} from './pages/forms';
import {filesPage, dataframesPage, entitiesPage, spacesPage, moleculesPage, bridgePage}
  from './pages/platform';
import {msaWorkbenchPage} from './pages/msa-workbench';

/** Repo-relative root of this package's sources; the source panel keys its bundled text by it. */
export const SRC_ROOT = 'packages/U2Demo/src';

export interface SourceRef {
  /** Repo-relative path in datagrok-ai/public, e.g. 'packages/U2Demo/src/pages/inputs.ts'. */
  file: string;
  /** Top-level function to slice out; omitted = the whole file. */
  symbol?: string;
}

export interface DemoContext {
  /** The last demo action, reported on the shell status bar. */
  status: Signal<string>;
  /** Opens another sub-demo — what the Overview's pointers are. */
  navigate: (id: string) => void;
}

export interface DemoGroup {
  id: string;
  label: string;
  children: DemoLeaf[];
}

export interface DemoLeaf {
  id: string;
  /** The id's noun, short enough not to truncate in the nav pane; what the page covers is the
   * description, which the source panel shows and a tree row can hang off `title`. */
  label: string;
  description: string;
  group: DemoGroup;
  build: (ctx: DemoContext) => HTMLElement;
  source: SourceRef;
}

type LeafSpec = Omit<DemoLeaf, 'group'>;

function group(id: string, label: string, leaves: LeafSpec[]): DemoGroup {
  const g: DemoGroup = {id, label, children: []};
  g.children = leaves.map((leaf) => ({...leaf, group: g}));
  return g;
}

export const DEMO_TREE: DemoGroup[] = [
  group('start', 'Start', [
    {id: 'overview', label: 'Overview',
      description: 'What u2 is, what each area covers, and the three sub-demos to read first',
      build: (ctx) => overviewPage(ctx),
      source: {file: `${SRC_ROOT}/pages/overview.ts`, symbol: 'overviewPage'}},
  ]),
  group('inputs', 'Inputs', [
    {id: 'all-inputs', label: 'All inputs',
      description: 'Every input type, Dart ui.input and u2 propertyForm side by side',
      build: () => convergencePage(),
      source: {file: `${SRC_ROOT}/convergence.ts`, symbol: 'convergencePage'}},
    {id: 'basic-inputs', label: 'Basic inputs',
      description: 'Signals, binding and validation across the basic inputs',
      build: () => basicInputsPage(),
      source: {file: `${SRC_ROOT}/pages/inputs.ts`, symbol: 'basicInputsPage'}},
    {id: 'range-slider', label: 'Range slider',
      description: 'Two-handle numeric range over a pair of signals, with minRange and step',
      build: () => rangeSliderPage(),
      source: {file: `${SRC_ROOT}/pages/inputs.ts`, symbol: 'rangeSliderPage'}},
    {id: 'multi-select', label: 'Multi-select',
      description: 'Popup multi-select and segmented button groups over a set of values',
      build: () => multiSelectPage(),
      source: {file: `${SRC_ROOT}/pages/inputs.ts`, symbol: 'multiSelectPage'}},
    {id: 'async', label: 'Async',
      description: 'Combobox over local and async items, AsyncView loading states',
      build: () => asyncPage(),
      source: {file: `${SRC_ROOT}/pages/inputs.ts`, symbol: 'asyncPage'}},
  ]),
  group('containers', 'Containers', [
    {id: 'layout', label: 'Layout',
      description: 'Splitter, accordion, breadcrumbs, panel-local toolbar',
      build: () => containersPage(),
      source: {file: `${SRC_ROOT}/pages/containers.ts`, symbol: 'containersPage'}},
    {id: 'popups', label: 'Popups',
      description: 'Menu, context menu, tooltip, dialog',
      build: (ctx) => popupsPage(ctx.status),
      source: {file: `${SRC_ROOT}/pages/containers.ts`, symbol: 'popupsPage'}},
  ]),
  group('collections', 'Collections', [
    {id: 'lists', label: 'Lists',
      description: 'VirtualList over 100,000 items',
      build: () => listsPage(),
      source: {file: `${SRC_ROOT}/pages/collections.ts`, symbol: 'listsPage'}},
    {id: 'trees', label: 'Trees',
      description: 'VirtualTree with lazy branches and expandPath',
      build: () => treesPage(),
      source: {file: `${SRC_ROOT}/pages/collections.ts`, symbol: 'treesPage'}},
  ]),
  group('display', 'Display', [
    {id: 'cards', label: 'Cards',
      description: 'Card surfaces and StatCard KPIs, clickable and selectable',
      build: () => cardsPage(),
      source: {file: `${SRC_ROOT}/pages/display.ts`, symbol: 'cardsPage'}},
    {id: 'feedback', label: 'Feedback',
      description: 'Progress bars, notification balloons and the guided tour',
      build: () => feedbackPage(),
      source: {file: `${SRC_ROOT}/pages/display.ts`, symbol: 'feedbackPage'}},
    {id: 'tables', label: 'Tables',
      description: 'BasicTable over small data and the virtualized VirtualGrid',
      build: () => tablesPage(),
      source: {file: `${SRC_ROOT}/pages/display.ts`, symbol: 'tablesPage'}},
    {id: 'sections', label: 'Sections & wizard',
      description: 'A collapsible section and the stepped wizard over gated content',
      build: () => sectionsPage(),
      source: {file: `${SRC_ROOT}/pages/display.ts`, symbol: 'sectionsPage'}},
    {id: 'messaging', label: 'Message input',
      description: 'The compose box: @-mentions as atomic chips, value as Datagrok markup',
      build: () => messagingPage(),
      source: {file: `${SRC_ROOT}/pages/display.ts`, symbol: 'messagingPage'}},
  ]),
  group('forms', 'Forms', [
    {id: 'form', label: 'Form',
      description: 'Form layout, sections and the summary checkbox',
      build: () => formPage(),
      source: {file: `${SRC_ROOT}/pages/forms.ts`, symbol: 'formPage'}},
    {id: 'property-grid', label: 'Property grid',
      description: 'PropertyGrid over a plain object',
      build: () => propertyGridPage(),
      source: {file: `${SRC_ROOT}/pages/forms.ts`, symbol: 'propertyGridPage'}},
    {id: 'objectform', label: 'Object form',
      description: 'objectForm over an object with typed properties',
      build: () => objectFormPage(),
      source: {file: `${SRC_ROOT}/pages/forms.ts`, symbol: 'objectFormPage'}},
    {id: 'funcs', label: 'Functions',
      description: 'Forms built from Func metadata: the platform form and u2 funcForm',
      build: () => funcsPage(),
      source: {file: `${SRC_ROOT}/funcs.ts`, symbol: 'funcsPage'}},
    {id: 'func-convergence', label: 'FuncCalls',
      description: 'FuncCall editors: DG.InputForm.forFuncCall and u2 funcForm, side by side',
      build: () => funcConvergencePage(),
      source: {file: `${SRC_ROOT}/func-convergence.ts`, symbol: 'funcConvergencePage'}},
    {id: 'func-history', label: 'Run history',
      description: 'FunctionInput over a funcForm with the standard Run button and history icon — runs save to the server and come back from the popup',
      build: () => funcHistoryPage(),
      source: {file: `${SRC_ROOT}/func-history.ts`, symbol: 'funcHistoryPage'}},
  ]),
  group('platform', 'Platform', [
    {id: 'dataframes', label: 'Dataframes',
      description: 'Table, column and columns pickers bound to live dataframes',
      build: () => dataframesPage(),
      source: {file: `${SRC_ROOT}/pages/platform.ts`, symbol: 'dataframesPage'}},
    {id: 'files', label: 'Files',
      description: 'File picker and file browsing',
      build: () => filesPage(),
      source: {file: `${SRC_ROOT}/pages/platform.ts`, symbol: 'filesPage'}},
    {id: 'entities', label: 'Entities',
      description: 'User and group pickers, entity chips and handler renderers',
      build: () => entitiesPage(),
      source: {file: `${SRC_ROOT}/pages/platform.ts`, symbol: 'entitiesPage'}},
    {id: 'spaces', label: 'Spaces',
      description: 'Space, project and dashboard pickers with a preview',
      build: () => spacesPage(),
      source: {file: `${SRC_ROOT}/pages/platform.ts`, symbol: 'spacesPage'}},
    {id: 'molecules', label: 'Molecules',
      description: 'moleculeInput, structure chips and a structure typeahead',
      build: () => moleculesPage(),
      source: {file: `${SRC_ROOT}/pages/platform.ts`, symbol: 'moleculesPage'}},
    {id: 'bridge', label: 'Bridge',
      description: 'DartInput bridge and the leak detector',
      build: () => bridgePage(),
      source: {file: `${SRC_ROOT}/pages/platform.ts`, symbol: 'bridgePage'}},
  ]),
  group('automation', 'Automation', [
    {id: 'msa-workbench', label: 'MSA workbench',
      description: 'A small app the behavioral tests in libraries/bdd drive: toolbar, form, dialog, tabs, list',
      build: () => msaWorkbenchPage(),
      source: {file: `${SRC_ROOT}/pages/msa-workbench.ts`, symbol: 'msaWorkbenchPage'}},
  ]),
];

export function leafById(id: string): DemoLeaf | undefined {
  for (const g of DEMO_TREE) {
    const leaf = g.children.find((l) => l.id === id);
    if (leaf)
      return leaf;
  }
  return undefined;
}

/** Resolves the tail of `/apps/U2demo/U2Demo/...` — leaf ids are globally unique, so the last
 * segment decides: `/forms/funcs` and `/funcs` both land on Functions. */
export function leafForPath(path: string | undefined): DemoLeaf | undefined {
  const segments = (path ?? '').split('/').filter((s) => s !== '');
  return segments.length === 0 ? undefined : leafById(segments[segments.length - 1]);
}

export function pathOf(leaf: DemoLeaf): string {
  return `/${leaf.group.id}/${leaf.id}`;
}

/** The navigation → source-panel seam: the panel installs the hook, so neither demo.ts nor
 * ribbon.ts ever imports it. */
export type PageShownHook = (leaf: DemoLeaf, open?: boolean) => void;
let hook: PageShownHook | undefined;
export function setPageShownHook(h: PageShownHook): void { hook = h; }
export function pageShown(leaf: DemoLeaf, open = false): void { hook?.(leaf, open); }

/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FuncInfo, getRegisteredFuncs, isWorkflowFunc, VIEWER_NODE_TYPES} from '../rete/node-factory';
import {tid, setTid} from '../utils/test-ids';
import {getFilesBrowser} from '../utils/files-browser-tree';
import {categorizeBySignature, domainCategory} from '../types/type-map';

/** Whether a function's output is a widget (→ the Widgets pane, not a category). */
export function funcOutputsWidget(f: FuncInfo): boolean {
  try {
    return f.func.outputs.some((o) => String(o.propertyType) === 'widget');
  } catch {
    return false;
  }
}

export type GroupByMode = 'category' | 'role' | 'tags' | 'package';

/** The task-oriented categories, in the order a scientist meets a pipeline:
 *  get data in → combine → reshape → derive columns → compute values →
 *  visualize → everything else. Used both to bucket and to order the tree. */
export const FUNC_CATEGORIES = [
  'Data Sources',
  'Workflows',
  'Combine Tables',
  'Transform Tables',
  'Column Operations',
  'Compute Values',
  'Visualize',
  'Cheminformatics',
  'Bioinformatics',
  'Other',
] as const;
export type FuncCategory = (typeof FUNC_CATEGORIES)[number];

/** Bucket a function by **what it does**, derived from its input/output
 *  signature (and role for viewers) — delegates to the shared
 *  `categorizeBySignature` (the same logic that colors node title bars). The
 *  key distinctions, validated against the live catalog:
 *  - **Data Sources** produce a table from *no* table input (OpenFile, DB
 *    queries, generators). A join is NOT a data source — it consumes tables.
 *  - **Combine Tables** take ≥2 tables (join / union / link / react).
 *  - **Transform Tables** take exactly one table and return a table or mutate
 *    it in place (filter, aggregate, pivot, sort, …).
 *  - **Column Operations** emit a column / column list (AddNewColumn, descriptors).
 *  - **Compute Values** emit a scalar (statistics, math, text, predicates).
 *  - **Visualize** emit a viewer / view / widget / graphics.
 *  - **Cheminformatics** / **Bioinformatics** — grouped by source package
 *    (domain wins over the signature category), so a scientist finds all
 *    chem/bio steps together regardless of what they do.
 *  - **Workflows** — saved flows (`DG.Script`, language `flow`); checked
 *    before everything else since their signature reads as a data source.
 *  - **Other** — expression/filter builders and dynamic helpers. */
export function categorizeFunc(func: DG.Func, role: string | null, packageName?: string): FuncCategory {
  if (isWorkflowFunc(func)) return 'Workflows';
  const inputTypes = func.inputs.map((i) => String(i.propertyType));
  const domain = domainCategory(packageName, inputTypes);
  if (domain) return domain;
  return categorizeBySignature(
    inputTypes,
    func.outputs.map((o) => String(o.propertyType)),
    role) as FuncCategory;
}

/** For each domain section, the flagship package whose functions lead the list
 *  (the rest of that domain's packages follow). */
export const DOMAIN_PRIMARY_PACKAGE: Record<string, string> = {
  Cheminformatics: 'Chem',
  Bioinformatics: 'Bio',
};

/** The "most-used" operations a scientist reaches for first inside a domain
 *  section — matched as lowercased substrings against the function name (lower
 *  index = higher up). Everything unlisted sorts after these, alphabetically —
 *  so descriptors / properties / similarity lead each domain list. */
export const DOMAIN_PRIORITY_KEYWORDS: Record<string, string[]> = {
  Cheminformatics: [
    'descriptor', 'propert', 'fingerprint', 'similar', 'substructure',
    'diversit', 'cluster', 'scaffold', 'r-group', 'rgroup', 'mcs',
    'activity', 'toxic', 'risk', 'elemental', 'mutate', 'curate', 'convert',
  ],
  Bioinformatics: [
    'descriptor', 'propert', 'composition', 'similar', 'sequence',
    'align', 'msa', 'logo', 'atomic', 'translate', 'helm', 'convert',
  ],
};

/** Order a domain section's items: flagship-package functions first (Chem under
 *  Cheminformatics, Bio under Bioinformatics), then the rest; within each group
 *  the curated "most-used" operations (descriptors, properties, …) lead, then
 *  alphabetical by package, then by display name. Pure — returns a new array. */
export function orderDomainSection(items: FuncInfo[], category: string): FuncInfo[] {
  const primary = DOMAIN_PRIMARY_PACKAGE[category];
  const keywords = DOMAIN_PRIORITY_KEYWORDS[category] ?? [];
  const priorityRank = (info: FuncInfo): number => {
    const hay = `${info.name} ${info.func.name}`.toLowerCase();
    const i = keywords.findIndex((k) => hay.includes(k));
    return i === -1 ? keywords.length : i;
  };
  return [...items].sort((a, b) => {
    const ap = a.packageName === primary ? 0 : 1;
    const bp = b.packageName === primary ? 0 : 1;
    if (ap !== bp) return ap - bp;
    const ar = priorityRank(a);
    const br = priorityRank(b);
    if (ar !== br) return ar - br;
    const pkg = a.packageName.localeCompare(b.packageName);
    if (pkg !== 0) return pkg;
    return a.name.localeCompare(b.name);
  });
}

/** MIME-ish key used to carry a node type name through HTML5 drag/drop.
 *  The canvas drop handler reads this and adds the corresponding node at
 *  the drop point. */
export const FF_DRAG_MIME = 'application/x-funcflow-node';

export interface FunctionBrowserCallbacks {
  onFunctionDoubleClick: (funcInfo: FuncInfo) => void;
  onBuiltinNodeDoubleClick: (nodeTypeName: string) => void;
  /** Double-click (or Enter) on a file in the Files tree → add an OpenFile node. */
  onFileDoubleClick: (file: DG.FileInfo) => void;
}

/** Connection display name for a query func, used to group queries in the
 *  Queries pane (`friendlyName` falls back to `name`). */
export function queryConnectionName(f: FuncInfo): string {
  if (!(f.func instanceof DG.DataQuery)) return '';
  try {
    const conn = (f.func as DG.DataQuery).connection as DG.DataConnection | null;
    return (conn?.friendlyName ?? conn?.name ?? 'Other') || 'Other';
  } catch {
    return 'Other';
  }
}

/** Persisted UI state for the browser (the group-by mode). Section
 *  expand/collapse state is NOT kept here — the platform accordion
 *  (`ui.accordion(key)`) persists it itself in `localStorage['Accordion:<key>']`,
 *  keyed by pane name. */
interface BrowserState {
  groupBy: GroupByMode;
}
const LS_KEY = 'funcflow.browser.v1';
/** Persistence key of the toolbox accordion; the nested per-connection
 *  accordion inside the Queries pane uses `<key>.queries`. */
export const TOOLBOX_ACCORDION_KEY = 'funcflow.toolbox';
const VALID_MODES: GroupByMode[] = ['category', 'role', 'tags', 'package'];

/** Case- AND whitespace-insensitive substring match — so "openfile" matches
 *  "Open File" and "table output" matches "TableOutput". */
export function nameMatchesQuery(text: string, query: string): boolean {
  if (!query) return true;
  const t = (text || '').toLowerCase();
  const q = query.toLowerCase();
  return t.includes(q) || t.replace(/\s+/g, '').includes(q.replace(/\s+/g, ''));
}

/** Whether a function matches a search query. Matches the display name AND the
 *  raw function name / friendlyName — the toolbox shows the friendly name (e.g.
 *  "Open File") but users search by the name they know (e.g. "OpenFile"), and
 *  vice versa — all whitespace-insensitive; plus description, tags, role, pkg. */
export function funcMatchesSearch(f: FuncInfo, query: string): boolean {
  if (!query) return true;
  const q = query.toLowerCase();
  return nameMatchesQuery(f.name, q) ||
    nameMatchesQuery(f.func.name || '', q) ||
    nameMatchesQuery(f.func.friendlyName || '', q) ||
    (f.func.description || '').toLowerCase().includes(q) ||
    f.tags.some((t) => t.toLowerCase().includes(q)) ||
    (f.role || '').toLowerCase().includes(q) ||
    f.packageName.toLowerCase().includes(q);
}

/** Left sidebar: searchable, groupable function catalog */
export class FunctionBrowser {
  root: HTMLElement;
  private searchInput!: HTMLInputElement;
  private groupBySelect!: HTMLSelectElement;
  private treeContainer!: HTMLElement;
  private groupBy: GroupByMode = 'category';
  /** The toolbox accordion, rebuilt on every render. A keyed `ui.accordion`
   *  restores pane states from localStorage at construction and persists them
   *  on header clicks — no manual bookkeeping. Exposed for tests. */
  accordion: DG.Accordion | null = null;
  /** The nested per-connection accordion inside the Queries pane (built lazily
   *  when that pane first expands). Exposed for tests. */
  queriesAccordion: DG.Accordion | null = null;
  private callbacks: FunctionBrowserCallbacks;
  /** The Files tree (KNIME-style file browser). Built lazily once and reused
   *  across renders so its expanded state and scroll position survive a search. */
  private filesTreeRoot?: HTMLElement;

  constructor(callbacks: FunctionBrowserCallbacks) {
    this.callbacks = callbacks;
    this.loadState();
    this.root = this.buildUI();
  }

  // ---------- persistence (localStorage) ----------

  private loadState(): void {
    try {
      const raw = localStorage.getItem(LS_KEY);
      if (!raw) return;
      const s = JSON.parse(raw) as Partial<BrowserState>;
      if (s.groupBy && VALID_MODES.includes(s.groupBy)) this.groupBy = s.groupBy;
    } catch {/* corrupt/blocked storage — fall back to defaults */}
  }

  private saveState(): void {
    try {
      const state: BrowserState = {groupBy: this.groupBy};
      localStorage.setItem(LS_KEY, JSON.stringify(state));
    } catch {/* storage blocked/full — non-fatal */}
  }

  private buildUI(): HTMLElement {
    // Search bar
    this.searchInput = document.createElement('input');
    this.searchInput.type = 'text';
    this.searchInput.placeholder = 'Search functions...';
    this.searchInput.className = 'funcflow-search-input';
    setTid(this.searchInput, 'browser-search');
    this.searchInput.addEventListener('input', () => this.render());

    // Clear (✕) affordance at the right edge of the search box — shown only
    // while there's text to clear.
    const clearBtn = document.createElement('span');
    clearBtn.className = 'funcflow-search-clear';
    clearBtn.innerHTML = '&times;';
    ui.tooltip.bind(clearBtn, 'Clear search');
    setTid(clearBtn, 'browser-search-clear');
    const syncClear = (): void => {clearBtn.style.display = this.searchInput.value ? 'flex' : 'none';};
    clearBtn.addEventListener('click', () => {
      this.searchInput.value = '';
      this.searchInput.dispatchEvent(new Event('input', {bubbles: true}));
      this.searchInput.focus();
      syncClear();
    });
    this.searchInput.addEventListener('input', syncClear);
    syncClear();
    const searchWrap = ui.div([this.searchInput, clearBtn], 'funcflow-search-wrap');

    // Group by selector
    this.groupBySelect = document.createElement('select');
    this.groupBySelect.className = 'funcflow-groupby-select';
    setTid(this.groupBySelect, 'browser-groupby');
    const groupByLabels: Record<GroupByMode, string> = {
      category: 'what it does',
      role: 'role',
      tags: 'tags',
      package: 'package',
    };
    for (const opt of Object.keys(groupByLabels) as GroupByMode[]) {
      const option = document.createElement('option');
      option.value = opt;
      option.textContent = `Group by: ${groupByLabels[opt]}`;
      this.groupBySelect.appendChild(option);
    }
    this.groupBySelect.value = this.groupBy; // restore persisted mode
    this.groupBySelect.addEventListener('change', () => {
      this.groupBy = this.groupBySelect.value as GroupByMode;
      this.saveState();
      this.render();
    });

    this.treeContainer = setTid(ui.div([], 'funcflow-tree-container'), 'browser-tree');

    const container = ui.divV([
      searchWrap,
      this.groupBySelect,
      this.treeContainer,
    ], 'funcflow-browser');

    return setTid(container, 'browser');
  }

  /** Domain sections floated to the very top of the function list (right after
   *  the Queries pane) — the science a chemist/biologist reaches for first. */
  private static readonly DOMAIN_CATEGORIES = ['Cheminformatics', 'Bioinformatics'];

  render(): void {
    this.treeContainer.innerHTML = '';
    // The platform accordion: a keyed one restores/persists pane states itself.
    // While a search is active every matching section is forced open — build an
    // UNKEYED accordion then, so the forced states don't overwrite the user's.
    const acc = ui.accordion(this.searchInput.value ? null : TOOLBOX_ACCORDION_KEY);
    this.accordion = acc;
    this.queriesAccordion = null;

    // DG function nodes — queries, widget-producers, and (already-excluded)
    // viewer functions don't appear here; they live in their dedicated panes.
    const funcs = this.filterBySearch(getRegisteredFuncs()
      .filter((f) => !(f.func instanceof DG.DataQuery) && !funcOutputsWidget(f)));
    const grouped = this.groupFunctions(funcs);
    const renderCategory = (category: string): void => {
      const raw = grouped[category];
      if (!raw || raw.length === 0) return;
      const items = FunctionBrowser.DOMAIN_CATEGORIES.includes(category) ?
        orderDomainSection(raw, category) : raw;
      this.addSection(acc, category, () => {
        const content = ui.div([], 'funcflow-section-content');
        for (const info of items) content.appendChild(this.createFuncItem(info));
        return content;
      }, {count: items.length});
    };

    // Section order: Files, Queries, the Cheminformatics / Bioinformatics
    // domain sections (the science a chemist/biologist reaches for first),
    // the task categories (Data Sources → … → Compute Values), Viewers,
    // Widgets, the building-block built-ins (Inputs / Outputs / Constants /
    // Utilities), then the Other catch-all, and Debug last.
    this.renderFilesSection(acc);
    this.renderQueriesSection(acc);
    for (const domain of FunctionBrowser.DOMAIN_CATEGORIES) renderCategory(domain);

    const sortedKeys = this.orderGroupKeys(Object.keys(grouped))
      .filter((k) => !FunctionBrowser.DOMAIN_CATEGORIES.includes(k));
    for (const category of sortedKeys.filter((k) => k !== 'Other')) renderCategory(category);

    this.renderViewersSection(acc);
    this.renderWidgetsSection(acc);
    this.renderBuiltinNodes(acc, ['Inputs', 'Outputs', 'Constants', 'Utilities']);
    if (sortedKeys.includes('Other')) renderCategory('Other');
    this.renderBuiltinNodes(acc, ['Debug']);

    acc.end();
    this.treeContainer.appendChild(acc.root);
  }

  /** Adds one toolbox section as an accordion pane. Content builds lazily on
   *  first expand; on a keyed accordion the remembered state wins over the
   *  `expanded` default. While a search is active the pane is forced open (the
   *  accordion is unkeyed then, so nothing persists). Pane headers keep the
   *  `data-section` attribute and `ff-browser-section-<title>` test id the
   *  guide and tests address sections by. */
  private addSection(acc: DG.Accordion, title: string, getContent: () => HTMLElement, opts: {
    expanded?: boolean; count?: number; tooltip?: string; tidParts?: Array<string | number>;
  } = {}): DG.AccordionPane {
    const pane = opts.count !== undefined ?
      acc.addCountPane(title, getContent, () => opts.count!, opts.expanded ?? false, null, false) :
      acc.addPane(title, getContent, opts.expanded ?? false, null, false);
    if (this.searchInput.value) pane.expanded = true;
    const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
    if (header) {
      header.dataset.section = title;
      setTid(header, 'browser-section', title);
      if (opts.tooltip) ui.tooltip.bind(header, opts.tooltip);
    }
    if (opts.tidParts) setTid(pane.root, ...opts.tidParts);
    return pane;
  }

  /** In "what it does" mode, order groups by the curated task sequence
   *  (Data Sources first); in every other mode fall back to alphabetical. */
  private orderGroupKeys(keys: string[]): string[] {
    if (this.groupBy !== 'category') return keys.sort();
    const idx = (k: string): number => {
      const i = (FUNC_CATEGORIES as readonly string[]).indexOf(k);
      return i === -1 ? FUNC_CATEGORIES.length : i;
    };
    return keys.sort((a, b) => idx(a) - idx(b) || a.localeCompare(b));
  }

  /** The KNIME-style Files browser — built once, reused across renders so its
   *  expanded folders and scroll position survive a search keystroke. */
  private getFilesTreeRoot(): HTMLElement {
    if (!this.filesTreeRoot) {
      const tree = getFilesBrowser(
        () => {},                                                  // selection: no-op
        (f) => this.callbacks.onFileDoubleClick(f),                // dbl-click / Enter → OpenFile node
        'funcflow.files.v1');
      this.filesTreeRoot = tree.root;
    }
    return this.filesTreeRoot;
  }

  private renderFilesSection(acc: DG.Accordion): void {
    this.addSection(acc, 'Files',
      () => ui.div([this.getFilesTreeRoot()], 'funcflow-section-content'), {
        expanded: true,
        tooltip: 'Browse data files. Double-click or drag a file onto the canvas to load it.',
        tidParts: ['browser-files'],
      });
  }

  private renderQueriesSection(acc: DG.Accordion): void {
    // Group every registered query by its connection (friendlyName ?? name).
    const queries = getRegisteredFuncs().filter((f) => f.func instanceof DG.DataQuery);
    const query = this.searchInput.value.toLowerCase();
    const matching = query ? queries.filter((f) => funcMatchesSearch(f, query)) : queries;
    if (matching.length === 0) return;

    const byConn = new Map<string, FuncInfo[]>();
    for (const f of matching) {
      const conn = queryConnectionName(f);
      let arr = byConn.get(conn);
      if (!arr) byConn.set(conn, arr = []);
      arr.push(f);
    }

    this.addSection(acc, 'Queries', () => this.buildQueriesContent(byConn), {
      count: matching.length,
      tooltip: 'Database queries, grouped by data connection. Double-click or drag to add.',
      tidParts: ['browser-queries'],
    });
  }

  /** The Queries pane content: a nested accordion with one pane per connection
   *  (its own persistence key — pane names are the connection names). */
  private buildQueriesContent(byConn: Map<string, FuncInfo[]>): HTMLElement {
    const hasSearch = !!this.searchInput.value;
    const inner = ui.accordion(hasSearch ? null : `${TOOLBOX_ACCORDION_KEY}.queries`);
    this.queriesAccordion = inner;
    for (const conn of [...byConn.keys()].sort((a, b) => a.localeCompare(b))) {
      const items = byConn.get(conn)!;
      items.sort((a, b) => a.name.localeCompare(b.name));
      const pane = inner.addCountPane(conn, () => {
        const subContent = ui.div([], 'funcflow-section-content');
        for (const info of items)
          subContent.appendChild(this.createFuncItem(info));
        return subContent;
      }, () => items.length, false, null, false);
      if (hasSearch) pane.expanded = true;
      pane.root.dataset.queryConn = conn;
      setTid(pane.root, 'browser-query-conn', conn);
      const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
      if (header) ui.tooltip.bind(header, `Queries from the “${conn}” connection`);
    }
    inner.end();
    return ui.div([inner.root], 'funcflow-section-content funcflow-subsections');
  }

  /** The Viewers pane — manually-built viewer nodes (core charts first, then
   *  discovered package viewers). Each adds a `Viewers/<label>` node. */
  private renderViewersSection(acc: DG.Accordion): void {
    const query = this.searchInput.value.toLowerCase();
    const types = query ? VIEWER_NODE_TYPES.filter((v) => nameMatchesQuery(v.label, query)) : VIEWER_NODE_TYPES;
    if (types.length === 0) return;

    this.addSection(acc, 'Viewers', () => {
      const content = ui.div([], 'funcflow-section-content');
      for (const vt of types) {
        const item = ui.div([], 'funcflow-func-item');
        item.textContent = vt.label;
        item.dataset.testid = tid('browser-item', vt.nodeTypeName);
        item.dataset.nodeTypeName = vt.nodeTypeName;
        ui.tooltip.bind(item, `Add a ${vt.label} viewer. Wire a table in, then run. Double-click or drag.`);
        item.addEventListener('dblclick', () => this.callbacks.onBuiltinNodeDoubleClick(vt.nodeTypeName));
        this.makeItemDraggable(item, vt.nodeTypeName);
        content.appendChild(item);
      }
      return content;
    }, {
      count: types.length,
      tooltip: 'Charts and viewers. Wire a table into one and run to see it in the preview panel.',
      tidParts: ['browser-viewers'],
    });
  }

  /** The Widgets pane — functions that produce a widget (info panels, search
   *  widgets, …), kept out of the categories. */
  private renderWidgetsSection(acc: DG.Accordion): void {
    const widgets = getRegisteredFuncs().filter(funcOutputsWidget);
    const query = this.searchInput.value.toLowerCase();
    const matching = query ? widgets.filter((f) => funcMatchesSearch(f, query)) : widgets;
    if (matching.length === 0) return;
    matching.sort((a, b) => a.name.localeCompare(b.name));

    this.addSection(acc, 'Widgets', () => {
      const content = ui.div([], 'funcflow-section-content');
      for (const info of matching)
        content.appendChild(this.createFuncItem(info));
      return content;
    }, {
      count: matching.length,
      tooltip: 'Functions that produce a widget. Double-click or drag to add.',
      tidParts: ['browser-widgets'],
    });
  }

  /** Renders the built-in node sections named in `titles` (lets the caller
   *  interleave them with other panes — Debug goes last in the toolbox). */
  private renderBuiltinNodes(acc: DG.Accordion, titles: string[]): void {
    const query = this.searchInput.value.toLowerCase();

    // Input nodes
    const inputNodes = [
      {name: 'Table Input', type: 'Inputs/Table Input', desc: 'Dataframe input parameter'},
      {name: 'Column Input', type: 'Inputs/Column Input', desc: 'Single column from a table'},
      {name: 'Column List Input', type: 'Inputs/Column List Input', desc: 'Multiple columns from a table'},
      {name: 'String Input', type: 'Inputs/String Input', desc: 'Text input parameter'},
      {name: 'Number Input', type: 'Inputs/Number Input', desc: 'Floating-point number input'},
      {name: 'Int Input', type: 'Inputs/Int Input', desc: 'Integer number input'},
      {name: 'Boolean Input', type: 'Inputs/Boolean Input', desc: 'True/false toggle input'},
      {name: 'DateTime Input', type: 'Inputs/DateTime Input', desc: 'Date and time input'},
      {name: 'File Input', type: 'Inputs/File Input', desc: 'File upload input'},
      {name: 'Map Input', type: 'Inputs/Map Input', desc: 'Key-value map input'},
      {name: 'Dynamic Input', type: 'Inputs/Dynamic Input', desc: 'Dynamically-typed input'},
      {name: 'String List Input', type: 'Inputs/String List Input', desc: 'List of strings input'},
      {name: 'Blob Input', type: 'Inputs/Blob Input', desc: 'Binary data input'},
    ];
    const outputNodes = [
      {name: 'Table Output', type: 'Outputs/Table Output', desc: 'Marks a dataframe as script output'},
      {name: 'Value Output', type: 'Outputs/Value Output', desc: 'Marks a value as script output (configurable type)'},
    ];
    const utilityNodes = [
      {name: 'Select Column', type: 'Utilities/Select Column', desc: 'Gets a column from a table by name'},
      {name: 'Select Columns', type: 'Utilities/Select Columns', desc: 'Gets multiple columns by names (comma-separated)'},
      {name: 'Select Table', type: 'Utilities/Select Table', desc: 'Gets an open table by name via grok.shell.tableByName()'},
      {name: 'Add Table View', type: 'Utilities/Add Table View', desc: 'Opens a table in a new view via grok.shell.addTableView()'},
      {name: 'Log', type: 'Utilities/Log', desc: 'Logs a value to the browser console'},
      {name: 'Info', type: 'Utilities/Info', desc: 'Shows an info balloon via grok.shell.info()'},
      {name: 'Warning', type: 'Utilities/Warning', desc: 'Shows a warning balloon via grok.shell.warning()'},
      {name: 'ToString', type: 'Utilities/ToString', desc: 'Converts any value to a string'},
      {name: 'FromJSON', type: 'Utilities/FromJSON', desc: 'Parses a JSON string into an object'},
      {name: 'ToJSON', type: 'Utilities/ToJSON', desc: 'Serializes a value to a JSON string'},
    ];
    const constantNodes = [
      {name: 'String', type: 'Constants/String', desc: 'A constant text value'},
      {name: 'Int', type: 'Constants/Int', desc: 'A constant integer value'},
      {name: 'Double', type: 'Constants/Double', desc: 'A constant floating-point value'},
      {name: 'Boolean', type: 'Constants/Boolean', desc: 'A constant true/false value'},
      {name: 'List', type: 'Constants/List', desc: 'A constant list of comma-separated values'},
    ];
    const debugNodes = [
      {name: 'Breakpoint', type: 'Debug/Breakpoint', desc: 'Pauses execution in debug mode until Continue is clicked'},
    ];

    // All built-in sections start collapsed: they are building-blocks a scientist
    // reaches for deliberately, not the first thing to scan. The data-function
    // categories below (Data Sources first) are what's expanded-ready instead.
    const sections: {title: string; nodes: {name: string; type: string; desc?: string}[]; tip?: string}[] = [
      {title: 'Inputs', nodes: inputNodes, tip: 'Script input parameters (become //input: lines)'},
      {title: 'Outputs', nodes: outputNodes, tip: 'Script output parameters (become //output: lines)'},
      {title: 'Constants', nodes: constantNodes, tip: 'Constant literal values'},
      {title: 'Utilities', nodes: utilityNodes, tip: 'Helper operations (logging, type conversion, etc.)'},
      {title: 'Debug', nodes: debugNodes, tip: 'Debugging and execution control nodes'},
    ];

    for (const section of sections.filter((s) => titles.includes(s.title))) {
      const filtered = query ?
        section.nodes.filter((n) => nameMatchesQuery(n.name, query)) :
        section.nodes;
      if (filtered.length === 0) continue;

      this.addSection(acc, section.title, () => {
        const content = ui.div([], 'funcflow-section-content');
        for (const node of filtered) {
          const item = ui.div([], 'funcflow-func-item');
          item.textContent = node.name;
          item.dataset.testid = tid('browser-item', node.type);
          item.dataset.nodeTypeName = node.type; // e.g. "Inputs/Table Input"
          const tip = node.desc ?
            `${node.desc}. Double-click or drag to add` :
            `Double-click or drag to add ${node.name}`;
          ui.tooltip.bind(item, tip);
          item.addEventListener('dblclick', () => {
            this.callbacks.onBuiltinNodeDoubleClick(node.type);
          });
          this.makeItemDraggable(item, node.type);
          content.appendChild(item);
        }
        return content;
      }, {tooltip: section.tip, count: filtered.length});
    }
  }

  /** One draggable, double-clickable catalog row for a DG function. Shared by
   *  the category sections and the Queries pane. */
  private createFuncItem(info: FuncInfo): HTMLElement {
    const item = ui.div([], 'funcflow-func-item');
    item.textContent = info.name;
    item.dataset.testid = tid('browser-item', info.nodeTypeName);
    item.dataset.nodeTypeName = info.nodeTypeName;     // registered factory name
    item.dataset.func = info.func.name;                // the DG function it adds
    if (info.packageName) item.dataset.package = info.packageName;
    let tip = info.func.description || info.name;
    if (info.packageName)
      tip += ` (${info.packageName})`;
    ui.tooltip.bind(item, tip);
    item.addEventListener('dblclick', () => {
      this.callbacks.onFunctionDoubleClick(info);
    });
    this.makeItemDraggable(item, info.nodeTypeName);
    return item;
  }

  /** Wire HTML5 drag on a browser item so dropping it on the canvas creates
   *  the matching node at the drop point. */
  private makeItemDraggable(item: HTMLElement, typeName: string): void {
    item.draggable = true;
    item.addEventListener('dragstart', (ev) => {
      if (!ev.dataTransfer) return;
      ev.dataTransfer.setData(FF_DRAG_MIME, typeName);
      ev.dataTransfer.setData('text/plain', typeName);
      ev.dataTransfer.effectAllowed = 'copy';
      item.classList.add('funcflow-func-item-dragging');
    });
    item.addEventListener('dragend', () => {
      item.classList.remove('funcflow-func-item-dragging');
    });
  }

  private filterBySearch(funcs: FuncInfo[]): FuncInfo[] {
    const query = this.searchInput.value.toLowerCase();
    if (!query) return funcs;
    return funcs.filter((f) => funcMatchesSearch(f, query));
  }

  private groupFunctions(funcs: FuncInfo[]): Record<string, FuncInfo[]> {
    const groups: Record<string, FuncInfo[]> = {};

    for (const f of funcs) {
      let keys: string[];
      // Saved flows get their own 'Workflows' section in EVERY grouping — their
      // role/tags/package say nothing useful (a script from whatever package
      // the user happened to save it in).
      if (isWorkflowFunc(f.func))
        keys = ['Workflows'];
      else {
        switch (this.groupBy) {
        case 'role':
          keys = [f.role || 'Uncategorized'];
          break;
        case 'tags':
          keys = f.tags.length > 0 ? f.tags : ['Untagged'];
          break;
        case 'package':
          keys = [f.packageName || 'Core'];
          break;
        case 'category':
          keys = [categorizeFunc(f.func, f.role, f.packageName)];
          break;
        default:
          keys = ['Other'];
        }
      }
      for (const key of keys) {
        if (!groups[key]) groups[key] = [];
        groups[key].push(f);
      }
    }

    return groups;
  }
}

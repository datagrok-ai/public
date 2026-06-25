/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FuncInfo, getRegisteredFuncs} from '../rete/node-factory';

export type GroupByMode = 'category' | 'role' | 'tags' | 'package';

/** The task-oriented categories, in the order a scientist meets a pipeline:
 *  get data in → combine → reshape → derive columns → compute values →
 *  visualize → everything else. Used both to bucket and to order the tree. */
export const FUNC_CATEGORIES = [
  'Data Sources',
  'Combine Tables',
  'Transform Tables',
  'Column Operations',
  'Compute Values',
  'Visualize',
  'Other',
] as const;
export type FuncCategory = (typeof FUNC_CATEGORIES)[number];

const VIS_TYPES = ['viewer', 'view', 'widget', 'graphics'];
const SCALAR_TYPES = ['string', 'int', 'double', 'bool', 'datetime', 'num', 'bigint', 'qnum'];
const COL_TYPES = ['column', 'column_list'];

/** Bucket a function by **what it does**, derived from its input/output
 *  signature (and role for viewers). The key distinctions, validated against
 *  the live catalog (see docs/func-catalog-snapshot.md):
 *  - **Data Sources** produce a table from *no* table input (OpenFile, DB
 *    queries, generators). A join is NOT a data source — it consumes tables.
 *  - **Combine Tables** take ≥2 tables (join / union / link / react).
 *  - **Transform Tables** take exactly one table and return a table or mutate
 *    it in place (filter, aggregate, pivot, sort, …).
 *  - **Column Operations** emit a column / column list (AddNewColumn, descriptors).
 *  - **Compute Values** emit a scalar (statistics, math, text, predicates).
 *  - **Visualize** emit a viewer / view / widget / graphics.
 *  - **Other** — expression/filter builders and dynamic helpers. */
export function categorizeFunc(func: DG.Func, role: string | null): FuncCategory {
  const outs = func.outputs.map((o) => String(o.propertyType));
  const ins = func.inputs.map((i) => String(i.propertyType));
  const has = (arr: string[], set: string[]): boolean => arr.some((t) => set.includes(t));

  const dfIn = ins.filter((t) => t === 'dataframe').length;
  const outDf = outs.includes('dataframe');
  const noOut = outs.length === 0;
  const roleHasViewer = !!role && role.split(',').some((r) => r.trim() === 'viewer');

  if (has(outs, VIS_TYPES) || roleHasViewer) return 'Visualize';
  if (dfIn >= 2) return 'Combine Tables';                          // join / union / link
  if (outDf && dfIn === 0) return 'Data Sources';                  // produce a table from non-table inputs
  if ((outDf && dfIn === 1) || (noOut && dfIn >= 1)) return 'Transform Tables';
  if (has(outs, COL_TYPES)) return 'Column Operations';
  if (has(outs, SCALAR_TYPES)) return 'Compute Values';
  return 'Other';
}

/** MIME-ish key used to carry a node type name through HTML5 drag/drop.
 *  The canvas drop handler reads this and adds the corresponding node at
 *  the drop point. */
export const FF_DRAG_MIME = 'application/x-funcflow-node';

export interface FunctionBrowserCallbacks {
  onFunctionDoubleClick: (funcInfo: FuncInfo) => void;
  onBuiltinNodeDoubleClick: (nodeTypeName: string) => void;
}

/** Persisted UI state for the browser (group mode + which sections the user
 *  has expanded), so it reopens exactly as the user left it. */
interface BrowserState {
  groupBy: GroupByMode;
  expanded: Record<string, boolean>;
}
const LS_KEY = 'funcflow.browser.v1';
const VALID_MODES: GroupByMode[] = ['category', 'role', 'tags', 'package'];

/** Left sidebar: searchable, groupable function catalog */
export class FunctionBrowser {
  root: HTMLElement;
  private searchInput!: HTMLInputElement;
  private groupBySelect!: HTMLSelectElement;
  private treeContainer!: HTMLElement;
  private groupBy: GroupByMode = 'category';
  /** Expanded-section memory, keyed `b:<title>` (built-ins) / `<mode>:<group>`. */
  private expanded: Record<string, boolean> = {};
  private callbacks: FunctionBrowserCallbacks;

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
      if (s.expanded && typeof s.expanded === 'object') this.expanded = s.expanded;
    } catch {/* corrupt/blocked storage — fall back to defaults */}
  }

  private saveState(): void {
    try {
      const state: BrowserState = {groupBy: this.groupBy, expanded: this.expanded};
      localStorage.setItem(LS_KEY, JSON.stringify(state));
    } catch {/* storage blocked/full — non-fatal */}
  }

  /** Whether a section should render expanded: forced open while searching,
   *  otherwise remembered (default collapsed). */
  private isExpanded(key: string, hasSearch: boolean): boolean {
    if (hasSearch) return true;
    return this.expanded[key] === true;
  }

  private setExpanded(key: string, value: boolean): void {
    this.expanded[key] = value;
    this.saveState();
  }

  private buildUI(): HTMLElement {
    // Search bar
    this.searchInput = document.createElement('input');
    this.searchInput.type = 'text';
    this.searchInput.placeholder = 'Search functions...';
    this.searchInput.className = 'funcflow-search-input';
    this.searchInput.addEventListener('input', () => this.render());

    // Group by selector
    this.groupBySelect = document.createElement('select');
    this.groupBySelect.className = 'funcflow-groupby-select';
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

    this.treeContainer = ui.div([], 'funcflow-tree-container');

    const container = ui.divV([
      this.searchInput,
      this.groupBySelect,
      this.treeContainer,
    ], 'funcflow-browser');

    return container;
  }

  render(): void {
    this.treeContainer.innerHTML = '';

    // Add built-in nodes section first
    this.renderBuiltinNodes();

    // Add DG function nodes
    const funcs = this.filterBySearch(getRegisteredFuncs());
    const grouped = this.groupFunctions(funcs);

    const sortedKeys = this.orderGroupKeys(Object.keys(grouped));
    for (const category of sortedKeys) {
      const items = grouped[category];
      if (items.length === 0) continue;
      const section = this.createCollapsibleSection(category, items);
      this.treeContainer.appendChild(section);
    }
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

  private renderBuiltinNodes(): void {
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
    const comparisonNodes = [
      {name: 'Equals (==)', type: 'Comparisons/Equals (==)', desc: 'Tests if two values are equal'},
      {name: 'Not Equals (!=)', type: 'Comparisons/Not Equals (!=)', desc: 'Tests if two values are not equal'},
      {name: 'Greater Than (>)', type: 'Comparisons/Greater Than (>)', desc: 'Tests if left is greater than right'},
      {name: 'Greater Or Equal (>=)', type: 'Comparisons/Greater Or Equal (>=)', desc: 'Tests if left is greater than or equal to right'},
      {name: 'Less Than (<)', type: 'Comparisons/Less Than (<)', desc: 'Tests if left is less than right'},
      {name: 'Less Or Equal (<=)', type: 'Comparisons/Less Or Equal (<=)', desc: 'Tests if left is less than or equal to right'},
      {name: 'Contains', type: 'Comparisons/Contains', desc: 'Tests if a string contains a substring'},
      {name: 'Starts With', type: 'Comparisons/Starts With', desc: 'Tests if a string starts with a prefix'},
      {name: 'Ends With', type: 'Comparisons/Ends With', desc: 'Tests if a string ends with a suffix'},
      {name: 'Is Null', type: 'Comparisons/Is Null', desc: 'Tests if a value is null or undefined'},
    ];

    const debugNodes = [
      {name: 'Breakpoint', type: 'Debug/Breakpoint', desc: 'Pauses execution in debug mode until Continue is clicked'},
    ];

    // All built-in sections start collapsed: they are building-blocks a scientist
    // reaches for deliberately, not the first thing to scan. The data-function
    // categories below (Data Sources first) are what's expanded-ready instead.
    const sections: {title: string; nodes: {name: string; type: string; desc?: string}[]; collapsed?: boolean; tip?: string}[] = [
      {title: 'Inputs', nodes: inputNodes, collapsed: true, tip: 'Script input parameters (become //input: lines)'},
      {title: 'Outputs', nodes: outputNodes, collapsed: true, tip: 'Script output parameters (become //output: lines)'},
      {title: 'Constants', nodes: constantNodes, collapsed: true, tip: 'Constant literal values'},
      {title: 'Comparisons', nodes: comparisonNodes, collapsed: true, tip: 'Comparison and logical operators'},
      {title: 'Utilities', nodes: utilityNodes, collapsed: true, tip: 'Helper operations (logging, type conversion, etc.)'},
      {title: 'Debug', nodes: debugNodes, collapsed: true, tip: 'Debugging and execution control nodes'},
    ];

    for (const section of sections) {
      const filtered = query ?
        section.nodes.filter((n) => n.name.toLowerCase().includes(query)) :
        section.nodes;
      if (filtered.length === 0) continue;

      const sectionEl = this.createBuiltinSection(section.title, filtered, section.tip);
      this.treeContainer.appendChild(sectionEl);
    }
  }

  private createBuiltinSection(title: string, nodes: {name: string; type: string; desc?: string}[], tooltip?: string): HTMLElement {
    const header = ui.div([], 'funcflow-section-header');
    header.textContent = title;
    header.style.cursor = 'pointer';
    if (tooltip)
      ui.tooltip.bind(header, tooltip);

    const content = ui.div([], 'funcflow-section-content');
    for (const node of nodes) {
      const item = ui.div([], 'funcflow-func-item');
      item.textContent = node.name;
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

    const key = `b:${title}`;
    const hasSearch = !!this.searchInput.value;
    let collapsed = !this.isExpanded(key, hasSearch);
    if (collapsed) {
      content.style.display = 'none';
      header.classList.add('collapsed');
    }
    header.addEventListener('click', () => {
      collapsed = !collapsed;
      content.style.display = collapsed ? 'none' : 'block';
      header.classList.toggle('collapsed', collapsed);
      if (!hasSearch) this.setExpanded(key, !collapsed); // remember explicit user toggles
    });

    return ui.divV([header, content]);
  }

  private createCollapsibleSection(category: string, items: FuncInfo[]): HTMLElement {
    const header = ui.div([], 'funcflow-section-header');
    header.textContent = `${category} (${items.length})`;
    header.style.cursor = 'pointer';

    const content = ui.div([], 'funcflow-section-content');

    for (const info of items) {
      const item = ui.div([], 'funcflow-func-item');
      item.textContent = info.name;
      let tip = info.func.description || info.name;
      if (info.packageName)
        tip += ` (${info.packageName})`;
      ui.tooltip.bind(item, tip);
      item.addEventListener('dblclick', () => {
        this.callbacks.onFunctionDoubleClick(info);
      });
      this.makeItemDraggable(item, info.nodeTypeName);
      content.appendChild(item);
    }

    const key = `${this.groupBy}:${category}`;
    const hasSearch = !!this.searchInput.value;
    let collapsed = !this.isExpanded(key, hasSearch);
    content.style.display = collapsed ? 'none' : 'block';
    header.classList.toggle('collapsed', collapsed);
    header.addEventListener('click', () => {
      collapsed = !collapsed;
      content.style.display = collapsed ? 'none' : 'block';
      header.classList.toggle('collapsed', collapsed);
      if (!hasSearch) this.setExpanded(key, !collapsed);
    });

    return ui.divV([header, content]);
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
    return funcs.filter((f) => {
      return f.name.toLowerCase().includes(query) ||
        (f.func.description || '').toLowerCase().includes(query) ||
        f.tags.some((t) => t.toLowerCase().includes(query)) ||
        (f.role || '').toLowerCase().includes(query) ||
        f.packageName.toLowerCase().includes(query);
    });
  }

  private groupFunctions(funcs: FuncInfo[]): Record<string, FuncInfo[]> {
    const groups: Record<string, FuncInfo[]> = {};

    for (const f of funcs) {
      let keys: string[];
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
        keys = [categorizeFunc(f.func, f.role)];
        break;
      default:
        keys = ['Other'];
      }
      for (const key of keys) {
        if (!groups[key]) groups[key] = [];
        groups[key].push(f);
      }
    }

    return groups;
  }
}

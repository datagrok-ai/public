/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FuncInfo, getRegisteredFuncs} from '../rete/node-factory';

export type GroupByMode = 'role' | 'tags' | 'package' | 'output';

/** Bucket a function by what kind of value it produces.
 *  Visualization wins over Data Sources when both apply (a func that builds
 *  a viewer is more naturally found under "Visualization"). No-output funcs
 *  are treated as Transformations — they typically mutate their dataframe
 *  input in place (e.g. addNewColumn, fillNullValues). */
function categorizeByOutput(func: DG.Func): string {
  const outputs = func.outputs;
  if (outputs.length === 0) return 'Transformations';
  const types = outputs.map((o) => String(o.propertyType));
  const has = (t: string): boolean => types.includes(t);
  if (has('viewer') || has('view') || has('widget') || has('graphics')) return 'Visualization';
  if (has('dataframe')) return 'Data Sources';
  if (has('column') || has('column_list')) return 'Column Operations';
  if (has('string') || has('int') || has('double') || has('bool') ||
      has('datetime') || has('num')) return 'Compute';
  return 'Utilities';
}

/** MIME-ish key used to carry a node type name through HTML5 drag/drop.
 *  The canvas drop handler reads this and adds the corresponding node at
 *  the drop point. */
export const FF_DRAG_MIME = 'application/x-funcflow-node';

export interface FunctionBrowserCallbacks {
  onFunctionDoubleClick: (funcInfo: FuncInfo) => void;
  onBuiltinNodeDoubleClick: (nodeTypeName: string) => void;
}

/** Left sidebar: searchable, groupable function catalog */
export class FunctionBrowser {
  root: HTMLElement;
  private searchInput!: HTMLInputElement;
  private groupBySelect!: HTMLSelectElement;
  private treeContainer!: HTMLElement;
  private groupBy: GroupByMode = 'role';
  private callbacks: FunctionBrowserCallbacks;

  constructor(callbacks: FunctionBrowserCallbacks) {
    this.callbacks = callbacks;
    this.root = this.buildUI();
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
      role: 'role',
      tags: 'tags',
      package: 'package',
      output: 'output (data sources / transformations / …)',
    };
    for (const opt of Object.keys(groupByLabels) as GroupByMode[]) {
      const option = document.createElement('option');
      option.value = opt;
      option.textContent = `Group by: ${groupByLabels[opt]}`;
      this.groupBySelect.appendChild(option);
    }
    this.groupBySelect.addEventListener('change', () => {
      this.groupBy = this.groupBySelect.value as GroupByMode;
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

    const sortedKeys = Object.keys(grouped).sort();
    for (const category of sortedKeys) {
      const items = grouped[category];
      if (items.length === 0) continue;
      const section = this.createCollapsibleSection(category, items);
      this.treeContainer.appendChild(section);
    }
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

    const sections: {title: string; nodes: {name: string; type: string; desc?: string}[]; collapsed?: boolean; tip?: string}[] = [
      {title: 'Inputs', nodes: inputNodes, tip: 'Script input parameters (become //input: lines)'},
      {title: 'Outputs', nodes: outputNodes, tip: 'Script output parameters (become //output: lines)'},
      {title: 'Constants', nodes: constantNodes, tip: 'Constant literal values'},
      {title: 'Comparisons', nodes: comparisonNodes, collapsed: true, tip: 'Comparison and logical operators'},
      {title: 'Utilities', nodes: utilityNodes, tip: 'Helper operations (logging, type conversion, etc.)'},
      {title: 'Debug', nodes: debugNodes, tip: 'Debugging and execution control nodes'},
    ];

    for (const section of sections) {
      const filtered = query ?
        section.nodes.filter((n) => n.name.toLowerCase().includes(query)) :
        section.nodes;
      if (filtered.length === 0) continue;

      const sectionEl = this.createBuiltinSection(section.title, filtered, section.collapsed, section.tip);
      this.treeContainer.appendChild(sectionEl);
    }
  }

  private createBuiltinSection(title: string, nodes: {name: string; type: string; desc?: string}[], startCollapsed?: boolean, tooltip?: string): HTMLElement {
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

    const hasSearch = !!this.searchInput.value;
    let collapsed = hasSearch ? false : !!startCollapsed;
    if (collapsed) {
      content.style.display = 'none';
      header.classList.add('collapsed');
    }
    header.addEventListener('click', () => {
      collapsed = !collapsed;
      content.style.display = collapsed ? 'none' : 'block';
      header.classList.toggle('collapsed', collapsed);
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

    const hasSearch = !!this.searchInput.value;
    let collapsed = !hasSearch;
    if (collapsed) {
      content.style.display = 'none';
      header.classList.add('collapsed');
    } else {
      content.style.display = 'block';
    }
    header.addEventListener('click', () => {
      collapsed = !collapsed;
      content.style.display = collapsed ? 'none' : 'block';
      header.classList.toggle('collapsed', collapsed);
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
      case 'output':
        keys = [categorizeByOutput(f.func)];
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

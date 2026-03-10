/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FuncInfo, getRegisteredFuncs} from '../nodes/node-factory';
import {getRole, getTags} from '../utils/dart-proxy-utils';

export type GroupByMode = 'role' | 'tags' | 'package';

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
    for (const opt of ['role', 'tags', 'package']) {
      const option = document.createElement('option');
      option.value = opt;
      option.textContent = `Group by: ${opt}`;
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
      {name: 'Table Input', type: 'Inputs/Table Input'},
      {name: 'Column Input', type: 'Inputs/Column Input'},
      {name: 'Column List Input', type: 'Inputs/Column List Input'},
      {name: 'String Input', type: 'Inputs/String Input'},
      {name: 'Number Input', type: 'Inputs/Number Input'},
      {name: 'Int Input', type: 'Inputs/Int Input'},
      {name: 'Boolean Input', type: 'Inputs/Boolean Input'},
      {name: 'DateTime Input', type: 'Inputs/DateTime Input'},
      {name: 'File Input', type: 'Inputs/File Input'},
      {name: 'Map Input', type: 'Inputs/Map Input'},
      {name: 'Dynamic Input', type: 'Inputs/Dynamic Input'},
      {name: 'String List Input', type: 'Inputs/String List Input'},
      {name: 'Blob Input', type: 'Inputs/Blob Input'},
    ];
    const outputNodes = [
      {name: 'Table Output', type: 'Outputs/Table Output'},
      {name: 'Value Output', type: 'Outputs/Value Output'},
    ];
    const utilityNodes = [
      {name: 'Select Column', type: 'Utilities/Select Column'},
      {name: 'Select Columns', type: 'Utilities/Select Columns'},
      {name: 'Add Table View', type: 'Utilities/Add Table View'},
      {name: 'Log', type: 'Utilities/Log'},
      {name: 'Info', type: 'Utilities/Info'},
      {name: 'Warning', type: 'Utilities/Warning'},
      {name: 'ToString', type: 'Utilities/ToString'},
      {name: 'FromJSON', type: 'Utilities/FromJSON'},
      {name: 'ToJSON', type: 'Utilities/ToJSON'},
    ];
    const constantNodes = [
      {name: 'String', type: 'Constants/String'},
      {name: 'Int', type: 'Constants/Int'},
      {name: 'Double', type: 'Constants/Double'},
      {name: 'Boolean', type: 'Constants/Boolean'},
      {name: 'List', type: 'Constants/List'},
    ];
    const comparisonNodes = [
      {name: 'Equals (==)', type: 'Comparisons/Equals (==)'},
      {name: 'Not Equals (!=)', type: 'Comparisons/Not Equals (!=)'},
      {name: 'Greater Than (>)', type: 'Comparisons/Greater Than (>)'},
      {name: 'Greater Or Equal (>=)', type: 'Comparisons/Greater Or Equal (>=)'},
      {name: 'Less Than (<)', type: 'Comparisons/Less Than (<)'},
      {name: 'Less Or Equal (<=)', type: 'Comparisons/Less Or Equal (<=)'},
      {name: 'Contains', type: 'Comparisons/Contains'},
      {name: 'Starts With', type: 'Comparisons/Starts With'},
      {name: 'Ends With', type: 'Comparisons/Ends With'},
      {name: 'Is Null', type: 'Comparisons/Is Null'},
    ];

    const sections: {title: string; nodes: {name: string; type: string}[]; collapsed?: boolean; tip?: string}[] = [
      {title: 'Inputs', nodes: inputNodes, tip: 'Script input parameters (become //input: lines)'},
      {title: 'Outputs', nodes: outputNodes, tip: 'Script output parameters (become //output: lines)'},
      {title: 'Constants', nodes: constantNodes, tip: 'Constant literal values'},
      {title: 'Comparisons', nodes: comparisonNodes, collapsed: true, tip: 'Comparison and logical operators'},
      {title: 'Utilities', nodes: utilityNodes, tip: 'Helper operations (logging, type conversion, etc.)'},
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

  private createBuiltinSection(title: string, nodes: {name: string; type: string}[], startCollapsed?: boolean, tooltip?: string): HTMLElement {
    const header = ui.div([], 'funcflow-section-header');
    header.textContent = title;
    header.style.cursor = 'pointer';
    if (tooltip)
      ui.tooltip.bind(header, tooltip);

    const content = ui.div([], 'funcflow-section-content');
    for (const node of nodes) {
      const item = ui.div([], 'funcflow-func-item');
      item.textContent = node.name;
      ui.tooltip.bind(item, `Double-click to add ${node.name}`);
      item.addEventListener('dblclick', () => {
        this.callbacks.onBuiltinNodeDoubleClick(node.type);
      });
      content.appendChild(item);
    }

    let collapsed = !!startCollapsed;
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
    // Start collapsed for DG functions since there can be many
    content.style.display = 'none';

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
      content.appendChild(item);
    }

    let collapsed = true;
    header.classList.add('collapsed');
    header.addEventListener('click', () => {
      collapsed = !collapsed;
      content.style.display = collapsed ? 'none' : 'block';
      header.classList.toggle('collapsed', collapsed);
    });

    return ui.divV([header, content]);
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

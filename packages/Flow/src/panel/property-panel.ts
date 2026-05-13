/* eslint-disable max-len */
/** Property panel — renders into Datagrok's native context panel.
 *
 * Reads from a `FlowNode` selected on the canvas and edits its `properties`
 * and `inputValues` maps. Re-renders the node on the canvas whenever a
 * property changes so visible state (title, etc.) stays in sync. */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {NodeExecState} from '../execution/execution-state';
import {buildExecutionMeta} from '../execution/value-inspector';

const PROP_TOOLTIPS: Record<string, string> = {
  'Title': 'Display name shown on the node',
  'Param Name': 'Variable name used in the generated script',
  'Description': 'Annotation rendered under the node title; for input/output nodes, also embedded in the //input:/output: line',
  'Default': 'Default value when no input is provided',
  'Nullable': 'Allow null/empty values for this input',
  'SemType': 'Semantic type annotation (e.g. Molecule)',
  'SemType Filter': 'Only show columns matching this semantic type',
  'Type Filter': 'Only show columns matching this data type',
  'Show Slider': 'Display a slider control in the script run dialog',
  'Choices (comma-sep)': 'Restrict input to predefined values, separated by commas',
  'Output Type': 'Datagrok type for this output parameter',
  'Caption': 'Display label shown to users in the script run dialog',
  'Min': 'Minimum allowed value',
  'Max': 'Maximum allowed value',
};

const UTILITY_PROP_TOOLTIPS: Record<string, Record<string, string>> = {
  'Select Column': {columnName: 'Column name to extract from the table'},
  'Select Columns': {columnNames: 'Comma-separated column names to extract'},
  'Log': {label: 'Optional label prefix for the log message'},
  'String': {value: 'The constant string value to output'},
  'Int': {value: 'The constant integer value to output'},
  'Double': {value: 'The constant floating-point value to output'},
  'Boolean': {value: 'The constant true/false value to output'},
  'List': {value: 'Comma-separated list of values'},
};

const TYPE_FILTER_VALUES = ['', 'numerical', 'categorical', 'string', 'int', 'double', 'bool'];
const SEMTYPE_VALUES = ['', 'Molecule', 'Macromolecule'];
const OUTPUT_TYPE_VALUES = [
  'string', 'int', 'double', 'bool',
  'dataframe', 'column', 'column_list',
  'object', 'dynamic', 'list',
  'view', 'viewer', 'widget',
  'graphics', 'grid_cell_renderer', 'filter',
  'map', 'datetime', 'blob', 'funccall',
];

function buildFuncInputTooltip(param: DG.Property): string {
  const parts: string[] = [];
  const desc = param.description;
  if (desc) parts.push(desc);
  parts.push(`Type: ${param.propertyType}`);
  if (param.defaultValue !== undefined && param.defaultValue !== null && param.defaultValue !== '')
    parts.push(`Default: ${param.defaultValue}`);
  if (param.nullable) parts.push('Nullable');
  return parts.join(' | ');
}

export class PropertyPanel {
  root: HTMLElement;
  private contentDiv: HTMLElement;
  private flow: FlowEditor;

  constructor(flow: FlowEditor) {
    this.flow = flow;
    this.contentDiv = ui.div([], 'funcflow-property-content');
    this.root = ui.divV([this.contentDiv], 'funcflow-property-panel');
  }

  showNode(node: FlowNode, execState?: NodeExecState): void {
    this.contentDiv.innerHTML = '';

    const titleInput = this.createTextarea('Title', node.label, (v) => {
      node.label = v;
      void this.flow.updateNode(node.id);
    });
    const typeBadge = ui.div([], 'funcflow-type-badge');
    typeBadge.textContent = node.dgNodeType || 'function';
    const titleRow = ui.div([titleInput, typeBadge], 'funcflow-title-row');
    this.contentDiv.appendChild(titleRow);

    // Per-node description: rendered under the title in the canvas, and
    // embedded as the [description] suffix in //input:/output: lines.
    this.contentDiv.appendChild(this.createTextarea('Description', node.description, (v) => {
      node.description = v;
      void this.flow.updateNode(node.id);
    }));

    const acc = ui.accordion('funcflow-context-panel');

    if (node.dgFunc) this.addFuncNodePanes(acc, node);
    if (node.dgNodeType === 'input') this.addInputNodePane(acc, node);
    if (node.dgNodeType === 'output') this.addOutputNodePane(acc, node);
    if (node.dgNodeType === 'utility') this.addUtilityNodePane(acc, node);

    this.addConnectionsPane(acc, node);

    this.contentDiv.appendChild(acc.root);

    // Execution metadata — status / duration / per-output dims / error.
    // Rich previews (grid, sample, image) live in the bottom-docked panel
    // instead, to keep this panel narrow.
    if (execState) {
      const header = ui.div([], 'funcflow-prop-section-header');
      header.textContent = 'Execution';
      this.contentDiv.appendChild(header);
      this.contentDiv.appendChild(buildExecutionMeta(execState));
    }
  }

  clear(): void {
    this.contentDiv.innerHTML = '';
    this.contentDiv.appendChild(ui.divText('Select a node to view its properties'));
  }

  // ---------- panes ----------

  private addFuncNodePanes(acc: DG.Accordion, node: FlowNode): void {
    const func = node.dgFunc;
    if (!func) return;

    acc.addPane('Function', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      if (func.description)
        content.appendChild(ui.div([ui.label('Description'), ui.divText(func.description)], 'funcflow-prop-row'));
      content.appendChild(ui.div([ui.label('Full Name'), ui.divText(node.dgFuncName ?? func.name)], 'funcflow-prop-row'));
      if (node.dgRole)
        content.appendChild(ui.div([ui.label('Role'), ui.divText(node.dgRole)], 'funcflow-prop-row'));
      return content;
    }, true);

    if (func.inputs.length > 0) {
      acc.addPane('Input Parameters', () => {
        const content = ui.div([], 'funcflow-accordion-content');
        for (const inp of func.inputs) {
          const tip = buildFuncInputTooltip(inp);
          const isPrimitive = inp.name in node.inputValues;
          if (!isPrimitive) {
            const row = ui.div([ui.divText(`${inp.name}: ${inp.propertyType} (connected only)`)], 'funcflow-prop-row');
            ui.tooltip.bind(row, tip);
            content.appendChild(row);
            continue;
          }
          if (this.flow.isInputConnected(node.id, inp.name)) {
            const row = ui.div([ui.divText(`${inp.name}: connected`)], 'funcflow-prop-row');
            ui.tooltip.bind(row, tip);
            content.appendChild(row);
            continue;
          }
          switch (inp.propertyType) {
          case 'string':
            content.appendChild(this.createTextarea(inp.name, String(node.inputValues[inp.name] ?? ''),
              (v) => {node.inputValues[inp.name] = v;}, tip));
            break;
          case 'int':
            content.appendChild(this.createNumberInput(inp.name, Number(node.inputValues[inp.name] ?? 0),
              (v) => {node.inputValues[inp.name] = Math.round(v);}, 0, 1, tip));
            break;
          case 'double':
          case 'num':
            content.appendChild(this.createNumberInput(inp.name, Number(node.inputValues[inp.name] ?? 0),
              (v) => {node.inputValues[inp.name] = v;}, 3, 0.1, tip));
            break;
          case 'bool':
            content.appendChild(this.createToggle(inp.name, Boolean(node.inputValues[inp.name]),
              (v) => {node.inputValues[inp.name] = v;}, tip));
            break;
          }
        }
        return content;
      }, true);
    }
  }

  // eslint-disable-next-line complexity
  private addInputNodePane(acc: DG.Accordion, node: FlowNode): void {
    acc.addPane('Input Configuration', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      content.appendChild(this.createTextarea('Param Name', String(node.properties['paramName'] ?? ''), (v) => {
        node.properties['paramName'] = v;
      }));

      const outputType = node.dgOutputType;
      if (outputType === 'bool') {
        content.appendChild(this.createToggle('Default', Boolean(node.properties['defaultValue']), (v) => {
          node.properties['defaultValue'] = v;
        }));
      } else if (outputType === 'int') {
        content.appendChild(this.createNumberInput('Default', Number(node.properties['defaultValue'] ?? 0), (v) => {
          node.properties['defaultValue'] = Math.round(v);
        }, 0, 1));
      } else if (outputType === 'double') {
        content.appendChild(this.createNumberInput('Default', Number(node.properties['defaultValue'] ?? 0), (v) => {
          node.properties['defaultValue'] = v;
        }, 3, 0.1));
      } else if (node.properties['defaultValue'] !== undefined && outputType !== 'dataframe' &&
                 outputType !== 'file' && outputType !== 'map' && outputType !== 'blob') {
        content.appendChild(this.createTextarea('Default', String(node.properties['defaultValue'] ?? ''), (v) => {
          node.properties['defaultValue'] = v;
        }));
      }

      if (node.properties['nullable'] !== undefined)
        content.appendChild(this.createToggle('Nullable', Boolean(node.properties['nullable']), (v) => {node.properties['nullable'] = v;}));
      if (node.properties['caption'] !== undefined)
        content.appendChild(this.createTextarea('Caption', String(node.properties['caption'] ?? ''), (v) => {node.properties['caption'] = v;}));
      if (node.properties['typeFilter'] !== undefined)
        content.appendChild(this.createCombo('Type Filter', String(node.properties['typeFilter'] ?? ''), TYPE_FILTER_VALUES, (v) => {node.properties['typeFilter'] = v;}));
      if (node.properties['semTypeFilter'] !== undefined)
        content.appendChild(this.createTextarea('SemType Filter', String(node.properties['semTypeFilter'] ?? ''), (v) => {node.properties['semTypeFilter'] = v;}));
      if (node.properties['semType'] !== undefined)
        content.appendChild(this.createCombo('SemType', String(node.properties['semType'] ?? ''), SEMTYPE_VALUES, (v) => {node.properties['semType'] = v;}));
      if (node.properties['choices'] !== undefined)
        content.appendChild(this.createTextarea('Choices (comma-sep)', String(node.properties['choices'] ?? ''), (v) => {node.properties['choices'] = v;}));
      if (node.properties['min'] !== undefined)
        content.appendChild(this.createTextarea('Min', String(node.properties['min'] ?? ''), (v) => {node.properties['min'] = v;}));
      if (node.properties['max'] !== undefined)
        content.appendChild(this.createTextarea('Max', String(node.properties['max'] ?? ''), (v) => {node.properties['max'] = v;}));
      if (node.properties['showSlider'] !== undefined)
        content.appendChild(this.createToggle('Show Slider', Boolean(node.properties['showSlider']), (v) => {node.properties['showSlider'] = v;}));

      return content;
    }, true);
  }

  private addOutputNodePane(acc: DG.Accordion, node: FlowNode): void {
    acc.addPane('Output Configuration', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      content.appendChild(this.createTextarea('Param Name', String(node.properties['paramName'] ?? ''), (v) => {
        node.properties['paramName'] = v;
      }));
      if (node.properties['outputType'] !== undefined) {
        content.appendChild(this.createCombo('Output Type', String(node.properties['outputType'] ?? 'dynamic'),
          OUTPUT_TYPE_VALUES, (v) => {node.properties['outputType'] = v;}));
      }
      return content;
    }, true);
  }

  private addUtilityNodePane(acc: DG.Accordion, node: FlowNode): void {
    const props = Object.entries(node.properties).filter(([k]) => !k.startsWith('_'));
    if (props.length === 0) return;

    acc.addPane('Configuration', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      const nodeTips = UTILITY_PROP_TOOLTIPS[node.label] ?? {};
      for (const [key, val] of props) {
        const tip = nodeTips[key];
        if (typeof val === 'boolean')
          content.appendChild(this.createToggle(key, val, (v) => {node.properties[key] = v;}, tip));
        else if (typeof val === 'number') {
          const isInt = Number.isInteger(val);
          content.appendChild(this.createNumberInput(key, val,
            (v) => {node.properties[key] = isInt ? Math.round(v) : v;},
            isInt ? 0 : 3, isInt ? 1 : 0.1, tip));
        } else
          content.appendChild(this.createTextarea(key, String(val ?? ''), (v) => {node.properties[key] = v;}, tip));
      }
      return content;
    }, true);
  }

  private addConnectionsPane(acc: DG.Accordion, node: FlowNode): void {
    acc.addPane('Connections', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      const ptCount = node.passthroughCount;

      const inputEntries = Object.entries(node.inputs) as Array<[string, {socket: {dgType: string}; label?: string} | undefined]>;
      const outputEntries = Object.entries(node.outputs) as Array<[string, {socket: {dgType: string}; label?: string} | undefined]>;

      if (inputEntries.length > 0) {
        content.appendChild(this.connGroupLabel('Inputs'));
        for (const [key, input] of inputEntries) {
          if (!input) continue;
          const connected = this.flow.isInputConnected(node.id, key);
          content.appendChild(this.buildConnRow('IN', key, input.socket.dgType, connected ? 'connected' : 'disconnected', connected));
        }
      }

      if (ptCount > 0) {
        content.appendChild(this.connSeparator());
        content.appendChild(this.connGroupLabel('Pass-through'));
        for (let i = 0; i < ptCount && i < outputEntries.length; i++) {
          const [key, out] = outputEntries[i];
          if (!out) continue;
          const baseName = key.endsWith('__pt') ? key.slice(0, -'__pt'.length) : key;
          const connected = this.flow.getConnections().some((c) => c.source === node.id && c.sourceOutput === key);
          content.appendChild(this.buildConnRow('PT', baseName, out.socket.dgType, connected ? 'connected' : 'disconnected', connected));
        }
      }

      if (outputEntries.length > ptCount) {
        content.appendChild(this.connSeparator());
        content.appendChild(this.connGroupLabel('Outputs'));
        for (let i = ptCount; i < outputEntries.length; i++) {
          const [key, out] = outputEntries[i];
          if (!out) continue;
          const connected = this.flow.getConnections().some((c) => c.source === node.id && c.sourceOutput === key);
          content.appendChild(this.buildConnRow('OUT', key, out.socket.dgType, connected ? 'connected' : 'disconnected', connected));
        }
      }
      return content;
    }, false);
  }

  private buildConnRow(dir: string, name: string, type: string, status: string, connected: boolean): HTMLElement {
    const dirSpan = ui.element('span');
    dirSpan.textContent = dir;
    dirSpan.className = 'funcflow-conn-dir';
    const detail = ui.element('span');
    detail.textContent = ` ${name} `;
    const typeSpan = ui.element('span');
    typeSpan.textContent = `(${type})`;
    typeSpan.className = 'funcflow-conn-type';
    const statusSpan = ui.element('span');
    statusSpan.textContent = ` — ${status}`;
    statusSpan.className = connected ? 'funcflow-conn-ok' : 'funcflow-conn-off';
    return ui.div([dirSpan, detail, typeSpan, statusSpan], 'funcflow-prop-row funcflow-conn-row');
  }

  private connSeparator(): HTMLElement {return ui.div([], 'funcflow-conn-separator');}
  private connGroupLabel(text: string): HTMLElement {
    const label = ui.div([], 'funcflow-conn-group-label');
    label.textContent = text;
    return label;
  }

  // ---------- editor helpers ----------

  private labelWithTooltip(text: string, explicitTip?: string): HTMLElement {
    const lbl = ui.label(text);
    const tip = explicitTip ?? PROP_TOOLTIPS[text];
    if (tip) ui.tooltip.bind(lbl, tip);
    return lbl;
  }

  private createTextarea(label: string, value: string, onChange: (v: string) => void, inputTooltip?: string): HTMLElement {
    const textarea = document.createElement('textarea');
    textarea.value = value;
    textarea.className = 'funcflow-prop-textarea';
    textarea.rows = 1;
    textarea.addEventListener('input', () => {
      textarea.style.height = 'auto';
      textarea.style.height = textarea.scrollHeight + 'px';
      onChange(textarea.value);
    });
    setTimeout(() => {
      textarea.style.height = 'auto';
      textarea.style.height = textarea.scrollHeight + 'px';
    }, 0);
    if (inputTooltip) ui.tooltip.bind(textarea, inputTooltip);
    return ui.div([this.labelWithTooltip(label, inputTooltip), textarea], 'funcflow-prop-row');
  }

  private createNumberInput(label: string, value: number, onChange: (v: number) => void, decimals: number, step: number, inputTooltip?: string): HTMLElement {
    const input = document.createElement('input');
    input.type = 'number';
    input.value = decimals === 0 ? String(Math.round(value)) : value.toFixed(decimals);
    input.step = String(step);
    input.className = 'funcflow-prop-input';
    input.addEventListener('change', () => {
      const parsed = parseFloat(input.value);
      if (!isNaN(parsed)) onChange(parsed);
    });
    if (inputTooltip) ui.tooltip.bind(input, inputTooltip);
    return ui.div([this.labelWithTooltip(label, inputTooltip), input], 'funcflow-prop-row');
  }

  private createToggle(label: string, value: boolean, onChange: (v: boolean) => void, inputTooltip?: string): HTMLElement {
    const input = document.createElement('input');
    input.type = 'checkbox';
    input.checked = value;
    input.className = 'funcflow-prop-checkbox';
    input.addEventListener('change', () => onChange(input.checked));
    if (inputTooltip) ui.tooltip.bind(input, inputTooltip);
    const lbl = this.labelWithTooltip(label, inputTooltip);
    return ui.div([input, lbl], 'funcflow-prop-row funcflow-prop-toggle-row');
  }

  private createCombo(label: string, value: string, options: string[], onChange: (v: string) => void, inputTooltip?: string): HTMLElement {
    const select = document.createElement('select');
    select.className = 'funcflow-prop-input';
    for (const opt of options) {
      const optEl = document.createElement('option');
      optEl.value = opt;
      optEl.textContent = opt || '(none)';
      if (opt === value) optEl.selected = true;
      select.appendChild(optEl);
    }
    select.addEventListener('change', () => onChange(select.value));
    if (inputTooltip) ui.tooltip.bind(select, inputTooltip);
    return ui.div([this.labelWithTooltip(label, inputTooltip), select], 'funcflow-prop-row');
  }
}

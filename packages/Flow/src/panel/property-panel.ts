import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LGraphNode} from 'litegraph.js';
import {NodeExecState} from '../execution/execution-state';
import {buildValuePanel} from '../execution/value-inspector';

/** Common shape for our custom node properties */
interface FuncFlowNode extends LGraphNode {
  dgNodeType?: string;
  dgFunc?: DG.Func;
  dgFuncName?: string;
  dgRole?: string | null;
  dgOutputType?: string;
}

/** Tooltips for property labels that aren't self-explanatory */
const PROP_TOOLTIPS: Record<string, string> = {
  'Title': 'Display name shown on the node',
  'Param Name': 'Variable name used in the generated script',
  'Description': 'Description shown in the script run dialog',
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

/** Tooltips for utility node properties keyed by node title → property name */
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

/** Build a rich tooltip for a DG.Func input parameter */
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

/** Known combo values for specific properties */
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

/** Property panel: renders into the Datagrok native context panel via grok.shell.o */
export class PropertyPanel {
  root: HTMLElement;
  private contentDiv: HTMLElement;

  constructor() {
    this.contentDiv = ui.div([], 'funcflow-property-content');
    this.root = ui.divV([
      this.contentDiv,
    ], 'funcflow-property-panel');
  }

  showNode(node: LGraphNode): void {
    this.contentDiv.innerHTML = '';
    const n = node as FuncFlowNode;

    // Title
    const titleInput = this.createTextarea('Title', node.title || '', (v) => {
      node.title = v;
    });
    this.contentDiv.appendChild(titleInput);

    // Node type label
    const nodeType = n.dgNodeType || 'function';
    this.contentDiv.appendChild(
      ui.div([ui.label('Type'), ui.divText(nodeType)], 'funcflow-prop-row'),
    );

    // Type-specific editors
    if (n.dgFunc)
      this.showFuncNodeProps(n);
    if (n.dgNodeType === 'input')
      this.showInputNodeProps(node);
    if (n.dgNodeType === 'output')
      this.showOutputNodeProps(node);
    if (n.dgNodeType === 'utility')
      this.showUtilityNodeProps(node);

    this.showConnectionsSummary(node);
  }

  private showFuncNodeProps(n: FuncFlowNode): void {
    const func = n.dgFunc;
    if (!func) return;

    // Function info section
    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Function')], 'funcflow-prop-section-header'));

    if (func.description)
      section.appendChild(ui.div([ui.label('Description'), ui.divText(func.description)], 'funcflow-prop-row'));

    const qualName = n.dgFuncName || func.name;
    section.appendChild(ui.div([ui.label('Full Name'), ui.divText(qualName)], 'funcflow-prop-row'));

    if (n.dgRole)
      section.appendChild(ui.div([ui.label('Role'), ui.divText(n.dgRole)], 'funcflow-prop-row'));

    this.contentDiv.appendChild(section);

    // Editable input parameters section
    if (func.inputs.length > 0) {
      const inputSection = ui.div([], 'funcflow-prop-section');
      inputSection.appendChild(ui.div([ui.label('Input Parameters')], 'funcflow-prop-section-header'));

      for (let i = 0; i < func.inputs.length; i++) {
        const inp = func.inputs[i];
        const propKey = `_input_${inp.name}`;
        const tip = buildFuncInputTooltip(inp);

        // Only show editors for primitive types that have stored values
        if (n.properties[propKey] === undefined) {
          const row = ui.div(
            [ui.divText(`${inp.name}: ${inp.propertyType} (connected only)`)],
            'funcflow-prop-row',
          );
          ui.tooltip.bind(row, tip);
          inputSection.appendChild(row);
          continue;
        }

        const connectedAtSlot = n.isInputConnected(i);
        if (connectedAtSlot) {
          const row = ui.div(
            [ui.divText(`${inp.name}: connected`)], 'funcflow-prop-row',
          );
          ui.tooltip.bind(row, tip);
          inputSection.appendChild(row);
          continue;
        }

        switch (inp.propertyType) {
        case 'string':
          inputSection.appendChild(this.createTextarea(
            inp.name, String(n.properties[propKey] || ''),
            (v) => {
              n.properties[propKey] = v;
            }, tip,
          ));
          break;
        case 'int':
          inputSection.appendChild(this.createNumberInput(
            inp.name, Number(n.properties[propKey] || 0),
            (v) => {
              n.properties[propKey] = Math.round(v);
            }, 0, 1, tip,
          ));
          break;
        case 'double':
        case 'num':
          inputSection.appendChild(this.createNumberInput(
            inp.name, Number(n.properties[propKey] || 0),
            (v) => {
              n.properties[propKey] = v;
            }, 3, 0.1, tip,
          ));
          break;
        case 'bool':
          inputSection.appendChild(this.createToggle(
            inp.name, Boolean(n.properties[propKey]),
            (v) => {
              n.properties[propKey] = v;
            }, tip,
          ));
          break;
        }
      }

      this.contentDiv.appendChild(inputSection);
    }
  }

  private showInputNodeProps(node: LGraphNode): void {
    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Input Configuration')], 'funcflow-prop-section-header'));

    // paramName
    section.appendChild(this.createTextarea('Param Name', node.properties['paramName'] || '', (v) => {
      node.properties['paramName'] = v;
    }));

    // description
    section.appendChild(this.createTextarea('Description', node.properties['description'] || '', (v) => {
      node.properties['description'] = v;
    }));

    const outputType = (node as any).dgOutputType;

    // defaultValue - varies by type
    if (outputType === 'bool') {
      section.appendChild(this.createToggle('Default', Boolean(node.properties['defaultValue']), (v) => {
        node.properties['defaultValue'] = v;
      }));
    } else if (outputType === 'int') {
      section.appendChild(this.createNumberInput('Default', Number(node.properties['defaultValue'] || 0), (v) => {
        node.properties['defaultValue'] = Math.round(v);
      }, 0, 1));
    } else if (outputType === 'double') {
      section.appendChild(this.createNumberInput('Default', Number(node.properties['defaultValue'] || 0), (v) => {
        node.properties['defaultValue'] = v;
      }, 3, 0.1));
    } else if (node.properties['defaultValue'] !== undefined && outputType !== 'dataframe' &&
               outputType !== 'file' && outputType !== 'map' && outputType !== 'blob') {
      section.appendChild(this.createTextarea('Default', String(node.properties['defaultValue'] || ''), (v) => {
        node.properties['defaultValue'] = v;
      }));
    }

    // nullable
    if (node.properties['nullable'] !== undefined) {
      section.appendChild(this.createToggle('Nullable', Boolean(node.properties['nullable']), (v) => {
        node.properties['nullable'] = v;
      }));
    }

    // caption
    if (node.properties['caption'] !== undefined) {
      section.appendChild(this.createTextarea('Caption', String(node.properties['caption'] || ''), (v) => {
        node.properties['caption'] = v;
      }));
    }

    // typeFilter (Column, Column List)
    if (node.properties['typeFilter'] !== undefined) {
      section.appendChild(this.createCombo('Type Filter', String(node.properties['typeFilter'] || ''),
        TYPE_FILTER_VALUES, (v) => {
          node.properties['typeFilter'] = v;
        }));
    }

    // semTypeFilter (Column, Column List)
    if (node.properties['semTypeFilter'] !== undefined) {
      section.appendChild(this.createTextarea('SemType Filter', String(node.properties['semTypeFilter'] || ''), (v) => {
        node.properties['semTypeFilter'] = v;
      }));
    }

    // semType (String Input)
    if (node.properties['semType'] !== undefined) {
      section.appendChild(this.createCombo('SemType', String(node.properties['semType'] || ''),
        SEMTYPE_VALUES, (v) => {
          node.properties['semType'] = v;
        }));
    }

    // choices (String Input)
    if (node.properties['choices'] !== undefined) {
      section.appendChild(this.createTextarea('Choices (comma-sep)', String(node.properties['choices'] || ''), (v) => {
        node.properties['choices'] = v;
      }));
    }

    // min, max (Number, Int)
    if (node.properties['min'] !== undefined) {
      section.appendChild(this.createTextarea('Min', String(node.properties['min'] || ''), (v) => {
        node.properties['min'] = v;
      }));
    }
    if (node.properties['max'] !== undefined) {
      section.appendChild(this.createTextarea('Max', String(node.properties['max'] || ''), (v) => {
        node.properties['max'] = v;
      }));
    }

    // showSlider (Number, Int)
    if (node.properties['showSlider'] !== undefined) {
      section.appendChild(this.createToggle('Show Slider', Boolean(node.properties['showSlider']), (v) => {
        node.properties['showSlider'] = v;
      }));
    }

    this.contentDiv.appendChild(section);
  }

  private showOutputNodeProps(node: LGraphNode): void {
    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Output Configuration')], 'funcflow-prop-section-header'));

    // paramName
    section.appendChild(this.createTextarea('Param Name', node.properties['paramName'] || '', (v) => {
      node.properties['paramName'] = v;
    }));

    // outputType (Value Output only)
    if (node.properties['outputType'] !== undefined) {
      section.appendChild(this.createCombo('Output Type', String(node.properties['outputType'] || 'dynamic'),
        OUTPUT_TYPE_VALUES, (v) => {
          node.properties['outputType'] = v;
        }));
    }

    this.contentDiv.appendChild(section);
  }

  private showUtilityNodeProps(node: LGraphNode): void {
    const props = Object.entries(node.properties).filter(([key]) => !key.startsWith('_'));
    if (props.length === 0) return;

    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Configuration')], 'funcflow-prop-section-header'));
    const nodeTips = UTILITY_PROP_TOOLTIPS[node.title] || {};

    for (const [key, val] of props) {
      const tip = nodeTips[key];
      if (typeof val === 'boolean') {
        section.appendChild(this.createToggle(key, val, (v) => {
          node.properties[key] = v;
        }, tip));
      } else if (typeof val === 'number') {
        const isInt = Number.isInteger(val);
        section.appendChild(this.createNumberInput(key, val, (v) => {
          node.properties[key] = isInt ? Math.round(v) : v;
          // Sync widget if it exists (e.g., ConstString)
          if (node.widgets) {
            const w = node.widgets.find((w) => w.options?.property === key);
            if (w) w.value = node.properties[key];
          }
        }, isInt ? 0 : 3, isInt ? 1 : 0.1, tip));
      } else {
        section.appendChild(this.createTextarea(key, String(val), (v) => {
          node.properties[key] = v;
          // Sync widget if it exists (e.g., ConstString)
          if (node.widgets) {
            const w = node.widgets.find((w) => w.options?.property === key);
            if (w) w.value = v;
          }
        }, tip));
      }
    }

    this.contentDiv.appendChild(section);
  }

  private showConnectionsSummary(node: LGraphNode): void {
    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Connections')], 'funcflow-prop-section-header'));

    const ptCount = (node as any).properties?.['_passthroughCount'] ?? 0;

    // --- Inputs ---
    if (node.inputs && node.inputs.length > 0) {
      section.appendChild(this.connGroupLabel('Inputs'));
      for (let i = 0; i < node.inputs.length; i++) {
        const inp = node.inputs[i];
        const connected = node.isInputConnected(i);
        section.appendChild(this.buildConnRow(
          'IN', inp.name, inp.type as string, connected ? 'connected' : 'disconnected', connected));
      }
    }

    // --- Pass-through outputs ---
    if (ptCount > 0 && node.outputs) {
      section.appendChild(this.connSeparator());
      section.appendChild(this.connGroupLabel('Pass-through'));
      for (let i = 0; i < ptCount && i < node.outputs.length; i++) {
        const out = node.outputs[i];
        const connected = node.isOutputConnected(i);
        const linkCount = out.links ? out.links.length : 0;
        const status = connected ? `${linkCount} link(s)` : 'disconnected';
        section.appendChild(this.buildConnRow('PT', out.name, out.type as string, status, connected));
      }
    }

    // --- Real outputs ---
    if (node.outputs && node.outputs.length > ptCount) {
      section.appendChild(this.connSeparator());
      section.appendChild(this.connGroupLabel('Outputs'));
      for (let i = ptCount; i < node.outputs.length; i++) {
        const out = node.outputs[i];
        const connected = node.isOutputConnected(i);
        const linkCount = out.links ? out.links.length : 0;
        const status = connected ? `${linkCount} link(s)` : 'disconnected';
        section.appendChild(this.buildConnRow('OUT', out.name, out.type as string, status, connected));
      }
    }

    this.contentDiv.appendChild(section);
  }

  private buildConnRow(
    dir: string, name: string, type: string, status: string, connected: boolean,
  ): HTMLElement {
    const dirSpan = ui.element('span');
    dirSpan.textContent = dir;
    dirSpan.className = 'funcflow-conn-dir';
    const detail = ui.element('span');
    detail.textContent = ` ${name} `;
    const typeSpan = ui.element('span');
    typeSpan.textContent = `(${type})`;
    typeSpan.className = 'funcflow-conn-type';
    const statusSpan = ui.element('span');
    statusSpan.textContent = ` \u2014 ${status}`;
    statusSpan.className = connected ? 'funcflow-conn-ok' : 'funcflow-conn-off';
    return ui.div([dirSpan, detail, typeSpan, statusSpan], 'funcflow-prop-row funcflow-conn-row');
  }

  private connSeparator(): HTMLElement {
    return ui.div([], 'funcflow-conn-separator');
  }

  private connGroupLabel(text: string): HTMLElement {
    const label = ui.div([], 'funcflow-conn-group-label');
    label.textContent = text;
    return label;
  }

  // --- Editor helpers ---

  /** Create a label element, attaching a tooltip if one is defined in PROP_TOOLTIPS or explicitly given */
  private labelWithTooltip(text: string, explicitTip?: string): HTMLElement {
    const lbl = ui.label(text);
    const tip = explicitTip || PROP_TOOLTIPS[text];
    if (tip)
      ui.tooltip.bind(lbl, tip);
    return lbl;
  }

  private createTextarea(
    label: string, value: string, onChange: (v: string) => void, inputTooltip?: string,
  ): HTMLElement {
    const textarea = document.createElement('textarea');
    textarea.value = value;
    textarea.className = 'funcflow-prop-textarea';
    textarea.rows = 1;
    textarea.addEventListener('input', () => {
      // Auto-resize
      textarea.style.height = 'auto';
      textarea.style.height = textarea.scrollHeight + 'px';
      onChange(textarea.value);
    });
    // Initial auto-size
    setTimeout(() => {
      textarea.style.height = 'auto';
      textarea.style.height = textarea.scrollHeight + 'px';
    }, 0);
    if (inputTooltip) ui.tooltip.bind(textarea, inputTooltip);
    return ui.div([this.labelWithTooltip(label, inputTooltip), textarea], 'funcflow-prop-row');
  }

  private createNumberInput(
    label: string, value: number, onChange: (v: number) => void,
    decimals: number, step: number, inputTooltip?: string,
  ): HTMLElement {
    const input = document.createElement('input');
    input.type = 'number';
    input.value = decimals === 0 ? String(Math.round(value)) : value.toFixed(decimals);
    input.step = String(step);
    input.className = 'funcflow-prop-input';
    input.addEventListener('change', () => {
      const parsed = parseFloat(input.value);
      if (!isNaN(parsed))
        onChange(parsed);
    });
    if (inputTooltip) ui.tooltip.bind(input, inputTooltip);
    return ui.div([this.labelWithTooltip(label, inputTooltip), input], 'funcflow-prop-row');
  }

  private createToggle(
    label: string, value: boolean, onChange: (v: boolean) => void, inputTooltip?: string,
  ): HTMLElement {
    const input = document.createElement('input');
    input.type = 'checkbox';
    input.checked = value;
    input.className = 'funcflow-prop-checkbox';
    input.addEventListener('change', () => {
      onChange(input.checked);
    });
    if (inputTooltip) ui.tooltip.bind(input, inputTooltip);
    const lbl = this.labelWithTooltip(label, inputTooltip);
    const row = ui.div([input, lbl], 'funcflow-prop-row funcflow-prop-toggle-row');
    return row;
  }

  private createCombo(
    label: string, value: string, options: string[], onChange: (v: string) => void,
    inputTooltip?: string,
  ): HTMLElement {
    const select = document.createElement('select');
    select.className = 'funcflow-prop-input';
    for (const opt of options) {
      const optEl = document.createElement('option');
      optEl.value = opt;
      optEl.textContent = opt || '(none)';
      if (opt === value) optEl.selected = true;
      select.appendChild(optEl);
    }
    select.addEventListener('change', () => {
      onChange(select.value);
    });
    if (inputTooltip) ui.tooltip.bind(select, inputTooltip);
    return ui.div([this.labelWithTooltip(label, inputTooltip), select], 'funcflow-prop-row');
  }

  /** Shows node properties with optional execution state section appended */
  showNodeWithExecution(node: LGraphNode, execState?: NodeExecState): void {
    this.showNode(node);
    if (execState) {
      const separator = ui.div([], 'funcflow-prop-section-header');
      separator.textContent = 'Execution';
      this.contentDiv.appendChild(separator);
      this.contentDiv.appendChild(buildValuePanel(execState));
    }
  }

  clear(): void {
    this.contentDiv.innerHTML = '';
    this.contentDiv.appendChild(ui.divText('Select a node to view its properties'));
  }
}

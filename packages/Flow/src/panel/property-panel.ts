import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LGraphNode} from 'litegraph.js';

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
  'Param Name': 'Variable name used in the generated script',
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

/** Right sidebar: shows properties of the selected node */
export class PropertyPanel {
  root: HTMLElement;
  private contentDiv: HTMLElement;

  constructor() {
    this.contentDiv = ui.div([], 'funcflow-property-content');
    this.root = ui.divV([
      ui.div([ui.label('Properties')], 'funcflow-panel-header'),
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

        // Only show editors for primitive types that have stored values
        if (n.properties[propKey] === undefined) {
          inputSection.appendChild(
            ui.div([ui.divText(`${inp.name}: ${inp.propertyType} (connected only)`)], 'funcflow-prop-row'),
          );
          continue;
        }

        const connectedAtSlot = n.isInputConnected(i);
        if (connectedAtSlot) {
          inputSection.appendChild(
            ui.div([ui.divText(`${inp.name}: connected`)], 'funcflow-prop-row'),
          );
          continue;
        }

        switch (inp.propertyType) {
        case 'string':
          inputSection.appendChild(this.createTextarea(inp.name, String(n.properties[propKey] || ''), (v) => {
            n.properties[propKey] = v;
          }));
          break;
        case 'int':
          inputSection.appendChild(this.createNumberInput(inp.name, Number(n.properties[propKey] || 0), (v) => {
            n.properties[propKey] = Math.round(v);
          }, 0, 1));
          break;
        case 'double':
        case 'num':
          inputSection.appendChild(this.createNumberInput(inp.name, Number(n.properties[propKey] || 0), (v) => {
            n.properties[propKey] = v;
          }, 3, 0.1));
          break;
        case 'bool':
          inputSection.appendChild(this.createToggle(inp.name, Boolean(n.properties[propKey]), (v) => {
            n.properties[propKey] = v;
          }));
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

    for (const [key, val] of props) {
      if (typeof val === 'boolean') {
        section.appendChild(this.createToggle(key, val, (v) => {
          node.properties[key] = v;
        }));
      } else if (typeof val === 'number') {
        const isInt = Number.isInteger(val);
        section.appendChild(this.createNumberInput(key, val, (v) => {
          node.properties[key] = isInt ? Math.round(v) : v;
          // Sync widget if it exists (e.g., ConstString)
          if (node.widgets) {
            const w = node.widgets.find((w) => w.options?.property === key);
            if (w) w.value = node.properties[key];
          }
        }, isInt ? 0 : 3, isInt ? 1 : 0.1));
      } else {
        section.appendChild(this.createTextarea(key, String(val), (v) => {
          node.properties[key] = v;
          // Sync widget if it exists (e.g., ConstString)
          if (node.widgets) {
            const w = node.widgets.find((w) => w.options?.property === key);
            if (w) w.value = v;
          }
        }));
      }
    }

    this.contentDiv.appendChild(section);
  }

  private showConnectionsSummary(node: LGraphNode): void {
    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Connections')], 'funcflow-prop-section-header'));

    if (node.inputs) {
      for (let i = 0; i < node.inputs.length; i++) {
        const inp = node.inputs[i];
        const connected = node.isInputConnected(i);
        const status = connected ? 'connected' : 'disconnected';
        section.appendChild(
          ui.div([ui.divText(`IN: ${inp.name} (${inp.type}) - ${status}`)], 'funcflow-prop-row'),
        );
      }
    }

    if (node.outputs) {
      for (let i = 0; i < node.outputs.length; i++) {
        const out = node.outputs[i];
        const connected = node.isOutputConnected(i);
        const linkCount = out.links ? out.links.length : 0;
        const status = connected ? `${linkCount} connection(s)` : 'disconnected';
        section.appendChild(
          ui.div([ui.divText(`OUT: ${out.name} (${out.type}) - ${status}`)], 'funcflow-prop-row'),
        );
      }
    }

    this.contentDiv.appendChild(section);
  }

  // --- Editor helpers ---

  /** Create a label element, attaching a tooltip if one is defined in PROP_TOOLTIPS */
  private labelWithTooltip(text: string): HTMLElement {
    const lbl = ui.label(text);
    const tip = PROP_TOOLTIPS[text];
    if (tip)
      ui.tooltip.bind(lbl, tip);
    return lbl;
  }

  private createTextarea(label: string, value: string, onChange: (v: string) => void): HTMLElement {
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
    return ui.div([this.labelWithTooltip(label), textarea], 'funcflow-prop-row');
  }

  private createNumberInput(
    label: string, value: number, onChange: (v: number) => void,
    decimals: number, step: number,
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
    return ui.div([this.labelWithTooltip(label), input], 'funcflow-prop-row');
  }

  private createToggle(label: string, value: boolean, onChange: (v: boolean) => void): HTMLElement {
    const input = document.createElement('input');
    input.type = 'checkbox';
    input.checked = value;
    input.className = 'funcflow-prop-checkbox';
    input.addEventListener('change', () => {
      onChange(input.checked);
    });
    const lbl = this.labelWithTooltip(label);
    const row = ui.div([input, lbl], 'funcflow-prop-row funcflow-prop-toggle-row');
    return row;
  }

  private createCombo(
    label: string, value: string, options: string[], onChange: (v: string) => void,
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
    return ui.div([this.labelWithTooltip(label), select], 'funcflow-prop-row');
  }

  clear(): void {
    this.contentDiv.innerHTML = '';
    this.contentDiv.appendChild(ui.divText('Select a node to view its properties'));
  }
}

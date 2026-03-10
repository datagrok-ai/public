import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LGraphNode} from 'litegraph.js';

/** Common shape for our custom node properties */
interface FuncFlowNode extends LGraphNode {
  dgNodeType?: string;
  dgFunc?: DG.Func;
  dgFuncName?: string;
  dgRole?: string | null;
}

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

    const titleInput = document.createElement('input');
    titleInput.type = 'text';
    titleInput.value = node.title || '';
    titleInput.className = 'funcflow-prop-input';
    titleInput.addEventListener('change', () => {
      node.title = titleInput.value;
    });
    this.contentDiv.appendChild(ui.div([ui.label('Title'), titleInput], 'funcflow-prop-row'));

    const nodeType = n.dgNodeType || 'function';
    this.contentDiv.appendChild(
      ui.div([ui.label('Type'), ui.divText(nodeType)], 'funcflow-prop-row'),
    );

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

    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Function')], 'funcflow-prop-section-header'));

    if (func.description)
      section.appendChild(ui.div([ui.label('Description'), ui.divText(func.description)], 'funcflow-prop-row'));

    const qualName = n.dgFuncName || func.name;
    section.appendChild(ui.div([ui.label('Full Name'), ui.divText(qualName)], 'funcflow-prop-row'));

    if (n.dgRole)
      section.appendChild(ui.div([ui.label('Role'), ui.divText(n.dgRole)], 'funcflow-prop-row'));

    if (func.inputs.length > 0) {
      const inputsDiv = ui.div([], 'funcflow-prop-section');
      inputsDiv.appendChild(ui.div([ui.label('Inputs')], 'funcflow-prop-section-header'));
      for (const inp of func.inputs)
        inputsDiv.appendChild(ui.div([ui.divText(`${inp.name}: ${inp.propertyType}`)], 'funcflow-prop-row'));
      section.appendChild(inputsDiv);
    }

    if (func.outputs.length > 0) {
      const outputsDiv = ui.div([], 'funcflow-prop-section');
      outputsDiv.appendChild(ui.div([ui.label('Outputs')], 'funcflow-prop-section-header'));
      for (const out of func.outputs)
        outputsDiv.appendChild(ui.div([ui.divText(`${out.name}: ${out.propertyType}`)], 'funcflow-prop-row'));
      section.appendChild(outputsDiv);
    }

    this.contentDiv.appendChild(section);
  }

  private showInputNodeProps(node: LGraphNode): void {
    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Input Configuration')], 'funcflow-prop-section-header'));

    const paramInput = document.createElement('input');
    paramInput.type = 'text';
    paramInput.value = node.properties['paramName'] || '';
    paramInput.className = 'funcflow-prop-input';
    paramInput.addEventListener('change', () => {
      node.properties['paramName'] = paramInput.value;
      if (node.widgets) {
        const w = node.widgets.find((w) => w.name === 'Param Name');
        if (w) w.value = paramInput.value;
      }
    });
    section.appendChild(ui.div([ui.label('Param Name'), paramInput], 'funcflow-prop-row'));

    const descInput = document.createElement('input');
    descInput.type = 'text';
    descInput.value = node.properties['description'] || '';
    descInput.className = 'funcflow-prop-input';
    descInput.addEventListener('change', () => {
      node.properties['description'] = descInput.value;
    });
    section.appendChild(ui.div([ui.label('Description'), descInput], 'funcflow-prop-row'));

    this.contentDiv.appendChild(section);
  }

  private showOutputNodeProps(node: LGraphNode): void {
    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Output Configuration')], 'funcflow-prop-section-header'));

    const paramInput = document.createElement('input');
    paramInput.type = 'text';
    paramInput.value = node.properties['paramName'] || '';
    paramInput.className = 'funcflow-prop-input';
    paramInput.addEventListener('change', () => {
      node.properties['paramName'] = paramInput.value;
      if (node.widgets) {
        const w = node.widgets.find((w) => w.name === 'Param Name');
        if (w) w.value = paramInput.value;
      }
    });
    section.appendChild(ui.div([ui.label('Param Name'), paramInput], 'funcflow-prop-row'));

    this.contentDiv.appendChild(section);
  }

  private showUtilityNodeProps(node: LGraphNode): void {
    const section = ui.div([], 'funcflow-prop-section');
    section.appendChild(ui.div([ui.label('Configuration')], 'funcflow-prop-section-header'));

    for (const [key, val] of Object.entries(node.properties)) {
      if (key.startsWith('_')) continue;
      const input = document.createElement('input');
      input.type = 'text';
      input.value = String(val);
      input.className = 'funcflow-prop-input';
      input.addEventListener('change', () => {
        node.properties[key] = input.value;
      });
      section.appendChild(ui.div([ui.label(key), input], 'funcflow-prop-row'));
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

  clear(): void {
    this.contentDiv.innerHTML = '';
    this.contentDiv.appendChild(ui.divText('Select a node to view its properties'));
  }
}

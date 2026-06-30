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
import {constLabel} from '../rete/nodes/utility-nodes';
import {NodeExecState} from '../execution/execution-state';
import {buildExecutionMeta} from '../execution/value-inspector';
import {setTid} from '../utils/test-ids';
import {ColumnPickRequest} from './column-picker';

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
  'Select Table': {tableName: 'Name of an open table (resolved via grok.shell.tableByName)'},
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

/** A property's `choices` as a non-empty string list, or `[]` when it has none
 *  (or the Dart proxy access fails). Empty entries are dropped — the nullable
 *  empty option is added separately by `stringChoiceOptions`. */
export function propertyChoices(param: DG.Property): string[] {
  try {
    const choices = (param as unknown as {choices?: unknown}).choices;
    if (!Array.isArray(choices)) return [];
    return choices.map((c) => String(c)).filter((c) => c.length > 0);
  } catch {
    return [];
  }
}

/** Build the option list for a string input that declares `choices`, or `null`
 *  when there are none (caller renders a free-text field instead). A nullable
 *  property gets a leading empty option; the current value is preserved as an
 *  option even when it isn't among the declared choices, so imported values are
 *  never silently dropped. */
export function stringChoiceOptions(choices: string[], nullable: boolean, current: string): string[] | null {
  if (choices.length === 0) return null;
  let options = [...choices];
  if (nullable) options = ['', ...options];
  if (current !== '' && !options.includes(current)) options = [current, ...options];
  return options;
}

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

  /** Set by the view: opens a column / columns picker dialog for a func-node
   *  column input, seeded by the upstream table (running the flow up to that
   *  point if needed). When unset, the picker icon is not rendered. */
  onPickColumns?: (req: ColumnPickRequest) => void;

  constructor(flow: FlowEditor) {
    this.flow = flow;
    this.contentDiv = setTid(ui.div([], 'funcflow-property-content'), 'property-content');
    this.root = setTid(ui.divV([this.contentDiv], 'funcflow-property-panel'), 'property-panel');
  }

  showNode(node: FlowNode, execState?: NodeExecState): void {
    this.contentDiv.innerHTML = '';

    // Coerce: labels can be derived from non-string values (constant nodes
    // title themselves after their value) and DG string inputs throw on
    // anything but a string.
    const titleInput = this.createTextarea('Title', String(node.label ?? ''), (v) => {
      node.label = String(v ?? '');
      void this.flow.updateNode(node.id);
    });
    const typeBadge = setTid(ui.div([], 'funcflow-type-badge'), 'property-type-badge');
    typeBadge.textContent = node.dgNodeType || 'function';
    const titleRow = setTid(ui.div([titleInput, typeBadge], 'funcflow-title-row'), 'property-title-row');
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
    if (node.properties['viewerType']) this.addViewerNodePane(acc, node);
    else if (node.dgNodeType === 'utility') this.addUtilityNodePane(acc, node);

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
      const dataframeParams = func.inputs.filter((p) => String(p.propertyType) === 'dataframe').map((p) => p.name);
      acc.addPane('Input Parameters', () => {
        const content = ui.div([], 'funcflow-accordion-content');
        for (const inp of func.inputs) {
          const tip = buildFuncInputTooltip(inp);
          const isEditable = inp.name in node.inputValues;
          if (!isEditable) {
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
          case 'string': {
            // When the property declares `choices`, render a combo (with a
            // leading empty option when nullable) instead of a free-text field.
            const current = String(node.inputValues[inp.name] ?? '');
            const options = stringChoiceOptions(propertyChoices(inp), Boolean(inp.nullable), current);
            if (options) {
              content.appendChild(this.createCombo(inp.name, current, options,
                (v) => {node.inputValues[inp.name] = v;}, tip));
            } else {
              content.appendChild(this.createTextarea(inp.name, current,
                (v) => {node.inputValues[inp.name] = v;}, tip));
            }
            break;
          }
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
          case 'column':
            content.appendChild(this.createColumnRow(node, inp.name, false, dataframeParams, tip));
            break;
          case 'column_list':
            content.appendChild(this.createColumnRow(node, inp.name, true, dataframeParams, tip));
            break;
          case 'string_list':
            // Comma-separated; the compiler trims and turns it into an array.
            // (`list<string>` arrives here as `string_list` — DG normalizes it.)
            content.appendChild(this.createTextarea(inp.name, String(node.inputValues[inp.name] ?? ''),
              (v) => {node.inputValues[inp.name] = v;}, `${tip} | Comma-separated list of strings`));
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

  /** Viewer node: the curated, exposed look options (column choices, title).
   *  Everything else is edited live via the preview's "Edit settings" button. */
  private addViewerNodePane(acc: DG.Accordion, node: FlowNode): void {
    if (!node.properties['viewerLook'] || typeof node.properties['viewerLook'] !== 'object')
      node.properties['viewerLook'] = {};
    const look = node.properties['viewerLook'] as Record<string, unknown>;
    const specs = (node.properties['viewerOptionSpecs'] as Array<{key: string; label: string; kind: string}>) ?? [];

    acc.addPane('Viewer', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      content.appendChild(ui.div([ui.label('Type'),
        ui.divText(String(node.properties['viewerType']))], 'funcflow-prop-row'));
      for (const o of specs) {
        // A column option can also be wired in (a column socket). When connected,
        // the wired column wins — show it as connected rather than an editor.
        if (o.kind === 'column' && this.flow.isInputConnected(node.id, o.key)) {
          const row = ui.div([ui.divText(`${o.label}: connected`)], 'funcflow-prop-row');
          ui.tooltip.bind(row, 'A column is wired into this option; its name is used.');
          content.appendChild(row);
          continue;
        }
        const tip = o.kind === 'column' ? 'Column name (or wire a column into the socket)' : undefined;
        content.appendChild(this.createTextarea(o.label, String(look[o.key] ?? ''), (v) => {
          const s = String(v ?? '').trim();
          if (s) look[o.key] = s;
          else delete look[o.key];
          void this.flow.updateNode(node.id);
        }, tip));
      }
      const note = ui.divText('Run the flow, then click the viewer in the preview panel and use ' +
        '“Edit settings” to change every other setting — your changes are saved on the node.');
      note.style.cssText = 'font-size:11px;color:#888;margin-top:6px;line-height:1.4;';
      content.appendChild(note);
      return content;
    }, true);
  }

  private addUtilityNodePane(acc: DG.Accordion, node: FlowNode): void {
    const props = Object.entries(node.properties).filter(([k]) => !k.startsWith('_'));
    if (props.length === 0) return;

    // Stable kind from the registered type — labels are user-editable.
    const kind = node.dgTypeName?.split('/').pop() ?? node.label;
    const isConstant = node.dgTypeName?.startsWith('Constants/') === true;
    const retitle = (value: unknown): void => {
      node.label = constLabel(kind, value);
      void this.flow.updateNode(node.id);
    };

    acc.addPane('Configuration', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      const nodeTips = UTILITY_PROP_TOOLTIPS[kind] ?? {};
      for (const [key, val] of props) {
        const tip = nodeTips[key];
        const isConstValue = isConstant && key === 'value';
        if (typeof val === 'boolean') {
          content.appendChild(this.createToggle(key, val, (v) => {
            node.properties[key] = v;
            if (isConstValue) retitle(v);
          }, tip));
        } else if (typeof val === 'number') {
          const isInt = Number.isInteger(val);
          content.appendChild(this.createNumberInput(key, val,
            (v) => {
              node.properties[key] = isInt ? Math.round(v) : v;
              if (isConstValue) retitle(node.properties[key]);
            },
            isInt ? 0 : 3, isInt ? 1 : 0.1, tip));
        } else {
          content.appendChild(this.createTextarea(key, String(val ?? ''), (v) => {
            node.properties[key] = v;
            if (isConstValue) retitle(v);
          }, tip));
        }
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

  private buildTextareaEl(value: string, onChange: (v: string) => void, inputTooltip?: string): HTMLTextAreaElement {
    const textarea = document.createElement('textarea');
    textarea.value = value;
    textarea.className = 'funcflow-prop-textarea';
    textarea.rows = 1;
    const autosize = (): void => {
      textarea.style.height = 'auto';
      textarea.style.height = textarea.scrollHeight + 'px';
    };
    textarea.addEventListener('input', () => {
      autosize();
      onChange(textarea.value);
    });
    setTimeout(autosize, 0);
    if (inputTooltip) ui.tooltip.bind(textarea, inputTooltip);
    return textarea;
  }

  /** Stamp an input row with its test-id + a human-findable `data-param` (the
   *  input/param name), so a specific field is addressable in the context panel. */
  private propRow(el: HTMLElement, label: string): HTMLElement {
    setTid(el, 'prop-input', label);
    el.dataset.param = label;
    return el;
  }

  private createTextarea(label: string, value: string, onChange: (v: string) => void, inputTooltip?: string): HTMLElement {
    return this.propRow(ui.div([this.labelWithTooltip(label, inputTooltip),
      this.buildTextareaEl(value, onChange, inputTooltip)], 'funcflow-prop-row'), label);
  }

  /** A column / column-list input laid out side by side with its table picker:
   *  the column-name field (≈70%) and, when the func has 2+ dataframe inputs, a
   *  table chooser (≈30%) writing the node's `columnTables` association. With a
   *  single dataframe input the field spans full width (the table is implicit). */
  private createColumnRow(
    node: FlowNode, paramName: string, isList: boolean, dataframeParams: string[], tip: string,
  ): HTMLElement {
    const colTip = isList ?
      `${tip} | Comma-separated column names` :
      `${tip} | Column name (compiled to table.col(...))`;
    const textarea = this.buildTextareaEl(String(node.inputValues[paramName] ?? ''),
      (v) => {node.inputValues[paramName] = v;}, colTip);

    // Which dataframe input this column resolves against: fixed for a single
    // dataframe input, or chosen per-row for multi-table funcs (JoinTables:
    // keys1→table1, keys2→table2). Read live so the picker uses the current pick.
    let getTableParam = (): string => dataframeParams[0];

    const cells: HTMLElement[] = [ui.div([textarea], 'funcflow-col-input-cell')];
    if (dataframeParams.length >= 2) {
      if (!node.properties['columnTables']) node.properties['columnTables'] = {};
      const associations = node.properties['columnTables'] as Record<string, string>;
      const current = associations[paramName] ?? dataframeParams[0];
      const select = this.buildSelectEl(current, dataframeParams,
        (v) => {associations[paramName] = v;}, 'Which table input this column refers to');
      getTableParam = (): string => select.value || dataframeParams[0];
      cells.push(ui.div([select], 'funcflow-col-table-cell'));
    }

    // Column chooser — opens a dialog seeded by the upstream table (running the
    // flow up to that point if needed) so users pick from a real column list.
    if (this.onPickColumns) {
      const pickBtn = ui.iconFA('list', () => {
        this.onPickColumns!({
          nodeId: node.id, paramName, isList,
          tableParam: getTableParam(),
          current: String(node.inputValues[paramName] ?? ''),
          apply: (value: string) => {
            textarea.value = value;
            node.inputValues[paramName] = value;
            textarea.style.height = 'auto';
            textarea.style.height = textarea.scrollHeight + 'px';
          },
        });
      }, isList ? 'Choose columns from the connected table' : 'Choose a column from the connected table');
      setTid(pickBtn, 'prop-pick-columns', paramName);
      pickBtn.classList.add('funcflow-col-pick');
      cells.push(ui.div([pickBtn], 'funcflow-col-pick-cell'));
    }

    return this.propRow(ui.div([this.labelWithTooltip(paramName, colTip), ui.div(cells, 'funcflow-col-grid')],
      'funcflow-prop-row'), paramName);
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
    return this.propRow(ui.div([this.labelWithTooltip(label, inputTooltip), input], 'funcflow-prop-row'), label);
  }

  private createToggle(label: string, value: boolean, onChange: (v: boolean) => void, inputTooltip?: string): HTMLElement {
    const input = document.createElement('input');
    input.type = 'checkbox';
    input.checked = value;
    input.className = 'funcflow-prop-checkbox';
    input.addEventListener('change', () => onChange(input.checked));
    if (inputTooltip) ui.tooltip.bind(input, inputTooltip);
    const lbl = this.labelWithTooltip(label, inputTooltip);
    return this.propRow(ui.div([input, lbl], 'funcflow-prop-row funcflow-prop-toggle-row'), label);
  }

  private buildSelectEl(value: string, options: string[], onChange: (v: string) => void, inputTooltip?: string): HTMLSelectElement {
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
    return select;
  }

  private createCombo(label: string, value: string, options: string[], onChange: (v: string) => void, inputTooltip?: string): HTMLElement {
    return this.propRow(ui.div([this.labelWithTooltip(label, inputTooltip),
      this.buildSelectEl(value, options, onChange, inputTooltip)], 'funcflow-prop-row'), label);
  }
}

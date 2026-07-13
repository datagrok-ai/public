/* eslint-disable max-len */
/** Property panel — renders into Datagrok's native context panel.
 *
 * Reads from a `FlowNode` selected on the canvas and edits its `properties`
 * and `inputValues` maps. Re-renders the node on the canvas whenever a
 * property changes so visible state (title, etc.) stays in sync. */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FlowEditor} from '../rete/flow-editor';
import {FlowNode, missingRequiredInputs, missingRequiredProps, isExecKey, EXEC_IN_KEY, EXEC_OUT_KEY} from '../rete/scheme';
import {constLabel} from '../rete/nodes/utility-nodes';
import {NodeExecState} from '../execution/execution-state';
import {buildExecutionMeta} from '../execution/value-inspector';
import {setTid} from '../utils/test-ids';
import {getParamDescription, getParamDisplayName, getFuncDisplayName, getTags} from '../utils/dart-proxy-utils';
import {propertyNameToFriendly} from '../utils/naming';
import {shouldUseFunctionEditor} from '../utils/func-editor-utils';
import {ColumnPickRequest} from './column-picker';

const PROP_TOOLTIPS: Record<string, string> = {
  'Title': 'Display name shown on the node',
  'Param Name': 'Variable name used in the generated script',
  'Description': 'What this node does — starts as the function\'s own description; edit to override. ' +
    'Rendered under the node title; for input/output nodes, also embedded in the //input:/output: line',
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

/** Primitive func-input types edited with a native Datagrok input (built from
 *  the property via `ui.input.forProperty`, which honours the declared type,
 *  numeric range, choices, and nullability). Structural (column/column_list)
 *  and `string_list` inputs are handled separately. */
const PRIMITIVE_INPUT_TYPES: ReadonlySet<string> = new Set([
  'string', 'int', 'double', 'num', 'qnum', 'datetime', 'bool',
]);

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
  const desc = getParamDescription(param) || param.description;
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
  /** The node the panel currently renders — the target of change reports. */
  private currentNode: FlowNode | null = null;

  /** Set by the view: opens a column / columns picker dialog for a func-node
   *  column input, seeded by the upstream table (running the flow up to that
   *  point if needed). When unset, the picker icon is not rendered. */
  onPickColumns?: (req: ColumnPickRequest) => void;

  /** Set by the view: opens the function's own editor dialog (functions with
   *  `editor:` meta or on the explicit allowlist) seeded with the node's real
   *  upstream tables, then writes the edited values back into `inputValues`.
   *  When unset, the editor icon is not rendered. */
  onEditFuncParams?: (node: FlowNode) => void;

  constructor(flow: FlowEditor) {
    this.flow = flow;
    this.contentDiv = setTid(ui.div([], 'funcflow-property-content'), 'property-content');
    this.root = setTid(ui.divV([this.contentDiv], 'funcflow-property-panel'), 'property-panel');
  }

  showNode(node: FlowNode, execState?: NodeExecState): void {
    this.contentDiv.innerHTML = '';
    this.currentNode = node;

    // Coerce: labels can be derived from non-string values (constant nodes
    // title themselves after their value) and DG string inputs throw on
    // anything but a string. Title is cosmetic — it never changes what the
    // flow computes, so it must not invalidate run results.
    const titleInput = this.createTextarea('Title', String(node.label ?? ''), (v) => {
      node.label = String(v ?? '');
      void this.flow.updateNode(node.id);
    }, undefined, true);
    const typeBadge = setTid(ui.div([], 'funcflow-type-badge'), 'property-type-badge');
    typeBadge.textContent = node.dgNodeType || 'function';
    const titleRow = setTid(ui.div([titleInput, typeBadge], 'funcflow-title-row'), 'property-title-row');

    // One header block (shared padding) so Title, chips, and Description line up.
    const header = setTid(ui.div([titleRow], 'funcflow-panel-header'), 'property-header');

    if (node.dgFunc) header.appendChild(this.buildFuncChips(node));

    // Per-node description: rendered under the title in the canvas, and
    // embedded as the [description] suffix in //input:/output: lines.
    // Cosmetic like the title — annotations don't change computed values.
    // For func nodes it starts as the function's own description; an edit
    // stores the override on the node (the function text stays the fallback).
    let funcDesc = '';
    try {
      funcDesc = node.dgFunc?.description ?? '';
    } catch {/* Dart proxy access can throw */}
    const descSeed = node.description?.trim() ? node.description : funcDesc;
    header.appendChild(this.createTextarea('Description', descSeed, (v) => {
      node.description = v;
      void this.flow.updateNode(node.id);
    }, undefined, true));
    this.contentDiv.appendChild(header);

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
    this.currentNode = null;
    this.contentDiv.innerHTML = '';
    this.contentDiv.appendChild(ui.divText('Select a node to view its properties'));
  }

  /** Report a (non-cosmetic) parameter edit on the shown node — routed to the
   *  editor so run results downstream of the node get invalidated (and autorun,
   *  when enabled, reruns the affected slice). Every semantic editor helper
   *  below funnels its change through a {@link changeReporter}, never here
   *  directly. */
  private paramsChanged(): void {
    if (this.currentNode) this.flow.notifyNodeParamsChanged(this.currentNode.id);
  }

  /** A change reporter for ONE editor: report a parameter edit only when the
   *  value ACTUALLY differs from the last seen one.
   *
   *  Why: creating/initializing a Datagrok input (`ui.input.forProperty`,
   *  `initInputValue` setting `stringValue`, …) can fire `onValueChanged`
   *  immediately — but not for every input type, so "skip the first event" is
   *  wrong too. Without this guard, merely clicking a node (which rebuilds the
   *  panel) would count as an edit — invalidating results and, with autorun
   *  on, rerunning the flow on every selection. */
  private changeReporter(initial: unknown): (v: unknown) => void {
    let last = initial;
    return (v: unknown): void => {
      if (PropertyPanel.sameValue(last, v)) return;
      last = v;
      this.paramsChanged();
    };
  }

  /** Loose equality across the forms an editor and the stored value can take:
   *  scalars compare by string form (`5` vs `'5'`, `null`/`undefined` vs `''`),
   *  arrays/objects (list inputs) by JSON. */
  static sameValue(a: unknown, b: unknown): boolean {
    if (a === b) return true;
    const isObj = (x: unknown): boolean => typeof x === 'object' && x !== null;
    if (isObj(a) || isObj(b)) {
      try {
        return JSON.stringify(a) === JSON.stringify(b);
      } catch {
        return false;
      }
    }
    return String(a ?? '') === String(b ?? '');
  }

  // ---------- panes ----------

  /** Compact chips replacing the old Function pane: full name, package, roles,
   *  tags — one wrapping row instead of a label+value row each. The function
   *  description lives in the header Description input now. */
  private buildFuncChips(node: FlowNode): HTMLElement {
    const chips = setTid(ui.div([], 'funcflow-chips'), 'prop-func-chips');
    const add = (text: string, tip: string, cls?: string, tid?: string): void => {
      const chip = ui.div([], 'funcflow-chip' + (cls ? ` ${cls}` : ''));
      chip.textContent = text;
      ui.tooltip.bind(chip, tip);
      if (tid) setTid(chip, tid);
      chips.appendChild(chip);
    };
    const fullName = node.dgFuncName ?? node.dgFunc?.name ?? '';
    if (fullName) add(fullName, 'Full function name', 'funcflow-chip-muted', 'prop-func-fullname');
    // Package disambiguates a vague function name (e.g. which "Descriptors").
    if (node.dgPackageName) add(node.dgPackageName, 'Package', undefined, 'prop-func-package');
    const roles = (node.dgRole ?? '').split(',').map((s) => s.trim()).filter(Boolean);
    for (const r of roles) add(r, 'Role');
    const tags = node.dgFunc ? getTags(node.dgFunc) : [];
    for (const t of tags.filter((t) => !roles.some((r) => r.toLowerCase() === t.toLowerCase())))
      add(`#${t}`, 'Tag');
    return chips;
  }

  private addFuncNodePanes(acc: DG.Accordion, node: FlowNode): void {
    const func = node.dgFunc;
    if (!func) return;

    if (func.inputs.length > 0) {
      // The pane is titled with the function itself — it IS the function's
      // parameter form (chips above carry package/role/tags).
      let paneTitle = '';
      try {
        paneTitle = getFuncDisplayName(func);
      } catch {/* Dart proxy access can throw */}
      if (!paneTitle) paneTitle = 'Parameters';
      const dataframeParams = func.inputs.filter((p) => String(p.propertyType) === 'dataframe').map((p) => p.name);
      const pane = acc.addPane(paneTitle, () => {
        const content = ui.div([], 'funcflow-accordion-content ui-form');
        for (const inp of func.inputs) {
          const tip = buildFuncInputTooltip(inp);
          // Display label — the property's caption when declared, else its name.
          // Purely cosmetic: the slot key / `inputValues` stay keyed by name.
          const label = getParamDisplayName(inp);
          const isEditable = inp.name in node.inputValues;
          if (!isEditable) {
            const row = ui.div([ui.divText(`${label}: ${inp.propertyType} (connected only)`)], 'funcflow-prop-row');
            ui.tooltip.bind(row, tip);
            content.appendChild(row);
            continue;
          }
          if (this.flow.isInputConnected(node.id, inp.name)) {
            const row = ui.div([ui.divText(`${label}: connected`)], 'funcflow-prop-row');
            ui.tooltip.bind(row, tip);
            content.appendChild(row);
            continue;
          }
          const pt = String(inp.propertyType);
          if (pt === 'column' || pt === 'column_list')
            content.appendChild(this.createColumnRow(node, inp.name, pt === 'column_list', dataframeParams, tip, label));
          else if (pt === 'string_list') {
            // Comma-separated string (native DG input, to match the primitives);
            // the compiler trims and turns it into an array. (`list<string>`
            // arrives here as `string_list` — DG normalizes it.)
            content.appendChild(this.createStringInput(inp.name, String(node.inputValues[inp.name] ?? ''),
              (v) => {node.inputValues[inp.name] = v;}, `${tip} | Comma-separated list of strings`, label));
          } else if (pt === 'list' || PRIMITIVE_INPUT_TYPES.has(pt)) {
            // A native Datagrok input built straight from the property — it
            // honours the declared type, numeric range, choices, and nullability
            // (choice-bearing strings render as a combo automatically). A `list`
            // property (incl. list<string>: propertyType 'list', subtype
            // 'string') gets DG's List input; its value is a JS array.
            content.appendChild(this.createPropertyInput(inp, node, tip));
          }
        }
        return content;
      }, true);
      this.decorateEditorHeader(pane, node, func);
    }
  }

  /** Functions with their own custom editor (an `editor:` meta, or the explicit
   *  allowlist — e.g. AddNewColumn) get a small "Open editor" button in the
   *  parameters pane header (the pane titled with the function name) that opens
   *  that editor seeded with the node's real upstream tables. Rendered only
   *  when the view wired `onEditFuncParams`. */
  private decorateEditorHeader(pane: DG.AccordionPane, node: FlowNode, func: DG.Func): void {
    if (!this.onEditFuncParams) return;
    let hasEditor = false;
    try {
      hasEditor = shouldUseFunctionEditor(func);
    } catch {/* Dart proxy access can throw — treat as no editor */}
    if (!hasEditor) return;
    const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
    if (!header) return;
    const btn = document.createElement('button');
    btn.textContent = 'Open editor';
    btn.classList.add('funcflow-func-editor-btn');
    ui.tooltip.bind(btn, 'Edit the parameters in the function’s own dialog (needs all table inputs connected)');
    setTid(btn, 'prop-func-editor');
    btn.onclick = (ev): void => {
      ev.stopPropagation(); // don't toggle the pane
      this.onEditFuncParams!(node);
    };
    header.appendChild(btn);
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
    // The table the viewer plots — its column options pick from it.
    const tableParam = this.dataframeInputKeys(node)[0];
    const setLook = (key: string) => (v: unknown): void => {
      const s = String(v ?? '').trim();
      if (s) look[key] = s;
      else delete look[key];
      void this.flow.updateNode(node.id);
    };

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
        // Column option → name field + picker dialog (seeded by the wired table);
        // non-column option (e.g. Title) → plain text.
        if (o.kind === 'column') {
          content.appendChild(this.createColumnFieldRow({
            nodeId: node.id, label: o.label, isList: false,
            tip: 'Column name (or wire a column into the socket)',
            getValue: () => String(look[o.key] ?? ''),
            setValue: setLook(o.key),
            tableParam,
          }));
          continue;
        }
        content.appendChild(this.createTextarea(o.label, String(look[o.key] ?? ''), setLook(o.key)));
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

    // A 'table' (dataframe) input means column-valued properties can be picked.
    const tableParam = this.dataframeInputKeys(node)[0];

    acc.addPane('Configuration', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      const nodeTips = UTILITY_PROP_TOOLTIPS[kind] ?? {};
      for (const [key, val] of props) {
        const tip = nodeTips[key];
        // Display caption mirrors core's fallback humanization; `key` stays
        // the row identity (data-param / storage).
        const caption = propertyNameToFriendly(key);
        const isConstValue = isConstant && key === 'value';
        // Column-valued props (Select Column → columnName, Select Columns →
        // columnNames) get the picker dialog when there's a table to pick from.
        if (tableParam && (key === 'columnName' || key === 'columnNames')) {
          content.appendChild(this.createColumnFieldRow({
            nodeId: node.id, label: key, caption, isList: key === 'columnNames',
            tip: tip ?? (key === 'columnNames' ? 'Comma-separated column names' : 'Column name'),
            getValue: () => String(node.properties[key] ?? ''),
            setValue: (v) => {node.properties[key] = v;},
            tableParam,
          }));
          continue;
        }
        if (typeof val === 'boolean') {
          content.appendChild(this.createToggle(key, val, (v) => {
            node.properties[key] = v;
            if (isConstValue) retitle(v);
          }, tip, caption));
        } else if (typeof val === 'number') {
          const isInt = Number.isInteger(val);
          content.appendChild(this.createNumberInput(key, val,
            (v) => {
              node.properties[key] = isInt ? Math.round(v) : v;
              if (isConstValue) retitle(node.properties[key]);
            },
            isInt ? 0 : 3, isInt ? 1 : 0.1, tip, caption));
        } else {
          content.appendChild(this.createTextarea(key, String(val ?? ''), (v) => {
            node.properties[key] = v;
            if (isConstValue) retitle(v);
          }, tip, false, caption));
        }
      }
      return content;
    }, true);
  }

  /** What's actually wired (disconnected slots are noise — the sockets on the
   *  canvas already show them), plus what's still MISSING: required inputs
   *  neither connected nor filled and required properties left empty. The pane
   *  opens expanded when something is missing. */
  private addConnectionsPane(acc: DG.Accordion, node: FlowNode): void {
    const isConnected = (key: string): boolean => this.flow.isInputConnected(node.id, key);
    const missingInputs = missingRequiredInputs(node, isConnected);
    const missingProps = missingRequiredProps(node);
    const hasMissing = missingInputs.length + missingProps.length > 0;

    acc.addPane('Connections', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      const ptCount = node.passthroughCount;

      const inputEntries = Object.entries(node.inputs) as Array<[string, {socket: {dgType: string}; label?: string} | undefined]>;
      const outputEntries = Object.entries(node.outputs) as Array<[string, {socket: {dgType: string}; label?: string} | undefined]>;

      if (hasMissing) {
        content.appendChild(this.connGroupLabel('Missing'));
        for (const label of missingInputs)
          content.appendChild(this.buildMissingRow(label, 'required — connect or set a value'));
        for (const key of missingProps)
          content.appendChild(this.buildMissingRow(propertyNameToFriendly(key), 'required value not set', key));
      }

      const conns = this.flow.getConnections();
      let anyConnected = false;
      const addGroup = (label: string, rows: HTMLElement[]): void => {
        if (rows.length === 0) return;
        if (anyConnected || hasMissing) content.appendChild(this.connSeparator());
        content.appendChild(this.connGroupLabel(label));
        rows.forEach((r) => content.appendChild(r));
        anyConnected = true;
      };
      const targetsOf = (key: string): string[] => conns
        .filter((c) => c.source === node.id && c.sourceOutput === key)
        .map((c) => this.endpointText(String(c.target), 'input', String(c.targetInput)));

      addGroup('Inputs', inputEntries
        .filter(([key, input]) => input && !isExecKey(key) && isConnected(key))
        .map(([key, input]) => {
          const src = this.flow.getInputSource(node.id, key);
          return this.buildConnRow('IN', input!.label ?? key, input!.socket.dgType,
            '←', src ? [this.endpointText(src.node.id, 'output', src.outputKey)] : [], key);
        }));

      addGroup('Pass-through', outputEntries.slice(0, ptCount)
        .filter(([key, out]) => out && targetsOf(key).length > 0)
        .map(([key, out]) => this.buildConnRow('PT',
          propertyNameToFriendly(key.endsWith('__pt') ? key.slice(0, -'__pt'.length) : key),
          out!.socket.dgType, '→', targetsOf(key), key)));

      addGroup('Outputs', outputEntries.slice(ptCount)
        .filter(([key, out]) => out && !isExecKey(key) && targetsOf(key).length > 0)
        .map(([key, out]) => this.buildConnRow('OUT', out!.label ?? key, out!.socket.dgType, '→', targetsOf(key), key)));

      // Order edges (exec ports) carry no data — show them as plain run-order
      // facts instead of IN/OUT rows with a raw `__exec_*` key.
      const nodeLabel = (id: string): string => String(this.flow.getNodeById(id)?.label ?? '?');
      addGroup('Run order', [
        ...conns.filter((c) => c.target === node.id && String(c.targetInput) === EXEC_IN_KEY)
          .map((c) => this.buildOrderRow('after', nodeLabel(String(c.source)))),
        ...conns.filter((c) => c.source === node.id && String(c.sourceOutput) === EXEC_OUT_KEY)
          .map((c) => this.buildOrderRow('before', nodeLabel(String(c.target)))),
      ]);

      if (!anyConnected && !hasMissing)
        content.appendChild(ui.divText('Nothing connected yet', 'funcflow-conn-empty'));
      return content;
    }, hasMissing);
  }

  /** "Node title · slot label" for the far end of a connection. A pass-through
   *  source renders as its humanized base input name (its literal label is
   *  just `→`). */
  private endpointText(nodeId: string, side: 'input' | 'output', key: string): string {
    const n = this.flow.getNodeById(nodeId);
    const name = String(n?.label ?? '?');
    let slot = propertyNameToFriendly(key.endsWith('__pt') ? key.slice(0, -'__pt'.length) : key);
    const ports = (side === 'input' ? n?.inputs : n?.outputs) as
      Record<string, {label?: string} | undefined> | undefined;
    const lbl = ports?.[key]?.label;
    if (lbl && lbl !== '→') slot = lbl;
    return `${name} · ${slot}`;
  }

  private buildConnRow(dir: string, name: string, type: string, arrow: string, ends: string[], key: string): HTMLElement {
    const dirSpan = ui.element('span');
    dirSpan.textContent = dir;
    dirSpan.className = 'funcflow-conn-dir';
    const detail = ui.element('span');
    detail.textContent = ` ${name} `;
    const typeSpan = ui.element('span');
    typeSpan.textContent = `(${type})`;
    typeSpan.className = 'funcflow-conn-type';
    const children = [dirSpan, detail, typeSpan];
    if (ends.length > 0) {
      const arrowSpan = ui.element('span');
      arrowSpan.textContent = ` ${arrow} `;
      arrowSpan.className = 'funcflow-conn-arrow';
      const endSpan = ui.element('span');
      endSpan.textContent = ends.join(', ');
      endSpan.className = 'funcflow-conn-endpoint';
      children.push(arrowSpan, endSpan);
    }
    const row = ui.div(children, 'funcflow-prop-row funcflow-conn-row');
    row.dataset.conn = key;
    return row;
  }

  private buildOrderRow(kind: 'after' | 'before', otherLabel: string): HTMLElement {
    const detail = ui.element('span');
    detail.textContent = `runs ${kind} `;
    detail.className = 'funcflow-conn-type';
    const endSpan = ui.element('span');
    endSpan.textContent = otherLabel;
    endSpan.className = 'funcflow-conn-endpoint';
    const row = ui.div([detail, endSpan], 'funcflow-prop-row funcflow-conn-row');
    row.dataset.conn = kind === 'after' ? EXEC_IN_KEY : EXEC_OUT_KEY;
    return row;
  }

  private buildMissingRow(label: string, why: string, key = label): HTMLElement {
    const warn = ui.element('span');
    warn.textContent = '⚠ ';
    const detail = ui.element('span');
    detail.textContent = `${label} `;
    const whySpan = ui.element('span');
    whySpan.textContent = `— ${why}`;
    whySpan.className = 'funcflow-conn-type';
    const row = ui.div([warn, detail, whySpan], 'funcflow-prop-row funcflow-conn-row funcflow-conn-missing');
    row.dataset.missing = key;
    return row;
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

  /** `caption` is display-only; `label` stays the row identity (data-param). */
  private createTextarea(
    label: string, value: string, onChange: (v: string) => void, inputTooltip?: string, cosmetic = false,
    caption?: string,
  ): HTMLElement {
    const report = this.changeReporter(value);
    const apply = cosmetic ? onChange : (v: string): void => {
      onChange(v);
      report(v);
    };
    return this.propRow(ui.div([this.labelWithTooltip(caption ?? label, inputTooltip),
      this.buildTextareaEl(value, apply, inputTooltip)], 'funcflow-prop-row'), label);
  }

  /** Initialize a DG input's editor from a stored value via the `stringValue`
   *  setter — `ui.input.forProperty` (and the `value` init option) does not
   *  reliably load the editor itself. Guarded: a value the editor can't parse
   *  just leaves it at its own blank/default state. */
  private static initInputValue(input: DG.InputBase, v: unknown): void {
    try {
      if (v !== undefined && v !== null && String(v) !== '') input.stringValue = String(v);
    } catch {/* leave the editor as-is */}
  }

  /** A native Datagrok input for a primitive func parameter, built from the
   *  property so it honours the declared type, numeric range, choices, and
   *  nullability — replacing the bespoke combo/number/toggle wiring. Edits are
   *  written back to the node's `inputValues`; the caption comes from the
   *  property. The row (not the DG root) carries the `data-param` / test-id so
   *  the field stays addressable by name. */
  private createPropertyInput(param: DG.Property, node: FlowNode, inputTooltip: string): HTMLElement {
    const report = this.changeReporter(node.inputValues[param.name]);
    const input = ui.input.forProperty(param, null, {
      tooltipText: inputTooltip,
      onValueChanged: (v) => {
        node.inputValues[param.name] = v;
        report(v);
      },
    });
    PropertyPanel.initInputValue(input, node.inputValues[param.name]);
    return this.propRow(ui.div([input.root], 'funcflow-prop-row funcflow-dg-row'), param.name);
  }

  /** A native Datagrok single-line string input (used where there's no
   *  `DG.Property` to drive `forProperty` — e.g. a `string_list` slot edited as
   *  a comma-separated string), so its styling matches the primitive inputs. */
  private createStringInput(
    label: string, value: string, onChange: (v: string) => void, inputTooltip?: string, caption?: string,
  ): HTMLElement {
    const report = this.changeReporter(value);
    const input = ui.input.string(caption ?? label, {
      tooltipText: inputTooltip,
      onValueChanged: (v) => {
        onChange(String(v ?? ''));
        report(String(v ?? ''));
      },
    });
    PropertyPanel.initInputValue(input, value);
    // Display caption may differ from the identity; keep the row keyed by name.
    return this.propRow(ui.div([input.root], 'funcflow-prop-row funcflow-dg-row'), label);
  }

  /** A column / column-list input laid out side by side with its table picker:
   *  the column-name field (≈70%) and, when the func has 2+ dataframe inputs, a
   *  table chooser (≈30%) writing the node's `columnTables` association. With a
   *  single dataframe input the field spans full width (the table is implicit). */
  /** Input keys on a node whose socket is a dataframe — the tables a column
   *  field can pick from. (Func nodes may have several; viewers / Select Column
   *  utilities have one named 'table'.) */
  private dataframeInputKeys(node: FlowNode): string[] {
    return (Object.entries(node.inputs) as Array<[string, {socket: {dgType: string}} | undefined]>)
      .filter(([, inp]) => inp?.socket.dgType === 'dataframe')
      .map(([k]) => k);
  }

  /** A column / column-list field laid out with its picker: the name field, an
   *  optional table-chooser (multi-table funcs), and the picker icon that opens
   *  a dialog seeded by the chosen table input. Storage is fully delegated via
   *  `getValue`/`setValue`, so this serves **any** node with a column field —
   *  func inputs (`inputValues`), viewer look options, Select Column utilities.
   *  `label` is the identity — the picker's `data-param` / test-id key; the
   *  visible caption is `caption ?? label` (display only). */
  private createColumnFieldRow(opts: {
    nodeId: string;
    label: string;
    isList: boolean;
    tip: string;
    getValue: () => string;
    setValue: (v: string) => void;
    /** Display caption when it differs from the identity `label` (optional). */
    caption?: string;
    /** Single dataframe input the columns come from (most nodes). */
    tableParam?: string;
    /** Per-row table chooser for multi-table funcs (JoinTables keys1→table1, …). */
    tableSelect?: {options: string[]; get: () => string; set: (v: string) => void};
  }): HTMLElement {
    // A column / column-list value is a plain (comma-separated) name string, so
    // it uses a native Datagrok string input with its own caption — not
    // `forProperty`, which would build a column picker bound to a live table we
    // don't have here. The table chooser (multi-table funcs) and the picker icon
    // are appended *inside* the input via `addOptions` (trailing controls).
    const report = this.changeReporter(opts.getValue());
    const nameInput = ui.input.string(opts.caption ?? opts.label, {
      tooltipText: opts.tip,
      onValueChanged: (v) => {
        opts.setValue(String(v ?? ''));
        report(String(v ?? ''));
      },
    });
    PropertyPanel.initInputValue(nameInput, opts.getValue());
    nameInput.input.style.minWidth = '70px';

    let getTableParam = (): string => opts.tableParam ?? '';
    if (opts.tableSelect) {
      const ts = opts.tableSelect;
      const select = this.buildSelectEl(ts.get(), ts.options, ts.set, 'Which table input this column refers to');
      select.classList.add('funcflow-col-table-select');
      getTableParam = (): string => select.value || ts.options[0];
      nameInput.addOptions(select);
    }

    // Column chooser — opens a dialog seeded by the upstream table (running the
    // flow up to that point if needed) so users pick from a real column list.
    if (this.onPickColumns && (opts.tableSelect || opts.tableParam)) {
      const pickBtn = ui.iconFA('list', () => {
        this.onPickColumns!({
          nodeId: opts.nodeId, paramName: opts.label, isList: opts.isList,
          tableParam: getTableParam(),
          current: opts.getValue(),
          anchor: pickBtn,
          apply: (value: string) => {
            nameInput.value = value; // fires onValueChanged → opts.setValue
            opts.setValue(value);
          },
        });
      }, opts.isList ? 'Choose columns from the connected table' : 'Choose a column from the connected table');
      setTid(pickBtn, 'prop-pick-columns', opts.label);
      pickBtn.classList.add('funcflow-col-pick');
      nameInput.addOptions(pickBtn);
    }

    return this.propRow(ui.div([nameInput.root], 'funcflow-prop-row funcflow-dg-row'), opts.label);
  }

  private createColumnRow(
    node: FlowNode, paramName: string, isList: boolean, dataframeParams: string[], tip: string, caption?: string,
  ): HTMLElement {
    const colTip = isList ?
      `${tip} | Comma-separated column names` :
      `${tip} | Column name (compiled to table.col(...))`;

    // Single dataframe input → fixed table; multi-table funcs → a per-row select
    // writing the node's `columnTables` association (keys1→table1, keys2→table2).
    let tableParam: string | undefined = dataframeParams[0];
    let tableSelect: {options: string[]; get: () => string; set: (v: string) => void} | undefined;
    if (dataframeParams.length >= 2) {
      if (!node.properties['columnTables']) node.properties['columnTables'] = {};
      const associations = node.properties['columnTables'] as Record<string, string>;
      tableParam = undefined;
      tableSelect = {
        options: dataframeParams,
        get: () => associations[paramName] ?? dataframeParams[0],
        set: (v) => {associations[paramName] = v;},
      };
    }

    return this.createColumnFieldRow({
      nodeId: node.id, label: paramName, caption, isList, tip: colTip,
      getValue: () => String(node.inputValues[paramName] ?? ''),
      setValue: (v) => {node.inputValues[paramName] = v;},
      tableParam, tableSelect,
    });
  }

  private createNumberInput(label: string, value: number, onChange: (v: number) => void, decimals: number, step: number, inputTooltip?: string, caption?: string): HTMLElement {
    const input = document.createElement('input');
    input.type = 'number';
    input.value = decimals === 0 ? String(Math.round(value)) : value.toFixed(decimals);
    input.step = String(step);
    input.className = 'funcflow-prop-input';
    const report = this.changeReporter(value);
    input.addEventListener('change', () => {
      const parsed = parseFloat(input.value);
      if (!isNaN(parsed)) {
        onChange(parsed);
        report(parsed);
      }
    });
    if (inputTooltip) ui.tooltip.bind(input, inputTooltip);
    return this.propRow(ui.div([this.labelWithTooltip(caption ?? label, inputTooltip), input], 'funcflow-prop-row'), label);
  }

  private createToggle(label: string, value: boolean, onChange: (v: boolean) => void, inputTooltip?: string, caption?: string): HTMLElement {
    const input = document.createElement('input');
    input.type = 'checkbox';
    input.checked = value;
    input.className = 'funcflow-prop-checkbox';
    const report = this.changeReporter(value);
    input.addEventListener('change', () => {
      onChange(input.checked);
      report(input.checked);
    });
    if (inputTooltip) ui.tooltip.bind(input, inputTooltip);
    const lbl = this.labelWithTooltip(caption ?? label, inputTooltip);
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
    const report = this.changeReporter(value);
    select.addEventListener('change', () => {
      onChange(select.value);
      report(select.value);
    });
    if (inputTooltip) ui.tooltip.bind(select, inputTooltip);
    return select;
  }

  private createCombo(label: string, value: string, options: string[], onChange: (v: string) => void, inputTooltip?: string): HTMLElement {
    return this.propRow(ui.div([this.labelWithTooltip(label, inputTooltip),
      this.buildSelectEl(value, options, onChange, inputTooltip)], 'funcflow-prop-row'), label);
  }
}

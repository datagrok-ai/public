/* eslint-disable max-len */
/** Property panel — renders into Datagrok's native context panel. */

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
import {shouldUseFunctionEditor, hasEditorShortcut} from '../utils/func-editor-utils';
import {
  hiddenInputsOf, customEditorFor, CustomInputEditorFactory, effectiveFuncInputs, inputHiddenByCondition,
  inputCaptionOf,
} from '../utils/func-input-overrides';
import {buildInputValueEditor} from '../utils/input-values';
import {ColumnPickRequest} from './column-picker';
import { processChoiceInput } from './choice-input-processor';

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

export function propertyChoices(param: DG.Property): string[] {
  try {
    const choices: unknown = param.choices;
    if (!Array.isArray(choices)) return [];
    return choices.map((c) => String(c)).filter((c) => c.length > 0);
  } catch {
    return [];
  }
}

/** Choice options for a string input; nullable adds a leading empty option, and the current value is kept even when undeclared. */
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
  private currentNode: FlowNode | null = null;
  private currentExecState?: NodeExecState;
  private execSection: HTMLElement | null = null;

  /** Set by the view: opens a column / columns picker dialog for a func-node column input; unset → no picker icon. */
  onPickColumns?: (req: ColumnPickRequest) => void;

  /** Set by the view: opens the function's own editor dialog and writes edits back into `inputValues`; unset → no editor icon. */
  onEditFuncParams?: (node: FlowNode) => void;

  constructor(flow: FlowEditor) {
    this.flow = flow;
    this.contentDiv = setTid(ui.div([], 'funcflow-property-content'), 'property-content');
    this.root = setTid(ui.divV([this.contentDiv], 'funcflow-property-panel'), 'property-panel');
  }

  showNode(node: FlowNode, execState?: NodeExecState): void {
    // Release watchers/custom editors BEFORE the DOM they feed is thrown away.
    this.disposeEditors();
    this.contentDiv.innerHTML = '';
    this.inputWatchers.clear();
    this.currentNode = node;
    this.currentExecState = execState;

    const titleInput = this.createTextarea('Title', String(node.label ?? ''), (v) => {
      node.label = String(v ?? '');
      void this.flow.updateNode(node.id);
    }, undefined, true);
    const titleRow = setTid(ui.div([titleInput], 'funcflow-title-row'), 'property-title-row');

    const header = setTid(ui.div([titleRow], 'funcflow-panel-header'), 'property-header');

    header.appendChild(this.buildFuncChips(node));

    let funcDesc = '';
    try {
      funcDesc = node.dgFunc?.description ?? '';
    } catch {/* Dart proxy access can throw */}
    const descSeed = node.description?.trim() ? node.description : funcDesc;
    header.appendChild(this.createTextAreaRow('Description', descSeed, (v) => {
      node.description = v;
      void this.flow.updateNode(node.id);
    }, true, 'Describe what this step does…'));
    this.contentDiv.appendChild(header);

    const acc = ui.accordion('funcflow-context-panel');

    if (node.dgFunc) this.addFuncNodePanes(acc, node);
    if (node.dgNodeType === 'input') this.addInputNodePane(acc, node);
    if (node.dgNodeType === 'output') this.addOutputNodePane(acc, node);
    if (node.properties['viewerType']) this.addViewerNodePane(acc, node);
    else if (node.dgNodeType === 'utility') this.addUtilityNodePane(acc, node);

    this.addConnectionsPane(acc, node);

    this.contentDiv.appendChild(acc.root);

    // Own container: a run finishing while the user types refreshes just this section (the full rebuild is focus-guarded).
    this.execSection = ui.div([]);
    this.renderExecSection();
    this.contentDiv.appendChild(this.execSection);
  }

  private renderExecSection(): void {
    if (!this.execSection) return;
    this.execSection.innerHTML = '';
    if (!this.currentExecState) return;
    const header = ui.div([], 'funcflow-prop-section-header');
    header.textContent = 'Execution';
    this.execSection.appendChild(header);
    this.execSection.appendChild(buildExecutionMeta(this.currentExecState));
  }

  clear(): void {
    this.disposeEditors();
    this.currentNode = null;
    this.currentExecState = undefined;
    this.execSection = null;
    this.contentDiv.innerHTML = '';
    this.contentDiv.appendChild(ui.divText('Select a node to view its properties'));
  }

  get shownNodeId(): string | null {
    return this.currentNode?.id ?? null;
  }

  updateExecState(nodeId: string, execState?: NodeExecState): void {
    if (this.currentNode?.id !== nodeId || this.currentExecState === execState) return;
    this.currentExecState = execState;
    // The Execution section holds no inputs, so it re-renders even when the full rebuild below is focus-guarded.
    this.renderExecSection();
    this.refreshShownNode();
  }

  /** Rebuild the panel for the shown node — skipped while the user is typing here (a rebuild would steal focus). */
  refreshShownNode(): void {
    if (!this.currentNode) return;
    if (this.root.contains(document.activeElement) && document.activeElement !== document.body) return;
    // Never rebuild under an open modal — a dialog launched from this panel holds live objects the rebuild destroys.
    if (DG.Dialog.getOpenDialogs().length > 0) return;
    if (!this.flow.getNodes().some((n) => n.id === this.currentNode!.id)) {
      this.clear();
      return;
    }
    this.showNode(this.currentNode, this.currentExecState);
  }

  private paramsChanged(): void {
    if (this.currentNode) {
      this.flow.notifyNodeParamsChanged(this.currentNode.id);
      this.syncMissingRows(this.currentNode);
    }
  }

  /** Per-editor guard: report only a REAL change — creating/initializing a DG input can fire `onValueChanged` immediately, so an unguarded report would turn every panel rebuild into an edit. */
  private changeReporter(initial: unknown): (v: unknown) => void {
    let last = initial;
    return (v: unknown): void => {
      if (PropertyPanel.sameValue(last, v)) return;
      last = v;
      this.paramsChanged();
    };
  }

  /** Loose equality: scalars by string form (`5` ≡ `'5'`, `null` ≡ `''`), arrays/objects by JSON. */
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

  private buildFuncChips(node: FlowNode): HTMLElement {
    const chips = setTid(ui.div([], 'funcflow-chips'), 'prop-func-chips');
    const add = (text: string, tip: string, cls?: string, tid?: string): void => {
      if (chips.childElementCount > 0) chips.appendChild(ui.div([], 'funcflow-chip-sep'));
      const chip = ui.div([], 'funcflow-chip' + (cls ? ` ${cls}` : ''));
      chip.textContent = text;
      ui.tooltip.bind(chip, tip);
      if (tid) setTid(chip, tid);
      chips.appendChild(chip);
    };
    const kindWords: Record<string, string> = {
      func: 'Function', utility: 'Utility', input: 'Input', output: 'Output',
    };
    add(kindWords[node.dgNodeType ?? ''] ?? 'Function', 'Node kind', 'funcflow-chip-kind', 'property-type-badge');
    // Package disambiguates a vague function name (e.g. which "Descriptors") —
    const fullName = node.dgFuncName ?? node.dgFunc?.name ?? '';
    const pkg = node.dgPackageName ?? '';
    const nameChip = pkg && fullName.toLowerCase().startsWith(`${pkg.toLowerCase()}:`) ?
      fullName.slice(pkg.length + 1) : fullName;
    if (nameChip) add(nameChip, 'Full function name', 'funcflow-chip-muted', 'prop-func-fullname');
    if (pkg) add(pkg, 'Package', undefined, 'prop-func-package');
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

    // The parameters the NODE exposes — the wrapper's when one is registered, so the form matches the node's own sockets.
    const funcInputs = effectiveFuncInputs(func);
    if (funcInputs.length > 0) {
      let paneTitle = '';
      try {
        paneTitle = getFuncDisplayName(func);
      } catch {/* Dart proxy access can throw */}
      if (!paneTitle) paneTitle = 'Parameters';
      const dataframeParams = funcInputs.filter((p) => String(p.propertyType) === 'dataframe').map((p) => p.name);
      const hidden = hiddenInputsOf(func);
      const pane = acc.addPane(paneTitle, () => {
        const content = ui.div([], 'funcflow-accordion-content ui-form');
        for (const inp of funcInputs) {
          if (hidden.has(inp.name)) continue;
          if (inputHiddenByCondition(func, inp.name, node)) continue;
          const tip = buildFuncInputTooltip(inp);
          const label = inputCaptionOf(func, inp.name) ?? getParamDisplayName(inp);
          const custom = customEditorFor(func, inp.name);
          const isEditable = custom !== null || inp.name in node.inputValues;
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
          if (custom) {
            content.appendChild(this.createCustomInputRow(custom, inp, node, tip));
            continue;
          }
          const pt = String(inp.propertyType);
          if (pt === 'column' || pt === 'column_list')
            content.appendChild(this.createColumnRow(node, inp.name, pt === 'column_list', dataframeParams, tip, label));
          else if (pt === 'string_list') {
            content.appendChild(this.createStringInput(inp.name, String(node.inputValues[inp.name] ?? ''),
              (v) => {node.inputValues[inp.name] = v;}, `${tip} | Comma-separated list of strings`, label));
          } else if (pt === 'list' || PRIMITIVE_INPUT_TYPES.has(pt)) {
            content.appendChild(this.createPropertyInput(inp, node, tip));
          }
        }
        return content;
      }, true);
      this.decorateEditorHeader(pane, node, func);
    }
  }

  private decorateEditorHeader(pane: DG.AccordionPane, node: FlowNode, func: DG.Func): void {
    if (!this.onEditFuncParams) return;
    let hasEditor = false;
    try {
      hasEditor = shouldUseFunctionEditor(func);
    } catch {/* Dart proxy access can throw — treat as no editor */}
    if (!hasEditor) return;
    const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
    if (!header) return;
    const btn = ui.button('Open editor', () => this.onEditFuncParams!(node));
    btn.classList.add('funcflow-func-editor-btn');
    ui.tooltip.bind(btn, 'Edit the parameters in the function’s own dialog (needs all table inputs connected)');
    setTid(btn, 'prop-func-editor');
    // The click must not bubble into the accordion header (it would toggle the pane).
    btn.addEventListener('click', (ev) => ev.stopPropagation());
    header.appendChild(btn);
  }

  // eslint-disable-next-line complexity
  private addInputNodePane(acc: DG.Accordion, node: FlowNode): void {
    acc.addPane('Input Configuration', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      content.appendChild(this.createTextarea('Param Name', String(node.properties['paramName'] ?? ''), (v) => {
        node.properties['paramName'] = v;
      }));

      const valueEditor = buildInputValueEditor(node, () => this.paramsChanged());
      if (valueEditor)
        content.appendChild(this.propRow(ui.div([valueEditor.root], 'funcflow-prop-row funcflow-dg-row'), 'Value'));

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
        content.appendChild(this.createTextarea('Choices (comma-sep)', String(node.properties['choices'] ?? ''), (v) => {node.properties['choices'] = v;},
          'Comma-separated list of allowed values — the run dialog and the value editor show them as a dropdown', false, 'Choices'));
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
        void this.flow.updateNode(node.id);
      }));
      if (node.properties['outputType'] !== undefined) {
        content.appendChild(this.createCombo('Output Type', String(node.properties['outputType'] ?? 'dynamic'),
          OUTPUT_TYPE_VALUES, (v) => {
            node.properties['outputType'] = v;
            void this.flow.updateNode(node.id);
          }));
      }
      return content;
    }, true);
  }

  private addViewerNodePane(acc: DG.Accordion, node: FlowNode): void {
    if (!node.properties['viewerLook'] || typeof node.properties['viewerLook'] !== 'object')
      node.properties['viewerLook'] = {};
    const look = node.properties['viewerLook'] as Record<string, unknown>;
    const specs = (node.properties['viewerOptionSpecs'] as Array<{key: string; label: string; kind: string}>) ?? [];
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
        if (o.kind === 'column' && this.flow.isInputConnected(node.id, o.key)) {
          const row = ui.div([ui.divText(`${o.label}: connected`)], 'funcflow-prop-row');
          ui.tooltip.bind(row, 'A column is wired into this option; its name is used.');
          content.appendChild(row);
          continue;
        }
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

    const tableParam = this.dataframeInputKeys(node)[0];

    acc.addPane('Configuration', () => {
      const content = ui.div([], 'funcflow-accordion-content');
      const nodeTips = UTILITY_PROP_TOOLTIPS[kind] ?? {};
      for (const [key, val] of props) {
        const tip = nodeTips[key];
        const caption = propertyNameToFriendly(key);
        const isConstValue = isConstant && key === 'value';
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

  /** "Node title · slot label" for the far end; a pass-through label's trailing arrow is canvas affordance and is stripped. */
  private endpointText(nodeId: string, side: 'input' | 'output', key: string): string {
    const n = this.flow.getNodeById(nodeId);
    const name = String(n?.label ?? '?');
    let slot = propertyNameToFriendly(key.endsWith('__pt') ? key.slice(0, -'__pt'.length) : key);
    const ports = (side === 'input' ? n?.inputs : n?.outputs) as
      Record<string, {label?: string} | undefined> | undefined;
    const lbl = ports?.[key]?.label?.replace(/\s*→$/, '');
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

  /** Surgically drop resolved "missing" rows after an edit — no rebuild, so typing keeps focus. */
  private syncMissingRows(node: FlowNode): void {
    const rows = Array.from(this.contentDiv.querySelectorAll<HTMLElement>('.funcflow-conn-missing'));
    if (rows.length === 0) return;
    const isConnected = (key: string): boolean => this.flow.isInputConnected(node.id, key);
    const still = new Set<string>([
      ...missingRequiredInputs(node, isConnected),
      ...missingRequiredProps(node),
    ]);
    let removed = false;
    for (const row of rows) {
      if (!still.has(row.dataset.missing ?? '')) {
        row.remove();
        removed = true;
      }
    }
    if (removed && still.size === 0) {
      for (const lbl of Array.from(this.contentDiv.querySelectorAll<HTMLElement>('.funcflow-conn-group-label'))) {
        if (lbl.textContent === 'Missing') {
          if (lbl.nextElementSibling?.classList.contains('funcflow-conn-separator'))
            lbl.nextElementSibling.remove();
          lbl.remove();
        }
      }
    }
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

  private static propTip(label: string, caption?: string, explicitTip?: string): string | undefined {
    return explicitTip ?? PROP_TOOLTIPS[caption ?? label] ?? PROP_TOOLTIPS[label];
  }

  private propRow(el: HTMLElement, label: string): HTMLElement {
    setTid(el, 'prop-input', label);
    el.dataset.param = label;
    return el;
  }

  /** Every editor in this panel is a real Datagrok input (`ui.input.*`), never a hand-rolled element. */
  private dgRow(label: string, input: DG.InputBase): HTMLElement {
    return this.propRow(ui.div([input.root], 'funcflow-prop-row funcflow-dg-row'), label);
  }

  /** Single-line DG string row; `cosmetic` edits (Title, Description) never report an invalidating change. */
  private createTextarea(
    label: string, value: string, onChange: (v: string) => void, inputTooltip?: string, cosmetic = false,
    caption?: string,
  ): HTMLElement {
    const report = this.changeReporter(value);
    const input = ui.input.string(caption ?? label, {
      tooltipText: PropertyPanel.propTip(label, caption, inputTooltip),
      onValueChanged: (v) => {
        const s = String(v ?? '');
        onChange(s);
        if (!cosmetic) report(s);
      },
    });
    PropertyPanel.initInputValue(input, value);
    return this.dgRow(label, input);
  }

  private createTextAreaRow(
    label: string, value: string, onChange: (v: string) => void, cosmetic = false, placeholder?: string,
  ): HTMLElement {
    const report = this.changeReporter(value);
    const input = ui.input.textArea(label, {
      tooltipText: PropertyPanel.propTip(label),
      onValueChanged: (v) => {
        const s = String(v ?? '');
        onChange(s);
        if (!cosmetic) report(s);
      },
    });
    if (placeholder) {
      try {
        (input.input as HTMLTextAreaElement).placeholder = placeholder;
      } catch {/* editor element not a textarea on this build */}
    }
    PropertyPanel.initInputValue(input, value);
    return this.dgRow(label, input);
  }

  /** Initialize via the guarded `stringValue` setter — `ui.input.forProperty` / the `value` init option don't reliably load the editor. */
  private static initInputValue(input: DG.InputBase, v: unknown, setStringValue = true): void {
    try {
      if (v !== undefined && v !== null && String(v) !== '') setStringValue ? (input.stringValue = String(v)) : (input.value = v);
    } catch {/* leave the editor as-is */}
  }

  private createPropertyInput(param: DG.Property, node: FlowNode, inputTooltip: string): HTMLElement {
    const report = this.changeReporter(node.inputValues[param.name]);
    const input = ui.input.forProperty(param, null, {
      tooltipText: inputTooltip,
      onValueChanged: (v) => {
        node.inputValues[param.name] = v;
        report(v);
        this.notifyInputChanged(param.name, v);
      },
    });
    if (param.choices && input instanceof DG.ChoiceInput && node.dgFunc) {
      PropertyPanel.initInputValue(input, node.inputValues[param.name], false);
      processChoiceInput(input, node.dgFunc, param);
    } else {
       PropertyPanel.initInputValue(input, node.inputValues[param.name]);
    }
    this.addEditorShortcut(input, node, param.name);
    return this.propRow(ui.div([input.root], 'funcflow-prop-row funcflow-dg-row'), param.name);
  }

  private addEditorShortcut(input: DG.InputBase, node: FlowNode, paramName: string): void {
    if (!this.onEditFuncParams || !node.dgFunc || !hasEditorShortcut(node.dgFunc, paramName)) return;
    const pencil = ui.iconFA('pencil', () => this.onEditFuncParams!(node),
      'Edit in the function’s own dialog (needs all table inputs connected)');
    pencil.classList.add('funcflow-input-editor-pencil');
    setTid(pencil, `prop-input-editor-${paramName}`);
    input.addOptions(pencil);
  }

  private createCustomInputRow(
    factory: CustomInputEditorFactory, param: DG.Property, node: FlowNode, tip: string,
  ): HTMLElement {
    const report = this.changeReporter(node.inputValues[param.name]);
    const ed = factory(param, {
      inputValue: (name) => node.inputValues[name],
      // Captured columns/tables only — resolving an uncomputed table would mean running the flow while a panel renders.
      columns: (tableParam) => this.upstreamColumns(node, tableParam),
      table: (tableParam) => this.upstreamTable(node, tableParam),
      isConnected: (tableParam) => this.flow.isInputConnected(node.id, tableParam),
      produceTable: (tableParam) => this.produceUpstreamTable(node, tableParam),
      watch: (name, cb) => this.watchInput(name, cb),
      node,
    });
    this.editorDisposers.push(ed);
    ed.onChanged = (v): void => {
      if (ed.isValid && !ed.isValid()) return;
      node.inputValues[param.name] = v;
      report(v);
      this.notifyInputChanged(param.name, v);
    };
    ed.setValue(node.inputValues[param.name]);
    ui.tooltip.bind(ed.element, tip);
    return this.propRow(ui.div([ed.element], 'funcflow-prop-row funcflow-dg-row'), param.name);
  }

  /** Live per-input subscriptions for editors depending on a sibling parameter — refreshShownNode skips itself while focus is inside the panel, so a rebuild can't deliver these. */
  private readonly inputWatchers = new Map<string, Array<(v: unknown) => void>>();

  private watchInput(name: string, cb: (v: unknown) => void): void {
    const list = this.inputWatchers.get(name);
    if (list) list.push(cb);
    else this.inputWatchers.set(name, [cb]);
  }

  private notifyInputChanged(name: string, value: unknown): void {
    for (const cb of this.inputWatchers.get(name) ?? []) {
      try {
        cb(value);
      } catch (e) {
        console.error(`Flow: input watcher for "${name}" failed`, e);
      }
    }
  }

  /** Wired by the view to `ExecutionController.cloneForNode`; unset in headless editors. */
  getUpstreamTable?: (sourceNodeId: string) => DG.DataFrame | null;

  /** Wired by the view to `ExecutionController.produceTableForNode`; only ever called from an explicit user action, never from a render. */
  runUpstreamNode?: (sourceNodeId: string) => Promise<DG.DataFrame | null>;

  private upstreamTable(node: FlowNode, tableParam: string): DG.DataFrame | null {
    if (!this.getUpstreamTable) return null;
    const src = this.flow.getInputSource(node.id, tableParam);
    return src ? this.getUpstreamTable(src.node.id) : null;
  }

  private upstreamColumns(node: FlowNode, tableParam: string): DG.Column[] | null {
    const table = this.upstreamTable(node, tableParam);
    return table ? Array.from(table.columns) : null;
  }

  private async produceUpstreamTable(node: FlowNode, tableParam: string): Promise<DG.DataFrame | null> {
    const src = this.flow.getInputSource(node.id, tableParam);
    if (!src || !this.runUpstreamNode) return null;
    return this.runUpstreamNode(src.node.id);
  }

  private readonly editorDisposers: Array<{detach?: () => void}> = [];

  private disposeEditors(): void {
    for (const ed of this.editorDisposers.splice(0)) {
      try {
        ed.detach?.();
      } catch (e) {
        console.error('Flow: custom editor cleanup failed', e);
      }
    }
  }

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
    return this.propRow(ui.div([input.root], 'funcflow-prop-row funcflow-dg-row'), label);
  }

  private dataframeInputKeys(node: FlowNode): string[] {
    return (Object.entries(node.inputs) as Array<[string, {socket: {dgType: string}} | undefined]>)
      .filter(([, inp]) => inp?.socket.dgType === 'dataframe')
      .map(([k]) => k);
  }

  /** Column / column-list field with its table chooser and picker icon; storage is delegated via `getValue`/`setValue`, so it serves any node with a column field. */
  private createColumnFieldRow(opts: {
    nodeId: string;
    label: string;
    isList: boolean;
    tip: string;
    getValue: () => string;
    setValue: (v: string) => void;
    caption?: string;
    tableParam?: string;
    tableSelect?: {options: string[]; labels?: string[]; get: () => string; set: (v: string) => void};
  }): HTMLElement {
    // ui.input.string, not `forProperty` — forProperty would build a column picker bound to a live table we don't have here.
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
      const items = ts.labels ?? ts.options;
      const toParam = (label: unknown): string => {
        const i = items.indexOf(String(label ?? ''));
        return ts.options[i >= 0 ? i : 0];
      };
      const toLabel = (param: string): string => items[Math.max(0, ts.options.indexOf(param))];
      const reportTable = this.changeReporter(ts.get());
      const tableChoice = ui.input.choice('', {
        items,
        tooltipText: 'Which table input this column refers to',
        onValueChanged: (v) => {
          const s = toParam(v);
          ts.set(s);
          reportTable(s);
        },
      });
      PropertyPanel.initInputValue(tableChoice, toLabel(ts.get()), false);
      tableChoice.root.classList.add('funcflow-col-table-select');
      getTableParam = (): string => toParam(tableChoice.value);
      nameInput.addOptions(tableChoice.root);
    }

    if (this.onPickColumns && (opts.tableSelect || opts.tableParam)) {
      const pickBtn = ui.iconFA('list', () => {
        this.onPickColumns!({
          nodeId: opts.nodeId, paramName: opts.label, isList: opts.isList,
          tableParam: getTableParam(),
          current: opts.getValue(),
          anchor: pickBtn,
          apply: (value: string) => {
            nameInput.value = value;
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

    let tableParam: string | undefined = dataframeParams[0];
    let tableSelect: {options: string[]; labels?: string[]; get: () => string; set: (v: string) => void} | undefined;
    if (dataframeParams.length >= 2) {
      if (!node.properties['columnTables']) node.properties['columnTables'] = {};
      const associations = node.properties['columnTables'] as Record<string, string>;
      tableParam = undefined;
      tableSelect = {
        options: dataframeParams,
        labels: dataframeParams.map((p) =>
          (node.dgFunc ? inputCaptionOf(node.dgFunc, p) : null) ?? p),
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
    const report = this.changeReporter(value);
    const opts = {
      tooltipText: PropertyPanel.propTip(label, caption, inputTooltip),
      onValueChanged: (v: number | null) => {
        if (v == null || isNaN(v)) return; // mid-edit blank — keep the stored value
        onChange(v);
        report(v);
      },
    };
    const input = decimals === 0 ? ui.input.int(caption ?? label, opts) : ui.input.float(caption ?? label, opts);
    PropertyPanel.initInputValue(input, value, false);
    return this.dgRow(label, input);
  }

  private createToggle(label: string, value: boolean, onChange: (v: boolean) => void, inputTooltip?: string, caption?: string): HTMLElement {
    const report = this.changeReporter(value);
    const input = ui.input.bool(caption ?? label, {
      tooltipText: PropertyPanel.propTip(label, caption, inputTooltip),
      onValueChanged: (v) => {
        onChange(Boolean(v));
        report(Boolean(v));
      },
    });
    PropertyPanel.initInputValue(input, Boolean(value), false);
    return this.dgRow(label, input);
  }

  private createCombo(label: string, value: string, options: string[], onChange: (v: string) => void, inputTooltip?: string): HTMLElement {
    const report = this.changeReporter(value);
    const input = ui.input.choice(label, {
      items: options,
      tooltipText: PropertyPanel.propTip(label, undefined, inputTooltip),
      onValueChanged: (v) => {
        const s = String(v ?? '');
        onChange(s);
        report(s);
      },
    });
    PropertyPanel.initInputValue(input, value, false);
    return this.dgRow(label, input);
  }
}

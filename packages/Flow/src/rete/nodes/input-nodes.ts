/** Input nodes — emit `//input:` lines. Qualifiers live in `node.properties`; each node
 *  carries an inline value editor mirrored in the side panel (see `utils/input-values.ts`). */

import {ClassicPreset} from 'rete';
import * as DG from 'datagrok-api/dg';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {categoricalColor, CAT, dgTypeToSlotType, isStringListType} from '../../types/type-map';
import {InputValueControl} from './input-value-control';
import {isInputOptional, getParamDefault, getParamDescription} from '../../utils/dart-proxy-utils';
import {isLiteralChoiceList, isChoicesRefString, isChoiceQueryRef} from '../../utils/choice-refs';
import {propertyNameToFriendly} from '../../utils/naming';

const COLOR_INPUT = categoricalColor(CAT.green);

abstract class InputBase extends FlowNode {
  constructor(label: string, paramName: string, dgType: string, slotName = 'value', extraProps: Record<string, any> = {}) {
    super(label);
    this.dgNodeType = 'input';
    this.dgOutputType = dgType;
    this.properties = {paramName, defaultValue: '', ...extraProps};
    (this as unknown as {color: string}).color = COLOR_INPUT;
    this.addOutput(slotName, new ClassicPreset.Output(getSocket(dgType), slotName));
    // Built lazily on first render so subclass dgOutputType overrides (Blob) apply.
    this.addControl('value', new InputValueControl(this));
  }
}

export class TableInputNode extends InputBase {
  constructor() { super('Table Input', 'df', 'dataframe', 'table'); }
}

export class ColumnInputNode extends InputBase {
  constructor() {
    super('Column Input', 'col', 'column', 'column', {typeFilter: '', semTypeFilter: ''});
  }
}

export class ColumnListInputNode extends InputBase {
  constructor() {
    super('Column List Input', 'cols', 'column_list', 'columns', {typeFilter: '', semTypeFilter: ''});
  }
}

export class StringInputNode extends InputBase {
  constructor() {
    super('String Input', 'text', 'string', 'value',
      {nullable: false, caption: '', choices: '', semType: ''});
  }
}

/** A String Input pre-tagged `semType: Molecule` — the panel's value editor becomes
 *  Chem's compact molecule input, while the NODE body embeds a real inplace sketcher
 *  (`buildInlineSketcherEditor` in utils/input-values.ts); emits an ordinary string input line. */
export class MoleculeInputNode extends InputBase {
  constructor() {
    super('Sketcher Input', 'molecule', 'string', 'molecule',
      {nullable: false, caption: '', choices: '', semType: 'Molecule'});
  }
}

/** The macromolecule counterpart of {@link MoleculeInputNode} — `semType: Macromolecule`
 *  routes the value editor to Helm's. */
export class HelmInputNode extends InputBase {
  constructor() {
    super('Helm Input', 'sequence', 'string', 'sequence',
      {nullable: false, caption: '', choices: '', semType: 'Macromolecule'});
  }
}

export class NumberInputNode extends InputBase {
  constructor() {
    super('Number Input', 'value', 'double', 'value',
      {nullable: false, caption: '', min: '', max: '', showSlider: false});
    this.properties['defaultValue'] = 0;
  }
}

export class IntInputNode extends InputBase {
  constructor() {
    super('Int Input', 'n', 'int', 'value',
      {nullable: false, caption: '', min: '', max: '', showSlider: false});
    this.properties['defaultValue'] = 0;
  }
}

export class BooleanInputNode extends InputBase {
  constructor() {
    super('Boolean Input', 'flag', 'bool', 'value', {nullable: false, caption: ''});
    this.properties['defaultValue'] = false;
  }
}

export class DateTimeInputNode extends InputBase {
  constructor() {
    super('DateTime Input', 'date', 'datetime', 'value', {nullable: false, caption: ''});
  }
}

export class FileInputNode extends InputBase {
  constructor() {
    super('File Input', 'file', 'file', 'value', {nullable: false, caption: ''});
  }
}

export class MapInputNode extends InputBase {
  constructor() {
    super('Map Input', 'params', 'map', 'value', {nullable: false, caption: ''});
  }
}

export class DynamicInputNode extends InputBase {
  constructor() {
    super('Dynamic Input', 'value', 'dynamic', 'value', {nullable: false, caption: ''});
  }
}

export class StringListInputNode extends InputBase {
  constructor() { super('String List Input', 'items', 'string_list', 'value'); }
}

export class BlobInputNode extends InputBase {
  constructor() {
    super('Blob Input', 'data', 'byte_array', 'value', {nullable: false, caption: ''});
    // The DG annotation type is `blob`, even though the slot type maps to `byte_array`.
    this.dgOutputType = 'blob';
  }
}

export const PROPERTY_INPUT_TYPE = 'Inputs/Property Input';

/** The universal input: starts `dynamic` (connectable to anything) and MIMICS the
 *  function input it gets wired into — `adoptInputProperty` copies the target
 *  parameter's `DG.Property` (type, semType, choices, bounds, default) onto the
 *  node, and everything downstream (value editor, panel rows, header emission,
 *  sketcher/helm routing) follows from those properties like on any input node. */
export class PropertyInputNode extends InputBase {
  constructor() {
    super('Property Input', 'value', 'dynamic', 'value',
      {nullable: false, caption: '', propertyType: 'dynamic'});
  }
}

/** Qualifier property keys per mimicked type — must match the dedicated input
 *  nodes exactly: the panel's rows and the emitted header qualifiers both read
 *  key existence. Types not listed carry the plain {nullable, caption} pair. */
const QUALIFIERS_BY_TYPE: Record<string, readonly string[]> = {
  string: ['nullable', 'caption', 'choices', 'semType'],
  double: ['nullable', 'caption', 'min', 'max', 'showSlider'],
  int: ['nullable', 'caption', 'min', 'max', 'showSlider'],
  column: ['typeFilter', 'semTypeFilter'],
  column_list: ['typeFilter', 'semTypeFilter'],
  dataframe: [],
  string_list: [],
};
const ALL_QUALIFIER_KEYS = ['nullable', 'caption', 'choices', 'semType', 'min', 'max',
  'showSlider', 'typeFilter', 'semTypeFilter'] as const;

/** Re-derive a Property Input's shape (output type, socket, qualifier keys) from
 *  `properties['propertyType']`. Idempotent; must run after adoption, a manual
 *  panel type change, and every path that assigns properties post-construction
 *  (flow load, duplicate/paste). A no-op for every other node type. */
export function applyPropertyInputShape(node: FlowNode): void {
  if (node.dgTypeName !== PROPERTY_INPUT_TYPE) return;
  const t = String(node.properties['propertyType'] ?? '') || 'dynamic';
  node.dgOutputType = t;
  const out = node.outputs['value'];
  if (out) out.socket = getSocket(dgTypeToSlotType(t));
  const keep = QUALIFIERS_BY_TYPE[t] ?? ['nullable', 'caption'];
  for (const k of ALL_QUALIFIER_KEYS) {
    if (!keep.includes(k)) delete node.properties[k];
    else if (node.properties[k] === undefined)
      node.properties[k] = k === 'nullable' || k === 'showSlider' ? false : '';
  }
  // Rides along with `choices` only (the connection a query-reference resolves on).
  if (!('choices' in node.properties)) delete node.properties['choicesConnection'];
}

/** Mimic a function input: copy its `DG.Property` onto the node. The param name and
 *  title follow the mimicked parameter only while still auto-derived — a user rename
 *  survives re-adoption; a configured value is never overwritten by the declared default.
 *  `func` (the input's owner) is only needed for query-reference choices — a
 *  `query("…")` list resolves on the owning query's connection, captured here. */
export function adoptInputProperty(node: FlowNode, prop: DG.Property, func?: DG.Func): void {
  if (node.dgTypeName !== PROPERTY_INPUT_TYPE) return;
  const read = <T, >(f: () => T): T | undefined => {
    try {
      return f() ?? undefined;
    } catch {
      return undefined; // Dart proxy reads can throw
    }
  };
  let t = String(read(() => prop.propertyType) ?? 'dynamic');
  if (isStringListType(t)) t = 'string_list';
  else if (t === 'num') t = 'double'; // no `num` editor/resolver — mimic as double
  const prevType = String(node.properties['propertyType'] ?? '');
  node.properties['propertyType'] = t;
  applyPropertyInputShape(node);

  const name = String(read(() => prop.name) ?? 'value');
  const set = (key: string, v: unknown): void => {
    if (key in node.properties) node.properties[key] = v;
  };
  set('nullable', isInputOptional(prop));
  const caption = String(read(() => prop.caption) ?? '');
  set('caption', caption !== '' && caption !== name ? caption : '');
  const semType = String(read(() => prop.semType) ?? '');
  set('semType', semType);
  set('semTypeFilter', semType);
  set('typeFilter', String(read(() => prop.columnTypeFilter) ?? ''));
  const choices = read(() => prop.choices);
  delete node.properties['choicesConnection'];
  // The platform can hand the reference JSON-escaped (`uery(\"SELECT …\"`) — unescape
  // BEFORE detecting, or an escaped query reference reads as a literal list.
  const rawRef = Array.isArray(choices) && choices.length === 1 ?
    String(choices[0]).replaceAll('\\"', '"') : null;
  if (rawRef != null && isChoicesRefString(rawRef)) {
    // A REFERENCE (func call / query) — stored verbatim, resolved live by the value
    // editors; a query reference runs on the owning query's connection.
    const raw = rawRef;
    set('choices', raw);
    if ('choices' in node.properties && isChoiceQueryRef(raw)) {
      const connId = read(() => (func instanceof DG.DataQuery ? func.connection?.id : undefined));
      if (connId) node.properties['choicesConnection'] = String(connId);
    }
  } else {
    set('choices', Array.isArray(choices) && isLiteralChoiceList(choices) ?
      choices.map((c) => String(c)).join(',') : '');
  }
  const min = read(() => prop.min);
  const max = read(() => prop.max);
  set('min', typeof min === 'number' && isFinite(min) ? String(min) : '');
  set('max', typeof max === 'number' && isFinite(max) ? String(max) : '');
  set('showSlider', read(() => prop.showSlider) === true);

  // Keep a configured value only across a same-type re-adoption — a value typed for
  // the previous type is meaningless (a SMILES riding into an int would emit NaN).
  const cur = node.properties['defaultValue'];
  if (t !== prevType || cur === undefined || cur === null || String(cur) === '') {
    node.properties['defaultValue'] = (getParamDefault(prop) as string | number | boolean | undefined) ?? '';
    node.transientValue = undefined;
  }

  const prev = String(node.properties['mimicParam'] ?? '');
  const autoTitle = (p: string): string => `${propertyNameToFriendly(p)} Input`;
  if (String(node.properties['paramName'] ?? '') === (prev || 'value'))
    node.properties['paramName'] = name;
  if (node.label === (prev ? autoTitle(prev) : 'Property Input')) node.label = autoTitle(name);
  if (!node.description) node.description = getParamDescription(prop);
  node.properties['mimicParam'] = name;
}

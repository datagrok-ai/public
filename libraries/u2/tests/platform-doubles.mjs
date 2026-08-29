/* Getter-backed doubles of the platform entities u2 reads — the one place the headless suite fakes
   a `DG.Property`, a `DG.Func`, a `DG.FileInfo`, a `DG.DataFrame` or an entity row. On the real
   entity every field is a PROTOTYPE getter over the `dart` handle, so a spread or `Object.keys` of
   one yields nothing but `dart` — the property plain-object fakes lacked, which is how two
   defects survived 600 tests. Every double keeps its whole state on `dart` and answers it through
   getters; a test that needs to poke the state behind the handle goes through `dart` too. (The
   real `DG.DataFrame` also keeps `columns, rows, filter, temp, tags` as own fields, data-frame.ts:48-53;
   the double is deliberately stricter, since no u2 path reads those — and so do `DG.Widget` and
   `DG.Viewer`, whose doubles below keep root, subs, look and events behind the handle.)
   tests/dg-stub.mjs and tests/platform-stub.mjs serve these same classes to the modules under test
   as `datagrok-api/dg` and `datagrok-api/grok`; dg-stub also installs the kill-walk globals over
   {@link platform}. */

import {Control, dfBindings, signal} from 'datagrok-api/u2core';

/** Prototype getters over the handle — what every field of a real entity is. */
function getters(cls, ...keys) {
  for (const key of keys)
    Object.defineProperty(cls.prototype, key, {get() { return this.dart[key]; }});
}

/** Prototype accessors over the handle, for the fields the platform lets a caller assign. */
function fields(cls, ...keys) {
  for (const key of keys)
    Object.defineProperty(cls.prototype, key, {get() { return this.dart[key]; }, set(x) { this.dart[key] = x; }});
}

/** The rxjs shape u2 subscribes to, plus what a test needs to drive it: `fire` and `count`. */
export class Stream {
  constructor() { this.subs = new Set(); }

  subscribe(next) {
    this.subs.add(next);
    return {unsubscribe: () => this.subs.delete(next)};
  }

  fire(x) {
    for (const next of [...this.subs])
      next(x);
  }

  get count() { return this.subs.size; }
}

export class Entity {
  constructor(name, data = {}) { this.dart = {name, ...data}; }

  get nqName() {
    return this.dart.nqName ??
      (this.dart.namespace == null ? this.dart.name : `${this.dart.namespace}:${this.dart.name}`);
  }
}
getters(Entity, 'id', 'name', 'friendlyName');

export class User extends Entity {}
getters(User, 'email');

export class Package extends Entity {}

/** The function hierarchy the source backends discriminate on. `registry` is what `find` answers
 * from; `run` is what a prepared call executes. */
export class Func extends Entity {
  static registry = [];

  constructor(name, data = {}) {
    super(name, {package: null, inputs: [], outputs: [], options: {}, run: async () => ({}), ...data});
  }

  /** `meta` matches the `options` record a `meta.role:` annotation surfaces as. */
  static find(query = {}) {
    return Func.registry.filter((f) =>
      (query.name === undefined || f.name === query.name) &&
      (query.package === undefined || f.package?.name === query.package) &&
      Object.entries(query.meta ?? {}).every(([key, value]) => f.options?.[key] === value));
  }

  prepare(params) { return new FuncCall(this, params); }
}
getters(Func, 'package', 'description', 'inputs', 'outputs', 'options', 'tags');

export class DataQuery extends Func {}

export class Script extends Func {}

/** `limit` lives on the shared entity and is read while the call EXECUTES — the whole point of the
 * design-time cap. */
export class TableQuery extends DataQuery {
  get limit() { return this.dart.limit; }
  set limit(x) { this.dart.limit = x; }
}

/** What a function declares its inputs and outputs as, and what the widget host mints out of a u2
 * IProperty: `create`'s arguments plus the options record `fromOptions` refines it with. */
export class Property {
  constructor(name, propertyType, options = {}) { this.dart = {name, propertyType, ...options}; }

  static create(name, type, get, set, defaultValue) { return new Property(name, type, {get, set, defaultValue}); }

  fromOptions(options) {
    Object.assign(this.dart, options);
    return this;
  }

  get type() { return this.dart.propertyType; }

  /** Unset answers null, as the Dart field does (prop_gen/lib/property.dart:63). */
  get choices() { return this.dart.choices ?? null; }

  /** A `list` property's element type (`string`, `map`) — what `grok_Property_Get_PropertySubType` reads. */
  get propertySubType() { return this.dart.subType ?? null; }
}
getters(Property, 'name', 'propertyType', 'semType', 'description', 'caption', 'friendlyName',
  'nullable', 'defaultValue', 'category', 'userEditable', 'min', 'max', 'step', 'inputType', 'editor',
  'format', 'units', 'showSlider', 'showPlusMinus');
fields(Property, 'get', 'set');

/** A function's declared parameter (ddt FuncParam): its `options` is the LIVE map the platform
 * answers for a tag-less FuncParam — a write through it reads back. */
export class FuncParam extends Property {
  constructor(name, propertyType, options = {}) {
    super(name, propertyType, options);
    this.dart.options ??= {};
  }

  get options() { return this.dart.options; }

  get editor() { return this.dart.options['editor'] ?? null; }

  get isOptional() { return this.dart.isOptional === true; }

  /** The `info ?? field` union DG.Property answers (`grok_api.dart:900`) — annotation `columns`
   * or `type` keys. The parent-table link is NOT a wrapper getter: tests set
   * `dart.parentTableParamName` directly and u2 reads it through `grok_Property_Get`. */
  get columnTypeFilter() {
    return this.dart.options['columns'] ?? this.dart.options['type'] ?? null;
  }
}

export class FuncCallParam {
  constructor(property, call) {
    this.dart = {property, call, value: undefined, onChanged: new Stream()};
  }

  get name() { return this.dart.property.name; }

  get property() { return this.dart.property; }

  /** Optional-default substitution, as func_call_param.dart:78-88. */
  get value() {
    const v = this.dart.value;
    return v == null && this.dart.property.isOptional ? this.dart.property.defaultValue ?? null : v;
  }

  get onChanged() { return this.dart.onChanged; }
}

/** What `Func.prepare` answers: params over the declared inputs (a plain-record input is adopted
 * into a FuncParam), `setParamValue` the one write path, `call()` the back-compat run — it still
 * resolves to an object whose `.outputs` are the run result. */
export class FuncCall {
  constructor(func, params = {}) {
    this.dart = {func, outputs: null, params: (func?.inputs ?? []).map((p) =>
      new FuncCallParam(p instanceof Property ? p : new FuncParam(p.name, p.propertyType ?? p.type, p), this))};
    for (const [k, v] of Object.entries(params))
      this.setParamValue(k, v);
  }

  get func() { return this.dart.func; }

  get inputParams() {
    const params = this.dart.params;
    const map = Object.fromEntries(params.map((p) => [p.name, p]));
    map.values = () => params.slice();
    return map;
  }

  get inputs() { return Object.fromEntries(this.dart.params.map((p) => [p.name, p.value])); }

  /** Same-value suppression at the source (func_call_param.dart:194) — EXCEPT a List value: the
   * real setter always fires for List/ColFilterCall values, a same-reference write included. */
  setParamValue(name, value) {
    const p = this.dart.params.find((x) => x.name === name);
    if (p.dart.value === value && !Array.isArray(value))
      return;
    p.dart.value = value;
    p.dart.onChanged.fire(p);
  }

  /** The `options['choices']` command, answered in the propagate shape: an array becomes
   * stringified items with identity `values` (the interop stringification, P-W2-1b) and no
   * lookup; a provider answering `{values, lookup}` passes through. */
  async evalParamChoices(name) {
    const p = this._findParam(name);
    const {result, dependsOn} =
      await this._evalCommand(p.property.options?.['choices'], name, () => p.value);
    if (Array.isArray(result)) {
      const items = result.map(String);
      return {items, values: Object.fromEntries(items.map((k) => [k, k])), lookup: null, dependsOn};
    }
    return {items: Object.keys(result.values), values: result.values,
      lookup: result.lookup ?? null, dependsOn};
  }

  /** The `options['suggestions']` command over the TYPED text — the fixed contract: a single
   * provider input binds `text`, never the stored param value. An array answer carries no
   * tooltips; an `{items, tooltips}` answer passes through. */
  async evalParamSuggestions(name, text) {
    const p = this._findParam(name);
    const {result} = await this._evalCommand(p.property.options?.['suggestions'], name, () => text);
    return Array.isArray(result) ? {items: result.map(String), tooltips: {}} :
      {items: result.items.map(String), tooltips: result.tooltips ?? {}};
  }

  /** The `options['default']` command; sibling params never reach it in practice (P-W2-5). */
  async evalParamDefault(name) {
    const p = this._findParam(name);
    return (await this._evalCommand(p.property.options?.['default'], name, () => p.value)).result;
  }

  /** The `validators:` run behind `FuncCall_EvalParamValidators`: builtins first, else the
   * registry (a missing func rejects; one marked `{sync: false}` rejects with the Dart
   * must-be-sync shape); each runs SYNC over the param's CURRENT value — `true`/`null` omitted,
   * `false` the didn't-pass message, a string the message, a full `{message, isError, isHelper}`
   * object passed through. Delivery is never synchronous (the `_evalCommand` clock). */
  async evalParamValidators(name) {
    const p = this._findParam(name);
    const delay = this.dart.evalDelayMs;
    if (delay != null)
      await new Promise((resolve) => setTimeout(resolve, delay));
    else {
      await Promise.resolve();
      await Promise.resolve();
    }
    const results = [];
    for (const vname of p.property.options?.['validators'] ?? []) {
      let run = BUILTIN_VALIDATORS[vname];
      let label = vname;
      if (run === undefined) {
        const func = Func.find({name: vname})[0];
        if (func === undefined)
          throw new Error(`Validator "${vname}" not found`);
        if (func.dart.sync === false)
          throw new Error(`Validator '${vname}' is a Script (async). Validators must be sync - ` +
            'use a package.ts export or grok.functions.register().');
        label = func.friendlyName ?? vname;
        run = (v) => func.dart.run({[func.inputs[0]?.name ?? 'value']: v});
      }
      const r = run(p.value);
      if (r == null || r === true)
        continue;
      if (r === false)
        results.push({message: `Value didn't pass the ${label} check`, isError: true, isHelper: false});
      else if (typeof r === 'string')
        results.push({message: r, isError: true, isHelper: false});
      else
        results.push(r);
    }
    return results;
  }

  /** Resolves an options command — `Name` or `Name(args)`, the parenthesized args accepted and
   * ignored — against {@link Func.registry}: provider inputs bind to same-named host params (a
   * single provider input matching none binds `fallback` — the param's own value, or the typed
   * text for suggestions); `dependsOn` is the FOREIGN host params bound, never the param itself;
   * an unknown func name rejects. Delivery is never synchronous: two microtasks by default,
   * `dart.evalDelayMs` milliseconds when a test sets the knob. */
  async _evalCommand(command, paramName, fallback) {
    const delay = this.dart.evalDelayMs;
    if (delay != null)
      await new Promise((resolve) => setTimeout(resolve, delay));
    else {
      await Promise.resolve();
      await Promise.resolve();
    }
    const name = String(command ?? '').trim().replace(/\(.*\)$/s, '').trim();
    const func = Func.find({name})[0];
    if (func === undefined)
      throw new Error(`Function not found: ${command}`);
    const inputs = {};
    const dependsOn = [];
    let bound = 0;
    for (const input of func.inputs) {
      const host = this.dart.params.find((p) => p.name === input.name);
      if (host === undefined)
        continue;
      inputs[input.name] = host.value;
      bound++;
      if (input.name !== paramName)
        dependsOn.push(input.name);
    }
    if (bound === 0 && func.inputs.length === 1)
      inputs[func.inputs[0].name] = fallback();
    return {result: await func.dart.run(inputs), dependsOn};
  }

  _findParam(name) {
    const p = this.dart.params.find((x) => x.name === name);
    if (p === undefined)
      throw new Error(`Parameter not found: ${name}`);
    return p;
  }

  get outputs() { return this.dart.outputs; }

  async call() {
    this.dart.outputs = await this.dart.func.dart.run(this.inputs);
    return this;
  }
}

/** The Dart `namedValidators` builtins the doubles carry (vld:79-85) — the minimal `notEmpty`. */
const BUILTIN_VALIDATORS = {
  notEmpty: (v) => v != null && String(v) !== '' ? null : 'Value cannot be empty.',
};

/** Module-level knob for {@link scriptSyncStub} — `true` throws a generic failure, a string
 * throws that message. */
export const dart = {scriptSyncError: null};

/** A test-only micro-translation of the DOCUMENTED expression forms — the production path never
 * runs it (`grok_ScriptSync` is the platform's). NOT auto-installed: tests assign
 * `globalThis.grok_ScriptSync = scriptSyncStub` in setup and `delete` it in teardown. */
export function scriptSyncStub(expr, vars = {}) {
  if (dart.scriptSyncError)
    throw new Error(dart.scriptSyncError === true ? 'script failed' : String(dart.scriptSyncError));
  // == before !=: the != rewrite mints a `!==` whose tail `==` the other rule must not touch
  const js = String(expr)
    .replace(/==(?!=)/g, '===').replace(/!=(?!=)/g, '!==')
    .replace(/\band\b/g, '&&').replace(/\bor\b/g, '||').replace(/\bnot\b/g, '!');
  const names = Object.keys(vars);
  return new Function(...names, `'use strict'; return (${js});`)(...names.map((n) => vars[n]));
}

/** A file in a share carries the connection it lives on; one built out of local bytes does not. */
export class FileInfo extends Entity {
  constructor(fullPath, options = {}) {
    super(fullPath.substring(fullPath.lastIndexOf('/') + 1),
      {fullPath, data: null, connection: null, directory: false, ...options});
  }

  static fromBytes(path, data) { return new FileInfo(path, {data}); }

  /** Relative to the connection: what follows `Namespace:Share/`, nothing for the share itself; a
   * local path is its own. */
  get path() {
    const full = this.dart.fullPath;
    const colon = full.indexOf(':');
    const slash = full.indexOf('/', colon);
    return colon < 0 ? full : slash < 0 ? '' : full.substring(slash + 1);
  }

  get fileName() { return this.dart.name; }

  get extension() {
    const dot = this.fileName.lastIndexOf('.');
    return dot < 0 ? '' : this.fileName.substring(dot + 1);
  }

  get isFile() { return !this.dart.directory; }
  get isDirectory() { return this.dart.directory; }
}
getters(FileInfo, 'fullPath', 'data', 'connection');

/** A deployed js-api bundle passes the wrapped handle to the accessor (`table-info.ts:73`), so
 * `connection` throws for every FileInfo. */
export class UnreadableFileInfo extends FileInfo {
  get connection() { throw new RangeError('Invalid argument'); }
}

export class BitSet {
  constructor(length = 0) { this.dart = {length}; }
}
getters(BitSet, 'length');

// datetime counts as numerical, as in ddt (date_time_column.dart:16)
const NUMERICAL = new Set(['int', 'double', 'bigint', 'qnum', 'datetime']);

export class Column {
  constructor(name, type = 'string', semType = null) { this.dart = {name, type, semType, frame: null}; }

  get name() { return this.dart.name; }

  /** A rename fires the frame's `onColumnNameChanged` ALONE (`column.dart:68`), names in the args. */
  set name(x) {
    const oldName = this.dart.name;
    this.dart.name = x;
    this.dart.frame?.onColumnNameChanged.fire({args: {oldName, newName: x}});
  }

  get isNumerical() { return NUMERICAL.has(this.dart.type); }

  // string and bool are the categorical types (ddt string_column.dart / bool_column.dart)
  get isCategorical() { return this.dart.type === 'string' || this.dart.type === 'bool'; }

  /** Distinct values off the frame rows — read only under a maxCategories cap. */
  get categories() {
    const frame = this.dart.frame;
    return frame == null ? [] : [...new Set(frame.dart.rows.map((r) => r[this.dart.name]))];
  }
}
getters(Column, 'type', 'semType');

class ColumnList {
  constructor(frame, columns) {
    this.dart = {frame, columns: []};
    for (const c of columns)
      this._adopt(c instanceof Column ? c : new Column(c.name, c.type, c.semType));
  }

  names() { return this.dart.columns.map((c) => c.name); }
  toList() { return this.dart.columns.slice(); }
  byName(name) { return this.dart.columns.find((c) => c.name === name) ?? null; }
  contains(name) { return this.byName(name) !== null; }

  add(column) {
    this._adopt(column);
    this.dart.frame.onColumnsChanged.fire();
    return column;
  }

  /** An unknown name removes nothing and still notifies, as the platform does (column_list.dart:254). */
  remove(name) {
    const at = this.dart.columns.findIndex((c) => c.name === name);
    if (at >= 0)
      this.dart.columns.splice(at, 1);
    this.dart.frame.onColumnsChanged.fire();
  }

  _adopt(column) {
    column.dart.frame = this.dart.frame;
    this.dart.columns.push(column);
  }
}

const EVENTS = ['onCurrentRowChanged', 'onValuesChanged', 'onSelectionChanged', 'onFilterChanged',
  'onColumnsChanged', 'onColumnNameChanged'];

/** As much of a frame as u2 reads: cells, the current row, the column list and the events. The
 * rows are the records given, behind the handle (`dart.rows`); the handle also counts the reads,
 * the writes and the row-index writes, which is what the laziness and echo tests assert on. */
export class DataFrame {
  constructor(columns = [], rows = [], name = 'table') {
    this.dart = {name, rows: rows.map((r) => ({...r})), currentRowIdx: rows.length > 0 ? 0 : -1,
      selection: new BitSet(rows.length), filter: new BitSet(rows.length),
      reads: 0, writes: [], idxWrites: 0,
      events: Object.fromEntries(EVENTS.map((event) => [event, new Stream()]))};
    this.dart.columns = new ColumnList(this, columns);
  }

  static fromByteArray(bytes) {
    const frame = new DataFrame();
    frame.dart.bytes = bytes;
    return frame;
  }

  /** An empty frame of `rows` rows — what an unbound viewer starts over (VP-8). */
  static create(rows = 0) { return new DataFrame([], Array.from({length: rows}, () => ({}))); }

  get name() { return this.dart.name; }
  set name(x) { this.dart.name = x; }
  get rowCount() { return this.dart.rows.length; }
  get currentRowIdx() { return this.dart.currentRowIdx; }

  set currentRowIdx(v) {
    this.dart.idxWrites++;
    if (v === this.dart.currentRowIdx)
      return;
    this.dart.currentRowIdx = v;
    this.dart.events.onCurrentRowChanged.fire(v);
  }

  get(column, row) {
    this.dart.reads++;
    return this.dart.rows[row][column];
  }

  set(column, row, value) {
    this.dart.writes.push([column, row, value]);
    this.dart.rows[row][column] = value;
    this.dart.events.onValuesChanged.fire({column, row});
  }

  liveSubscriptions() { return EVENTS.reduce((n, event) => n + this.dart.events[event].count, 0); }
}
getters(DataFrame, 'columns', 'selection', 'filter');
for (const event of EVENTS)
  Object.defineProperty(DataFrame.prototype, event, {get() { return this.dart.events[event]; }});

/** The shell: the current object — which the platform rebuilds the property panel for on EVERY
 * assignment, so `dart.writes` is what a test counts — the balloons, and the open tables. The
 * tables are FRAMES on the handle (`dart.tables`); `tableNames` derives from them, and its
 * setter keeps platform-stub's `resetShell` name-list reset clearing the same backing store. */
export class Shell {
  constructor() {
    const dart = {o: null, writes: [], tables: [], t: null, tableAdded: new Stream(),
      tableRemoved: new Stream(), windows: {showContextPanel: false}};
    Object.defineProperty(dart, 'tableNames', {
      get: () => dart.tables.map((table) => table.name),
      set: (names) => dart.tables = dart.tables.filter((table) => names.includes(table.name)),
    });
    this.dart = dart;
  }

  get o() { return this.dart.o; }

  set o(x) {
    this.dart.o = x;
    this.dart.writes.push(x);
  }

  info() {}
  warning() {}
  error() {}

  tableByName(name) { return this.dart.tables.find((table) => table.name === name) ?? null; }

  /** The platform never opens two tables under one name; the second becomes `demog (2)`. */
  addTable(table) {
    let unique = table.name;
    for (let n = 2; this.dart.tableNames.includes(unique); n++)
      unique = `${table.name} (${n})`;
    table.name = unique;
    this.dart.tables = [...this.dart.tables, table];
    this.dart.tableAdded.fire();
    return table;
  }

  /** Takes the frame or its name; the current table never answers a closed frame. */
  closeTable(x) {
    const name = typeof x === 'string' ? x : x?.name;
    this.dart.tables = this.dart.tables.filter((table) => table.name !== name);
    if (this.dart.t !== null && this.dart.t.name === name)
      this.dart.t = null;
    this.dart.tableRemoved.fire();
  }
}
getters(Shell, 'windows', 'tableNames', 'tables');
fields(Shell, 't');

/** What the kill-walk globals of tests/dg-stub.mjs work over: the elements killed, the cleanups
 * registered and not yet run, and the widgets the platform knows of — every viewer built, every
 * widget that went through `toDart()`. */
export const platform = {
  kills: [], cleanups: [], widgets: new Set(),
  reset() {
    this.kills.length = 0;
    this.cleanups.length = 0;
    this.widgets.clear();
  },
};

/** An event a descriptor declares: `name` is the id `onEvent()` takes, `eventName` a label. */
export class EventType {
  constructor(dart) { this.dart = dart; }
}
getters(EventType, 'name', 'eventName', 'description');

/** The bag over a widget's declared properties (js-api widgets/base.ts:76): a key read or written
 * goes through the property's `get`/`set` with `target` as the receiver — the widget itself, or a
 * viewer's LOOK. A name nothing declares throws, as the platform's does. */
export class ObjectPropertyBag {
  constructor(source, target = source) {
    this.source = source;
    const own = ['source', 'getProperties', 'getProperty', 'hasProperty', 'get', 'set', 'setAll'];
    return new Proxy(this, {
      ownKeys: () => this.getProperties().map((p) => p.name),
      has: (_, name) => this.hasProperty(name),
      getOwnPropertyDescriptor: () => ({enumerable: true, configurable: true}),
      get: (bag, name) => typeof name !== 'string' || own.includes(name) ? bag[name] : bag.getProperty(name).get(target),
      set: (bag, name, value) => {
        bag.getProperty(name).set(target, value);
        return true;
      },
    });
  }

  getProperties() { return this.source.getProperties(); }

  getProperty(name) {
    const property = this.getProperties().find((p) => p.name === name);
    if (property === undefined)
      throw `Property not found: ${name}`;
    return property;
  }

  hasProperty(name) { return this.getProperties().some((p) => p.name === name); }

  /** Over `source`, never the look — so on a viewer these throw, as the platform's do. */
  get(name) { return this.getProperty(name).get(this.source); }
  set(name, value) { this.getProperty(name).set(this.source, value); }

  setAll(params) {
    for (const [name, value] of Object.entries(params))
      this.set(name, value);
  }
}

/** The table step every viewer answers first (js-api viewer.ts keeps it module-private). */
const TABLE_STEP = Object.freeze({name: 'table', type: 'dataframe', walkable: true});

/** The base widget contract (js-api widgets/base.ts:265): a real u2 `Control` by inheritance, as
 * the platform's is — so `_u2` and `_root` are own fields, as on the real Widget (`scope` is the
 * lazily-minting accessor, never an own field) — with the platform state behind the handle: the
 * root, the subscriptions `detach` cancels, the properties `addProperty` declares and the bag over
 * them. Ambient adoption is off; the lifecycle is wired at the first scope mint, as in js-api. */
export class Widget extends Control {
  constructor(root) {
    super(root);
    this.dart = {root, subs: [], temp: {}, _properties: [], _functions: [], isDetached: false};
  }

  static fromRoot(root) { return new Widget(root); }

  /** Registers a cleanup to run when `element` is killed (base.ts:310). */
  static registerCleanup(element, cleanup) { globalThis.grok_Widget_RegisterCleanup(element, cleanup); }

  get type() { return 'Unknown'; }
  get props() { return this.dart.props ??= new ObjectPropertyBag(this); }
  get ambientAdoption() { return false; }
  get propertyTier() { return true; }

  /** Registers the widget with the platform — what the kill-walk finds (`grok_Widget_Wrap`). */
  toDart() {
    platform.widgets.add(this);
    return this.dart;
  }

  sub(subscription) { this.subs.push(subscription); }
  getProperties() { return [...super.getProperties(), ...this._properties]; }
  getFunctions() { return this._functions; }

  /** The virtual notify site: a subclass override that never calls super is covered by the
   * setter site in {@link addProperty}, exactly as in js-api. */
  onPropertyChanged(property) { this._notifyPropertyChange(property); }

  _notifyPropertyChange(p) {
    for (const next of this._u2.propertyChangeListeners ?? [])
      next(p?.name ?? null);
  }

  _ownPropertyChanges() {
    return {
      subscribe: (next) => {
        const listeners = this._u2.propertyChangeListeners ??= new Set();
        listeners.add(next);
        return {unsubscribe: () => listeners.delete(next)};
      },
    };
  }

  propertyChanges() { return this._ownPropertyChanges(); }

  /** Declares a property over the field of that name; its setter notifies the tier, then
   * {@link onPropertyChanged} — the two-site seam of js-api base.ts. */
  addProperty(propertyName, propertyType, defaultValue = null, options = null) {
    const fieldName = options?.fieldName ?? propertyName;
    const p = Property.create(propertyName, propertyType, () => this[fieldName], null, defaultValue);
    p.set = (_, x) => {
      this[fieldName] = x;
      this._notifyPropertyChange(p);
      this.onPropertyChanged(p);
    };
    for (const [key, value] of Object.entries(options ?? {}))
      if (key !== 'fieldName')
        p.dart[key] = value;
    this._properties.push(p);
    return p.defaultValue;
  }

  onEvent(_eventId = null) { return {subscribe: () => ({unsubscribe() {}})}; }

  getWidgetStatus() {
    return {parts: {}, hitAreas: {}, shortcuts: {}, events: [], description: null, error: null};
  }

  detach() {
    this._u2.detaching = true;
    for (const s of this.subs)
      s.unsubscribe();
    this.isDetached = true;
    this._u2.scope?.dispose();
  }

  _scopeMinted() { this._wireLifecycle(); }

  _wireLifecycle() { this._wireJsLifecycle(); }

  _wireJsLifecycle() {
    this.scope.own(() => {
      if (!this.isDetached && !this._u2.detaching)
        this.detach();
    });
  }

  _wireDartLifecycle() {
    const root = this.root;
    globalThis.grok_Widget_RegisterCleanup?.(root, () => {
      this._u2.platformKilled = true;
      this.scope.dispose();
    });
    this.scope.own(() => {
      if (!this._u2.platformKilled)
        globalThis.grok_Widget_Kill?.(root);
    });
  }
}
fields(Widget, 'subs', 'temp', '_properties', '_functions', 'isDetached');

/** A widget implemented in Dart (widgets/base.ts): its properties' `get`/`set` take the handle,
 * and its lifecycle is the platform's. */
export class DartWidget extends Widget {
  constructor(dart) {
    super(dart.root ?? document.createElement('div'));
    this.dart = {type: 'DartWidget', properties: [], functions: [], ...this.dart, ...dart};
  }

  get type() { return this.dart.type; }
  get props() { return this.dart.props ??= new ObjectPropertyBag(this, this.dart); }
  get propertyTarget() { return this.dart; }
  _wireLifecycle() { this._wireDartLifecycle(); }
  getProperties() { return this.dart.properties; }
  getFunctions() { return this.dart.functions; }
}

/** The core property grid (js-api widgets/data-widgets.ts:177): every update is recorded on the
 * handle — `update` here, `grok_PropertyGrid_Update` with the table when a caller goes to the global. */
export class PropertyGrid extends DartWidget {
  constructor() {
    super({root: document.createElement('div'), updates: []});
    this.root.classList.add('d4-property-grid-widget');
  }

  update(src, props) { this.dart.updates.push({src, props, table: undefined}); }
}

/** A platform view over a root (js-api view.ts): what `appView` builds on. */
export class View {
  constructor(root) { this.dart = {root, name: '', ribbonPanels: [], statusBarPanels: [], toolbox: null}; }

  static fromRoot(root) { return new View(root); }

  setRibbonPanels(panels) { this.dart.ribbonPanels = panels; }
}
getters(View, 'root', 'ribbonPanels');
fields(View, 'name', 'path', 'toolbox', 'statusBarPanels');

/** A look's owner — the viewer whose `onPropertyValueChanged` a property write fires — and the
 * proof a receiver IS a look: any other receiver is the platform's NoSuchMethodError. */
const LOOKS = new WeakMap();

/** What the platform knows about a viewer type without instantiating it (viewer.ts:32). The
 * registry is the test's to fill, as `Func.registry` is. A descriptor property's `get`/`set`,
 * unless given, are defined over the LOOK; a write fires the owning viewer's
 * `onPropertyValueChanged` with the property — for a same value too. */
export class WidgetDescriptor {
  static registry = [];

  constructor(name, properties = [], data = {}) {
    this.dart = {name, description: '', synonyms: [], events: [], ...data, properties};
    for (const p of properties)
      WidgetDescriptor._overLook(p);
  }

  static getDescriptors() { return WidgetDescriptor.registry.slice(); }
  static getByName(name) { return WidgetDescriptor.registry.find((d) => d.name === name) ?? null; }

  get events() { return this.dart.events.map((e) => new EventType(e)); }

  static _overLook(p) {
    const name = p.name;
    const look = (target) => {
      if (!LOOKS.has(target))
        throw new Error(`NoSuchMethodError: method not found: '${name}' Receiver: Instance of '${target?.constructor?.name ?? target}'`);
      return target;
    };
    if (!('get' in p.dart))
      p.dart.get = (target) => look(target)[name];
    if (!('set' in p.dart))
      p.dart.set = (target, value) => {
        look(target)[name] = value;
        LOOKS.get(target).onPropertyValueChanged.fire({args: {property: p}});
      };
  }
}
getters(WidgetDescriptor, 'name', 'description', 'synonyms', 'properties');

export class ViewerMetaHelper {
  constructor(viewer) { this.dart = {viewer, formulaLines: {}, annotationRegions: {}}; }
}
getters(ViewerMetaHelper, 'formulaLines', 'annotationRegions');

/** A platform viewer (viewer.ts:82): the handle carries the type, the LOOK its properties are
 * defined over, the frame, one event stream per id and the descriptor it was built from. The
 * V-4 overrides are the real ones: `propertyTarget` is the look, `propertyChanges` the platform
 * property event, `bindStep('table')` the frame's binding surface, and the lifecycle is Dart's —
 * disposing the scope kills the viewer through the kill-walk, and vice versa. */
export class Viewer extends Widget {
  constructor(dart, root) {
    super(root ?? dart.root ?? Viewer._root());
    this.initDartObject(dart);
  }

  static _root() {
    const root = document.createElement('div');
    root.setAttribute('data-widget', 'true');
    return root;
  }

  initDartObject(dart) {
    this.dart = {look: {}, events: {}, detached: new Stream(), propertyReads: 0, ...this.dart, ...dart};
    LOOKS.set(this.dart.look, this);
    platform.widgets.add(this);
  }

  /** The viewer of `viewerType` from its registered descriptor: the look at the descriptor's
   * defaults, `options` applied over it. */
  static fromType(viewerType, table, options = null) {
    const descriptor = WidgetDescriptor.getByName(viewerType);
    if (descriptor === null)
      throw new Error(`Unknown viewer type: ${viewerType}`);
    const viewer = new (VIEWER_CLASSES[viewerType] ?? Viewer)({type: viewerType, descriptor, dataFrame: table});
    for (const p of descriptor.properties)
      viewer.dart.look[p.name] = p.defaultValue ?? null;
    viewer.setOptions(options ?? {});
    return viewer;
  }

  static grid(t, options) { return Viewer.fromType('Grid', t, options); }
  static filters(t, options) { return Viewer.fromType('Filters', t, options); }
  static form(t, options) { return Viewer.fromType('Form', t, options); }
  static scatterPlot(t, options) { return Viewer.fromType('Scatter plot', t, options); }

  get type() { return this.dart.type; }
  get root() { return this.dart.root; }
  get descriptor() { return this.dart.descriptor; }
  get table() { return this.dart.dataFrame; }
  get dataFrame() { return this.dart.dataFrame; }

  set dataFrame(t) {
    this.dart.dataFrame = t;
    this.onDataFrameChanged.fire();
  }

  /** The bag over the LOOK (viewer.ts:154) — why `viewer.props.x` works where `prop.get(viewer.dart)` throws. */
  get props() { return this.dart.props ??= new ObjectPropertyBag(this, this.dart.look); }

  /** A viewer's property lives on its look, which the `props` bag is over (viewer.ts:219). */
  get propertyTarget() { return this.dart.look; }

  propertyChanges() {
    return {
      subscribe: (next) => {
        const s = this.onPropertyValueChanged.subscribe((e) => next(e?.args?.property?.name ?? null));
        return {unsubscribe: () => s.unsubscribe()};
      },
    };
  }

  _wireLifecycle() {
    if (this instanceof JsViewer)
      this._wireJsLifecycle();
    else
      this._wireDartLifecycle();
  }

  /** The walkable `table` step (viewer.ts): the frame's binding surface, repointing with it. */
  bindStep(name) {
    if (name !== 'table')
      return super.bindStep(name);
    const u2 = this._u2;
    if (u2.table === undefined) {
      const frame = signal(this.dataFrame ?? undefined);
      const s = this.onDataFrameChanged.subscribe(() => frame.value = this.dataFrame ?? undefined);
      this.scope.own(() => s.unsubscribe());
      u2.table = dfBindings(frame, this.scope);
    }
    return u2.table;
  }

  bindProps() { return [TABLE_STEP, ...super.bindProps()]; }

  /** The descriptor's properties, a fresh list per call as the platform answers; counted on the handle. */
  getProperties() {
    this.dart.propertyReads++;
    return this.dart.descriptor.properties.slice();
  }

  getFunctions() { return []; }

  /** Declared keys go through their property, one event each; `type` and unknown keys are ignored. */
  setOptions(map) {
    for (const [name, value] of Object.entries(map))
      if (name !== 'type' && this.props.hasProperty(name))
        this.props[name] = value;
  }

  getOptions(includeDefaults = false) {
    const defaultOf = (name) => this.descriptor.properties.find((p) => p.name === name)?.defaultValue ?? null;
    const look = Object.fromEntries(Object.entries(this.dart.look)
      .filter(([name, value]) => includeDefaults || JSON.stringify(value) !== JSON.stringify(defaultOf(name))));
    return {type: this.type, look};
  }

  onEvent(eventId) { return this.dart.events[eventId] ??= new Stream(); }
  get onDataEvent() { return this.onEvent('d4-data-event'); }
  get onPropertyValueChanged() { return this.onEvent('d4-property-value-changed'); }
  get onDataFrameChanged() { return this.onEvent('d4-data-frame-changed'); }

  /** The Dart detach — the kill-walk's; JS `detach()` never reaches it. */
  get onDetached() { return this.dart.detached; }

  /** Getter-only, as the platform's: `viewer.meta = x` throws. */
  get meta() { return this.dart.meta ??= new ViewerMetaHelper(this); }

  liveSubscriptions() {
    return Object.values(this.dart.events).reduce((n, stream) => n + stream.count, this.dart.detached.count);
  }
}

export class Grid extends Viewer {
  static create(table) { return Viewer.fromType('Grid', table); }
}

export class FilterGroup extends Viewer {
  /** The one write that rebuilds the panel (P14): the entry of the same column and type is
   * replaced in place, a new one appended; the `filters` property fires. */
  updateOrAdd(state, _requestFilter) {
    const filters = [...(this.dart.look.filters ?? [])];
    const at = filters.findIndex((f) => f.column === state.column && f.type === state.type);
    if (at < 0)
      filters.push(state);
    else
      filters[at] = state;
    this.props.filters = filters;
  }
}

/** A viewer implemented in JS (viewer.ts:478) — JS-owned: its properties are the ones
 * `addProperty` declared, over the viewer itself, and `detach()` is the real thing — it also
 * disposes the scope and deliberately never unregisters. */
export class JsViewer extends Viewer {
  constructor() {
    super({type: 'JsViewer', descriptor: null, dataFrame: null}, document.createElement('div'));
  }

  get root() { return this.dart.root; }
  set root(r) { this.dart.root = r; }
  get props() { return this.dart.props ??= new ObjectPropertyBag(this); }
  get propertyTarget() { return this; }
  propertyChanges() { return this._ownPropertyChanges(); }
  getProperties() { return this._properties; }

  detach() {
    this._u2.detaching = true;
    for (const s of this.subs)
      s.unsubscribe();
    this.isDetached = true;
    this._u2.scope?.dispose();
  }
}

/** The js-api class a type wraps into: `Viewer.grid` and `heatMap` answer a `Grid`, `filters` a `FilterGroup`. */
const VIEWER_CLASSES = {'Grid': Grid, 'Heat map': Grid, 'Filters': FilterGroup};

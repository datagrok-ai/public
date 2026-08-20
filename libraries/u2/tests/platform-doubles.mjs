/* Getter-backed doubles of the platform entities u2 reads — the one place the headless suite fakes
   a `DG.Property`, a `DG.Func`, a `DG.FileInfo`, a `DG.DataFrame` or an entity row. On the real
   entity every field is a PROTOTYPE getter over the `dart` handle, so a spread or `Object.keys` of
   one yields nothing but `dart` — the property plain-object fakes lacked, which is how two P4
   defects survived 600 tests. Every double keeps its whole state on `dart` and answers it through
   getters; a test that needs to poke the state behind the handle goes through `dart` too. (The
   real `DG.DataFrame` also keeps `columns, rows, filter, temp, tags` as own fields, data-frame.ts:48-53;
   the double is deliberately stricter, since no u2 path reads those.)
   tests/dg-stub.mjs and tests/platform-stub.mjs serve these same classes to the modules under test
   as `datagrok-api/dg` and `datagrok-api/grok`. */

/** Prototype getters over the handle — what every field of a real entity is. */
function getters(cls, ...keys) {
  for (const key of keys)
    Object.defineProperty(cls.prototype, key, {get() { return this.dart[key]; }});
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

  prepare(params) { return {call: async () => ({outputs: await this.dart.run(params)})}; }
}
getters(Func, 'package', 'description', 'inputs', 'outputs', 'options');

export class DataQuery extends Func {}

export class Script extends Func {}

/** `limit` lives on the shared entity and is read while the call EXECUTES — the whole point of the
 * design-time cap. */
export class TableQuery extends DataQuery {
  get limit() { return this.dart.limit; }
  set limit(x) { this.dart.limit = x; }
}

/** What a function declares its inputs and outputs as, and what the widget host mints out of a u2
 * PropertyLike: `create`'s arguments plus the options record `fromOptions` refines it with. */
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
}
getters(Property, 'name', 'propertyType', 'semType', 'description', 'caption', 'friendlyName',
  'nullable', 'defaultValue', 'category', 'get', 'set');

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
 * assignment, so `dart.writes` is what a test counts — the balloons, and the open tables. */
export class Shell {
  constructor() {
    this.dart = {o: null, writes: [], tableNames: [], tableAdded: new Stream(),
      tableRemoved: new Stream(), windows: {showContextPanel: false}};
  }

  get o() { return this.dart.o; }

  set o(x) {
    this.dart.o = x;
    this.dart.writes.push(x);
  }

  info() {}
  warning() {}
  error() {}

  /** The platform never opens two tables under one name; the second becomes `demog (2)`. */
  addTable(table) {
    let unique = table.name;
    for (let n = 2; this.dart.tableNames.includes(unique); n++)
      unique = `${table.name} (${n})`;
    table.name = unique;
    this.dart.tableNames = [...this.dart.tableNames, unique];
    this.dart.tableAdded.fire();
    return table;
  }

  closeTable(name) {
    this.dart.tableNames = this.dart.tableNames.filter((n) => n !== name);
    this.dart.tableRemoved.fire();
  }
}
getters(Shell, 'windows', 'tableNames');

/* Module hooks that serve minimal `datagrok-api/grok` and `datagrok-api/dg` modules to node, so
   the platform-bound pickers and column selectors (src/dg/pickers.ts, src/dg/columns.ts) load
   headlessly. Same idea as tests/dg-stub.mjs, which serves the input bridge's slice of the API;
   this one serves the shell (open tables + their events) and the enums the column modules read.

   Nothing here fakes a DataFrame — tests build their own, since a picker only ever asks a table
   for `columns` and `onColumnsChanged`. Register with
   `register('./platform-stub.mjs', import.meta.url)` before importing the modules under test, then
   import 'datagrok-api/grok' to drive the shell (`shell.tableNames`, `emitTableAdded()`). */

const GROK = `
export const shell = {
  tableNames: [],
  error: () => {},
  /** Adds a table to the workspace and returns it, under the name it ended up with. */
  addTable(table) {
    table.name = uniqueName(table.name ?? 'table');
    openTable(table.name);
    return table;
  },
};

export const data = {
  /** A DataFrame as a picker sees it — a name is all it ever reads. */
  parseCsv: (csv) => ({name: 'table', csv}),
};

const listeners = {added: [], removed: []};

function stream(key) {
  return {
    subscribe(fn) {
      listeners[key].push(fn);
      return {unsubscribe() {
        const i = listeners[key].indexOf(fn);
        if (i >= 0)
          listeners[key].splice(i, 1);
      }};
    },
  };
}

export const events = {onTableAdded: stream('added'), onTableRemoved: stream('removed')};

/** The platform never opens two tables under one name; the second becomes \`demog (2)\`. */
function uniqueName(name) {
  let unique = name;
  for (let n = 2; shell.tableNames.includes(unique); n++)
    unique = \`\${name} (\${n})\`;
  return unique;
}

/** Opens or closes tables the way the shell would: the name list first, then the event. */
export function openTable(name) {
  shell.tableNames = [...shell.tableNames, name];
  for (const fn of listeners.added.slice())
    fn();
}

export function closeTable(name) {
  shell.tableNames = shell.tableNames.filter((n) => n !== name);
  for (const fn of listeners.removed.slice())
    fn();
}

export function resetShell() {
  shell.tableNames = [];
  shell.error = () => {};
  data.parseCsv = (csv) => ({name: 'table', csv});
  listeners.added.length = 0;
  listeners.removed.length = 0;
}

export function liveSubscriptions() {
  return listeners.added.length + listeners.removed.length;
}
`;

/* The enum members the column modules read, copied from js-api's const.ts, plus the one DataFrame
   entry point the table import uses. */
const DG = `
/** Only what a picker reads off a table it just opened. \`fromByteArray\` takes the d42 blob the
 * platform writes (\`data-frame.ts:75\`), so the bytes are carried, not parsed. */
export class DataFrame {
  constructor(name, bytes) {
    this.name = name;
    this.bytes = bytes;
  }

  static fromByteArray(bytes) { return new DataFrame('table', bytes); }
}

export const COLUMN_TYPE = {
  STRING: 'string', INT: 'int', FLOAT: 'double', BOOL: 'bool', BYTE_ARRAY: 'byte_array',
  DATE_TIME: 'datetime', BIG_INT: 'bigint', QNUM: 'qnum', DATA_FRAME: 'dataframe',
  OBJECT: 'object',
};

export const AGG = {
  KEY: 'key', PIVOT: 'pivot', FIRST: 'first', TOTAL_COUNT: 'count', VALUE_COUNT: 'values',
  UNIQUE_COUNT: 'unique', MISSING_VALUE_COUNT: 'nulls', MIN: 'min', MAX: 'max', SUM: 'sum',
  MED: 'med', AVG: 'avg', STDEV: 'stdev', VARIANCE: 'variance', SKEW: 'skew', KURT: 'kurt',
  Q1: 'q1', Q2: 'q2', Q3: 'q3', SELECTED_ROWS_COUNT: '#selected',
};
`;

const URLS = {'datagrok-api/grok': 'u2-test:datagrok-api-grok', 'datagrok-api/dg': 'u2-test:datagrok-api-dg'};
const SOURCES = {'u2-test:datagrok-api-grok': GROK, 'u2-test:datagrok-api-dg': DG};

export async function resolve(specifier, context, next) {
  const url = URLS[specifier];
  if (url)
    return {url, format: 'module', shortCircuit: true};
  return next(specifier, context);
}

export async function load(url, context, next) {
  const source = SOURCES[url];
  if (source)
    return {source, format: 'module', shortCircuit: true};
  return next(url, context);
}

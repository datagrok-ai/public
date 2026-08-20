/* Module hooks that serve minimal `datagrok-api/grok` and `datagrok-api/dg` modules to node, so
   the platform-bound pickers and column selectors (src/dg/pickers.ts, src/dg/columns.ts) load
   headlessly. Same idea as tests/dg-stub.mjs, which serves the input bridge's slice of the API;
   this one serves the shell (open tables + their events) and the enums the column modules read.

   The shell and the frames are the getter-backed doubles of tests/platform-doubles.mjs. Register
   with `register('./platform-stub.mjs', import.meta.url)` before importing the modules under test,
   then import 'datagrok-api/grok' to drive the shell (`openTable()`, `closeTable()`). */

const DOUBLES = new URL('./platform-doubles.mjs', import.meta.url).href;

const GROK = `
import {DataFrame, Shell} from '${DOUBLES}';

export const shell = new Shell();

/** A DataFrame as a picker sees it — a name is all it ever reads. */
export const data = {parseCsv: () => new DataFrame()};

export const events = {onTableAdded: shell.dart.tableAdded, onTableRemoved: shell.dart.tableRemoved};

/** Opens or closes tables the way the shell would: the name list first, then the event. */
export const openTable = (name) => shell.addTable(new DataFrame([], [], name));

export const closeTable = (name) => shell.closeTable(name);

export function resetShell() {
  shell.dart.tableNames = [];
  delete shell.error;
  data.parseCsv = () => new DataFrame();
  events.onTableAdded.subs.clear();
  events.onTableRemoved.subs.clear();
}

export function liveSubscriptions() {
  return events.onTableAdded.count + events.onTableRemoved.count;
}
`;

/* The enum members the column modules read, copied from js-api's const.ts, plus the DataFrame
   the table import builds. */
const DG = `
export {Column, DataFrame} from '${DOUBLES}';

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

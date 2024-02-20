import * as DG from 'datagrok-api/dg';

// custom checkers interfaces API for adding custom checks
export interface CheckerError {
  path?: string[],
  msg: string,
}

export type NestedChecker = (path: string[], actual: any, expected: any) => void
export type ExpectPredicate = (actual: any, expected: any) => boolean;
export type ExpectChecker =
  (actual: any, expected: any, state: Readonly<ExpectRunnerState>, checkNested: NestedChecker) => CheckerError[] | void;

export interface CheckerItem {
  name: string;
  predicate: ExpectPredicate;
  checker: ExpectChecker;
}


// default checkers sequence, order matters!
export const defaulCheckersSeq: CheckerItem[] = [
  {name: 'Null', predicate: nullPredicate, checker: nullChecker},
  {name: 'Boolean', predicate: boolPredicate, checker: eqChecker},
  {name: 'String', predicate: stringPredicate, checker: eqChecker},
  {name: 'Integer', predicate: integerPredicate, checker: eqChecker},
  {name: 'Float', predicate: floatPredicate, checker: floatChecker},
  {name: 'Map', predicate: mapPredicate, checker: mapChecker},
  {name: 'Set', predicate: setPredicate, checker: setChecker},
  {name: 'ArrayBuffer', predicate: arrayBufferPredicate, checker: arrayBufferChecker},
  {name: 'Array', predicate: arrayPredicate, checker: arrayChecker},
  {name: 'DataFrame', predicate: dataframePredicate, checker: dataframeChecker},
  {name: 'Object', predicate: objectPredicate, checker: objectChecker},
];


// custom error class, just in case
export class ExpectError extends Error {}

// expectDeepEqual options
export interface ExpectDeepEqualOptions {
  floatTolerance?: number,
  maxErrorsReport?: number,
  maxDepth?: number,
  prefix?: string,
  checkersSeq?: CheckerItem[],
}

// public API
export function expectDeepEqual(actual: any, expected: any, options?: ExpectDeepEqualOptions): void {
  const runner = new ExpectRunner(options);
  runner.run(actual, expected);
}

// implementation details
export interface ExpectRunnerState {
  options: Required<ExpectDeepEqualOptions>;
  errors: CheckerError[];
  currentPath: string[];
}


function nullPredicate(actual: any, expected: any) {
  return actual == null && expected == null;
}

function nullChecker() {
  return;
}


function integerPredicate(actual: any, expected: any) {
  return (Number.isInteger(actual) && Number.isInteger(expected));
}

function stringPredicate(actual: any, expected: any) {
  return (typeof actual === 'string' && typeof expected === 'string');
}

function boolPredicate(actual: any, expected: any) {
  return (typeof actual === 'boolean' && typeof expected === 'boolean');
}

function eqChecker(actual: any, expected: any) {
  const areEqual = actual === expected;
  if (!areEqual) {
    const msg = `Expected ${expected}, got ${actual}`;
    return [{msg}];
  }
}


function floatPredicate(actual: any, expected: any) {
  return (typeof actual === 'number' && typeof expected === 'number');
}

function floatChecker(actual: any, expected: any, state: Readonly<ExpectRunnerState>): CheckerError[] | undefined {
  if ((actual === Number.POSITIVE_INFINITY && expected === Number.POSITIVE_INFINITY) ||
    (actual === Number.NEGATIVE_INFINITY && expected === Number.NEGATIVE_INFINITY) ||
    (isNaN(actual) && isNaN(expected)))
    return;
  const areEqual = Math.abs(actual - expected) < state.options.floatTolerance;
  if (!areEqual) {
    const msg = `Expected ${expected}, got ${actual} (tolerance = ${state.options.floatTolerance})`;
    return [{msg}];
  }
}


function mapPredicate(actual: any, expected: any) {
  return (actual instanceof Map && expected instanceof Map);
}

function mapChecker(actual: Map<any, any>, expected: Map<any, any>,
  _state: Readonly<ExpectRunnerState>, checkDeep: NestedChecker
): CheckerError[] {
  const errors: CheckerError[] = [];
  if (actual.size !== expected.size) {
    const err = {msg: `Maps are of different size: actual map size is ${actual.size} ` +
      `and expected map size is ${expected.size}`};
    errors.push(err);
  }
  for (const k of expected.keys()) {
    const aval = actual.get(k);
    const exval = expected.get(k);
    checkDeep([String(k)], aval, exval)
  }
  return errors;
}


function setPredicate(actual: any, expected: any) {
  return (actual instanceof Set && expected instanceof Set);
}

function setChecker(actual: Set<any>, expected: Set<any>,
  _state: Readonly<ExpectRunnerState>, checkDeep: NestedChecker
): CheckerError[] {
  const errors: CheckerError[] = [];
  if (actual.size !== expected.size) {
    const err = {msg: `Sets are of different size: actual set size is ${actual.size} ` +
      `and expected set size is ${expected.size}`};
    errors.push(err);
  }
  for (const exval of expected.values()) {
    const aval = actual.has(exval);
    checkDeep([String(exval)], aval, true);
  }
  return errors;
}


function arrayPredicate(actual: any, expected: any) {
  return (actual instanceof Array && expected instanceof Array);
}

function arrayChecker(actual: any[], expected: any[],
  _state: Readonly<ExpectRunnerState>, checkDeep: NestedChecker): CheckerError[] {
  const errors: CheckerError[] = [];
  if (actual.length !== expected.length) {
    const err = {msg: `Arrays are of different length: actual array length is ${actual.length} ` +
      `and expected array length is ${expected.length}`};
    errors.push(err);
  }
  for (let i = 0; i < actual.length; i++)
    checkDeep([String(i)], actual[i], expected[i]);

  return errors;
}


function dataframePredicate(actual: any, expected: any) {
  return (actual instanceof DG.DataFrame && expected instanceof DG.DataFrame);
}

function dataframeChecker(actual: DG.DataFrame, expected: DG.DataFrame,
  _state: Readonly<ExpectRunnerState>, checkDeep: NestedChecker): CheckerError[] {
  const errors: CheckerError[] = [];
  if (actual.rowCount !== expected.rowCount) {
    const err = {msg: `Dataframes has different row count: actual row count is ${actual.rowCount} ` +
      `and expected row count is ${expected.rowCount}`};
    errors.push(err);
  }
  for (const column of expected.columns) {
    const actualColumn = actual.columns.byName(column.name);
    if (!actualColumn) {
      const err = {msg: `Column ${column.name} not found`};
      errors.push(err);
      continue;
    }
    if (actualColumn.type !== column.type) {
      const err = {msg: `Column ${column.name} type expected ${column.type} got ${actualColumn.type}`};
      errors.push(err);
      continue;
    }
    for (let i = 0; i < expected.rowCount; i++) {
      const actualValue = actualColumn.get(i);
      const expectedValue = column.get(i);
      checkDeep([column.name, String(i)], actualValue, expectedValue);
    }
  }
  return errors;
}


function arrayBufferPredicate(actual: any, expected: any) {
  return (actual instanceof ArrayBuffer && expected instanceof ArrayBuffer);
}

function arrayBufferChecker(actual: ArrayBuffer, expected: ArrayBuffer,
  state: Readonly<ExpectRunnerState>, checkDeep: NestedChecker) {
  return arrayChecker(Array.from(new Uint8Array(actual)), Array.from(new Uint8Array(expected)), state, checkDeep);
}


function objectPredicate(actual: any, expected: any) {
  return (typeof actual === 'object' && typeof expected === 'object');
}

function objectChecker(actual: Record<any, any>, expected: Record<any, any>,
  _state: Readonly<ExpectRunnerState>, checkDeep: NestedChecker): CheckerError[] {
  const errors: CheckerError[] = [];
  for (const [expectedKey, expectedValue] of Object.entries(expected)) {
    const actualValue = actual[expectedKey];
    checkDeep([expectedKey], actualValue, expectedValue);
  }
  return errors;
}

const defaultOptions: Required<ExpectDeepEqualOptions> = {
  floatTolerance: 0.001,
  maxErrorsReport: 5,
  maxDepth: 128,
  checkersSeq: defaulCheckersSeq,
  prefix: '',
};


class ExpectRunner {
  private options: Required<ExpectDeepEqualOptions>;
  private state: ExpectRunnerState;

  constructor(ops: ExpectDeepEqualOptions = {}) {
    this.options = {...defaultOptions, ...ops};
    this.state = {
      errors: [],
      currentPath: [],
      options: this.options,
    };
  }

  public run(actual: any, expected: any): void {
    this.matchChecker(actual, expected);
    this.checkErrorsLimit(0);
  }

  private checkDeep(path: string[], actual: any, expected: any) {
    const pathLength = path.length;
    this.state.currentPath.push(...path);
    if (this.state.currentPath.length > this.options.maxDepth)
      throw new ExpectError(`Max depth ${this.options.maxDepth} is exceeded, potentially a circular structure`);

    this.matchChecker(actual, expected);
    this.state.currentPath.splice(-pathLength);
  }

  private matchChecker(actual: any, expected: any) {
    let found = false;
    for (const {predicate, checker} of this.options.checkersSeq) {
      if (predicate(actual, expected)) {
        found = true;
        const errors = checker(actual, expected, this.state, this.checkDeep.bind(this));
        if (errors)
          errors.map((e) => this.addError(e));

        break;
      }
    }
    if (!found) {
      const err = {msg: `Type mismatch or no checker defined: expected ${expected}, got ${actual}`};
      this.addError(err);
    }
  }

  private addError(err: CheckerError): void {
    const fullPath = [...(this.state.currentPath ?? []), ...(err.path ?? [])];
    const ferr: CheckerError = {
      path: fullPath,
      msg: err.msg,
    };
    this.state.errors.push(ferr);
    this.checkErrorsLimit();
  }

  private checkErrorsLimit(limit = this.options.maxErrorsReport) {
    if (this.state.errors.length > limit) {
      const report = this.makeReport();
      throw new ExpectError(report);
    }
  }

  private makeReport() {
    const msg = this.state.errors.map((e) => {
      const prefix = (e.path ?? []).map((a) => `['${String(a)}']`).join('.');
      const msg = prefix ? `${prefix}: ${e.msg}` : e.msg;
      return this.options.prefix ? `${this.options.prefix}: ${msg}` : msg;
    }).slice(0, this.options.maxErrorsReport).join('\n');
    return msg;
  }
}

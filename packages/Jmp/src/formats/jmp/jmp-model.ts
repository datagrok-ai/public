// The design matrix of a JMP Fit Model (Standard Least Squares), built from its Effects() the way JMP builds
// it, so that a least-squares fit of it predicts what JMP's fit predicts. See JMP_PLUGIN.md §11 (Phase 3).

import {JslNode, columnNames, findArg, isCall, jslKey, numberValue} from './jsl';
import {JMP_INT_MISSING, JmpColumn, modelingType} from './jmp-types';

type Call = Extract<JslNode, {type: 'call'}>;

/** A column's values as the model reads them: numbers (NaN missing), or level codes (-1 missing). */
export type ModelColumn =
  | {name: string, values: Float64Array, coding?: [number, number]}
  | {name: string, levels: string[], codes: Int32Array};

/** The low and high of a column's Coding property (`Coding({-1, 0.38})`), which a designed experiment sets:
 * JMP fits such a factor scaled to -1..1. */
export function codingOf(property: JslNode | undefined): [number, number] | undefined {
  const list = property?.type === 'call' ? property.args[0] : undefined;
  const [lo, hi] = list?.type === 'list' ? list.items.map(numberValue) : [];
  return lo !== undefined && hi !== undefined && hi !== lo ? [lo, hi] : undefined;
}

/** A table column as the model reads it: a Nominal or Ordinal one, text or numbers, by its levels (numbers in
 * ascending order), a Continuous one by its numbers and Coding. */
export function modelColumn(c: JmpColumn): ModelColumn {
  const d = c.data;
  const name = c.info.name;
  if (d.type === 'string')
    return {name, levels: d.categories, codes: d.indexes.map((k) => d.categories[k] === '' ? -1 : k)};
  const values = d.type === 'float' ? d.values : d.type === 'int' ?
    Float64Array.from(d.values, (v) => v === JMP_INT_MISSING ? NaN : v) : new Float64Array(0);
  if (modelingType(c) === 'continuous') {
    return {name, values,
      coding: codingOf(c.properties?.properties.find((p) => p.name === 'Coding')?.node)};
  }
  const levels = [...new Set(values.filter((v) => !Number.isNaN(v)))].sort((a, b) => a - b);
  return {name, levels: levels.map(String), codes: Int32Array.from(values, (v) => levels.indexOf(v))};
}

/** Why a model can't be fit as plain least squares. */
export class UnfittableModel extends Error {}

export interface ModelDesign {
  /** Names of the design's columns: "Intercept", then one or more per effect. */
  terms: string[];
  /** Row-major, `rows × terms.length`; a row with a missing value in an effect is NaN throughout. */
  x: Float64Array;
  /** The response, transformed as the model names it (`Y(Log(:y))`). */
  y: Float64Array;
  /** The rows JMP fits on: every value present, inside Where(), in the training set of a Validation(). */
  fit: Uint8Array;
}

/** JMP's transforms of a column inside a model, as in `Y(Log(:y))` or `Effects(Sqrt(:dose))`. */
const TRANSFORMS: Record<string, (v: number) => number> = {
  log: Math.log, log10: Math.log10, sqrt: Math.sqrt, exp: Math.exp,
};

/** The design of `call`, a Fit Model, for the response `y` (one of its Y), over `rows` rows of the table. */
export function modelDesign(call: Call, y: JslNode, rows: number, column: (name: string) => ModelColumn):
  ModelDesign {
  if (findArg(call, 'Random Effects'))
    throw new UnfittableModel('Random effects need a mixed model, which isn\'t converted yet');
  for (const f of ['Weight', 'Freq']) {
    if (findArg(call, f))
      throw new UnfittableModel(`The ${f} role isn't applied in a refit`);
  }
  const fit = new Uint8Array(rows).fill(1);
  const response = numeric(y, column, rows);
  for (let i = 0; i < rows; i++)
    fit[i] &= Number(!Number.isNaN(response[i]));
  const where = findArg(call, 'Where');
  if (where)
    keepWhere(where.args[0], column, fit);
  const validation = findArg(call, 'Validation');
  if (validation)
    keepTraining(columnNames(validation), column, fit);

  // JMP centers the continuous factors of crossed terms on their means, unless told not to or they're coded
  const center = numberValue(findArg(call, 'Center Polynomials')?.args[0]) !== 0;
  const effects = (findArg(call, 'Effects')?.args ?? []).map((e) => factors(e, column, rows));
  for (const f of effects.flat())
    f.values?.forEach((v, i) => fit[i] &= Number(!Number.isNaN(v)));
  for (const f of effects.flat())
    f.codes?.forEach((c, i) => fit[i] &= Number(c >= 0));
  const means = new Map<string, number>();
  for (const f of effects.flat().filter((f) => f.values && !f.coded)) {
    let sum = 0;
    let n = 0;
    f.values!.forEach((v, i) => {
      if (fit[i]) {
        sum += v;
        n++;
      }
    });
    means.set(f.name, sum / n);
  }

  const intercept = numberValue(findArg(call, 'No Intercept')?.args[0]) !== 1;
  const terms: string[] = intercept ? ['Intercept'] : [];
  const columns: ((row: number) => number)[] = intercept ? [() => 1] : [];
  for (const effect of effects) {
    const crossed = effect.length > 1;
    // A crossed term is the product of its factors' codes: one column per combination of their columns
    let parts: {name: string, value: (row: number) => number}[] = [{name: '', value: () => 1}];
    for (const f of effect) {
      const own = f.values ? [{name: f.name, value: (r: number) =>
        f.values![r] - (crossed && center && !f.coded ? means.get(f.name)! : 0)}] : effectCodes(f);
      parts = parts.flatMap((p) => own.map((o) => ({name: p.name ? `${p.name}*${o.name}` : o.name,
        value: (r: number) => p.value(r) * o.value(r)})));
    }
    for (const p of parts) {
      terms.push(p.name);
      columns.push(p.value);
    }
  }

  const x = new Float64Array(rows * terms.length);
  for (let r = 0; r < rows; r++) {
    const present = effects.flat().every((f) => f.values ? !Number.isNaN(f.values[r]) : f.codes![r] >= 0);
    for (let j = 0; j < terms.length; j++)
      x[r * terms.length + j] = present ? columns[j](r) : NaN;
  }
  return {terms, x, y: response, fit};
}

interface Factor {
  name: string;
  values?: Float64Array;
  coded?: boolean;
  levels?: string[];
  codes?: Int32Array;
}

/** JMP's effect coding of a nominal factor: a column per level but the last, which is -1 in all of them. */
function effectCodes(f: Factor): {name: string, value: (row: number) => number}[] {
  const last = f.levels!.length - 1;
  return f.levels!.slice(0, last).map((level, k) => ({name: `${f.name}[${level}]`,
    value: (r: number) => f.codes![r] === k ? 1 : f.codes![r] === last ? -1 : 0}));
}

/** The factors of one effect: `:a`, `:a * :b`, `Log(:a)`, `:a & RS` (a response-surface mark). */
function factors(e: JslNode, column: (name: string) => ModelColumn, rows: number): Factor[] {
  if (e.type === 'binary' && e.op === '*')
    return [...factors(e.left, column, rows), ...factors(e.right, column, rows)];
  if (e.type === 'binary' && e.op === '&') {
    if (e.right.type === 'name' && jslKey(e.right.name) === 'random')
      throw new UnfittableModel('Random effects need a mixed model, which isn\'t converted yet');
    return factors(e.left, column, rows);
  }
  if (e.type === 'index')
    throw new UnfittableModel('Nested effects aren\'t converted yet');
  if (e.type === 'column') {
    const c = column(e.name);
    if ('levels' in c)
      return [{name: c.name, levels: usedLevels(c), codes: recode(c)}];
    if (!c.coding)
      return [{name: c.name, values: c.values}];
    const [lo, hi] = c.coding;
    return [{name: c.name, coded: true, values: c.values.map((v) => (v - (hi + lo) / 2) / ((hi - lo) / 2))}];
  }
  if (e.type === 'call' && TRANSFORMS[jslKey(e.name)])
    return [{name: `${e.name}(${referenced(e)})`, values: numeric(e, column, rows)}];
  throw new UnfittableModel('An effect that isn\'t a column, a product of columns or a transform of one');
}

function referenced(e: Call): string {
  return e.args[0]?.type === 'column' ? e.args[0].name : '';
}

/** The levels a factor takes: a level no row holds has no column, as in JMP. */
function usedLevels(c: {levels: string[], codes: Int32Array}): string[] {
  const used = new Set(c.codes);
  return c.levels.filter((_, k) => used.has(k));
}

function recode(c: {levels: string[], codes: Int32Array}): Int32Array {
  const used = usedLevels(c);
  const to = c.levels.map((l) => used.indexOf(l));
  return c.codes.map((k) => k < 0 ? -1 : to[k]);
}

/** A numeric column, or a transform of one, as the model reads it. */
function numeric(e: JslNode, column: (name: string) => ModelColumn, rows: number): Float64Array {
  if (e.type === 'column') {
    const c = column(e.name);
    if ('levels' in c)
      throw new UnfittableModel(`The response "${c.name}" is categorical, which least squares doesn't fit`);
    return c.values;
  }
  if (e.type === 'call' && TRANSFORMS[jslKey(e.name)] && e.args.length === 1) {
    const f = TRANSFORMS[jslKey(e.name)];
    return numeric(e.args[0], column, rows).map((v) => f(v));
  }
  throw new UnfittableModel('A response that isn\'t a column or a transform of one');
}

/** Where(:col == value), the only form saved scripts use: other rows aren't fit. */
function keepWhere(cond: JslNode | undefined, column: (name: string) => ModelColumn, fit: Uint8Array): void {
  if (!(cond?.type === 'binary' && cond.op === '==' && cond.left.type === 'column' &&
    (cond.right.type === 'string' || cond.right.type === 'number')))
    throw new UnfittableModel('A Where() other than a column equal to a value isn\'t applied in a refit');
  const c = column(cond.left.name);
  const value = cond.right.value;
  if ('levels' in c) {
    const k = c.levels.indexOf(String(value));
    c.codes.forEach((code, i) => fit[i] &= Number(code === k));
  } else
    c.values.forEach((v, i) => fit[i] &= Number(v === Number(value)));
}

/** A Validation() column: JMP fits on its training rows, those holding 0. */
function keepTraining(names: string[], column: (name: string) => ModelColumn, fit: Uint8Array): void {
  const c = names.length === 1 ? column(names[0]) : undefined;
  if (!c || 'levels' in c)
    throw new UnfittableModel('A Validation() column that isn\'t numeric isn\'t applied in a refit');
  c.values.forEach((v, i) => fit[i] &= Number(v === 0));
}

/** The responses of a Fit Model: its Y's columns and transforms of columns. */
export function modelResponses(call: Call): JslNode[] {
  return findArg(call, 'Y')?.args.filter((a) => a.type === 'column' || isCall(a, 'Log') || isCall(a, 'Log10') ||
    isCall(a, 'Sqrt') || isCall(a, 'Exp')) ?? [];
}

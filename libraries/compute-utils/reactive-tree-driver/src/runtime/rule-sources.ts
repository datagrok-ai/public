import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {IControllerBase} from '../RuntimeControllers';
import {RuleExpr, RuleSource} from '../config/PipelineConfiguration';
import {evaluate, RuleContext, usedAliases} from './rule-expressions';

export type ValidatorVerdict = {message: string, isError: boolean, isHelper: boolean};

let warnedNoEval = false;

// the annotation's own validators, resolved and run by the platform (1.28+)
function annotationVerdicts(
  controller: IControllerBase, spec: {input: string, call?: string},
): ValidatorVerdict[] | Promise<ValidatorVerdict[]> {
  if (!spec.call || !controller.hasCall(spec.call))
    return [];
  const call = controller.getFirst<DG.FuncCall | undefined>(spec.call);
  if (!call)
    return [];
  if (typeof (call as any).evalParamValidators !== 'function') {
    if (!warnedNoEval) {
      warnedNoEval = true;
      console.warn('RTD: FuncCall.evalParamValidators is not available on this platform, named validators are skipped');
    }
    return [];
  }
  const io = controller.getMatchedPositions(spec.input)[0]?.ioName;
  if (!io)
    return [];
  return call.evalParamValidators(io).catch((e) => [failed('', e)]);
}

function failed(name: string, e: unknown): ValidatorVerdict {
  const detail = e instanceof Error ? e.message : String(e);
  return {message: `Couldn't validate (${name ? name + ': ' : ''}${detail})`, isError: false, isHelper: false};
}

// the documented validator contract: one argument, null/true valid, string is the
// message, false is a generic failure naming the function
async function callValidator(name: string, value: any): Promise<ValidatorVerdict | undefined> {
  const func = DG.Func.byName(name);
  if (!func)
    return failed(name, 'function not found');
  const param = func.inputs[0]?.name;
  if (!param)
    return failed(name, 'validator takes no input');
  try {
    const result = await grok.functions.call(name, {[param]: value});
    if (result == null || result === true)
      return undefined;
    if (result === false)
      return {message: `Validation failed: ${name}`, isError: true, isHelper: false};
    return {message: String(result), isError: true, isHelper: false};
  } catch (e) {
    return failed(name, e);
  }
}

function namedVerdicts(value: any, names: string[]): Promise<ValidatorVerdict[]> {
  return Promise.all(names.map((name) => callValidator(name, value)))
    .then((verdicts) => verdicts.filter((verdict): verdict is ValidatorVerdict => !!verdict));
}

// synchronous whenever nothing is called (absent value, mock mode, old platform), so those
// links stay batchable and virtual-time tests see their results immediately
function resolveValidators(
  controller: IControllerBase, spec: ValidatorsSource,
): ValidatorVerdict[] | Promise<ValidatorVerdict[]> {
  if (spec.names) {
    const value = controller.getFirst(spec.input);
    if (value == null || value === '')
      return [];
    return namedVerdicts(value, spec.names);
  }
  return annotationVerdicts(controller, spec);
}

export type ChoicesVerdict = {
  items: string[], values: Record<string, any>, inList: boolean, row: Record<string, any> | null, rowErrors: string[],
};
type ChoicesResult = Awaited<ReturnType<DG.FuncCall['evalParamChoices']>>;
type ChoicesEntry = {deps: string[], values: any[], result: Promise<ChoicesResult>, landed?: ChoicesResult};

// one evaluation per call and io, shared by the links reading it and redone only when a param
// named in its `dependsOn` changes, as the function form does
const choicesCache = new WeakMap<DG.FuncCall, Map<string, ChoicesEntry>>();
let warnedNoChoices = false;

const notConverted = Symbol('notConverted');

// mirrors the form, which parses the cell's text with the input's editor; a cell it would empty is reported
function convertCell(cell: any, type: string): any {
  if (cell == null)
    return null;
  if (type === DG.TYPE.BOOL)
    return cell === true || cell === 'true';
  if (cell === '' && type !== DG.TYPE.STRING)
    return null;
  switch (type) {
  case DG.TYPE.INT:
  case DG.TYPE.FLOAT:
  case DG.TYPE.NUM: {
    const num = typeof cell === 'number' ? cell : typeof cell === 'string' && cell.trim() ? Number(cell) : NaN;
    return Number.isFinite(num) && (type !== DG.TYPE.INT || Number.isInteger(num)) ? num : notConverted;
  }
  case DG.TYPE.STRING:
    if (typeof cell === 'string')
      return cell;
    return typeof cell === 'number' || typeof cell === 'boolean' ? String(cell) : notConverted;
  case DG.TYPE.DATE_TIME: {
    // evalParamChoices leaves datetime cells as platform objects, which toJs turns into dayjs
    const date = typeof cell === 'string' || typeof cell === 'number' || cell instanceof Date ? dayjs(cell) :
      dayjs.isDayjs(cell) ? cell : DG.toJs(cell);
    return dayjs.isDayjs(date) && date.isValid() ? date : notConverted;
  }
  default:
    return cell;
  }
}

function convertRow(call: DG.FuncCall, cells: Record<string, any>) {
  const types = new Map(call.func.inputs.map((prop) => [prop.name.toLowerCase(), prop.propertyType as string]));
  const row: Record<string, any> = {};
  const rowErrors: string[] = [];
  for (const [column, cell] of Object.entries(cells)) {
    const type = types.get(column.toLowerCase());
    const value = type ? convertCell(cell, type) : cell;
    if (value === notConverted)
      rowErrors.push(`${column}: ${JSON.stringify(cell)} is not a valid ${type}`);
    else
      row[column] = value;
  }
  return {row, rowErrors};
}

function choicesVerdict(call: DG.FuncCall, r: ChoicesResult, value: any): ChoicesVerdict {
  const key = value == null || value === '' ? undefined : String(value);
  const cells = key === undefined ? undefined : r.lookup?.[key];
  const {row, rowErrors} = cells ? convertRow(call, cells) : {row: null, rowErrors: []};
  return {
    items: r.items, values: r.values,
    inList: key === undefined || r.items.includes(key),
    row, rowErrors,
  };
}

function resolveChoices(
  controller: IControllerBase, spec: ChoicesSource,
): ChoicesVerdict | Promise<ChoicesVerdict> | undefined {
  if (!spec.call || !controller.hasCall(spec.call))
    return undefined;
  const call = controller.getFirst<DG.FuncCall | undefined>(spec.call);
  const io = controller.getMatchedPositions(spec.input)[0]?.ioName;
  if (!call || !io)
    return undefined;
  if (typeof (call as any).evalParamChoices !== 'function') {
    if (!warnedNoChoices) {
      warnedNoChoices = true;
      console.warn('RTD: FuncCall.evalParamChoices is not available on this platform, annotation choices are skipped');
    }
    return undefined;
  }
  const value = controller.getFirst(spec.input);
  const byIo = choicesCache.get(call) ?? new Map<string, ChoicesEntry>();
  choicesCache.set(call, byIo);
  let entry = byIo.get(io);
  if (!entry || entry.deps.some((dep, idx) => call.inputs[dep] !== entry!.values[idx])) {
    const snapshot = Object.fromEntries(call.func.inputs.map((prop) => [prop.name, call.inputs[prop.name]]));
    const created: ChoicesEntry = {deps: [], values: [], result: call.evalParamChoices(io)};
    created.result.then((r) => {
      created.deps = r.dependsOn;
      created.values = r.dependsOn.map((dep) => snapshot[dep]);
      created.landed = r;
    }, () => byIo.delete(io));
    byIo.set(io, created);
    entry = created;
  }
  return entry.landed ? choicesVerdict(call, entry.landed, value) :
    entry.result.then((r) => choicesVerdict(call, r, value));
}

type ValidatorsSource = Extract<RuleSource, {validators: any}>['validators'];
type ChoicesSource = Extract<RuleSource, {choices: any}>['choices'];
type JsSource = Extract<RuleSource, {js: any}>['js'];
type FuncSource = Extract<RuleSource, {func: any}>['func'];
type QuerySource = Extract<RuleSource, {query: any}>['query'];
type TableSource = Extract<RuleSource, {table: any}>['table'];

function resolveTable(spec: TableSource) {
  if (spec instanceof DG.DataFrame)
    return spec;
  return typeof spec === 'string' ? DG.DataFrame.fromCsv(spec) : DG.DataFrame.fromCsv(spec.csv, spec.options);
}

function resolveFile(path: string) {
  return /^[a-z][a-z0-9+.-]*:\/\//i.test(path) ? grok.data.loadTable(path) : grok.data.files.openTable(path);
}

// a source that reads no input is loaded once per link
function isConstant(source: RuleSource) {
  if ('file' in source || 'table' in source)
    return true;
  const args = 'func' in source ? source.func.args : 'query' in source ? source.query.args : undefined;
  return args !== undefined && Object.values(args).every((expr) => !usedAliases(expr).length);
}

function resolveJs(controller: IControllerBase, spec: JsSource) {
  return spec.fn(...spec.args.map((alias) => controller.getFirst(alias)));
}

function evaluateArgs(args: Record<string, RuleExpr> | undefined, ctx: RuleContext) {
  return Object.fromEntries(Object.entries(args ?? {}).map(([param, expr]) => [param, evaluate(expr, ctx)]));
}

function resolveFunc(spec: FuncSource, ctx: RuleContext) {
  return grok.functions.call(spec.name, evaluateArgs(spec.args, ctx));
}

function grokType(value: any): string {
  if (typeof value === 'number')
    return Number.isInteger(value) ? DG.TYPE.INT : DG.TYPE.FLOAT;
  if (typeof value === 'boolean')
    return DG.TYPE.BOOL;
  if (value instanceof DG.DataFrame)
    return DG.TYPE.DATA_FRAME;
  if (value instanceof Date)
    return DG.TYPE.DATE_TIME;
  return DG.TYPE.STRING;
}

// an ad-hoc query binds @name only for parameters its text declares
async function resolveQuery(spec: QuerySource, ctx: RuleContext) {
  const args = evaluateArgs(spec.args, ctx);
  const declared = new Set([...spec.sql.matchAll(/^--input:\s*\w+\s+(\w+)/gm)].map((m) => m[1]));
  const header = Object.entries(args)
    .filter(([name]) => !declared.has(name))
    .map(([name, value]) => `--input: ${grokType(value)} ${name}`).join('\n');
  const sql = header ? `${header}\n${spec.sql}` : spec.sql;
  const connection: DG.DataConnection = await grok.functions.eval(spec.connection);
  return connection.query('adhoc', sql).apply(args);
}

const sourceKinds = ['validators', 'choices', 'js', 'func', 'query', 'file', 'table'];

/** Resolves the values a rule declares in `sources`; each alias becomes a context variable.
 *  Returns a plain object when nothing had to be awaited. */
export function resolveSources(
  controller: IControllerBase, sources: Record<string, RuleSource> | undefined, ctx: RuleContext = {$all: {}},
): Record<string, any> | Promise<Record<string, any>> {
  const resolved: Record<string, any> = {};
  const pending: Promise<void>[] = [];
  for (const [alias, source] of Object.entries(sources ?? {})) {
    if (!sourceKinds.some((kind) => kind in source))
      throw new Error(`Unknown rule source ${JSON.stringify(source)} for alias ${alias}`);
    const cache = controller.sourceCache;
    const constant = !!cache && isConstant(source);
    let value: any;
    if (constant && cache.has(alias))
      value = cache.get(alias);
    else {
      value = 'js' in source ? resolveJs(controller, source.js) :
        'func' in source ? resolveFunc(source.func, ctx) :
          'query' in source ? resolveQuery(source.query, ctx) :
            'file' in source ? resolveFile(source.file) :
              'table' in source ? resolveTable(source.table) :
                'choices' in source ? resolveChoices(controller, source.choices) :
                  resolveValidators(controller, source.validators);
      if (constant) {
        cache.set(alias, value);
        if (value instanceof Promise)
          value.catch(() => cache.delete(alias));
      }
    }
    if (value instanceof Promise)
      pending.push(value.then((result) => {resolved[alias] = result;}));
    else
      resolved[alias] = value;
  }
  return pending.length ? Promise.all(pending).then(() => resolved) : resolved;
}

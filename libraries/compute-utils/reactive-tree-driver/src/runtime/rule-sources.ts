import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {defer, forkJoin, from, isObservable, Observable, of} from 'rxjs';
import {catchError, concatMap, defaultIfEmpty, map, shareReplay, take, tap} from 'rxjs/operators';
import {IControllerBase} from '../RuntimeControllers';
import {RuleExpr, RuleSource, SourceKind, SourceSpec} from '../config/PipelineConfiguration';
import {evaluate, RuleContext, usedAliases} from './rule-expressions';

export type ValidatorVerdict = {message: string, isError: boolean, isHelper: boolean};

let warnedNoEval = false;

// the annotation's own validators, resolved and run by the platform (1.28+)
function annotationVerdicts(
  controller: IControllerBase, spec: {input: string, call?: string},
): Observable<ValidatorVerdict[]> {
  if (!spec.call || !controller.hasCall(spec.call))
    return of([]);
  const call = controller.getFirst<DG.FuncCall | undefined>(spec.call);
  if (!call)
    return of([]);
  if (typeof (call as any).evalParamValidators !== 'function') {
    if (!warnedNoEval) {
      warnedNoEval = true;
      console.warn('RTD: FuncCall.evalParamValidators is not available on this platform, named validators are skipped');
    }
    return of([]);
  }
  const io = controller.getMatchedPositions(spec.input)[0]?.ioName;
  if (!io)
    return of([]);
  return defer(() => call.evalParamValidators(io)).pipe(catchError((e) => of([failed('', e)])));
}

function failed(name: string, e: unknown): ValidatorVerdict {
  const detail = e instanceof Error ? e.message : String(e);
  return {message: `Couldn't validate (${name ? name + ': ' : ''}${detail})`, isError: false, isHelper: false};
}

// the documented validator contract: one argument, null/true valid, string is the
// message, false is a generic failure naming the function
function callValidator(name: string, value: any): Observable<ValidatorVerdict | undefined> {
  const func = DG.Func.byName(name);
  if (!func)
    return of(failed(name, 'function not found'));
  const param = func.inputs[0]?.name;
  if (!param)
    return of(failed(name, 'validator takes no input'));
  return defer(() => grok.functions.call(name, {[param]: value})).pipe(
    map((result): ValidatorVerdict | undefined => result == null || result === true ? undefined :
      {message: result === false ? `Validation failed: ${name}` : String(result), isError: true, isHelper: false}),
    catchError((e) => of(failed(name, e))),
  );
}

// synchronous whenever nothing is called (absent value, no names, mock mode, old platform), so those
// links stay batchable and virtual-time tests see their results immediately
function resolveValidators(controller: IControllerBase, spec: ValidatorsSource): Observable<ValidatorVerdict[]> {
  if (!spec.names)
    return annotationVerdicts(controller, spec);
  const value = controller.getFirst(spec.input);
  if (value == null || value === '' || !spec.names.length)
    return of([]);
  return forkJoin(spec.names.map((name) => callValidator(name, value))).pipe(
    map((verdicts) => verdicts.filter((verdict): verdict is ValidatorVerdict => !!verdict)),
  );
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
  const choices = evalChoices(call, io);
  return choices instanceof Promise ?
    choices.then((r) => choicesVerdict(call, r, value)) :
    choicesVerdict(call, choices, value);
}

// a pending evaluation does not know its dependencies yet, so a run waits for it and checks them then
function evalChoices(call: DG.FuncCall, io: string): ChoicesResult | Promise<ChoicesResult> {
  const byIo = choicesCache.get(call) ?? new Map<string, ChoicesEntry>();
  choicesCache.set(call, byIo);
  const entry = byIo.get(io);
  if (entry && !entry.landed)
    return entry.result.then(() => evalChoices(call, io));
  if (entry?.landed && entry.deps.every((dep, idx) => call.inputs[dep] === entry.values[idx]))
    return entry.landed;
  const snapshot = Object.fromEntries(call.func.inputs.map((prop) => [prop.name, call.inputs[prop.name]]));
  const created: ChoicesEntry = {deps: [], values: [], result: call.evalParamChoices(io)};
  created.result.then((r) => {
    created.deps = r.dependsOn;
    created.values = r.dependsOn.map((dep) => snapshot[dep]);
    created.landed = r;
  }, () => byIo.delete(io));
  byIo.set(io, created);
  return created.result;
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

// typed per kind, so a new RuleSource kind fails the build until it says which names it reads
const sourceNames: {[K in SourceKind]: (spec: SourceSpec<K>) => string[]} = {
  validators: (spec) => [spec.input],
  choices: (spec) => [spec.input],
  js: (spec) => spec.args,
  func: (spec) => Object.values(spec.args ?? {}).flatMap((arg) => usedAliases(arg)),
  query: (spec) => Object.values(spec.args ?? {}).flatMap((arg) => usedAliases(arg)),
  file: () => [],
  table: () => [],
};

/** The inputs, driver values and earlier sources a source reads. */
export function sourceReads(source: RuleSource): string[] {
  // expandRule has checked the source, so it has exactly one known kind
  const kind = Object.keys(source)[0] as SourceKind;
  return (sourceNames[kind] as (spec: unknown) => string[])((source as Record<SourceKind, unknown>)[kind]);
}

// a source that reads no input, directly or through the sources it reads, is loaded once per link
function isConstant(source: RuleSource, constants: Set<string>) {
  if ('file' in source || 'table' in source)
    return true;
  if (!('func' in source) && !('query' in source))
    return false;
  return sourceReads(source).every((name) => constants.has(name));
}

function resolveJs(controller: IControllerBase, spec: JsSource, ctx: RuleContext) {
  return spec.fn(...spec.args.map((alias) => Object.hasOwn(ctx, alias) ? ctx[alias] : controller.getFirst(alias)));
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

type SourceResolver = (controller: IControllerBase, spec: any, ctx: RuleContext) => any;

// typed per kind, so a new RuleSource kind fails the build until it gets a resolver
const resolvers: {[K in SourceKind]: (controller: IControllerBase, spec: SourceSpec<K>, ctx: RuleContext) => any} = {
  validators: (controller, spec) => resolveValidators(controller, spec),
  choices: (controller, spec) => resolveChoices(controller, spec),
  js: (controller, spec, ctx) => resolveJs(controller, spec, ctx),
  func: (_controller, spec, ctx) => resolveFunc(spec, ctx),
  query: (_controller, spec, ctx) => resolveQuery(spec, ctx),
  file: (_controller, spec) => resolveFile(spec),
  table: (_controller, spec) => resolveTable(spec),
};

// a resolver or a js fn returns a value, a promise or an observable; a source takes the first value
function toSource$(value: unknown): Observable<any> {
  const value$ = isObservable(value) ? value : value instanceof Promise ? from(value) : of(value);
  return value$.pipe(take(1), defaultIfEmpty<any>(undefined));
}

/** Resolves the values a rule declares in `sources`; each alias becomes a context variable.
 *  A source reads the sources declared before it and waits only for those it reads.
 *  Emits synchronously when nothing has to be awaited. */
export function resolveSources(
  controller: IControllerBase & {sourceCache?: Map<string, any>},
  sources: Record<string, RuleSource> | undefined, ctx: RuleContext = {$all: {}},
): Observable<Record<string, any>> {
  const values: Record<string, Observable<any>> = {};
  const constants = new Set<string>();
  const cache = controller.sourceCache;
  for (const [alias, source] of Object.entries(sources ?? {})) {
    const constant = isConstant(source, constants);
    if (constant)
      constants.add(alias);
    const cached = constant ? cache?.get(alias) : undefined;
    if (cached) {
      values[alias] = cached;
      continue;
    }
    // expandRule has checked the source, so it has exactly one known kind
    const kind = Object.keys(source)[0] as SourceKind;
    const spec = (source as Record<SourceKind, unknown>)[kind];
    const read = sourceReads(source).filter((name) => Object.hasOwn(values, name));
    const value$ = (read.length ? forkJoin(read.map((name) => values[name])) : of([])).pipe(
      concatMap((results) => toSource$((resolvers[kind] as SourceResolver)(controller, spec,
        {...ctx, ...Object.fromEntries(read.map((name, idx) => [name, results[idx]]))}))),
      tap({error: () => cache?.delete(alias)}),
      // a source read by several others runs once, and a loaded cached one replays synchronously
      shareReplay({bufferSize: 1, refCount: false}),
    );
    if (constant)
      cache?.set(alias, value$);
    values[alias] = value$;
  }
  return Object.keys(values).length ? forkJoin(values) : of({});
}

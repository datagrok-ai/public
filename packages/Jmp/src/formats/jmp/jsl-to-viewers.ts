// Turns a saved JMP analysis into Datagrok viewers (JMP_PLUGIN.md §11, Phase 3). Reads only the platform
// and its column roles; everything else in the script is ignored or named in `notes`. No platform imports,
// so the Node tests can run it.

import {
  JslNode, columnNames, findAll, findArg, isCall, jslKey, referencedColumns, statements,
} from './jsl';
import {JmpScript, isSourceScript} from './jmp-scripts';
import type {JmpColumn} from './jmp-types';

/** A viewer to add to the table view: a Datagrok viewer type (DG.VIEWER) and its properties. */
export interface JmpViewerSpec {
  type: string;
  options: Record<string, unknown>;
}

export interface JmpScriptConversion {
  name: string;
  text: string;
  /** The JMP platform the script launches, e.g. "Graph Builder"; undefined when it launches none. */
  platform?: string;
  viewers: JmpViewerSpec[];
  /** Parts of the analysis that didn't come across. */
  notes: string[];
  /** Why nothing came across; set exactly when `viewers` is empty. */
  reason?: string;
}

/** What the converter needs to know about a table column. */
export interface JmpColumnRef {
  name: string;
  numeric: boolean;
  /** For a saved prediction formula: the column it predicts (JMP's `Predicting` column property), when
   * the property names one, and the model terms of its formula ({@link modelTerms}). */
  predicts?: string;
  terms?: string[];
}

/** What the converter needs to know about the table's columns. */
export function columnRefs(columns: JmpColumn[]): JmpColumnRef[] {
  return columns.map((c) => {
    const property = (name: string) => c.properties?.properties.find((p) => p.name === name)?.node;
    const predicting = property('Predicting');
    const formula = property('Formula')?.args[0];
    return {name: c.info.name, numeric: c.kind !== 'string',
      ...(predicting && {predicts: columnNames(predicting)[0], terms: formula ? formulaTerms(formula) : []})};
  });
}

/** A model term as a key: the columns it involves, sorted, so `:A * :B`, `(:B - 1.5) * :A` and `:A & RS`
 * with `:B & RS` crossed all read "A*B". */
function termKey(columns: string[]): string {
  return [...new Set(columns.map(jslKey))].sort().join('*');
}

/** The terms of a saved prediction formula: the columns of each of its top-level summands, without the
 * intercept. */
function formulaTerms(formula: JslNode): string[] {
  const summands = (n: JslNode): JslNode[] =>
    n.type === 'binary' && (n.op === '+' || n.op === '-') ? [...summands(n.left), ...summands(n.right)] : [n];
  return [...new Set(summands(formula).map((t) => termKey(referencedColumns(t))).filter((k) => k !== ''))].sort();
}

/** The terms of a Fit Model, from its Effects(), in the form of {@link formulaTerms}. */
function modelTerms(call: Call): string[] {
  return [...new Set((findArg(call, 'Effects')?.args ?? []).map((e) => termKey(referencedColumns(e))))].sort();
}

type Call = Extract<JslNode, {type: 'call'}>;

const HISTOGRAM = 'Histogram';
const SCATTER = 'Scatter plot';
const BOX = 'Box plot';
const BAR = 'Bar chart';
const LINE = 'Line chart';
const PIE = 'Pie chart';
const PIVOT = 'Pivot table';
const TRELLIS = 'Trellis plot';

class Unconvertible extends Error {}

class Columns {
  private readonly byKey = new Map<string, JmpColumnRef>();
  private readonly predictions: JmpColumnRef[];
  /** @param transforms names of the columns the script makes itself, with Transform Column() */
  constructor(columns: JmpColumnRef[], private readonly transforms: string[]) {
    for (const c of columns)
      this.byKey.set(jslKey(c.name), c);
    this.predictions = columns.filter((c) => c.terms !== undefined);
  }

  /** The prediction formulas saved for `y`: those whose Predicting property names it, and, as JMP keeps the
   * name a response had when it was fit, those named "Pred Formula <y>…" whose Predicting names no column
   * of the table (the response was renamed since) or none at all. */
  predictionsOf(y: string): JmpColumnRef[] {
    const prefix = jslKey(`Pred Formula ${y}`);
    return this.predictions.filter((c) => c.predicts !== undefined && this.byKey.has(jslKey(c.predicts)) ?
      jslKey(c.predicts) === jslKey(y) : jslKey(c.name).startsWith(prefix));
  }

  /** The table's column for a name used in the script; JSL matches names ignoring case and spaces. */
  get(name: string): JmpColumnRef {
    const c = this.byKey.get(jslKey(name));
    if (!c && this.transforms.some((t) => jslKey(t) === jslKey(name)))
      throw new Unconvertible(`Column "${name}" is a transform the script computes, which isn't converted yet`);
    if (!c)
      throw new Unconvertible(`Column "${name}" isn't in the imported table`);
    return c;
  }
}

/** JMP statistic names, as in Tabulate's Statistics() and Graph Builder's Summary Statistic(). */
const STATISTICS: Record<string, string> = {
  mean: 'avg', sum: 'sum', n: 'count', min: 'min', max: 'max', median: 'med', stddev: 'stdev',
  variance: 'variance', nmissing: 'nulls', nunique: 'unique',
};

function statistic(name: string | undefined, notes: string[]): string | undefined {
  if (name === undefined)
    return undefined;
  const agg = STATISTICS[jslKey(name)];
  if (!agg)
    notes.push(`Statistic "${name}" shown as the mean`);
  return agg ?? 'avg';
}

export function convertScript(script: JmpScript, columns: JmpColumnRef[]): JmpScriptConversion {
  const result: JmpScriptConversion = {name: script.name, text: script.text, viewers: [], notes: []};
  if (!script.body) {
    result.reason = `The script couldn't be read: ${script.error ?? 'unknown error'}`;
    return result;
  }
  if (isSourceScript(script)) {
    result.reason = 'Re-imports the file the table came from; not an analysis';
    return result;
  }
  const launches = platformCalls(script.body);
  if (launches.length === 0) {
    result.reason = isProgram(script.body) ? 'A JSL program rather than a saved analysis' :
      'Doesn\'t launch an analysis';
    return result;
  }
  const cols = new Columns(columns,
    findAll(script.body, 'Transform Column').map(stringArg).filter((n) => n !== undefined));
  const reasons: string[] = [];
  for (const call of launches) {
    result.platform ??= call.name;
    const converter = PLATFORMS[jslKey(call.name)];
    if (!converter) {
      reasons.push(`JMP's ${call.name} isn't converted to Datagrok viewers yet`);
      continue;
    }
    try {
      const notes: string[] = [];
      const viewers = converter(call, cols, notes);
      for (const v of viewers)
        v.options.title ??= script.name;
      noteIgnored(call, notes);
      result.viewers.push(...viewers);
      result.notes.push(...notes);
    } catch (e) {
      if (!(e instanceof Unconvertible))
        throw e;
      reasons.push(e.message);
    }
  }
  if (result.viewers.length === 0)
    result.reason = [...new Set(reasons)].join('; ');
  else
    result.notes = [...new Set([...result.notes, ...reasons])];
  return result;
}

/** The platform launches in a script: `Bivariate(...)`, `dt << Bivariate(...)`, one per statement. */
function platformCalls(body: JslNode): Call[] {
  if (isProgram(body))
    return [];
  const out: Call[] = [];
  for (const s of statements(body)) {
    const node = topCall(s) ?? s;
    if (node.type === 'call' && !NOT_PLATFORMS.has(jslKey(node.name)))
      out.push(node);
    else if (node.type === 'list')
      out.push(...node.items.flatMap(platformCalls));
  }
  return out;
}

function isProgram(body: JslNode): boolean {
  return statements(body).some((s) => s.type === 'list' ? s.items.some(isProgram) :
    PROGRAM.has(jslKey(topCall(s)?.name ?? '')));
}

/** The call a statement makes: `Platform(...)`, `dt << Platform(...)` or `obj = dt << Platform(...)`. */
function topCall(s: JslNode): Call | undefined {
  let node = s;
  if (node.type === 'binary' && node.op === '=')
    node = node.right;
  while (node.type === 'binary' && node.op === '<<')
    node = node.right;
  return node.type === 'call' ? node : undefined;
}

/** Calls that set things up rather than launch an analysis. */
const NOT_PLATFORMS = new Set(['namesdefaulttohere', 'currentdatatable', 'datatable', 'open', 'clearglobals',
  'deletesymbols', 'print', 'show', 'write', 'wait', 'include', 'subset', 'roweditor']);

/** Calls of a hand-written JSL program, which a table can keep as a script: loops, dialogs, windows built
 * from parts. Such a script isn't a saved analysis, whatever platforms it launches along the way. */
const PROGRAM = new Set(['for', 'while', 'if', 'foreachrow', 'getcolumnnames', 'newwindow', 'columndialog',
  'dialog', 'function', 'try', 'eval', 'expr']);

/** Column names as notes and reasons write them: quoted, so tools can tell them from JSL's own words. */
function quoted(names: string[]): string {
  return names.map((n) => `"${n}"`).join(', ');
}

/** Options that change what a viewer shows and that we don't carry over. */
function noteIgnored(call: Call, notes: string[]): void {
  const by = findArg(call, 'By');
  if (by)
    notes.push(`Separate analyses by ${quoted(columnNames(by))} not reproduced`);
  if (findArg(call, 'Where'))
    notes.push('The Where() row filter isn\'t applied');
  if (findArg(call, 'Weight') || findArg(call, 'Freq'))
    notes.push('Weight and Freq roles aren\'t applied');
}

type Converter = (call: Call, cols: Columns, notes: string[]) => JmpViewerSpec[];

const PLATFORMS: Record<string, Converter> = {
  distribution: distribution,
  bivariate: bivariate,
  oneway: oneway,
  contingency: contingency,
  fitybyx: fitYByX,
  fitmodel: fitModel,
  fitgroup: fitGroup,
  graphbuilder: graphBuilder,
  tabulate: tabulate,
};

/** Distribution: a histogram per continuous column, a bar chart of counts per categorical one. A column
 * launched as Nominal or Ordinal Distribution is categorical even when it holds numbers. */
function distribution(call: Call, cols: Columns): JmpViewerSpec[] {
  const categorical = new Set([
    ...findAll(call, 'Nominal Distribution'), ...findAll(call, 'Ordinal Distribution'),
  ].flatMap((c) => columnNames(findArg(c, 'Column')).map((n) => cols.get(n).name)));
  const names = [
    ...columnNames(findArg(call, 'Column')), ...columnNames(findArg(call, 'Y')),
    ...findAll(call, 'Continuous Distribution').flatMap((c) => columnNames(findArg(c, 'Column'))),
    ...categorical,
  ].map((n) => cols.get(n).name);
  if (names.length === 0)
    throw new Unconvertible('Distribution names no columns');
  return [...new Set(names)].map((n) => {
    const c = cols.get(n);
    return c.numeric && !categorical.has(c.name) ?
      {type: HISTOGRAM, options: {valueColumnName: c.name, title: `Distribution of ${c.name}`}} :
      {type: BAR, options: {splitColumnName: c.name, valueAggrType: 'count', title: `Distribution of ${c.name}`}};
  });
}

/** Pairs of Y and X columns; JMP launches one analysis per pair. */
function pairs(call: Call, cols: Columns): {y: JmpColumnRef, x: JmpColumnRef}[] {
  const ys = columnNames(findArg(call, 'Y')).map((n) => cols.get(n));
  const xs = columnNames(findArg(call, 'X')).map((n) => cols.get(n));
  if (ys.length === 0 || xs.length === 0)
    throw new Unconvertible(`${call.name} needs both Y and X columns`);
  return ys.flatMap((y) => xs.map((x) => ({y, x})));
}

function bivariate(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  const linear = findAll(call, 'Fit Line').length > 0 || findAll(call, 'Fit Polynomial')
    .some((f) => f.args[0]?.type === 'number' && f.args[0].value === 1);
  for (const fit of ['Fit Polynomial', 'Fit Spline', 'Fit Each Value', 'Fit Orthogonal', 'Density Ellipse',
    'Nonpar Density', 'Fit Mean', 'Kernel Smoother']) {
    if (findAll(call, fit).some((f) => !(jslKey(fit) === 'fitpolynomial' && f.args[0]?.type === 'number' &&
      f.args[0].value === 1)))
      notes.push(`${fit} not shown`);
  }
  return pairs(call, cols).map(({x, y}) => ({type: SCATTER, options: {
    xColumnName: x.name, yColumnName: y.name, title: `${y.name} by ${x.name}`,
    ...(linear ? {showRegressionLine: true} : {}),
  }}));
}

function oneway(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  for (const test of ['Means', 'Means and Std Dev', 'Anova', 't Test', 'All Pairs', 'Each Pair', 'With Best',
    'With Control', 'Nonparametric', 'Unequal Variances', 'Wilcoxon Test', 'Compare Densities']) {
    if (findAll(call, test).length > 0)
      notes.push(`${test} results not shown`);
  }
  return pairs(call, cols).map(({x, y}) => {
    if (!y.numeric)
      throw new Unconvertible(`${call.name} needs a continuous Y column; "${y.name}" is categorical`);
    return {type: BOX, options: {
      valueColumnName: y.name, categoryColumnNames: [x.name], title: `${y.name} by ${x.name}`,
    }};
  });
}

/** Fit Model: JMP's Actual by Predicted plot for each Y, drawn from the prediction formula the fit saved to
 * the table, with the Y = X line a perfect fit would follow. The model isn't refit, so a Y without a saved
 * prediction formula has no plot. */
function fitModel(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  return fitModels([call], cols, notes);
}

/** Fit Group: the Fit Models it holds, each converted on its own. JMP separates them with `,` or `;`. */
function fitGroup(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  const models = call.args.flatMap(statements).filter((a): a is Call => isCall(a, 'Fit Model'));
  if (models.length === 0)
    throw new Unconvertible('Fit Group holds no Fit Model');
  return fitModels(models, cols, notes);
}

/** The Actual by Predicted plots of Fit Models; what didn't come across is told once for them all. */
function fitModels(models: Call[], cols: Columns, notes: string[]): JmpViewerSpec[] {
  const viewers: JmpViewerSpec[] = [];
  const reasons: string[] = [];
  const unsaved: string[] = [];
  const otherModel: string[] = [];
  let profiler = false;
  for (const call of models) {
    try {
      const personality = stringArg(findArg(call, 'Personality'));
      if (personality !== undefined && jslKey(personality) !== 'standardleastsquares')
        throw new Unconvertible(`Fit Model's ${personality} personality isn't converted yet`);
      const ys = columnNames(findArg(call, 'Y')).map((n) => cols.get(n));
      if (ys.length === 0)
        throw new Unconvertible('Fit Model names no Y column');
      const terms = modelTerms(call).join(',');
      for (const y of ys) {
        const saved = cols.predictionsOf(y.name);
        // A table often keeps the predictions of several fits of one response; only the one with this
        // model's terms is its Actual by Predicted.
        const predicted = saved.find((c) => c.terms!.join(',') === terms);
        if (!predicted) {
          (saved.length > 0 ? otherModel : unsaved).push(y.name);
          continue;
        }
        viewers.push({type: SCATTER, options: {
          xColumnName: predicted.name, yColumnName: y.name, title: `${y.name} Actual by Predicted`,
          formulaLines: JSON.stringify([{type: 'line', formula: 'Y = X', title: 'Actual = Predicted'}]),
        }});
      }
      profiler ||= findAll(call, 'Profiler').some((p) => p.args[0]?.type === 'number' && p.args[0].value === 1);
    } catch (e) {
      if (!(e instanceof Unconvertible))
        throw e;
      reasons.push(e.message);
    }
  }
  if (unsaved.length > 0) {
    reasons.push(`The ${unsaved.length === 1 ? 'fit' : 'fits'} of ${quoted(unsaved)} saved no prediction formula ` +
      'to the table, and the model isn\'t refit');
  }
  if (otherModel.length > 0) {
    reasons.push(`The prediction formulas saved for ${quoted(otherModel)} come from other models than ` +
      `${otherModel.length === 1 ? 'this one' : 'these'}, and the model isn't refit`);
  }
  if (viewers.length === 0)
    throw new Unconvertible([...new Set(reasons)].join('; '));
  notes.push(...new Set(reasons));
  notes.push('The model\'s reports (fit summary, estimates, effect tests) and its other plots not shown');
  if (profiler)
    notes.push('Profiler not shown');
  return viewers;
}

function contingency(call: Call, cols: Columns): JmpViewerSpec[] {
  return pairs(call, cols).map(({x, y}) => ({type: BAR, options: {
    splitColumnName: x.name, stackColumnName: y.name, valueAggrType: 'count', title: `${y.name} by ${x.name}`,
  }}));
}

/** Fit Y by X picks its analysis from the columns' types, as JMP does: bivariate, oneway, contingency or
 * (categorical Y by continuous X) logistic, which isn't converted. */
function fitYByX(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  const out: JmpViewerSpec[] = [];
  const failed: string[] = [];
  for (const {x, y} of pairs(call, cols)) {
    const one: Call = {type: 'call', name: call.name, args: [
      {type: 'call', name: 'Y', args: [{type: 'column', name: y.name}]},
      {type: 'call', name: 'X', args: [{type: 'column', name: x.name}]},
      ...call.args.filter((a) => !isCall(a, 'Y') && !isCall(a, 'X'))]};
    try {
      if (!y.numeric && x.numeric)
        throw new Unconvertible(`"${y.name}" by "${x.name}" is a logistic fit, which isn't converted yet`);
      out.push(...(!y.numeric ? contingency(one, cols) :
        x.numeric ? bivariate(one, cols, notes) : oneway(one, cols, notes)));
    } catch (e) {
      if (!(e instanceof Unconvertible))
        throw e;
      failed.push(e.message);
    }
  }
  if (out.length === 0)
    throw new Unconvertible(failed.join('; '));
  notes.push(...failed);
  return out;
}

/** Graph Builder: variables by role, and the elements drawn with them. The first element decides the
 * viewer; a Line of Fit over points becomes the scatter plot's regression line. Panels (Group X, Group Y,
 * Wrap) make it a trellis plot of that viewer. */
function graphBuilder(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  const vars = findArg(call, 'Variables');
  if (!vars)
    throw new Unconvertible('Graph Builder names no variables');
  const role = (r: string): JmpColumnRef[] =>
    vars.args.filter((a) => isCall(a, r)).flatMap((a) => columnNames(a as Call)).map((n) => cols.get(n));
  const inner = graphBuilderViewer(call, role, notes);
  const groupX = [...role('Group X'), ...role('Wrap')].map((c) => c.name);
  const groupY = role('Group Y').map((c) => c.name);
  if (role('Page').length > 0)
    notes.push(`Pages by ${quoted(role('Page').map((c) => c.name))} not reproduced`);
  if (groupX.length === 0 && groupY.length === 0)
    return [inner];
  return [{type: TRELLIS, options: {
    viewerType: inner.type, xColumnNames: groupX, yColumnNames: groupY, innerViewerLook: inner.options,
    ...(role('Wrap').length > 0 ? {useTiledView: true} : {}),
  }}];
}

function graphBuilderViewer(call: Call, role: (r: string) => JmpColumnRef[], notes: string[]): JmpViewerSpec {
  const xs = role('X');
  const ys = role('Y');
  const color = role('Color')[0] ?? role('Overlay')[0];
  const size = role('Size')[0];
  if (xs.length > 1 || ys.length > 1)
    notes.push('Only the first X and Y variables are shown');

  const elements = findAll(call, 'Elements').flatMap((e) => e.args.filter((a): a is Call =>
    a.type === 'call' && !isCall(a, 'Position')));
  if (findAll(call, 'Elements').length > 1)
    notes.push('Only the first panel is shown');
  const kinds = elements.map((e) => jslKey(e.name));
  const first = kinds[0] ?? 'points';
  const x = xs[0];
  const y = ys[0];
  const others = elements.slice(1).filter((e) => !['lineoffit', 'points', 'smoother'].includes(jslKey(e.name)));
  for (const e of others)
    notes.push(`${e.name} not shown`);

  switch (first) {
  case 'points':
  case 'smoother':
  case 'lineoffit': {
    if (!x || !y)
      throw new Unconvertible('A Graph Builder scatter needs X and Y variables');
    if (kinds.includes('smoother'))
      notes.push('Smoother not shown');
    // JMP fits a line per Overlay group; Datagrok can only split the regression by the color column
    const fit = kinds.includes('lineoffit');
    const overlay = role('Overlay')[0];
    const perGroup = fit && overlay !== undefined && overlay === color && !overlay.numeric;
    if (fit && overlay && !perGroup)
      notes.push(`Line of Fit by "${overlay.name}" shown as one line`);
    return {type: SCATTER, options: {
      xColumnName: x.name, yColumnName: y.name,
      ...(color ? {colorColumnName: color.name} : {}), ...(size ? {sizeColumnName: size.name} : {}),
      ...(fit ? {showRegressionLine: true} : {}), ...(perGroup ? {regressionPerCategory: true} : {}),
    }};
  }
  case 'line': {
    if (!x || !y)
      throw new Unconvertible('A Graph Builder line needs X and Y variables');
    return {type: LINE, options: {
      xColumnName: x.name, yColumnNames: [y.name],
      ...(color && !color.numeric ? {splitColumnName: color.name} : {}),
    }};
  }
  case 'bar': {
    const category = x ?? y;
    if (!category)
      throw new Unconvertible('A Graph Builder bar chart needs a variable');
    const value = x && y ? y : undefined;
    const stat = statistic(stringArg(findAll(elements[0], 'Summary Statistic')[0]), notes);
    return {type: BAR, options: {
      splitColumnName: category.name,
      ...(value ? {valueColumnName: value.name, valueAggrType: stat ?? 'avg'} : {valueAggrType: 'count'}),
      ...(color && !color.numeric && color !== category ? {stackColumnName: color.name} : {}),
    }};
  }
  case 'boxplot': {
    const value = y?.numeric ? y : x?.numeric ? x : undefined;
    if (!value)
      throw new Unconvertible('A Graph Builder box plot needs a continuous variable');
    const category = value === y ? x : y;
    return {type: BOX, options: {
      valueColumnName: value.name, ...(category ? {categoryColumnNames: [category.name]} : {}),
    }};
  }
  case 'histogram': {
    const value = x?.numeric ? x : y?.numeric ? y : undefined;
    if (!value)
      throw new Unconvertible('A Graph Builder histogram needs a continuous variable');
    return {type: HISTOGRAM, options: {valueColumnName: value.name, ...(color ? {colorColumnName: color.name} : {})}};
  }
  case 'pie': {
    const category = x ?? y;
    if (!category)
      throw new Unconvertible('A Graph Builder pie chart needs a variable');
    return {type: PIE, options: {categoryColumnName: category.name}};
  }
  default:
    throw new Unconvertible(`Graph Builder's ${elements[0].name} element isn't converted yet`);
  }
}

function stringArg(call: Call | undefined): string | undefined {
  const a = call?.args[0];
  return a?.type === 'string' ? a.value : a?.type === 'name' ? a.name : undefined;
}

/** Tabulate: grouping columns in the row table become the pivot table's rows, those in the column
 * table its columns; analysis columns are aggregated with the statistics named. */
function tabulate(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  const grouping = (table: string): string[] => findAll(call, table)
    .flatMap((t) => findAll(t, 'Grouping Columns')).flatMap(columnNames).map((n) => cols.get(n).name);
  const groupBy = grouping('Row Table');
  const pivot = grouping('Column Table');
  const analysis = [...new Set(findAll(call, 'Analysis Columns').flatMap(columnNames).map((n) => cols.get(n).name))];
  const stats = findAll(call, 'Statistics').flatMap((s) => s.args)
    .map((a) => a.type === 'name' ? a.name : a.type === 'string' ? a.value : undefined)
    .filter((s): s is string => s !== undefined)
    .map((s) => statistic(s, notes)!);
  if (groupBy.length === 0 && pivot.length === 0 && analysis.length === 0)
    throw new Unconvertible('Tabulate names no columns');
  const aggs = analysis.length === 0 ? [] : stats.length === 0 ? ['sum'] : [...new Set(stats)];
  return [{type: PIVOT, options: {
    groupByColumnNames: groupBy, pivotColumnNames: pivot,
    aggregateColumnNames: aggs.flatMap(() => analysis), aggregateAggTypes: aggs.flatMap((a) => analysis.map(() => a)),
  }}];
}

/** Converts every script; the "Source" import script and failed reads are kept, so the caller can list them. */
export function convertScripts(scripts: JmpScript[], columns: JmpColumnRef[]): JmpScriptConversion[] {
  return scripts.map((s) => convertScript(s, columns));
}

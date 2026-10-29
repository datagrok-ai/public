// Turns a saved JMP analysis into Datagrok viewers (JMP_PLUGIN.md §11, Phase 3). Reads only the platform
// and its column roles; everything else in the script is ignored or named in `notes`. No platform imports,
// so the Node tests can run it.

import {
  JslNode, columnNames, findAll, findArg, isCall, jslKey, numberValue, referencedColumns, statements,
} from './jsl';
import {JmpScript, isSourceScript} from './jmp-scripts';
import type {JmpColumn} from './jmp-types';

/** A viewer to add to the table view: a Datagrok viewer type (DG.VIEWER) and its properties. */
export interface JmpViewerSpec {
  type: string;
  options: Record<string, unknown>;
  /** What tells this viewer from the others the analysis opens, such as one panel's Y; added to its title. */
  subtitle?: string;
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
        v.options.title ??= v.subtitle ? `${script.name}: ${v.subtitle}` : script.name;
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
  'deletesymbols', 'print', 'show', 'write', 'wait', 'include', 'subset', 'roweditor', 'transformcolumn']);

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

/** Fit Group: the analyses it holds, each converted as if launched on its own, except that its Fit Models
 * are told about together. JMP separates them with `,` or `;`. */
function fitGroup(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  const members = call.args.flatMap(statements)
    .filter((a): a is Call => a.type === 'call' && !NOT_PLATFORMS.has(jslKey(a.name)) &&
      !isCall(a, 'Where') && !isCall(a, 'By'));
  if (members.length === 0)
    throw new Unconvertible('Fit Group holds no analysis');
  const models = members.filter((m) => isCall(m, 'Fit Model'));
  const viewers: JmpViewerSpec[] = [];
  const reasons: string[] = [];
  const convert = (f: (notes: string[]) => JmpViewerSpec[]): void => {
    const own: string[] = [];
    try {
      viewers.push(...f(own));
      notes.push(...own);
    } catch (e) {
      if (!(e instanceof Unconvertible))
        throw e;
      reasons.push(e.message);
    }
  };
  if (models.length > 0)
    convert((own) => fitModels(models, cols, own));
  for (const m of members.filter((m) => !isCall(m, 'Fit Model'))) {
    const converter = PLATFORMS[jslKey(m.name)];
    if (!converter) {
      reasons.push(`JMP's ${m.name} isn't converted to Datagrok viewers yet`);
      continue;
    }
    convert((own) => {
      const out = converter(m, cols, own);
      noteIgnored(m, own);
      return out;
    });
  }
  if (viewers.length === 0)
    throw new Unconvertible([...new Set(reasons)].join('; '));
  notes.push(...new Set(reasons));
  return viewers;
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
      // A response such as Log(:y) is fit on its transform, which no saved column holds
      const transforms = findArg(call, 'Y')?.args.filter((a): a is Call => a.type === 'call') ?? [];
      for (const t of transforms) {
        reasons.push(`The response ${t.name}(${quoted(referencedColumns(t))}) is a transform, ` +
          'which isn\'t converted yet');
      }
      if (ys.length === 0 && transforms.length === 0)
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

/** One panel of a Graph Builder: the X and Y variables of its place in the grid, and its elements. */
interface GraphPanel {
  xs: JmpColumnRef[];
  ys: JmpColumnRef[];
  elements: Call[];
}

/** Graph Builder: variables by role, and the elements drawn with them. Each panel of the grid, set by the X
 * and Y variables (`Y(:a), Y(:b)` stacks two; `Y(:b, Position(1))` overlays b on the first), becomes a
 * viewer. A panel's first element decides the viewer; a Line of Fit over points becomes the scatter plot's
 * regression line. Group X, Group Y and Wrap make each a trellis plot of that viewer. */
function graphBuilder(call: Call, cols: Columns, notes: string[]): JmpViewerSpec[] {
  const vars = findArg(call, 'Variables');
  if (!vars)
    throw new Unconvertible('Graph Builder names no variables');
  const roleCalls = (r: string): Call[] => vars.args.filter((a): a is Call => isCall(a, r));
  const role = (r: string): JmpColumnRef[] => roleCalls(r).flatMap((a) => columnNames(a)).map((n) => cols.get(n));
  const slots = (r: string): JmpColumnRef[][] => {
    const out: JmpColumnRef[][] = [];
    for (const a of roleCalls(r)) {
      const position = numberValue(findArg(a, 'Position')?.args[0]);
      const columns = columnNames(a).map((n) => cols.get(n));
      if (position !== undefined && out[position - 1])
        out[position - 1].push(...columns);
      else
        out.push(columns);
    }
    return out;
  };
  const xSlots = slots('X');
  const ySlots = slots('Y');
  const blocks = findAll(call, 'Elements');
  const panels: GraphPanel[] = (blocks.length > 0 ? blocks : [undefined]).map((b) => {
    const position = b && findArg(b, 'Position');
    const [i, j] = [0, 1].map((k) => numberValue(position?.args[k]) ?? 1);
    return {xs: xSlots[i - 1] ?? [], ys: ySlots[j - 1] ?? [],
      elements: b?.args.filter((a): a is Call => a.type === 'call' && !isCall(a, 'Position')) ?? []};
  });
  const report = graphBuilderReport(call);
  const groupX = [...role('Group X'), ...role('Wrap')].map((c) => c.name);
  // JMP stacks a graph per page, as the rows of a trellis plot do
  const groupY = [...role('Page'), ...role('Group Y')].map((c) => c.name);
  const applied = new Set<string>();
  const viewers = panels.map((panel): JmpViewerSpec => {
    const inner = graphBuilderViewer(panel, role, notes);
    for (const axis of applyAxes(inner, report.axes, cols))
      applied.add(axis);
    const subtitle = panels.length > 1 ? (panel.ys[0] ?? panel.xs[0])?.name : undefined;
    if (report.title)
      inner.options.title = subtitle ? `${report.title}: ${subtitle}` : report.title;
    if (groupX.length === 0 && groupY.length === 0)
      return {...inner, subtitle};
    return {type: TRELLIS, subtitle, options: {
      viewerType: inner.type, xColumnNames: groupX, yColumnNames: groupY, innerViewerLook: inner.options,
      ...(role('Wrap').length > 0 ? {useTiledView: true} : {}),
      ...(inner.options.title ? {title: inner.options.title} : {}),
    }};
  });
  const unapplied = [...report.axes].filter(([k, a]) => !applied.has(k) &&
    (a.min !== undefined || a.max !== undefined || a.log || a.refLines.length > 0)).map(([, a]) => a.column);
  if (unapplied.length > 0)
    notes.push(`Axis settings of ${quoted(unapplied)} not applied`);
  return viewers;
}

function graphBuilderViewer(panel: GraphPanel, role: (r: string) => JmpColumnRef[], notes: string[]):
  JmpViewerSpec {
  const {xs, ys, elements} = panel;
  const color = role('Color')[0] ?? role('Overlay')[0];
  const size = role('Size')[0];
  const kinds = elements.map((e) => jslKey(e.name));
  const first = kinds[0] ?? 'points';
  const x = xs[0];
  const y = ys[0];
  const others = elements.slice(1)
    .filter((e) => ![first, 'lineoffit', 'points', 'smoother'].includes(jslKey(e.name)));
  // JMP overlays every Y of a line chart on one axis, as Datagrok does; the other viewers take one
  const overlaid = [...xs.slice(1), ...first === 'line' ? [] : ys.slice(1)].map((c) => c.name);
  if (overlaid.length > 0)
    notes.push(`${quoted(overlaid)} overlaid on the same axis not shown`);
  for (const e of others) {
    // The box plot draws each median already; only the line through them is lost
    const medianLine = first === 'boxplot' && jslKey(e.name) === 'line' &&
      jslKey(stringArg(findAll(e, 'Summary Statistic')[0]) ?? '') === 'median';
    notes.push(medianLine ? 'The line joining the group medians isn\'t drawn' : `${e.name} not shown`);
  }

  switch (first) {
  case 'points':
  case 'smoother':
  case 'lineoffit': {
    if (!x || !y)
      throw new Unconvertible('A Graph Builder scatter needs X and Y variables');
    // JMP fits a line or a smoother per Overlay group; Datagrok can only split them by the color column
    const fit = kinds.includes('lineoffit');
    const smoother = kinds.includes('smoother');
    const overlay = role('Overlay')[0];
    const perGroup = overlay !== undefined && overlay === color && !overlay.numeric;
    if (fit && overlay && !perGroup)
      notes.push(`Line of Fit by "${overlay.name}" shown as one line`);
    // Datagrok has no smoothing spline; its moving average is the closest curve through the points
    if (smoother)
      notes.push('Smoother shown as a moving average, not JMP\'s smoothing spline');
    return {type: SCATTER, options: {
      xColumnName: x.name, yColumnName: y.name,
      ...(color ? {colorColumnName: color.name} : {}), ...(size ? {sizeColumnName: size.name} : {}),
      ...(fit ? {showRegressionLine: true} : {}), ...(fit && perGroup ? {regressionPerCategory: true} : {}),
      ...(smoother ? {showMovingAverageLine: true} : {}),
      ...(smoother && perGroup ? {movingAveragePerCategory: true} : {}),
    }};
  }
  case 'line': {
    if (!x || !y)
      throw new Unconvertible('A Graph Builder line needs X and Y variables');
    return {type: LINE, options: {
      xColumnName: x.name, yColumnNames: ys.map((c) => c.name),
      ...(color && !color.numeric ? {splitColumnName: color.name} : {}),
    }};
  }
  case 'bar': {
    // The bars stand on the categorical variable: on X they are vertical, on Y horizontal
    const category = x && !(x.numeric && y && !y.numeric) ? x : y;
    if (!category)
      throw new Unconvertible('A Graph Builder bar chart needs a variable');
    const other = category === x ? y : x;
    const value = other?.numeric ? other : undefined;
    const stat = statistic(stringArg(findAll(elements[0], 'Summary Statistic')[0]), notes);
    return {type: BAR, options: {
      splitColumnName: category.name, orientation: category === x ? 'vertical' : 'horizontal',
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

/** An axis as Graph Builder's report sets it: `Dispatch({}, "<column>", ScaleBox, {Min(0), Max(60), ...})`. */
interface JmpAxis {
  column: string;
  min?: number;
  max?: number;
  log: boolean;
  refLines: {value: number, label: string}[];
}

/** The report settings of a Graph Builder that have a Datagrok counterpart: its title and its axes, by
 * column name. Only the first panel's axes (`ScaleBox`, not `ScaleBox(2)`) are read, as only it is shown. */
function graphBuilderReport(call: Call): {title?: string, axes: Map<string, JmpAxis>} {
  let title: string | undefined;
  const axes = new Map<string, JmpAxis>();
  for (const d of findAll(call, 'Dispatch')) {
    const [, target, box, settings] = d.args;
    if (target?.type !== 'string' || box?.type !== 'name' || settings?.type !== 'list')
      continue;
    const items = settings.items.filter((i): i is Call => i.type === 'call');
    const setting = (name: string) => items.find((i) => isCall(i, name));
    const num = (c: Call | undefined) => numberValue(c?.args[0]);
    if (jslKey(box.name) === 'texteditbox' && jslKey(target.value) === 'graphtitle')
      title = stringArg(setting('Set Text')) || title;
    else if (jslKey(box.name) === 'scalebox') {
      axes.set(jslKey(target.value), {
        column: target.value, min: num(setting('Min')), max: num(setting('Max')),
        log: jslKey(stringArg(setting('Scale')) ?? '') === 'log',
        refLines: items.filter((i) => isCall(i, 'Add Ref Line') && num(i) !== undefined)
          .map((i) => ({value: num(i)!, label: i.args[3]?.type === 'string' ? i.args[3].value : ''})),
      });
    }
  }
  return {title, axes};
}

/** Viewer properties for an axis bound to a column, by viewer: [column property, min, max, axis type]. */
const AXES: Record<string, [string, string, string, string][]> = {
  [SCATTER]: [['xColumnName', 'xMin', 'xMax', 'xAxisType'], ['yColumnName', 'yMin', 'yMax', 'yAxisType']],
  [HISTOGRAM]: [['valueColumnName', 'valueMin', 'valueMax', 'xAxisType']],
  [BOX]: [['valueColumnName', 'valueMin', 'valueMax', 'axisType']],
  [BAR]: [['valueColumnName', 'valueMin', 'valueMax', 'axisType']],
};

/** Sets a viewer's axis ranges, log scales and reference lines from Graph Builder's axes; returns the axes
 * it used. */
function applyAxes(viewer: JmpViewerSpec, axes: Map<string, JmpAxis>, cols: Columns): string[] {
  const o = viewer.options;
  const lines: Record<string, unknown>[] = [];
  const applied = new Set<string>();
  for (const [prop, min, max, type] of AXES[viewer.type] ?? []) {
    const column = o[prop] as string | undefined;
    const axis = column !== undefined ? axes.get(jslKey(column)) : undefined;
    if (!axis)
      continue;
    applied.add(jslKey(column!));
    if (axis.min !== undefined)
      o[min] = axis.min;
    if (axis.max !== undefined)
      o[max] = axis.max;
    // A Datagrok scatter plot's range filters the table by default; JMP's axis only frames the plot
    if (viewer.type === SCATTER && (axis.min !== undefined || axis.max !== undefined))
      o.zoomAndFilter = 'no action';
    if (axis.log)
      o[type] = 'logarithmic';
    for (const r of axis.refLines) {
      lines.push({type: 'line', formula: `\${${cols.get(column!).name}} = ${r.value}`,
        ...(r.label ? {title: r.label} : {})});
    }
  }
  if (lines.length > 0)
    o.formulaLines = JSON.stringify(lines);
  return [...applied];
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

import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {JslNode, JslSyntaxError, columnNames, findArg, isCall, parseJsl, statements} from '../formats/jmp/jsl';
import {readJmp} from '../formats/jmp/jmp-parser';
import {JmpScript, readTableScripts} from '../formats/jmp/jmp-scripts';
import {JmpInfo} from '../formats/jmp/jmp-types';
import {JmpColumnRef, JmpScriptConversion, columnRefs, convertScript} from '../formats/jmp/jsl-to-viewers';
import {fixture} from './jmp-reader-tests';

// Scripts in the form JMP saves them; the column names follow JMP's Big Class sample.
const BIG_CLASS: JmpColumnRef[] = [
  {name: 'name', numeric: false}, {name: 'age', numeric: true}, {name: 'sex', numeric: false},
  {name: 'height', numeric: true}, {name: 'weight', numeric: true},
];

function convert(text: string, name = 'script'): JmpScriptConversion {
  return convertScript({name, text, body: parseJsl(text)}, BIG_CLASS);
}

function only(c: JmpScriptConversion): {type: string, options: Record<string, unknown>} {
  expect(c.viewers.length, 1, `viewers: ${JSON.stringify(c.viewers)}, reason: ${c.reason}`);
  return c.viewers[0];
}

async function testFile(name: string): Promise<Uint8Array> {
  return await _package.files.readAsBytes(`tests/${name}.jmp`);
}

category('JMP: JSL', () => {
  test('calls, columns and literals', async () => {
    const n = parseJsl('Bivariate(Y(:weight), X(:"height in"n), Fit Line({Line Color({212, 73, 88})}), ' +
      'Title("a\\!"b"))');
    expect(n.type, 'call');
    const call = n as Extract<JslNode, {type: 'call'}>;
    expect(call.name, 'Bivariate');
    expectArray(columnNames(findArg(call, 'Y')), ['weight']);
    expectArray(columnNames(findArg(call, 'X')), ['height in']);
    const title = findArg(call, 'Title')!.args[0];
    expect(title.type === 'string' && title.value, 'a"b');
  });

  test('names with spaces, any case', async () => {
    const n = parseJsl('fit y BY x( y( :Name("weight") ), x( :height ) )') as Extract<JslNode, {type: 'call'}>;
    expect(n.name, 'fit y BY x');
    expect(isCall(n, 'Fit Y by X'), true);
    expectArray(columnNames(findArg(n, 'Y')), ['weight']);
  });

  test('statements, sends and comments', async () => {
    const n = parseJsl('// set up\nNames Default To Here(1); dt = Current Data Table(); ' +
      '/* run */ dt << Distribution(Y(:age));');
    expect(statements(n).length, 3);
  });

  test('numbers, dates and operators', async () => {
    const n = parseJsl('{1, -2.5e3, .5, 01Jan2020, a + b * c, x[1, 2], !y}');
    expect(n.type === 'list' && n.items.length, 7);
  });

  test('raw strings and empty arguments', async () => {
    const n = parseJsl('f("\\[a "quoted" path]\\", , 2,)') as Extract<JslNode, {type: 'call'}>;
    expect(n.args.length, 2);
    expect(n.args[0].type === 'string' && n.args[0].value, 'a "quoted" path');
  });

  test('missing values, leading sends, names called', async () => {
    const n = parseJsl('Add Response(Minimize, "Time", ., ., .)') as Extract<JslNode, {type: 'call'}>;
    expect(n.args.length, 5);
    expect(n.args[2].type === 'number' && Number.isNaN(n.args[2].value), true);
    const g = parseJsl('Fit Group(Fit Model(Y(:y)), <<{Profiler(1)})') as Extract<JslNode, {type: 'call'}>;
    expect(g.args[1].type === 'unary' && g.args[1].op, '<<');
    const named = parseJsl('Name("1a fit")(Fit Model(Y(:y)))') as Extract<JslNode, {type: 'call'}>;
    expect(named.name, '1a fit');
    expect(isCall(named.args[0], 'Fit Model'), true);
    expect(isCall(parseJsl('"Value Ordering"({"Yes", "No"})'), 'Value Ordering'), true);

    const bare = parseJsl('Y(:HC Met255 %Oxidation, :ints * 2)') as Extract<JslNode, {type: 'call'}>;
    expect(bare.args[0].type === 'column' && bare.args[0].name, 'HC Met255 %Oxidation');
    expect(bare.args[1].type, 'binary');
  });

  test('names with #, %, \' and dots, columns given a format', async () => {
    const names = (text: string) => (parseJsl(text) as Extract<JslNode, {type: 'call'}>).args
      .map((a) => a.type === 'call' ? a.name : '');
    expectArray(names('Profiler(DAS28#0M(5.1984, Lock(0)), Peptides%CV(30, Hi Limit(30)))'),
      ['DAS28#0M', 'Peptides%CV']);
    expectArray(names('Multivariate(Spearman\'s ρ(1), Scatterplot Matrix(1))'),
      ['Spearman\'s ρ', 'Scatterplot Matrix']);
    const assign = parseJsl('Ratio5_25.5_50_mean = Numeric');
    expect(assign.type === 'binary' && assign.left.type === 'name' && assign.left.name, 'Ratio5_25.5_50_mean');
    const format = parseJsl('Set Format(Mean(:"Age (years)"n(10, 1), :Sqrt AEX(10, 2)))') as
      Extract<JslNode, {type: 'call'}>;
    const mean = format.args[0] as Extract<JslNode, {type: 'call'}>;
    expect(mean.args.length, 2);
    const first = mean.args[0];
    expect(first.type === 'index' && first.target.type === 'column' && first.target.name, 'Age (years)');
    expect(first.type === 'index' && first.args.length, 2);
    const aliases = parseJsl('Aliases(["en" => "Pred Formula", "ja" => "予測式"])') as Extract<JslNode, {type: 'call'}>;
    const map = aliases.args[0];
    expect(map.type === 'list' && map.items.length, 2);
    expect(map.type === 'list' && map.items[1].type === 'binary' && map.items[1].op, '=>');
  });

  test('statements in arguments, matrices, loops, any script', async () => {
    const glued = parseJsl('Trend(dt = Current Data Table(); L = dt << Get Column Names(continuous); )') as
      Extract<JslNode, {type: 'call'}>;
    expect(glued.args.length, 1);
    expect(statements(glued.args[0]).length, 2);
    const m = parseJsl('Formula([-0.5 1e-3 2; 3 .])') as Extract<JslNode, {type: 'call'}>;
    const rows = m.args[0].type === 'list' ? m.args[0].items : [];
    expect(rows.length, 2);
    expect(rows[0].type === 'list' && rows[0].items.map((n) => n.type === 'number' ? n.value : '').join(' '),
      '-0.5 0.001 2');
    expect(parseJsl('[1 2 3]').type === 'list', true);
    const loop = parseJsl('For(k = 1, k <= 3, k++, Oneway(Y(:y)))') as Extract<JslNode, {type: 'call'}>;
    expect(loop.args[2].type === 'unary' && loop.args[2].op, '++');
    expect(isCall(parseJsl('Contingency(Total %(1), Col %(1))'), 'Contingency'), true);
    expect((parseJsl('Open("x.xls", Arbeitsblätter("Daten"))') as Extract<JslNode, {type: 'call'}>).args[1].type,
      'call');
    const jp = parseJsl('Y(:（説明）年齢)') as Extract<JslNode, {type: 'call'}>;
    expect(jp.args[0].type === 'column' && jp.args[0].name, '（説明）年齢');
  });

  test('syntax errors', async () => {
    for (const bad of ['f(', 'f(a b', '"open', 'f(:)', 'f(1))', '@']) {
      let error: unknown;
      try {
        parseJsl(bad);
      } catch (e) {
        error = e;
      }
      expect(error instanceof JslSyntaxError, true, bad);
    }
  });

  test('deep nesting is refused, not a stack overflow', async () => {
    await expectExceptionAsync(async () => {
      parseJsl('f('.repeat(5000) + ')'.repeat(5000));
    }, (e) => e instanceof JslSyntaxError && /nested too deeply/.test(e.message));
  });

  // Truncated and mutated scripts either parse or fail with a JslSyntaxError.
  test('fuzz', async () => {
    const src = 'Graph Builder(Size(534, 456), Variables(X(:height), Y(:weight), Overlay(:sex)), ' +
      'Elements(Points(X, Y, Legend(9)), Line Of Fit(X, Y, Legend(11))), SendToReport(Dispatch({}, "x", ' +
      'ScaleBox, {Min(50), Max(70)})));';
    const chars = '(){}[],;:"\\! <<=a1.n';
    let seed = 777;
    const rnd = () => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    for (let k = 0; k < 2000; k++) {
      let s = src;
      if (k % 2 === 0)
        s = s.slice(0, Math.floor(rnd() * s.length));
      else {
        const p = Math.floor(rnd() * s.length);
        s = s.slice(0, p) + chars[Math.floor(rnd() * chars.length)] + s.slice(p + 1);
      }
      try {
        const body = parseJsl(s);
        convertScript({name: 'fuzz', text: s, body}, BIG_CLASS);
      } catch (e) {
        if (!(e instanceof JslSyntaxError))
          throw new Error(`mutation ${k} ${JSON.stringify(s)}: ${e instanceof Error ? e.name + ': ' + e.message : e}`);
      }
    }
  });
}, {node: true});

category('JMP: saved scripts', () => {
  test('the Source script of a JMP-made file', async () => {
    const t = await fixture('compact_UInt8');
    expect(t.scripts.length, 1);
    expect(t.scripts[0].name, 'Source');
    expect(t.scripts[0].text.startsWith('Open("/C:/'), true);
    expect(t.scripts[0].error === undefined, true, t.scripts[0].error);
    const c = convertScript(t.scripts[0], [{name: 'A', numeric: false}]);
    expect(c.viewers.length, 0);
    expect(/not an analysis/.test(c.reason!), true, c.reason);
  });

  test('a Graph Builder saved by JMP 15', async () => {
    const t = await readJmp(await _package.files.readAsBytes('tests/jmp15-graph-builder.jmp'));
    expect(t.scripts.map((sc) => sc.name).join('|'), 'Source|HC Met255 %Oxidation vs. Week_5C');
    expect(t.scripts[1].error === undefined, true, t.scripts[1].error);
    const c = convertScript(t.scripts[1], [{name: 'Condition', numeric: false}, {name: 'Temperature', numeric: false},
      {name: 'Week', numeric: true}, {name: 'HC Met255 %Oxidation', numeric: true}]);
    expect(c.viewers.length, 1);
    expect(c.viewers[0].type, 'Scatter plot');
    expect(c.viewers[0].options.yColumnName, 'HC Met255 %Oxidation');
    expect(c.viewers[0].options.colorColumnName, 'Condition');
    expect(c.viewers[0].options.showRegressionLine, true);
    expect(c.viewers[0].options.regressionPerCategory, true);
  });

  test('files without scripts', async () => {
    for (const name of ['example1', 'rowstate', 'compressed', 'date'])
      expect((await fixture(name)).scripts.length, 0, name);
  });

  test('a saved Graph Builder (synthetic file)', async () => {
    const t = await readJmp(await testFile('saved-bar-chart'));
    expect(t.scripts.length, 1);
    expect(t.scripts[0].name, 'Bar chart of A');
    const v = only(convertScript(t.scripts[0], [{name: 'A', numeric: false}]));
    expect(v.type, 'Bar chart');
    expect(v.options.splitColumnName, 'A');
    expect(v.options.valueAggrType, 'count');
  });

  test('a script whose name JSL can\'t write bare', async () => {
    const text = new TextEncoder().encode('Name("1a_logistic-regr")(Fit Model(Y(:sex), Effects(:height)))');
    const bytes = new Uint8Array([0x03, ...new Uint8Array(new Uint32Array([text.length]).buffer), ...text]);
    const scripts = readTableScripts(bytes, {propertiesRange: [0, bytes.length]} as JmpInfo);
    expect(scripts.length, 1);
    expect(scripts[0].name, '1a_logistic-regr');
    expect(scripts[0].error === undefined, true);
    expect(scripts[0].text, 'Fit Model(Y(:sex), Effects(:height))');
  });

  test('a script name JMP 15 wrote bare', async () => {
    const text = new TextEncoder().encode('HC Met255 %Oxidation vs. Week_5C(Graph Builder(Variables(X(:Week))))');
    const bytes = new Uint8Array([0x03, ...new Uint8Array(new Uint32Array([text.length]).buffer), ...text]);
    const scripts = readTableScripts(bytes, {propertiesRange: [0, bytes.length]} as JmpInfo);
    expect(scripts.length, 1);
    expect(scripts[0].name, 'HC Met255 %Oxidation vs. Week_5C');
    expect(scripts[0].error === undefined, true);
    expect(isCall(scripts[0].body, 'Graph Builder'), true);
  });

  test('New Script form', async () => {
    const text = 'New Script("Weight by height", Bivariate(Y(:weight), X(:height)))';
    const root = parseJsl(text) as Extract<JslNode, {type: 'call'}>;
    const script: JmpScript = {name: 'Weight by height', text, body: root.args[1]};
    expect(only(convertScript(script, BIG_CLASS)).type, 'Scatter plot');
  });
}, {node: true});

category('JMP: analyses', () => {
  test('Distribution', async () => {
    const c = convert('Distribution(Continuous Distribution(Column(:weight)), Nominal Distribution(Column(:age)), ' +
      'Nominal Distribution(Column(:sex)))');
    expectArray(c.viewers.map((v) => v.type), ['Histogram', 'Bar chart', 'Bar chart']);
    expect(c.viewers[0].options.valueColumnName, 'weight');
    expect(c.viewers[1].options.splitColumnName, 'age');
  });

  test('Distribution, short form', async () => {
    expectArray(convert('Distribution(Column(:height, :sex))').viewers.map((v) => v.type), ['Histogram', 'Bar chart']);
  });

  test('Bivariate with a line fit', async () => {
    const c = convert('Bivariate(Y(:weight), X(:height), Fit Line({Line Color({212, 73, 88})}), Fit Spline(0.1))');
    const v = only(c);
    expect(v.type, 'Scatter plot');
    expect(v.options.xColumnName, 'height');
    expect(v.options.yColumnName, 'weight');
    expect(v.options.showRegressionLine, true);
    expectArray(c.notes, ['Fit Spline not shown']);
  });

  test('Oneway, with By noted', async () => {
    const c = convert('Oneway(Y(:height), X(:sex), Means(1), Box Plots(1), By(:age))');
    const v = only(c);
    expect(v.type, 'Box plot');
    expect(v.options.valueColumnName, 'height');
    expectArray(v.options.categoryColumnNames as string[], ['sex']);
    expect(c.notes.some((n) => /Separate analyses by "age"/.test(n)), true);
  });

  test('Fit Y by X picks the analysis from the column types', async () => {
    const c = convert('Fit Y by X(Y(:weight, :sex), X(:height, :age))');
    // weight~height and weight~age are bivariate; sex by a continuous X is logistic, not converted
    expectArray(c.viewers.map((v) => v.type), ['Scatter plot', 'Scatter plot']);
    expect(c.notes.filter((n) => /logistic/.test(n)).length, 2);
    expect(only(convert('Fit Y by X(Y(:sex), X(:name))')).type, 'Bar chart');
    expect(only(convert('Fit Y by X(Y(:height), X(:sex))')).type, 'Box plot');
  });

  test('Graph Builder elements', async () => {
    const gb = (vars: string, elements: string) => only(convert(
      `Graph Builder(Size(534, 456), Show Control Panel(0), Variables(${vars}), Elements(${elements}))`));
    let v = gb('X(:height), Y(:weight), Overlay(:sex)', 'Points(X, Y, Legend(9)), Line Of Fit(X, Y, Legend(11))');
    expect(v.type, 'Scatter plot');
    expect(v.options.colorColumnName, 'sex');
    expect(v.options.showRegressionLine, true);
    expect(v.options.regressionPerCategory, true);
    v = gb('X(:height), Y(:weight)', 'Points(X, Y), Line Of Fit(X, Y)');
    expect(v.options.regressionPerCategory === undefined, true);
    const byAge = convert('Graph Builder(Variables(X(:height), Y(:weight), Overlay(:age)), ' +
      'Elements(Points(X, Y), Line Of Fit(X, Y)))');
    expect(byAge.viewers[0].options.regressionPerCategory === undefined, true);
    expectArray(byAge.notes, ['Line of Fit by "age" shown as one line']);
    v = gb('X(:age), Y(:height)', 'Line(X, Y, Legend(3))');
    expect(v.type, 'Line chart');
    expectArray(v.options.yColumnNames as string[], ['height']);
    v = gb('X(:sex), Y(:height)', 'Bar(X, Y, Legend(4), Summary Statistic("Median"))');
    expect(v.type, 'Bar chart');
    expect(v.options.valueAggrType, 'med');
    expect(v.options.orientation, 'vertical');
    v = gb('X(:sex)', 'Bar(X, Legend(4))');
    expect(v.options.valueAggrType, 'count');
    v = gb('X(:height), Y(:sex)', 'Bar(X, Y)');
    expect(v.options.splitColumnName, 'sex');
    expect(v.options.valueColumnName, 'height');
    expect(v.options.orientation, 'horizontal');
    v = gb('X(:sex), Y(:height)', 'Box Plot(X, Y, Legend(5))');
    expect(v.type, 'Box plot');
    expect(v.options.valueColumnName, 'height');
    v = gb('X(:weight)', 'Histogram(X, Legend(6))');
    expect(v.type, 'Histogram');
    v = gb('X(:sex)', 'Pie(X, Legend(7))');
    expect(v.type, 'Pie chart');
  });

  test('Graph Builder panels make a trellis plot', async () => {
    let v = only(convert('Graph Builder(Variables(X(:height), Y(:weight), Group X(:sex), Group Y(:age)), ' +
      'Elements(Points(X, Y)))', 'Panels'));
    expect(v.type, 'Trellis plot');
    expect(v.options.viewerType, 'Scatter plot');
    expectArray(v.options.xColumnNames as string[], ['sex']);
    expectArray(v.options.yColumnNames as string[], ['age']);
    expect((v.options.innerViewerLook as Record<string, unknown>).xColumnName, 'height');
    expect(v.options.title, 'Panels');
    v = only(convert('Graph Builder(Variables(X(:sex), Y(:height), Wrap(:age)), Elements(Bar(X, Y)))'));
    expect(v.options.viewerType, 'Bar chart');
    expectArray(v.options.xColumnNames as string[], ['age']);
    expect(v.options.useTiledView, true);
    v = only(convert('Graph Builder(Variables(X(:height), Y(:weight), Page(:sex), Group Y(:age)), ' +
      'Elements(Points(X, Y)))'));
    expectArray(v.options.yColumnNames as string[], ['sex', 'age']);
  });

  test('Graph Builder: a viewer per panel, overlays named', async () => {
    const c = convert('Graph Builder(Variables(X(:sex), Y(:height), Y(:weight)), ' +
      'Elements(Position(1, 1), Bar(X, Y), Bar(X, Y)), Elements(Position(1, 2), Box Plot(X, Y)), ' +
      'SendToReport(Dispatch({}, "height", ScaleBox, {Add Ref Line(60, "Solid", "Red", "", 1)}), ' +
      'Dispatch({}, "weight", ScaleBox, {Min(0)})))', 'Sizes');
    expectArray(c.viewers.map((v) => `${v.type}: ${v.options.title}`),
      ['Bar chart: Sizes: height', 'Box plot: Sizes: weight']);
    expect(c.viewers[0].options.valueColumnName, 'height');
    expect(JSON.parse(c.viewers[0].options.formulaLines as string)[0].formula, '${height} = 60');
    expect(c.viewers[1].options.valueMin, 0);
    expectArray(c.notes, []);
    const overlay = convert('Graph Builder(Variables(X(:age), Y(:height), Y(:weight, Position(1))), ' +
      'Elements(Points(X, Y(1), Y(2))))');
    expectArray(overlay.notes, ['"weight" overlaid on the same axis not shown']);
  });

  test('Graph Builder: what doesn\'t come across is named', async () => {
    const c = convert('Graph Builder(Variables(X(:height), Y(:weight)), Elements(Points(X, Y), Smoother(X, Y)))');
    only(c);
    expect(c.notes.some((n) => /Smoother shown as a moving average/.test(n)), true);
    expect(c.viewers[0].options.showMovingAverageLine, true);
    expect(/Contour/.test(convert('Graph Builder(Variables(X(:height), Y(:weight)), Elements(Contour(X, Y)))').reason!),
      true);
  });

  test('Graph Builder: title, axis ranges, log scales and reference lines', async () => {
    const c = convert('Graph Builder(Variables(X(:height), Y(:weight), Overlay(:sex)), ' +
      'Elements(Points(X, Y), Smoother(X, Y)), SendToReport(' +
      'Dispatch({}, "height", ScaleBox, {Format("Fixed Dec", 12, 0), Min(-50), Max(75), Inc(5), Minor Ticks(1)}), ' +
      'Dispatch({}, "weight", ScaleBox, {Scale("Log"), Min(60), Max(200), Add Ref Line(100, "Solid", "Red", ' +
      '"Limit", 1)}), Dispatch({}, "height", ScaleBox(2), {Min(0)}), ' +
      'Dispatch({}, "graph title", TextEditBox, {Set Text("Weight by height")}), ' +
      'Dispatch({}, "Y title", TextEditBox, {Set Text("kg")})))', 'Growth');
    const o = only(c).options;
    expect(o.title, 'Weight by height');
    expect(o.xMin, -50);
    expect(o.xMax, 75);
    expect(o.zoomAndFilter, 'no action');
    expect(o.xAxisType === undefined, true);
    expect(o.yAxisType, 'logarithmic');
    expect(o.yMin, 60);
    expectArray(JSON.parse(o.formulaLines as string).map((l: {formula: string, title: string}) =>
      `${l.title}: ${l.formula}`), ['Limit: ${weight} = 100']);
    expect(o.movingAveragePerCategory, true);
    const box = only(convert('Graph Builder(Variables(X(:sex), Y(:height)), Elements(Box Plot(X, Y)), ' +
      'SendToReport(Dispatch({}, "height", ScaleBox, {Min(0), Add Ref Line(60, "Solid", "Red", "", 1)}), ' +
      'Dispatch({}, "graph title", TextEditBox, {Set Text("")})))', 'Heights'));
    expect(box.options.valueMin, 0);
    expect(box.options.title, 'Heights');
    const unapplied = convert('Graph Builder(Variables(X(:sex), Y(:height)), Elements(Box Plot(X, Y)), ' +
      'SendToReport(Dispatch({}, "Height (in)", ScaleBox, {Min(0)})))');
    expectArray(unapplied.notes, ['Axis settings of "Height (in)" not applied']);
  });

  test('Graph Builder: a line per Y', async () => {
    const c = convert('Graph Builder(Variables(X(:age), Y(:height), Y(:weight, Position(1))), ' +
      'Elements(Line(X, Y(1), Y(2))))');
    expectArray(only(c).options.yColumnNames as string[], ['height', 'weight']);
    expectArray(c.notes, []);
  });

  test('Tabulate', async () => {
    const v = only(convert('Tabulate(Show Control Panel(0), Add Table(' +
      'Column Table(Analysis Columns(:height, :weight), Statistics(Mean, Std Dev)), ' +
      'Row Table(Grouping Columns(:sex))))'));
    expect(v.type, 'Pivot table');
    expectArray(v.options.groupByColumnNames as string[], ['sex']);
    expectArray(v.options.aggregateColumnNames as string[], ['height', 'weight', 'height', 'weight']);
    expectArray(v.options.aggregateAggTypes as string[], ['avg', 'avg', 'stdev', 'stdev']);
  });

  test('launched from a data table, with column names in other forms', async () => {
    const v = only(convert('Names Default To Here(1); dt = Current Data Table(); ' +
      'dt << Graph Builder(Variables(X(:"SEX"n), Y(:Name("Weight"))), Elements(Box Plot(X, Y)));'));
    expect(v.options.valueColumnName, 'weight');
    expectArray(v.options.categoryColumnNames as string[], ['sex']);
  });

  test('viewers are titled after the script', async () => {
    expect(only(convert('Graph Builder(Variables(X(:height), Y(:weight)), Elements(Points(X, Y)))', 'Growth'))
      .options.title, 'Growth');
  });

  test('Fit Model: Actual by Predicted, from the saved prediction formula', async () => {
    const cols = [...BIG_CLASS, {name: 'Pred Formula weight', numeric: true, predicts: 'weight', terms: ['height']}];
    const fit = (text: string) => convertScript({name: 'Fit', text, body: parseJsl(text)}, cols);
    const one = fit('Fit Model(Y(:weight), Effects(:height), Personality("Standard Least Squares"), ' +
      'Run(:weight << {Profiler(1)}))');
    expect(one.viewers.length, 1);
    const o = one.viewers[0].options;
    expect(one.viewers[0].type, 'Scatter plot');
    expect(o.xColumnName, 'Pred Formula weight');
    expect(o.yColumnName, 'weight');
    expect(o.title, 'weight Actual by Predicted');
    expect(JSON.parse(o.formulaLines as string)[0].formula, 'Y = X');
    expect(one.notes.includes('Profiler not shown'), true);
    const group = fit('Fit Group(Fit Model(Y(:weight), Effects(:height)), Fit Model(Y(:height), Effects(:age)))');
    expect(group.platform, 'Fit Group');
    expect(group.viewers.length, 1);
    expect(group.notes.filter((n) => /reports/.test(n)).length, 1);
    expect(/The fit of "height" saved no prediction formula/.test(group.notes.join()), true);
    expect(fit('Fit Model(Y(:weight), Effects(:height, :age))').reason,
      'The prediction formulas saved for "weight" come from other models than this one, and the model isn\'t refit');
    expect(/Generalized Regression personality isn't converted/
      .test(fit('Fit Model(Y(:weight), Effects(:height), Personality("Generalized Regression"))').reason!), true);
  });

  test('not converted, with a reason', async () => {
    expect(/"weight" saved no prediction formula/
      .test(convert('Fit Model(Y(:weight), Effects(:height), Run)').reason!), true);
    expect(/Column "nope"/.test(convert('Bivariate(Y(:weight), X(:nope))').reason!), true);
    expect(/Doesn't launch an analysis/.test(convert('Names Default To Here(1)').reason!), true);
    const program = convert('Names Default To Here(1); cols = Current Data Table() << Get Column Names(Numeric); ' +
      'For(i = 1, i <= N Items(cols), i++, Distribution(Y(cols[i])))');
    expect(program.platform === undefined, true);
    expect(program.reason, 'A JSL program rather than a saved analysis');
    const unread = convertScript({name: 'broken', text: 'f(', error: 'Expected ")"'}, BIG_CLASS);
    expect(/couldn't be read/.test(unread.reason!), true);
  });
}, {node: true});

// Files researchers published openly (CREDITS.md §5): what JMP 10–18 saved, in English and Japanese.
category('JMP: public files', () => {
  const read = async (name: string) => readJmp(await _package.files.readAsBytes(`public/${name}.jmp`));

  test('JMP 13, in Japanese: Oneway analyses', async () => {
    const t = await read('dryad-t76hdr7xq-hadmsc');
    expect(t.info.jmpVersion, '13.2.1');
    const cols = t.columns.map((c) => ({name: c.info.name, numeric: c.kind !== 'string'}));
    const analyses = t.scripts.filter((s) => s.name !== 'Source').map((s) => convertScript(s, cols));
    expect(analyses.length, 3);
    expect(analyses[0].name, 'No.によるsynteninの二変量の関係');
    for (const a of analyses) {
      expect(a.platform, 'Oneway');
      expect(a.viewers[0]?.type, 'Box plot');
    }
    expect(analyses[2].viewers[0].options.valueColumnName, 'Tsg101');
    const tac = await read('dryad-t76hdr7xq-tac2wks');
    expect(tac.columns.map((c) => c.info.name).includes('LVPWｄ'), true);
    expect(tac.scripts.filter((s) => s.name !== 'Source' && !s.error).length, 4);
  });

  test('JMP 11 and 12: formulas and notes', async () => {
    const keiser = await read('dryad-60nn7-keiser-2016');
    expect(keiser.info.jmpVersion, '12.1.0');
    expect(keiser.info.nrows, 76);
    const diff = keiser.columns.find((c) => c.info.name === 'Difference in boldness between pair')!;
    expect(diff.properties!.properties[0].text,
      'Formula(:Exposed spider boldness value - :Susceptible spider boldness value)');
    const condition = keiser.columns.find((c) => c.info.name === 'Susceptible spider body condition')!;
    expect(condition.properties!.properties[0].text, 'Notes("Bivariate Save Residuals Linear Fit")');
    const acoustic = await read('dryad-6b340-acoustic');
    expect(acoustic.info.jmpVersion, '11.0.0');
    expect(acoustic.info.nrows, 161);
    const groups = acoustic.columns.find((c) => c.info.name === 'MaxFreq Groups')!;
    expect(groups.properties!.properties[0].node !== undefined, true);
    expect(groups.properties!.properties[0].text, 'Formula(If(:MaxFreq_mean > 1900, 1, 2))');
  });

  test('JMP 10: hidden columns, value labels and colors', async () => {
    const t = await read('dataverse-jylsxj-rlrva');
    expect(t.info.jmpVersion, '10.0.0');
    expect(t.info.hiddenIndex.length, 21);
    expect(t.info.visibleIndex.length + t.info.hiddenIndex.length, t.info.ncols);
    expect(t.columns[t.info.hiddenIndex[0]].info.name, 'Source Table');
    const labels = t.columns.find((c) => c.info.name === 'Labels')!.properties!;
    expectArray(labels.valueLabels!.keys, [0, 1, 2, 3, 4]);
    expectArray(labels.valueLabels!.labels, ['', 'A', 'V', 'R', 'D']);
    expect(labels.properties.some((p) => p.name === 'Value Colors' && p.node !== undefined), true);
  });

  test('JMP 11: value labels and a Graph Builder', async () => {
    const t = await read('dataverse-hqaihu-figure4');
    expect(t.info.jmpVersion, '11.2.0');
    expectArray(t.columns.find((c) => c.info.name === 'Adopted')!.properties!.valueLabels!.labels, ['', 'A']);
    const cols = t.columns.map((c) => ({name: c.info.name, numeric: c.kind !== 'string'}));
    const figure = convertScript(t.scripts.find((s) => s.name === 'Figure 4')!, cols);
    expect(figure.platform, 'Graph Builder');
    expect(figure.viewers.length, 1);
  });

  test('JMP 15 and 18: p-value and date display formats', async () => {
    const das = await read('figshare-15044331-das28');
    const p = das.columns.find((c) => c.info.name === 'Prob>|t|')!;
    expect(p.kind, 'float');
    expect(p.descriptor[3], 0x61);
    expect(p.data.type === 'float' && Array.from(p.data.values).every((v) => isNaN(v) || (v >= 0 && v <= 1)), true);
    expectArray(das.columns.find((c) => c.info.name === 'Source Table')!.properties!.valueLabels!.labels.slice(0, 2),
      ['FIGS', 'Module 1']);
    const rat = await read('dataverse-jvinzu-rat');
    expect(rat.info.jmpVersion, '18.1.0');
    const date = rat.columns.find((c) => c.info.name === 'Date')!;
    expect(date.kind, 'date');
    expect(date.descriptor[3], 0x90);
  });

  test('JMP 15: a designed experiment, with its DOE scripts', async () => {
    const t = await read('zenodo-6149318-ivt-rna-dsd');
    expect(t.info.jmpVersion, '15.2.1');
    expect(t.info.nrows, 18);
    const props = (name: string) => t.columns.find((c) => c.info.name === name)!.properties!.properties;
    expectArray(props('DMSO').map((p) => p.text), ['Coding({0, 10})', 'Design Role(DesignRole(Continuous))',
      'Factor Changes(Easy)']);
    expect(props('Block')[0].text, 'RunsPerBlock(9)');
    expect(props('RNA Yield (ng/uL)')[0].text, 'Response Limits({Goal(Maximize), Importance(.), Show Limits(0)})');
    const cols = t.columns.map((c) => ({name: c.info.name, numeric: c.kind !== 'string'}));
    const analyses = t.scripts.filter((s) => s.name !== 'Source').map((s) => convertScript(s, cols));
    expectArray(analyses.map((a) => a.platform!), ['Fit Definitive Screening', 'DOE', 'DOE']);
    expect(analyses.every((a) => a.viewers.length === 0 && /isn't converted/.test(a.reason!)), true);
  });

  test('JMP 18: spec limits on predicted quality attributes, and their fits', async () => {
    const t = await read('mdpi-bioengineering-11-01089-wp1');
    expect(t.info.jmpVersion, '18.1.0');
    expect(t.info.nrows, 21);
    const spec = (name: string) => t.columns.find((c) => c.info.name === name)!.properties!.properties
      .find((p) => p.name === 'Spec Limits')!.text;
    expect(spec('Pred Formula CQA2 - WP1'), 'Spec Limits({USL(-0.51833593600000005), Show Limits(1)})');
    expect(spec('Pred Formula CQA3 - WP1'),
      'Spec Limits({LSL(-2.3376879380000002), USL(1.8659177870000001), Show Limits(1)})');
    const fits = convertScript(t.scripts.find((s) => s.name === 'Fit OLS')!, columnRefs(t.columns));
    expect(fits.platform, 'Fit Group');
    expectArray(fits.viewers.map((v) => `${v.options.yColumnName} by ${v.options.xColumnName}`),
      ['CQA1 - WP1 by Pred Formula CQA1 - WP1', 'CQA2 - WP1 by Pred Formula CQA2 - WP1',
        'CQA3 - WP1 by Pred Formula CQA3 - WP1']);
    const design = await read('mdpi-bioengineering-11-01089-wp1-3');
    const convert3 = (name: string) => convertScript(design.scripts.find((s) => s.name === name)!,
      columnRefs(design.columns));
    expect(convert3('Evaluate Design').platform, 'DOE');
    expect(convert3('Fit Least Squares').viewers.length, 3);
    expect(convert3('Generalized Regression').reason,
      'Fit Model\'s Generalized Regression personality isn\'t converted yet');
  });

  test('JMP 19: each fit gets only the prediction saved from its own model', async () => {
    const t = await read('btpr-42-e88508-phase-dependent-doe');
    expect(t.info.jmpVersion, '19.0.0');
    expect(t.info.nrows, 78);
    const convert = (name: string) => convertScript(t.scripts.find((s) => s.name === name)!, columnRefs(t.columns));
    // The saved formulas were fit without the outliers; APG and Oligo were named Acidic_PG1 and Man5 then.
    const plots = (name: string) => convert(name).viewers
      .map((v) => `${v.options.yColumnName} by ${v.options.xColumnName}`);
    expectArray(plots('Fit Group: WITHOUT outliers'),
      ['Viability by Pred Formula Viability_pd', 'Titer by Pred Formula Titer_pd', 'APG by Pred Formula APG_pd',
        'BPG by Pred Formula BPG_pd', 'LMW by Pred Formula LMW_pd', 'HMW by Pred Formula HMW_pd',
        'Oligo by Pred Formula Oligo_pd']);
    expect(convert('Fit Group: WITH outliers').reason,
      'Column "Log[Lactate]" is a transform the script computes, which isn\'t converted yet; ' +
      'The prediction formulas saved for "Viability", "Titer", "LMW", "HMW", "APG", "BPG", "Oligo" come from ' +
      'other models than these, and the model isn\'t refit');
  });
}, {node: true});

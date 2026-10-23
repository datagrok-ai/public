import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {JslNode, JslSyntaxError, columnNames, findArg, isCall, parseJsl, statements} from '../formats/jmp/jsl';
import {readJmp} from '../formats/jmp/jmp-parser';
import {JmpScript} from '../formats/jmp/jmp-scripts';
import {JmpColumnRef, JmpScriptConversion, convertScript} from '../formats/jmp/jsl-to-viewers';
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
    expect(c.notes.some((n) => /Separate analyses by age/.test(n)), true);
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
    v = gb('X(:age), Y(:height)', 'Line(X, Y, Legend(3))');
    expect(v.type, 'Line chart');
    expectArray(v.options.yColumnNames as string[], ['height']);
    v = gb('X(:sex), Y(:height)', 'Bar(X, Y, Legend(4), Summary Statistic("Median"))');
    expect(v.type, 'Bar chart');
    expect(v.options.valueAggrType, 'med');
    v = gb('X(:sex)', 'Bar(X, Legend(4))');
    expect(v.options.valueAggrType, 'count');
    v = gb('X(:sex), Y(:height)', 'Box Plot(X, Y, Legend(5))');
    expect(v.type, 'Box plot');
    expect(v.options.valueColumnName, 'height');
    v = gb('X(:weight)', 'Histogram(X, Legend(6))');
    expect(v.type, 'Histogram');
    v = gb('X(:sex)', 'Pie(X, Legend(7))');
    expect(v.type, 'Pie chart');
  });

  test('Graph Builder: what doesn\'t come across is named', async () => {
    const c = convert('Graph Builder(Variables(X(:height), Y(:weight), Group X(:sex), Wrap(:age)), ' +
      'Elements(Points(X, Y), Smoother(X, Y)))');
    only(c);
    expect(c.notes.some((n) => /Panels by sex/.test(n)), true);
    expect(c.notes.some((n) => /Panels by age/.test(n)), true);
    expect(c.notes.some((n) => /Smoother/.test(n)), true);
    expect(/Contour/.test(convert('Graph Builder(Variables(X(:height), Y(:weight)), Elements(Contour(X, Y)))').reason!),
      true);
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

  test('not converted, with a reason', async () => {
    expect(/Fit Model isn't converted/.test(convert('Fit Model(Y(:weight), Effects(:height), Run)').reason!), true);
    expect(/Column "nope"/.test(convert('Bivariate(Y(:weight), X(:nope))').reason!), true);
    expect(/Doesn't launch an analysis/.test(convert('Names Default To Here(1)').reason!), true);
    const unread = convertScript({name: 'broken', text: 'f(', error: 'Expected ")"'}, BIG_CLASS);
    expect(/couldn't be read/.test(unread.reason!), true);
  });
}, {node: true});

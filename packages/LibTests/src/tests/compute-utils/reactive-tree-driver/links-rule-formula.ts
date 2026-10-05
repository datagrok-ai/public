import {category, test, before, expect} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineLinkConfigurationInput} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {
  compileCheckFormulas, compileEffect, compileExpression, compileRuleFormulas, compileSource, compileValue, formulaOps,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/rule-formula';
import {annotationRules} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/rule-expansion';
import {expandChecks} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/checks';
import {evaluate} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-expressions';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrowsAsync} from '../../../test-utils';

const twoSteps = (links: PipelineLinkConfigurationInput<string | string[]>[]): PipelineConfiguration => ({
  id: 'pipeline1',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {id: 'step2', nqName: 'LibTests:TestMul2'},
  ],
  links,
});

const symbolicOps = new Set(['!', '!!', '==', '!=', '===', '!==', '>', '>=', '<', '<=', '+', '-', '*', '/', '%']);

// JSON Logic wraps a single argument in an array itself, and formulas always produce the array
function normalize(value: any): any {
  if (Array.isArray(value))
    return value.map(normalize);
  if (value == null || typeof value !== 'object')
    return value;
  const keys = Object.keys(value);
  const op = keys[0];
  if (keys.length === 1 && op !== 'literal' && (symbolicOps.has(op) || formulaOps.has(op)))
    return {[op]: (Array.isArray(value[op]) ? value[op] : [value[op]]).map(normalize)};
  return Object.fromEntries(keys.map((key) => [key, normalize(value[key])]));
}

const same = (actual: any, expected: any, prefix: string) =>
  expectDeepEqual(normalize(actual), normalize(expected), {prefix});

const present = {'!': {missing: ['value']}};

const expressionCases: [string, any][] = [
  ['gt(m, 0)', {'>': [{var: 'm'}, 0]}],
  ['eq(mode, "advanced")', {'==': [{var: 'mode'}, 'advanced']}],
  ['ne(m, 0)', {'!=': [{var: 'm'}, 0]}],
  ['same(a, b)', {'===': [{var: 'a'}, {var: 'b'}]}],
  ['notSame(a, b)', {'!==': [{var: 'a'}, {var: 'b'}]}],
  ['and(gte(a, 1), lte(a, 2, 3), lt(a1, b1))', {and: [{'>=': [{var: 'a'}, 1]}, {'<=': [{var: 'a'}, 2, 3]},
    {'<': [{var: 'a1'}, {var: 'b1'}]}]}],
  ['add(1, sub(m), mul(m, 2), div(m, 2), mod(m, 2))', {'+': [1, {'-': [{var: 'm'}]}, {'*': [{var: 'm'}, 2]},
    {'/': [{var: 'm'}, 2]}, {'%': [{var: 'm'}, 2]}]}],
  ['not(v)', {'!': {var: 'v'}}],
  ['bool(row)', {'!!': {var: 'row'}}],
  ['not(missing(table))', {'!': {missing: ['table']}}],
  ['missing("a", b)', {missing: ['a', 'b']}],
  ['missing_some(1, [m, $all.m])', {missing_some: [1, ['m', '$all.m']]}],
  ['bool(var(dflt, 1))', {'!!': [{var: ['dflt', 1]}]}],
  ['len($all.m)', {len: {var: '$all.m'}}],
  ['col.name', {var: 'col.name'}],
  ['$all.m.0', {var: '$all.m.0'}],
  ['and(not(missing(column)), not(in(column, columns)))',
    {and: [{'!': {missing: ['column']}}, {'!': {in: [{var: 'column'}, {var: 'columns'}]}}]}],
  ['in(key, column(presets, "preset"))', {in: [{var: 'key'}, {column: [{var: 'presets'}, 'preset']}]}],
  ['if(gt(m, 0), ["x", "y"], ["z"])', {if: [{'>': [{var: 'm'}, 0]}, ['x', 'y'], ['z']]}],
  ['cat("Must exceed ambient ", amb)', {cat: ['Must exceed ambient ', {var: 'amb'}]}],
  ['bool(columnsMissing(df, [["smiles", "Molecule"], "activity"]))',
    {'!!': {columnsMissing: [{var: 'df'}, [['smiles', 'Molecule'], 'activity']]}}],
  ['columns(df, "numerical")', {columns: [{var: 'df'}, 'numerical']}],
  ['row(df, "s", missing)', {row: [{var: 'df'}, 's', {var: 'missing'}]}],
  ['regex(value, "^[0-9]{4}$", "i")', {regex: [{var: 'value'}, '^[0-9]{4}$', 'i']}],
  ['eq(script(`type == "ICE"`), false)', {'==': [{script: 'type == "ICE"'}, false]}],
  ['scriptVerdict(`startsWith(value, "12")`)', {scriptVerdict: 'startsWith(value, "12")'}],
  ['gt(nulls(value), 0)', {'>': [{nulls: {var: 'value'}}, 0]}],
  ['map(verdicts, message)', {map: [{var: 'verdicts'}, {var: 'message'}]}],
  ['filter(list, gt($it, 0))', {filter: [{var: 'list'}, {'>': [{var: ''}, 0]}]}],
  ['map(lists, filter(xs, gt($it, 0)))', {map: [{var: 'lists'}, {filter: [{var: 'xs'}, {'>': [{var: ''}, 0]}]}]}],
  ['[.5, 1., -.5, 0.25e1]', [0.5, 1, -0.5, 2.5]],
  ['regex(code, `^\\d{4}$`)', {regex: [{var: 'code'}, '^\\d{4}$']}],
  ['"tab\\t, quote \\", backslash \\\\"', 'tab\t, quote ", backslash \\'],
  ['reduce(list, add(current, accumulator), 0)',
    {reduce: [{var: 'list'}, {'+': [{var: 'current'}, {var: 'accumulator'}]}, 0]}],
  ['[true, false, null, -1.5e2, "a\\"b\\n", `raw \\n "q"`]', [true, false, null, -150, 'a"b\n', 'raw \\n "q"']],
  ['obj(foo: 1)', {literal: {foo: 1}}],
  ['obj(var: "not an alias")', {literal: {var: 'not an alias'}}],
  ['if(gt(m, 0), obj(var: "kept"), "no")', {if: [{'>': [{var: 'm'}, 0]}, {literal: {var: 'kept'}}, 'no']}],
  ['[obj(a: 1), obj(b: [1, obj(c: "d")])]', [{literal: {a: 1}}, {literal: {b: [1, {c: 'd'}]}}]],
  ['  and( a ,\n b )  ', {and: [{var: 'a'}, {var: 'b'}]}],
  ['=gt(m, 0)', {'>': [{var: 'm'}, 0]}],
];

const effectCases: [string, any][] = [
  ['hide(t1)', {effect: 'hide', targets: 't1'}],
  ['show([tol, col])', {effect: 'show', targets: ['tol', 'col']}],
  ['show(u, when: gt(m, 0))', {effect: 'show', targets: 'u', when: {'>': [{var: 'm'}, 0]}}],
  ['items(col, columns(df, "numerical"))',
    {effect: 'items', targets: 'col', items: {columns: [{var: 'df'}, 'numerical']}}],
  ['meta(t, units: "K", twice: mul(m, 2))',
    {effect: 'meta', targets: 't', meta: {units: 'K', twice: {'*': [{var: 'm'}, 2]}}}],
  ['meta(t, cfg: obj(foo: 1), when: m)',
    {effect: 'meta', targets: 't', meta: {cfg: {literal: {foo: 1}}}, when: {var: 'm'}}],
  ['error(init, cat("Must exceed ambient ", amb))',
    {effect: 'error', targets: 'init', message: {cat: ['Must exceed ambient ', {var: 'amb'}]}}],
  ['warning(u, "negative", when: lt(a1, 0))',
    {effect: 'warning', targets: 'u', message: 'negative', when: {'<': [{var: 'a1'}, 0]}}],
  ['notification(t, "done")', {effect: 'notification', targets: 't', message: 'done'}],
  ['verdicts(t, v)', {effect: 'verdicts', targets: 't', source: 'v'}],
  ['set(init, amb, restriction: "restricted")',
    {effect: 'set', targets: 'init', value: {var: 'amb'}, restriction: 'restricted'}],
  ['set(t, 42)', {effect: 'set', targets: 't', value: 42}],
  ['clear(c, when: and(not(missing(column)), not(in(column, columns))))', {effect: 'clear', targets: 'c',
    when: {and: [{'!': {missing: ['column']}}, {'!': {in: [{var: 'column'}, {var: 'columns'}]}}]}}],
  ['clear(u, restriction: "info")', {effect: 'clear', targets: 'u', restriction: 'info'}],
  ['assign(row(presets, "preset", key), restriction: "restricted", when: in(key, column(presets, "preset")))',
    {effect: 'assign', values: {row: [{var: 'presets'}, 'preset', {var: 'key'}]}, restriction: 'restricted',
      when: {in: [{var: 'key'}, {column: [{var: 'presets'}, 'preset']}]}}],
  ['assign(row, targets: a, ignoreCase: true)',
    {effect: 'assign', values: {var: 'row'}, targets: 'a', ignoreCase: true}],
  ['=hide(t)', {effect: 'hide', targets: 't'}],
];

const sourceCases: [string, any][] = [
  ['file("System:AppData/Pkg/presets.csv")', {file: 'System:AppData/Pkg/presets.csv'}],
  ['table("preset,a\\nfast,1")', {table: 'preset,a\nfast,1'}],
  ['table("p;a\\nf;1,5", delimiter: ";", decimalSeparator: ",")',
    {table: {csv: 'p;a\nf;1,5', options: {delimiter: ';', decimalSeparator: ','}}}],
  ['func("LibTests:TestPresets")', {func: {name: 'LibTests:TestPresets'}}],
  ['func("LibTests:TestAdd2", a: x, b: 5)', {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'x'}, b: 5}}}],
  ['func("OpenFile", fullPath: "System:AppData/Pkg/presets.csv")',
    {func: {name: 'OpenFile', args: {fullPath: 'System:AppData/Pkg/presets.csv'}}}],
  ['query("System:Datagrok", `select login from users where login = @login`, login: login)',
    {query: {connection: 'System:Datagrok', sql: 'select login from users where login = @login',
      args: {login: {var: 'login'}}}}],
  ['validators(tol, names: ["Pkg:checkTolerance"])', {validators: {input: 'tol', names: ['Pkg:checkTolerance']}}],
  ['validators(tol)', {validators: {input: 'tol'}}],
  ['choices(region)', {choices: {input: 'region'}}],
];

const errorCases: [(text: string) => any, string, RegExp][] = [
  [compileEffect, 'gt(m, 0)', /Expected an effect call at the root, got an expression op gt at column 1/],
  [compileEffect, 'tol', /got a reference/],
  [compileEffect, '"too small"', /got a literal/],
  [compileEffect, 'file("x")', /got a source call file/],
  [compileSource, 'hide(t)', /Expected a source call at the root, got an effect call hide/],
  [compileExpression, 'hide(t)', /hide is an effect call and cannot be used inside an expression/],
  [compileExpression, 'and(m, choices(x))',
    /choices is a source call and cannot be used inside an expression at column 8/],
  [compileEffect, 'set(1)', /set needs value \(positional arguments: targets, value\)/],
  [compileEffect, 'hide(t, u)', /hide takes 1 positional argument\(s\) \(targets\)/],
  [compileEffect, 'hide(t, color: "red")', /hide has no option color \(expected when\) at column 16/],
  [compileEffect, 'set(t, 1, restriction: "strict")', /Unknown restriction strict at column 24/],
  [compileEffect, 'set(t, 1, restriction: r)', /restriction must be a string/],
  [compileEffect, 'assign(row, ignoreCase: "yes")', /ignoreCase must be true or false/],
  [compileEffect, 'hide(t.x)', /Targets are an output alias or a list of them/],
  [compileEffect, 'verdicts(t, "v")', /source is an alias/],
  [compileEffect, 'hide(t, when: gt(m, 0), when: m)', /Duplicate argument when/],
  [compileEffect, 'set(t, value: 1)', /set needs value/],
  [compileEffect, 'set(restriction: "info", t, 1)', /Positional argument after a named one/],
  [compileExpression, 'colums(df)', /Unknown operation colums at column 1/],
  [compileExpression, 'columns(df, kind: "numerical")', /columns takes positional arguments only at column 19/],
  [compileExpression, 'obj(n: len(df))', /obj values must be constants \(n\) at column 8/],
  [compileExpression, 'obj(1)', /obj takes named arguments only/],
  [compileExpression, 'missing(1)', /missing takes aliases/],
  [compileExpression, 'map(list, $it.name)', /Read the element's fields by name, name/],
  [compileExpression, 'gt($it, 0)', /\$it is the element inside map, filter, all, some and none at column 4/],
  [compileExpression, 'reduce(list, add($it, 1), 0)', /\$it is the element inside map/],
  [compileExpression, 'regex(code, "^\\d{4}$")', /Unknown escape \\d, write patterns .* in backticks at column 15/],
  [compileExpression, 'gt(m, .)', /Syntax error/],
  [compileExpression, 'gt(m, 0', /Syntax error at column 3/],
  [compileExpression, 'too small', /Syntax error at column 5/],
  [compileExpression, 'm > 0', /Syntax error/],
  [compileExpression, 'eq(mode, \'advanced\')', /Syntax error/],
  [compileSource, 'func(a: x)', /func needs name/],
  [compileSource, 'func(name)', /name must be a string/],
  [compileSource, 'query(`select 1`)', /query needs sql \(positional arguments: connection, sql\)/],
  [compileSource, 'file("x", cache: true)', /file has no option cache/],
  [compileSource, 'table("a", delimiter: d)', /delimiter must be a constant/],
  [compileSource, 'validators(tol, names: [n])', /names must be a constant/],
  [compileSource, 'validators(tol, names: "Pkg:f")', /names must be a list of function names/],
  [compileSource, 'validators(tol, call: c)', /validators has no option call \(expected names\)/],
  [compileSource, 'choices("region")', /input is an alias/],
  [compileExpression, 'var()', /var takes an alias and an optional default at column 1/],
  [compileExpression, 'var(a, 1, 2)', /var takes an alias and an optional default/],
  [compileExpression, 'missing_some(1, [a], c)', /missing_some takes a count and a list of aliases/],
  [compileExpression, 'map(xs, var($it))', /var takes aliases, read the element as \$it and its fields by name/],
  [compileExpression, 'missing($it.x)', /missing takes aliases, read the element as \$it/],
  [compileExpression, 'toString(1)', /Unknown operation toString/],
  [compileEffect, 'toString(t)', /Expected an effect call at the root, got an unknown call toString/],
  [compileSource, 'constructor("x")', /Expected a source call at the root, got an unknown call constructor/],
  [compileEffect, 'set(t, 1, restriction: "constructor")', /Unknown restriction constructor/],
  [compileEffect, 'hied(t)', /got an unknown call hied/],
  [compileSource, 'js(x)', /got an unknown call js/],
];

category('ComputeUtils: Driver rule formulas', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Expressions translate to JSON Logic', async () => {
    for (const [text, expected] of expressionCases)
      same(compileExpression(text), expected, text);
  });

  test('Expressions keep their meaning when evaluated', async () => {
    const ctx = {$all: {m: [2, 3]}, m: 2, list: [-1, 2, 3], s: 'abc', col: {name: 'x'}};
    const cases: [string, any][] = [
      ['add(m, mul(m, 3))', 8], ['if(gt(m, 1), "big", "small")', 'big'], ['len($all.m)', 2],
      ['filter(list, gt($it, 0))', [2, 3]], ['reduce(list, add(current, accumulator), 0)', 4],
      ['regex("1234", `^\\d{4}$`)', true], ['regex("dddd", `^\\d{4}$`)', false],
      ['col.name', 'x'], ['missing(m, gone)', ['gone']], ['obj(a: 1)', {a: 1}], ['lt(1, m, 3)', true],
      ['substr(s, 1)', 'bc'], ['merge([1], [2])', [1, 2]], ['bool(var(gone, 1))', true],
    ];
    for (const [text, expected] of cases)
      expectDeepEqual(evaluate(compileExpression(text), ctx), expected, {prefix: text});
  });

  test('Literals work inside element arguments', async () => {
    const ctx = {$all: {}, xs: [1, -2, 3]};
    expectDeepEqual(evaluate({map: [{var: 'xs'}, {literal: {a: 1}}]} as any, ctx), [{a: 1}, {a: 1}, {a: 1}],
      {prefix: 'object form'});
    expectDeepEqual(evaluate(compileExpression('map(xs, if(gt($it, 0), obj(pos: true), obj(pos: false)))'), ctx),
      [{pos: true}, {pos: false}, {pos: true}], {prefix: 'formula'});
    expectDeepEqual(evaluate(compileExpression('map(xs, [obj(n: 1), $it])'), ctx),
      [[{n: 1}, 1], [{n: 1}, -2], [{n: 1}, 3]], {prefix: 'nested in a list'});
    expectDeepEqual(evaluate(compileExpression('[obj(a: 1), if(true, obj(b: 2), 0)]'), ctx), [{a: 1}, {b: 2}],
      {prefix: 'outside elements'});
  });

  test('Effect calls translate to effect objects', async () => {
    for (const [text, expected] of effectCases)
      same(compileEffect(text), expected, text);
  });

  test('Source calls translate to source objects', async () => {
    for (const [text, expected] of sourceCases)
      same(compileSource(text), expected, text);
  });

  test('Every op a formula names is registered', async () => {
    const args: Record<string, any[]> = {
      columnsMissing: [null, []], column: [null, 'x'], row: [null, 'x', null], regex: ['a', 'a'],
      columnIs: [null, 'int'],
    };
    for (const op of formulaOps) {
      try {
        evaluate({[op]: args[op] ?? [[]]} as any, {$all: {}});
      } catch (e) {
        if (String(e).includes('Unrecognized operation'))
          throw new Error(`${op} is not a registered operation`);
      }
    }
  });

  test('Formula errors name the problem and the column', async () => {
    for (const [compile, text, pattern] of errorCases) {
      let message = '';
      try {
        compile(text);
      } catch (e) {
        message = e instanceof Error ? e.message : String(e);
      }
      expect(pattern.test(message), true, `${text}: got "${message}"`);
    }
  });

  test('Annotation rules have formula forms', async () => {
    const io = [
      {id: 'region', type: 'string', nullable: false, direction: 'input' as const, dynamicChoices: {propagate: true}},
      {id: 'n', type: 'int', nullable: false, direction: 'input' as const},
    ];
    const [choices, lookup] = annotationRules('LibTests:Fake', io);
    same(choices.sources!.region_choices, compileSource('choices(region)'), 'choices source');
    same(choices.effects, [
      compileEffect('items(region_target, region_choices.items, when: bool(region_choices))'),
      compileEffect('warning(region_target, region_choices.rowErrors, when: bool(region_choices))'),
    ], 'choices effects');
    same(lookup.effects, [compileEffect(
      'assign(region_choices.row, ignoreCase: true, restriction: "restricted", when: bool(region_choices.row))')],
    'lookup effects');
  });

  test('Check templates have formula forms', async () => {
    const params = (options: any) => expandChecks(options).map((check) => check.params);
    same(params({nullable: false})[0].when, compileExpression('missing(value)'), 'required');
    same(params({min: 0})[0].when, compileExpression('and(not(missing(value)), lt(value, 0))'), 'min');
    same(params({validator: '/^[0-9]{4}$/i'})[0].when,
      compileExpression('and(not(missing(value)), not(regex(value, "^[0-9]{4}$", "i")))'), 'regex');
    same(params({validator: 'bar > 3'})[0].when,
      compileExpression('and(not(missing(value)), bool(scriptVerdict("bar > 3")))'), 'expression');
    const visible = params({visible: 'type == "ICE"'})[0];
    same(visible.when, compileExpression('eq(script(`type == "ICE"`), false)'), 'visible when');
    same(visible.effects, [compileEffect('hide([$target])')], 'visible effect');
    const validators = params({validators: ['Chem:isSmiles']})[0];
    same(validators.sources!.$verdicts, compileSource('validators(value, names: ["Chem:isSmiles"])'), 'validators');
    same(validators.effects, [compileEffect('verdicts([$target], $verdicts)')], 'verdicts');
    same(params({choices: ['fast', 'exact']})[0].when,
      compileExpression('and(not(missing(value)), not(in(value, ["fast", "exact"])))'), 'choices');
    same(params({table: 'x'})[0].when, compileExpression(
      'and(not(missing(value)), and(not(missing($table)), not(in(value.name, columns($table)))))'), 'table');
    same(params({allowNulls: false})[0].when,
      compileExpression('and(not(missing(value)), gt(nulls(value), 0))'), 'nulls');
    same(compileExpression('not(missing(value))'), present, 'present');
  });

  test('Calls always emit argument arrays', async () => {
    expectDeepEqual(compileExpression('not(v)'), {'!': [{var: 'v'}]});
    expectDeepEqual(compileExpression('var(a, 1)'), {var: ['a', 1]});
    expectDeepEqual(compileExpression('len(df)'), {len: [{var: 'df'}]});
  });

  test('Names of object properties are ordinary names', async () => {
    same(compileEffect('meta(t, toString: 1, constructor: "c")'),
      {effect: 'meta', targets: 't', meta: {toString: 1, constructor: 'c'}}, 'meta keys');
    expectDeepEqual(compileExpression('obj(valueOf: 1, hasOwnProperty: 2)'),
      {literal: {valueOf: 1, hasOwnProperty: 2}});
    expectDeepEqual(compileExpression('[``, ""]'), ['', '']);
    same(compileSource('file(``)'), {file: ''}, 'empty backtick string');
  });

  test('Malformed objects pass through to the existing checks', async () => {
    const malformed: any = {
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      effects: [null, {effect: 'meta', targets: 't', meta: 'units'}],
      sources: {a: {func: {name: 'Pkg:F', args: ['x']}}, b: 7},
    };
    const compiled = compileRuleFormulas(malformed);
    expectDeepEqual(compiled.effects[0], null);
    expectDeepEqual((compiled.effects[1] as any).meta, 'units');
    expectDeepEqual(compiled.sources, malformed.sources);
    expectDeepEqual(compileRuleFormulas({...malformed, effects: {}, sources: 'x'} as any).sources, 'x');
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      sources: {v: {func: {name: 'LibTests:TestAdd2', args: [{var: 'x'}] as any}}}, effects: ['error(t, "m")'],
    }])), /source v args must map parameters to expressions/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a', effects: {} as any,
    }])), /effects list is empty/);
  });

  test('Strings compile by position', async () => {
    const rule = compileRuleFormulas({
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      when: '=gt(m, 0)',
      sources: {
        f: {func: {name: 'Pkg:F', args: {path: 'System:x.csv', n: '=add(m, 1)', nested: ['=a']}}},
        q: {query: {connection: 'Pkg:Conn', sql: 'select 1', args: {n: '=m'}}},
        g: 'func("Pkg:G")',
      },
      effects: [
        {effect: 'error', targets: 't', message: 'too small', when: 'gt(m, 1)'},
        {effect: 'warning', targets: 't', message: '=cat("bad ", m)'},
        {effect: 'items', targets: 't', items: ['=a', 'b']},
        {effect: 'meta', targets: 't', meta: {units: 'K', twice: '=mul(m, 2)', raw: {cat: ['=', {var: 'm'}]}}},
        {effect: 'set', targets: 't', value: '=m'},
        {effect: 'assign', values: '=row'},
        'hide(t)',
      ],
    });
    same(rule.when, {'>': [{var: 'm'}, 0]}, 'when');
    same(rule.sources, {
      f: {func: {name: 'Pkg:F', args: {path: 'System:x.csv', n: {'+': [{var: 'm'}, 1]}, nested: ['=a']}}},
      q: {query: {connection: 'Pkg:Conn', sql: 'select 1', args: {n: {var: 'm'}}}},
      g: {func: {name: 'Pkg:G'}},
    }, 'sources');
    same(rule.effects, [
      {effect: 'error', targets: 't', message: 'too small', when: {'>': [{var: 'm'}, 1]}},
      {effect: 'warning', targets: 't', message: {cat: ['bad ', {var: 'm'}]}},
      {effect: 'items', targets: 't', items: ['=a', 'b']},
      {effect: 'meta', targets: 't', meta: {units: 'K', twice: {'*': [{var: 'm'}, 2]}, raw: {cat: ['=', {var: 'm'}]}}},
      {effect: 'set', targets: 't', value: {var: 'm'}},
      {effect: 'assign', values: {var: 'row'}},
      {effect: 'hide', targets: 't'},
    ], 'effects');
    const check = compileCheckFormulas({
      id: 'c', type: 'check', io: 'step1/a', check: {min: 0}, when: 'gt(value, 1)', message: '=cat("min ", value)',
    });
    same(check.when, {'>': [{var: 'value'}, 1]}, 'check when');
    same(check.message, {cat: ['min ', {var: 'value'}]}, 'check message');
    const literal = compileCheckFormulas({id: 'c', type: 'check', io: 'step1/a', check: {min: 0}, message: 'Too low'});
    expectDeepEqual(literal.message, 'Too low');
    expectDeepEqual(compileValue('plain'), 'plain');
    expectDeepEqual(compileValue({var: 'x'} as any), {var: 'x'});
  });

  test('Errors carry the position in the config', async () => {
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a', effects: ['hide(t)', 'hied(t)'],
    }])), /Rule r: effects\[1\]: Expected an effect call at the root, got an unknown call hied/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a', sources: {v: 'hide(t)'}, effects: ['hide(t)'],
    }])), /Rule r: sources.v: Expected a source call/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a', when: 'gt(m, 0', effects: ['hide(t)'],
    }])), /Rule r: when: Syntax error at column 3/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      effects: [{effect: 'set', targets: 't', value: '=gt(m'}],
    }])), /Rule r: effects\[0\].value: Syntax error/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'c', type: 'check', io: 'step1/a', check: {min: 0}, when: 'gt(value 0)',
    }])), /Check c: when: Syntax error/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a', when: 'advanced', effects: ['hide(t)'],
    }])), /unknown input alias advanced/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a', effects: ['hide(zzz)'],
    }])), /unknown output alias zzz/);
  });

  test('Element arguments cannot read rule aliases', async () => {
    const rule = (extra: any) => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: ['xs:step1/a', 'limit:step1/b'], to: 't:step2/a', effects: ['hide(t)'], ...extra,
    }]));
    const message = /Rule r: limit inside map.* is a field of the element, not the alias limit; .* var\("limit"\)/;
    await expectThrowsAsync(() => rule({when: 'some(xs, gt($it, limit))'}), message);
    await expectThrowsAsync(() => rule({when: 'some(xs, some(ys, gt($it, limit.x)))'}), message);
    await expectThrowsAsync(() => rule({when: {some: [{var: 'xs'}, {'>': [{var: ''}, {var: 'limit'}]}]}}), message);
    await expectThrowsAsync(() => rule({effects: ['set(t, reduce(xs, add(current, limit), 0))']}), message);
    await expectThrowsAsync(() => rule({sources: {lim: 'func("Pkg:F")'}, when: 'map(xs, lim)'}),
      /lim inside map.* is a field of the element, not the alias lim/);
    const pconf = await rule({when: 'some(xs, and(gt($it, var("limit")), not(message)))'});
    same(pconf.links![0].params!.when, {some: [{var: 'xs'}, {and: [{'>': [{var: ''}, {var: ['limit']}]},
      {'!': [{var: 'message'}]}]}]}, 'var("limit") and other fields');
    await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: ['xs:step1/a', 'current:step1/b'], to: 't:step2/a',
      when: 'gt(reduce(xs, add(current, accumulator), 0), current)', effects: ['hide(t)'],
    }]));
  });

  test('Aliases named like keywords are read with var', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: ['it:step1/a', 'null:step1/b'], to: 't:step2/a',
      when: 'and(gt(it, 0), bool(var("null")), not(null))', effects: ['hide(t)'],
    }]));
    same(pconf.links![0].params!.when,
      {and: [{'>': [{var: 'it'}, 0]}, {'!!': [{var: ['null']}]}, {'!': [null]}]}, 'it and null aliases');
  });
  test('Formula and object rules process to the same links', async () => {
    const objects = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: ['m:step1/a', 'n:step1/b'], to: ['t1:step2/a', 't2:step2/b'],
      when: {'>': [{var: 'm'}, 0]},
      sources: {sum: {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'm'}, b: 5}}}},
      effects: [
        {effect: 'hide', targets: 't1', when: {'!': {var: 'n'}}},
        {effect: 'items', targets: ['t2'], items: ['x', 'y']},
        {effect: 'error', targets: 't1', message: {cat: ['bad ', {var: '$all.n'}]}},
        {effect: 'set', targets: 't2', value: {var: 'sum'}, restriction: 'restricted'},
      ],
    }, {
      id: 'c', type: 'check', io: 'step1/a', check: {min: 0}, when: {'>': [{var: 'value'}, 1]},
      message: {cat: ['min ', {var: 'value'}]},
    }]));
    const formulas = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: ['m:step1/a', 'n:step1/b'], to: ['t1:step2/a', 't2:step2/b'],
      when: 'gt(m, 0)',
      sources: {sum: 'func("LibTests:TestAdd2", a: m, b: 5)'},
      effects: [
        'hide(t1, when: not(n))',
        {effect: 'items', targets: ['t2'], items: ['x', 'y']},
        'error(t1, cat("bad ", $all.n))',
        'set(t2, sum, restriction: "restricted")',
      ],
    }, {
      id: 'c', type: 'check', io: 'step1/a', check: {min: 0}, when: 'gt(value, 1)', message: '=cat("min ", value)',
    }]));
    expectDeepEqual(formulas.links!.map((l) => l.id), objects.links!.map((l) => l.id));
    formulas.links!.forEach((link, idx) => same(link.params, objects.links![idx].params, link.id));
  });

  test('Docs rules match their object forms', async () => {
    const pairs: [any, any][] = [[{
      when: 'eq(mode, "advanced")',
      effects: ['show([tol, col])', 'items(col, columns(df, "numerical"))', 'clear(col)',
        'error(init, cat("Must exceed ambient ", amb))'],
    }, {
      when: {'==': [{var: 'mode'}, 'advanced']},
      effects: [
        {effect: 'show', targets: ['tol', 'col']},
        {effect: 'items', targets: 'col', items: {columns: [{var: 'df'}, 'numerical']}},
        {effect: 'clear', targets: 'col'},
        {effect: 'error', targets: 'init', message: {cat: ['Must exceed ambient ', {var: 'amb'}]}},
      ],
    }], [{
      sources: {v: 'validators(tol, names: ["Pkg:checkTolerance"])'},
      when: 'eq(mode, "advanced")',
      effects: ['show(t)', 'verdicts(t, v)', 'warning(t, "Custom tolerance is slower")',
        'error(init, cat("Must exceed ambient ", amb))', 'set(init, amb, restriction: "restricted")'],
    }, {
      sources: {v: {validators: {input: 'tol', names: ['Pkg:checkTolerance']}}},
      when: {'==': [{var: 'mode'}, 'advanced']},
      effects: [
        {effect: 'show', targets: 't'},
        {effect: 'verdicts', targets: 't', source: 'v'},
        {effect: 'warning', targets: 't', message: 'Custom tolerance is slower'},
        {effect: 'error', targets: 'init', message: {cat: ['Must exceed ambient ', {var: 'amb'}]}},
        {effect: 'set', targets: 'init', value: {var: 'amb'}, restriction: 'restricted'},
      ],
    }], [{
      when: 'gt(threshold, 0.5)',
      effects: ['warning(n, "Consider more iterations for a strict threshold")'],
    }, {
      when: {'>': [{var: 'threshold'}, 0.5]},
      effects: [{effect: 'warning', targets: 'n', message: 'Consider more iterations for a strict threshold'}],
    }], [{
      sources: {presets: 'file("System:AppData/Pkg/presets.csv")'},
      effects: ['items(k, column(presets, "preset"))',
        'assign(row(presets, "preset", key), restriction: "restricted", when: in(key, column(presets, "preset")))'],
    }, {
      sources: {presets: {file: 'System:AppData/Pkg/presets.csv'}},
      effects: [
        {effect: 'items', targets: 'k', items: {column: [{var: 'presets'}, 'preset']}},
        {effect: 'assign', values: {row: [{var: 'presets'}, 'preset', {var: 'key'}]}, restriction: 'restricted',
          when: {in: [{var: 'key'}, {column: [{var: 'presets'}, 'preset']}]}},
      ],
    }]];
    pairs.forEach(([formula, object], idx) => {
      const compiled = compileRuleFormulas({id: 'r', type: 'rule', from: [], to: [], ...formula});
      same({when: compiled.when, sources: compiled.sources, effects: compiled.effects},
        {when: object.when, sources: object.sources, effects: object.effects}, `docs rule ${idx}`);
    });
  });

  test('Formula rules run like object rules', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: ['t:step2/a', 'u:step2/b'],
      when: 'gt(m, 0)',
      effects: ['hide(t)', 'set(u, mul(m, 2), when: gt(m, 1))'],
    }]));
    testScheduler.run(({expectObservable, cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inStore = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as any;
      cold('-a').subscribe(() => inStore.setState('a', 1));
      cold('--a').subscribe(() => inStore.setState('a', 2));
      cold('---a').subscribe(() => inStore.setState('a', -1));
      expectObservable(outBridge.meta.a).toBe('abbc', {a: undefined, b: {hidden: true}, c: {hidden: false}});
      expectObservable(outBridge.getStateChanges('b')).toBe('a-b', {a: undefined, b: 4});
    });
  });

  test('Check formulas run like object checks', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'c', type: 'check', io: 'step1/a', check: {max: 10}, when: 'lt(value, 100)',
      message: '=cat("At most 10, got ", value)', debounce: 0,
    }]));
    const snapshots: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const snap = () => snapshots.push(node.validationInfo$.value);
      cold('-a').subscribe(() => node.getStateStore().setState('a', 12));
      cold('--a').subscribe(snap);
      cold('---a').subscribe(() => node.getStateStore().setState('a', 500));
      cold('----a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [
      {a: {errors: [{description: 'At most 10, got 12'}], warnings: [], notifications: []}},
      {},
    ]);
  });
});

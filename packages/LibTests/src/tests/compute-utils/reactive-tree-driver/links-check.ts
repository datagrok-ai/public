import {category, test, before, expect} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {expandChecks, parseRegexLiteral, validateCheckOptions} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/checks';
import {evaluate} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-expressions';
import {resolveSources} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-sources';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineLinkConfigurationInput} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import * as DG from 'datagrok-api/dg';
import {createTestScheduler, expectThrowsAsync} from '../../../test-utils';

const ANNOTATED = 'LibTests:TestAnnotatedInputs';
const NAMED = 'LibTests:TestNamedValidators';

const annotatedStep = (links: PipelineLinkConfigurationInput<string | string[]>[] = []): PipelineConfiguration => ({
  id: 'pipeline1',
  type: 'static',
  steps: [{id: 's', nqName: ANNOTATED}],
  links,
});

const twoSteps = (links: PipelineLinkConfigurationInput<string | string[]>[]): PipelineConfiguration => ({
  id: 'pipeline1',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {id: 'step2', nqName: 'LibTests:TestMul2'},
  ],
  links,
});

const errors = (...descriptions: string[]) =>
  ({errors: descriptions.map((description) => ({description})), warnings: [], notifications: []});
const warnings = (...descriptions: string[]) =>
  ({errors: [], warnings: descriptions.map((description) => ({description})), notifications: []});

function makeTable(withNulls = false) {
  const mol = DG.Column.fromList('string', 'mol', ['C', 'CC']);
  mol.semType = 'Molecule';
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('int', 'x', withNulls ? [1, null] : [1, 2]),
    DG.Column.fromList('string', 's', ['a', 'b']),
    mol,
  ]);
}

// the config-side twin of the TestAnnotatedInputs annotations
const annotationsAsChecks: PipelineLinkConfigurationInput<string | string[]>[] = [
  {id: 'c', type: 'check', io: 's/c', check: {nullable: false}},
  {id: 'v', type: 'check', io: 's/v', check: {nullable: false, min: 0, max: 10}},
  {id: 'code', type: 'check', io: 's/code', check: {nullable: false, validator: '/^[0-9]{4}$/i'}},
  {id: 'mode', type: 'check', io: 's/mode', check: {nullable: false, choices: ['fast', 'exact']}},
  {id: 'df', type: 'check', io: 's/df', check: {nullable: false}},
  {id: 'col', type: 'check', io: 's/col',
    check: {nullable: false, type: 'numerical', table: 's/df', allowNulls: false}},
  {id: 'mol', type: 'check', io: 's/mol', check: {nullable: false, semType: 'Molecule', table: 's/df'}},
];

category('ComputeUtils: Driver links check', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Annotations parse into checks', async () => {
    const pconf: any = await getProcessedConfig(annotatedStep());
    const io = Object.fromEntries(pconf.steps[0].io.map((item: any) => [item.id, item]));
    expect(io.a.nullable, true, 'nullable: true');
    expect(io.b.nullable, true, 'optional: true');
    expect(io.c.nullable, false, 'plain input');
    expect(io.res.nullable, false, 'output');
    expectDeepEqual(io.v.checks, {min: 0, max: 10});
    expectDeepEqual(io.code.checks, {validator: '/^[0-9]{4}$/i'});
    expectDeepEqual(io.mode.checks, {choices: ['fast', 'exact']});
    expectDeepEqual(io.col.checks, {type: 'numerical', table: 'df', allowNulls: false});
    expectDeepEqual(io.mol.checks, {semType: 'Molecule', table: 'df'});
    expect(io.c.checks === undefined, true);
    expect(io.df.checks === undefined, true);
    expect(io.a.checks === undefined, true);
  });

  test('Expand checks into validator params', async () => {
    const expanded = expandChecks({
      nullable: false, min: 0, max: 10, validator: '/^a/i', choices: ['x', 'y'],
      type: 'numerical', semType: 'Molecule', table: 'df', allowNulls: false,
    });
    expectDeepEqual(expanded.map((item) => item.key),
      ['required', 'min', 'max', 'validator', 'choices', 'type', 'semType', 'table', 'allowNulls']);
    expectDeepEqual(expanded.map((item) => item.needsTable),
      [false, false, false, false, false, false, false, true, false]);
    const present = {'!': {missing: ['value']}};
    expectDeepEqual(expanded[0].params, {
      when: {missing: ['value']},
      effects: [{effect: 'error', targets: ['target'], message: 'Missing value'}],
    });
    expectDeepEqual(expanded[1].params, {
      when: {and: [present, {'<': [{var: 'value'}, 0]}]},
      effects: [{effect: 'error', targets: ['target'], message: 'Must be at least 0'}],
    });
    expectDeepEqual(expanded[3].params.when, {and: [present, {'!': {regex: [{var: 'value'}, '^a', 'i']}}]});
    expectDeepEqual(expanded[4].params.when, {and: [present, {'!': {in: [{var: 'value'}, ['x', 'y']]}}]});
    expectDeepEqual(expanded[5].params.when, {and: [present, {'!': {columnIs: [{var: 'value'}, 'numerical']}}]});
    expectDeepEqual(expanded[7].params.when, {and: [present, {and: [
      {'!': {missing: ['table']}},
      {'!': {in: [{var: 'value.name'}, {columns: [{var: 'table'}]}]}},
    ]}]});
    expectDeepEqual(expanded[8].params.when, {and: [present, {'>': [{nulls: {var: 'value'}}, 0]}]});

    const custom = expandChecks({min: 1}, {when: {var: 'table'}, message: 'Too small', severity: 'warning'});
    expectDeepEqual(custom[0].params, {
      when: {and: [{var: 'table'}, {and: [present, {'<': [{var: 'value'}, 1]}]}]},
      effects: [{effect: 'warning', targets: ['target'], message: 'Too small'}],
    });
    expectDeepEqual(expandChecks({nullable: true}), []);
    expectDeepEqual(expandChecks({optional: true}), []);
    expectDeepEqual(expandChecks({optional: false}), expandChecks({nullable: false}));
    expectDeepEqual(parseRegexLiteral('/^[0-9]+$/im'), {pattern: '^[0-9]+$', flags: 'im'});
    expect(parseRegexLiteral('a > 1') === undefined, true);
  });

  test('Reject invalid checks', async () => {
    const bad = (check: any) =>
      expectThrowsAsync(() => getProcessedConfig(twoSteps([{id: 'bad', type: 'check', io: 'step1/a', check}])));
    await bad({foo: 1});
    await bad({nullable: true});
    await bad({optional: true});
    await bad({nullable: false, optional: true});
    await bad({validator: 'a > 1'});
    await bad({min: '1'});
    await bad({choices: 'a,b'});
    await bad({});
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([
      {id: 'bad', type: 'check', io: 'step1/a', check: {min: 1}, when: {var: 'other'}},
    ])));
    validateCheckOptions('ok', {max: 2});
  });

  test('Runtime operations', async () => {
    const ctx = (data: Record<string, any>) => ({all: {}, ...data});
    expect(evaluate({regex: [{var: 'v'}, '^[0-9]{4}$']}, ctx({v: '1234'})), true);
    expect(evaluate({regex: [{var: 'v'}, '^[a-z]+$', 'i']}, ctx({v: 'ABC'})), true);
    expect(evaluate({regex: [{var: 'v'}, '^[0-9]{4}$']}, ctx({v: '12'})), false);
    expect(evaluate({regex: [{var: 'v'}, '^[0-9]{4}$']}, ctx({v: 1234})), false);

    const table = makeTable(true);
    expect(evaluate({nulls: {var: 'c'}}, ctx({c: table.col('x')})), 1);
    expect(evaluate({nulls: {var: 'c'}}, ctx({c: table.col('s')})), 0);
    expect(evaluate({nulls: {var: 'c'}}, ctx({c: null})), 0);

    const dt = DG.Column.fromType(DG.TYPE.DATE_TIME, 'd', 2);
    const allKinds = ['numerical', 'numerical_no_datetime', 'categorical', 'datetime', 'categorical_or_datetime'];
    const kinds = (col: any) => allKinds.map((kind) => evaluate({columnIs: [{var: 'c'}, kind]}, ctx({c: col})));
    expectDeepEqual(kinds(table.col('x')), [true, true, false, false, false]);
    expectDeepEqual(kinds(table.col('s')), [false, false, true, false, true]);
    expectDeepEqual(kinds(dt), [true, false, false, true, true]);
    expectDeepEqual(kinds(null), [false, false, false, false, false]);
    expect(evaluate({columnIs: [{var: 'c'}, 'int']}, ctx({c: table.col('x')})), true);
    expect(evaluate({columnIs: [{var: 'c'}, 'Molecule']}, ctx({c: table.col('mol')})), true);
    expect(evaluate({columnIs: [{var: 'c'}, 'Molecule']}, ctx({c: table.col('s')})), false);
    expectDeepEqual(evaluate({columns: [{var: 't'}, 'categorical_or_datetime']}, ctx({t: table})), ['s', 'mol']);
  });

  test('Check links validate scalar values', async () => {
    const pconf = await getProcessedConfig(twoSteps([
      {id: 'range', type: 'check', io: 'step1/a', check: {min: 0, max: 10}, debounce: 0},
      {id: 'pick', type: 'check', io: 'step1/b', check: {choices: [1, 2]}, severity: 'warning', debounce: 0},
      {id: 'gated', type: 'check', io: 'step2/a', check: {max: 1}, when: {'>': [{var: 'value'}, 100]},
        message: 'Way too big', debounce: 0},
    ]));
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const n1 = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const n2 = tree.nodeTree.getNode([{idx: 1}]).getItem() as FuncCallNode;
      const snap = () => snapshots.push([n1.validationInfo$.value, n2.validationInfo$.value]);
      cold('-a').subscribe(() => {
        n1.getStateStore().setState('a', -1);
        n1.getStateStore().setState('b', 3);
        n2.getStateStore().setState('a', 50);
      });
      cold('--a').subscribe(snap);
      cold('---a').subscribe(() => {
        n1.getStateStore().setState('a', 11);
        n2.getStateStore().setState('a', 500);
      });
      cold('----a').subscribe(snap);
      cold('-----a').subscribe(() => {
        n1.getStateStore().setState('a', 5);
        n1.getStateStore().setState('b', 1);
        n2.getStateStore().setState('a', null);
      });
      cold('------a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [
      [{a: errors('Must be at least 0'), b: warnings('Must be one of: 1, 2')}, {}],
      [{a: errors('Must be at most 10'), b: warnings('Must be one of: 1, 2')}, {a: errors('Way too big')}],
      [{}, {}],
    ]);
  });

  test('Check links validate columns against the table', async () => {
    const pconf = await getProcessedConfig(annotatedStep([
      {id: 'col', type: 'check', io: 's/col', debounce: 0,
        check: {type: 'numerical', table: 's/df', allowNulls: false}},
    ]));
    const clean = makeTable();
    const dirty = makeTable(true);
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const store = node.getStateStore();
      const snap = () => snapshots.push(node.validationInfo$.value);
      cold('-a').subscribe(() => {
        store.setState('df', clean);
        store.setState('col', clean.col('x'));
      });
      cold('--a').subscribe(snap);
      cold('---a').subscribe(() => store.setState('col', clean.col('s')));
      cold('----a').subscribe(snap);
      cold('-----a').subscribe(() => store.setState('col', dirty.col('x')));
      cold('------a').subscribe(snap);
      cold('-------a').subscribe(() => store.setState('df', dirty));
      cold('--------a').subscribe(snap);
      cold('---------a').subscribe(() => store.setState('col', DG.Column.fromList('int', 'z', [1, 2])));
      cold('----------a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [
      {},
      {col: errors('Column must be numerical')},
      {col: errors('Column has missing values')},
      {col: errors('Column has missing values')},
      {col: errors('Column does not belong to the table')},
    ]);
  });

  // the initial link run happens only with defaultValidators, so parity is compared from the first change on
  function annotationScenario(pconf: any, defaultValidators: boolean, initialSnapshot = true) {
    const clean = makeTable();
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, defaultValidators});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const store = node.getStateStore();
      const snap = () => snapshots.push(node.validationInfo$.value);
      if (initialSnapshot)
        cold('-a').subscribe(snap);
      cold('--a').subscribe(() => {
        store.setState('c', 1);
        store.setState('v', 5);
        store.setState('code', '1234');
        store.setState('mode', 'fast');
        store.setState('df', clean);
        store.setState('col', clean.col('x'));
        store.setState('mol', clean.col('mol'));
      });
      cold('---a').subscribe(snap);
      cold('----a').subscribe(() => {
        store.setState('v', 20);
        store.setState('code', '12');
        store.setState('mode', 'slow');
        store.setState('mol', clean.col('s'));
      });
      cold('-----a').subscribe(snap);
      cold('------a').subscribe(() => {
        store.setState('v', 5);
        store.setState('code', 'ABCD');
        store.setState('mode', 'exact');
        store.setState('mol', clean.col('mol'));
      });
      cold('-------a').subscribe(snap);
    });
    return snapshots;
  }

  test('Annotation checks run as default validators', async () => {
    const pconf = await getProcessedConfig(annotatedStep());
    const snapshots = annotationScenario(pconf, true);
    const initial = Object.keys(snapshots[0]).sort();
    const watched = ['a', 'b', 'c', 'df', 'col', 'mol'];
    expectDeepEqual(initial.filter((key) => watched.includes(key)), ['c', 'col', 'df', 'mol']);
    for (const key of initial)
      expectDeepEqual(snapshots[0][key], errors('Missing value'));
    expectDeepEqual(snapshots.slice(1), [
      {},
      {v: errors('Must be at most 10'), code: errors('Must match /^[0-9]{4}$/i'),
        mode: errors('Must be one of: fast, exact'), mol: errors('Column must have semantic type Molecule')},
      {code: errors('Must match /^[0-9]{4}$/i')},
    ]);
  });

  test('Annotations parse named validators', async () => {
    const pconf: any = await getProcessedConfig({id: 'p', type: 'static', steps: [{id: 's', nqName: NAMED}]});
    const io = Object.fromEntries(pconf.steps[0].io.map((item: any) => [item.id, item]));
    expectDeepEqual(io.x.checks, {validators: ['LibTests:MockValidator']});
    expect(io.y.checks === undefined, true);
  });

  test('Expand the validators check', async () => {
    const [check] = expandChecks({validators: ['Pkg:f']});
    expect(check.key, 'validators');
    expect(check.needsCall, false);
    expectDeepEqual(check.params.sources, {verdicts: {validators: {input: 'value', names: ['Pkg:f']}}});
    expectDeepEqual(check.params.when, {'!': {missing: ['value']}});
    expectDeepEqual(check.params.effects, [{effect: 'verdicts', targets: ['target'], source: 'verdicts'}]);
    const [lowered] = expandChecks({validators: ['Pkg:f']}, {severity: 'notification'});
    expectDeepEqual(lowered.params.effects,
      [{effect: 'notification', targets: ['target'], message: {map: [{var: 'verdicts'}, {var: 'message'}]}}]);
    expect(expandChecks({min: 0})[0].needsCall, false);
  });

  test('Check links run the named validators directly', async () => {
    const pconf: any = await getProcessedConfig(twoSteps([
      {id: 'named', type: 'check', io: 'step1/a', check: {validators: ['LibTests:MockValidator']}},
    ]));
    const link = pconf.links[0];
    expect(link.id, 'named::validators');
    expectDeepEqual(link.from.map((io: any) => [io.name, io.flags ?? []]), [['value', []]]);
    expectDeepEqual(link.params.sources,
      {verdicts: {validators: {input: 'value', names: ['LibTests:MockValidator']}}});
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([
      {id: 'bad', type: 'check', io: 'step1/a', check: {validators: 'LibTests:MockValidator' as any}},
    ])));
  });

  test('Annotation validators are silent without a FuncCall', async () => {
    const pconf = await getProcessedConfig({id: 'p', type: 'static', steps: [{id: 's', nqName: NAMED}]});
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, defaultValidators: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      cold('-a').subscribe(() => {
        node.getStateStore().setState('x', 20);
        node.getStateStore().setState('y', 1);
      });
      cold('--a').subscribe(() => snapshots.push(node.validationInfo$.value));
      cold('---a').subscribe(() => node.getStateStore().setState('x', null));
      cold('----a').subscribe(() => snapshots.push(node.validationInfo$.value));
    });
    expectDeepEqual(snapshots, [{}, {x: errors('Missing value')}]);
  });

  test('Validator source resolves named and annotation validators alike', async () => {
    const fc = DG.Func.byName(NAMED).prepare({x: 20, y: 1});
    const controllerFor = (value: any, call?: any) => ({
      hasCall: (name: string) => name === 'call',
      getFirst: (name: string) => name === 'call' ? call : value,
      getMatchedPositions: () => [{path: [], position: 0, ioName: 'x'}],
    }) as any;
    const named = (...names: string[]) => ({verdicts: {validators: {input: 'value', names}}});
    const annotation = {verdicts: {validators: {input: 'value', call: 'call'}}};
    const tooBig = {verdicts: [{message: 'too big', isError: true, isHelper: false}]};
    expectDeepEqual(await resolveSources(controllerFor(20), named('LibTests:MockValidator')), tooBig);
    expectDeepEqual(await resolveSources(controllerFor(20, fc), annotation), tooBig);
    expectDeepEqual(await resolveSources(controllerFor(1), named('LibTests:MockValidator')), {verdicts: []});
    expectDeepEqual(await resolveSources(controllerFor(null), named('LibTests:MockValidator')), {verdicts: []});
    expectDeepEqual(await resolveSources(controllerFor(20), named('LibTests:MockValidatorBool')),
      {verdicts: [{message: 'Validation failed: LibTests:MockValidatorBool', isError: true, isHelper: false}]});
    expectDeepEqual(await resolveSources(controllerFor(1), named('LibTests:MockValidatorBool')), {verdicts: []});
    const thrown: any = await resolveSources(controllerFor(1), named('LibTests:MockValidatorThrow'));
    expect(thrown.verdicts.length, 1);
    expect(thrown.verdicts[0].isError, false);
    expect(thrown.verdicts[0].message.includes('boom'), true);
    const unknown: any = await resolveSources(controllerFor(1), named('LibTests:NoSuchValidator'));
    expect(unknown.verdicts[0].isError, false);
    expectDeepEqual(await resolveSources(controllerFor(20, undefined), annotation), {verdicts: []});
    expectDeepEqual(await resolveSources(controllerFor(20, {}), annotation), {verdicts: []});
    expectDeepEqual(resolveSources(controllerFor(20), undefined), {});
  });

  test('Rules accept sources, verdicts and array messages', async () => {
    const pconf: any = await getProcessedConfig(twoSteps([{
      id: 'named',
      type: 'rule',
      from: 'x:step1/a',
      to: 't:step1/a',
      sources: {v: {validators: {input: 'x'}}},
      debounce: 0,
      effects: [{effect: 'verdicts', targets: 't', source: 'v'}],
    }, {
      id: 'list',
      type: 'rule',
      from: 'b:step1/b',
      to: 'tb:step1/b',
      debounce: 0,
      when: {'>': [{var: 'b'}, 0]},
      effects: [{effect: 'warning', targets: 'tb', message: {if: [{'>': [{var: 'b'}, 1]}, ['first', 'second'], []]}}],
    }]));
    const namedLink = pconf.links.find((link: any) => link.id === 'named::validator');
    expectDeepEqual(namedLink.from.map((io: any) => [io.name, io.flags ?? []]),
      [['x', []], ['call', ['call', 'optional']]]);
    expectDeepEqual(namedLink.params.sources, {v: {validators: {input: 'x', call: 'call'}}});
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const store = node.getStateStore();
      cold('-a').subscribe(() => {
        store.setState('a', 20);
        store.setState('b', 2);
      });
      cold('--a').subscribe(() => snapshots.push(node.validationInfo$.value));
      cold('---a').subscribe(() => store.setState('b', 1));
      cold('----a').subscribe(() => snapshots.push(node.validationInfo$.value));
    });
    expectDeepEqual(snapshots, [{b: warnings('first', 'second')}, {}]);
    const badRule = (rule: any) => expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'bad', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      effects: [{effect: 'error', targets: 't', message: 'm'}], ...rule,
    }])));
    await badRule({sources: {x: {validators: {input: 'x'}}}});
    await badRule({sources: {v: {validators: {input: 'nope'}}}});
    await badRule({sources: {v: {validators: {input: 'x', names: 'Pkg:f'}}}});
    await badRule({sources: {v: {other: {}}}});
    await badRule({sources: {v: {validators: {input: 'x'}}},
      effects: [{effect: 'verdicts', targets: 't', source: 'w'}]});
    await badRule({from: ['x:step1/a', 'call:step1/b'], sources: {v: {validators: {input: 'x'}}}});
  });

  test('Check links give the same results as annotations', async () => {
    const annotated = annotationScenario(await getProcessedConfig(annotatedStep()), true, false);
    const configured = annotationScenario(await getProcessedConfig(annotatedStep(annotationsAsChecks)), false, false);
    const brief = (snapshot: any) => Object.fromEntries(Object.entries(snapshot).map(([io, res]: [string, any]) =>
      [io, [...res.errors, ...res.warnings].map((item: any) => item.description)]));
    expect(configured.length, annotated.length);
    for (let i = 0; i < annotated.length; i++) {
      const [c, a] = [brief(configured[i]), brief(annotated[i])];
      for (const io of new Set([...Object.keys(c), ...Object.keys(a)])) {
        if (JSON.stringify(c[io] ?? []) !== JSON.stringify(a[io] ?? [])) {
          throw new Error(`snapshot ${i} io ${io} configured ${JSON.stringify(c[io])} ` +
            `annotated ${JSON.stringify(a[io])}`);
        }
      }
    }
    expectDeepEqual(configured, annotated);
  });
});

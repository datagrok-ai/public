import * as DG from 'datagrok-api/dg';
import {category, test, before, expect} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {compileExpression, compileRuleFormulas, compileSource} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/rule-formula';
import {evaluate} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-expressions';
import {resolveSources} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-sources';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FuncCallInstancesBridge} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineLinkConfigurationInput} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';

// the examples of help/compute/workflows/rules-and-checks.mdx, with Pkg: written as LibTests:

type Links = PipelineLinkConfigurationInput<string | string[]>[];

const pk = (links: Links): PipelineConfiguration => ({
  id: 'pk',
  type: 'static',
  steps: [
    {id: 'data', nqName: 'LibTests:LoadProfile'},
    {id: 'model', nqName: 'LibTests:PkModel'},
  ],
  links,
});

const PROFILE = 'subject,time,parent\nS1,0.5,4.2\nS1,1,6.1\nS1,2,5.3\nS1,4,3.4\nS1,8,1.2';
const profile = (csv = PROFILE) => DG.DataFrame.fromCsv(csv);

const errors = (...descriptions: string[]) =>
  ({errors: descriptions.map((description) => ({description})), warnings: [], notifications: []});
const warnings = (...descriptions: string[]) =>
  ({errors: [], warnings: descriptions.map((description) => ({description})), notifications: []});

type Nodes = {data: FuncCallNode, model: FuncCallNode};
type Step = (nodes: Nodes) => void;

const nodeAt = (tree: StateTree, ...idx: number[]) =>
  tree.nodeTree.getNode(idx.map((i) => ({idx: i}))).getItem() as FuncCallNode;
const bridge = (node: FuncCallNode) => node.getStateStore() as FuncCallInstancesBridge;
const set = (node: FuncCallNode, values: Record<string, any>) => {
  for (const [io, value] of Object.entries(values))
    node.getStateStore().editState(io, value);
};
const meta = (node: FuncCallNode, io: string) => bridge(node).meta[io].value;
const restriction = (node: FuncCallNode, io: string) => bridge(node).inputRestrictions$.value[io];

category('ComputeUtils: Driver rules guide cases', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  // each step runs a second apart and is read half a second later, past the 250 ms debounce of rule messages
  function scenario(pconf: any, steps: Step[], read: (nodes: Nodes) => any) {
    const snapshots: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const nodes = {data: nodeAt(tree, 0), model: nodeAt(tree, 1)};
      steps.forEach((step, i) => {
        cold(`${i * 1000 + 1}ms a`).subscribe(() => step(nodes));
        cold(`${i * 1000 + 500}ms a`).subscribe(() => snapshots.push(read(nodes)));
      });
    });
    return snapshots;
  }

  const validations = ({model}: Nodes) => model.validationInfo$.value;

  test('Checks on the model inputs', async () => {
    const pconf = await getProcessedConfig(pk([
      {id: 'doseMin', type: 'check', io: 'model/dose', check: {min: 0}},
      {id: 'weightRange', type: 'check', io: 'model/weight', check: {min: 1, max: 300}},
      {id: 'compoundId', type: 'check', io: 'model/compound', check: {validator: '/^CMP-[0-9]{4}$/'}},
      {id: 'validatedMethods', type: 'check', io: 'model/method', check: {choices: ['LSODA']}, severity: 'warning'},
      {id: 'timeColumn', type: 'check', io: 'model/time',
        check: {table: 'model/profile', allowNulls: false}},
    ]));
    const gappy = profile('subject,time,parent\nS1,,4.2\nS1,1,6.1');
    const clean = profile();
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {dose: -1, weight: 0.5, compound: 'C-1', method: 'RK45', profile: gappy,
        time: gappy.col('time')}),
      ({model}) => set(model, {dose: 100, weight: 70, compound: 'CMP-0001', method: 'LSODA', profile: clean,
        time: clean.col('time')}),
      ({model}) => set(model, {time: DG.Column.fromList('double', 'hours', [1, 2])}),
    ], validations);
    expectDeepEqual(snapshots, [
      {
        dose: errors('Must be at least 0'),
        weight: errors('Must be at least 1'),
        compound: errors('Must match /^CMP-[0-9]{4}$/'),
        method: warnings('Must be one of: LSODA'),
        time: errors('Column has missing values'),
      },
      {},
      {time: errors('Column does not belong to the table')},
    ]);
  });

  test('Expression checks read other inputs', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'flipFlop',
      type: 'check',
      io: 'model/ka',
      check: {validator: 'value > clearance / volume'},
      vars: {clearance: 'model/clearance', volume: 'model/volume'},
      severity: 'warning',
      message: 'Absorption is slower than elimination (flip-flop kinetics)',
    }, {
      id: 'kaOral',
      type: 'check',
      io: 'model/ka',
      check: {visible: 'route == "oral"'},
      vars: {route: 'model/route'},
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {route: 'oral', clearance: 5, volume: 40, ka: 0.1}),
      ({model}) => set(model, {ka: 1}),
      ({model}) => set(model, {ka: 0.1, route: 'iv'}),
    ], (nodes) => [validations(nodes), meta(nodes.model, 'ka')]);
    expectDeepEqual(snapshots, [
      [{ka: warnings('Absorption is slower than elimination (flip-flop kinetics)')}, {hidden: false}],
      [{}, {hidden: false}],
      [{}, {hidden: true}],
    ]);
  });

  test('A check message is a formula over the value', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'doseCap',
      type: 'check',
      io: 'model/dose',
      check: {max: 2000},
      severity: 'warning',
      message: '=cat("Above the usual maximum of 2000 mg, got ", value)',
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {dose: 2500}),
      ({model}) => set(model, {dose: 100}),
    ], validations);
    expectDeepEqual(snapshots, [{dose: warnings('Above the usual maximum of 2000 mg, got 2500')}, {}]);
  });

  test('Show an input for oral dosing', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'oralOnly',
      type: 'rule',
      from: 'route:model/route',
      to: 'k:model/ka',
      when: 'eq(route, "oral")',
      effects: ['show(k)'],
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {route: 'oral'}),
      ({model}) => set(model, {route: 'iv'}),
    ], ({model}) => meta(model, 'ka'));
    expectDeepEqual(snapshots, [{hidden: false}, {hidden: true}]);
  });

  test('The formula examples evaluate as described', async () => {
    const ctx = {
      $all: {}, dose: 10, route: 'oral', ka: 1, doseUnit: 'mg/kg', weight: 70, compound: 'CMP-0001',
      subject: 'S1', profile: profile(),
    };
    const value = (formula: string) => evaluate(compileExpression(formula), ctx);
    expect(value('gt(dose, 0)'), true);
    expect(value('and(eq(route, "oral"), not(missing(ka)))'), true);
    expect(value('if(eq(doseUnit, "mg/kg"), mul(dose, weight), dose)'), 700);
    expect(value('cat("Dose ", dose, " ", doseUnit)'), 'Dose 10 mg/kg');
    expect(value('regex(compound, `^CMP-\\d{4}$`)'), true);
    expect(value('in(subject, column(profile, "subject"))'), true);
    expectDeepEqual(value('map(filter(column(profile, "time"), lt($it, 1)), cat("Early sample at ", $it))'),
      ['Early sample at 0.5']);
  });

  test('Object fields are text unless they start with =', async () => {
    const rule = compileRuleFormulas({
      id: 'weightMessages',
      type: 'rule',
      from: 'weight:model/weight',
      to: 'w:model/weight',
      effects: [
        {effect: 'error', targets: 'w', message: 'Body weight is required'},
        {effect: 'error', targets: 'w', message: '=cat("Weight ", weight, " kg is too low")'},
        {effect: 'error', targets: 'w', message: 'cat("Weight ", weight)'},
        {effect: 'error', targets: 'w', message: {var: 'weight'}},
      ],
    });
    const ctx = {$all: {}, weight: 0.5};
    expectDeepEqual((rule.effects as any[]).map((effect) => evaluate(effect.message, ctx)),
      ['Body weight is required', 'Weight 0.5 kg is too low', 'cat("Weight ", weight)', 0.5]);
  });

  test('An error while the condition holds', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'weightNeeded',
      type: 'rule',
      from: ['unit:model/doseUnit', 'weight:model/weight'],
      to: 'w:model/weight',
      when: 'and(eq(unit, "mg/kg"), missing(weight))',
      effects: ['error(w, "Body weight is required for a dose in mg/kg")'],
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {doseUnit: 'mg/kg'}),
      ({model}) => set(model, {weight: 70}),
      ({model}) => set(model, {weight: null, doseUnit: 'mg'}),
    ], validations);
    expectDeepEqual(snapshots, [{weight: errors('Body weight is required for a dose in mg/kg')}, {}, {}]);
  });

  test('One message per negative sampling time', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'sampleTimes',
      type: 'rule',
      from: 'profile:model/profile',
      to: 'p:model/profile',
      effects: ['error(p, map(filter(column(profile, "time"), lt($it, 0)), cat("Negative sampling time ", $it)))'],
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {profile: profile('subject,time,parent\nS1,-1,4.2\nS1,-0.5,5\nS1,1,6.1')}),
      ({model}) => set(model, {profile: profile()}),
    ], validations);
    expectDeepEqual(snapshots, [
      {profile: errors('Negative sampling time -1', 'Negative sampling time -0.5')},
      {},
    ]);
  });

  const subjects = (effect: string) => getProcessedConfig(pk([{
    id: 'subjects',
    type: 'rule',
    from: ['profile:model/profile', 'subject:model/subject'],
    to: 's:model/subject',
    sources: {ids: {js: {args: ['profile'], fn: (df?: DG.DataFrame) => df?.col('subject')?.categories ?? []}}},
    when: 'not(missing(profile))',
    effects: ['items(s, ids)', effect],
  }]));
  const twoSubjects = 'subject,time,parent\nS1,1,6.1\nS2,1,5.8';

  test('A subject list from the profile, with its own validation', async () => {
    const pconf = await subjects(
      'error(s, "Not a subject of the profile", when: and(not(missing(subject)), not(in(subject, ids))))');
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {profile: profile(twoSubjects), subject: 'S3'}),
      ({model}) => set(model, {subject: 'S2'}),
      ({model}) => set(model, {profile: profile()}),
    ], (nodes) => [meta(nodes.model, 'subject'), validations(nodes)]);
    expectDeepEqual(snapshots, [
      [{items: ['S1', 'S2']}, {subject: errors('Not a subject of the profile')}],
      [{items: ['S1', 'S2']}, {}],
      [{items: ['S1']}, {subject: errors('Not a subject of the profile')}],
    ]);
  });

  test('A stale subject is cleared instead', async () => {
    const pconf = await subjects('clear(s, when: and(not(missing(subject)), not(in(subject, ids))))');
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {profile: profile(twoSubjects), subject: 'S2'}),
      ({model}) => set(model, {profile: profile()}),
    ], ({model}) => [model.getStateStore().getState('subject'), meta(model, 'subject')]);
    expectDeepEqual(snapshots, [['S2', {items: ['S1', 'S2']}], [null, {items: ['S1']}]]);
  });

  test('A sensitivity analysis range per dose unit', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'doseRange',
      type: 'rule',
      from: 'unit:model/doseUnit',
      to: 'd:model/dose',
      effects: ['meta(d, rangeSA: if(eq(unit, "mg/kg"), obj(min: 1, max: 20), obj(min: 50, max: 2000)))'],
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {doseUnit: 'mg/kg'}),
      ({model}) => set(model, {doseUnit: 'mg'}),
    ], ({model}) => meta(model, 'dose'));
    expectDeepEqual(snapshots, [{rangeSA: {min: 1, max: 20}}, {rangeSA: {min: 50, max: 2000}}]);
  });

  test('A default tolerance per method, untracked and tracked', async () => {
    const rule = (effect: string) => getProcessedConfig(pk([{
      id: 'methodTolerance',
      type: 'rule',
      from: 'method:model/method',
      to: 't:model/tolerance',
      effects: [effect],
    }]));
    const steps: Step[] = [
      ({model}) => set(model, {method: 'RK45'}),
      ({model}) => set(model, {tolerance: 0.01}),
      ({model}) => set(model, {method: 'LSODA'}),
    ];
    const read = ({model}: Nodes) =>
      [model.getStateStore().getState('tolerance'), model.consistencyInfo$.value['tolerance']];
    const untracked = scenario(await rule('set(t, if(eq(method, "LSODA"), 0.000001, 0.0001))'), steps, read);
    expectDeepEqual(untracked, [[0.0001, undefined], [0.01, undefined], [0.000001, undefined]]);
    const tracked = scenario(
      await rule('set(t, if(eq(method, "LSODA"), 0.000001, 0.0001), restriction: "restricted")'), steps, read);
    expectDeepEqual(tracked, [
      [0.0001, {restriction: 'restricted', inconsistent: false, assignedValue: 0.0001}],
      [0.01, {restriction: 'restricted', inconsistent: true, assignedValue: 0.0001}],
      [0.000001, {restriction: 'restricted', inconsistent: false, assignedValue: 0.000001}],
    ]);
  });

  test('runOnInit applies the default at creation', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'methodTolerance',
      type: 'rule',
      runOnInit: true,
      from: 'method:model/method',
      to: 't:model/tolerance',
      effects: ['set(t, if(eq(method, "LSODA"), 0.000001, 0.0001))'],
    }]));
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      cold('-a').subscribe(() => values.push(nodeAt(tree, 1).getStateStore().getState('tolerance')));
    });
    // mock steps carry no annotation defaults, so the method is absent at init
    expectDeepEqual(values, [0.0001]);
  });

  test('Conditions per effect', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'doseChecks',
      type: 'rule',
      from: ['dose:model/dose', 'unit:model/doseUnit', 'weight:model/weight'],
      to: ['d:model/dose', 'w:model/weight'],
      when: 'eq(unit, "mg/kg")',
      effects: [
        'show(w)',
        'error(w, "Body weight is required for a dose in mg/kg", when: missing(weight))',
        'warning(d, cat("Total dose ", mul(dose, weight), " mg is above the usual maximum of 2000 mg"), ' +
          'when: gt(mul(dose, weight), 2000))',
      ],
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {doseUnit: 'mg/kg', dose: 30}),
      ({model}) => set(model, {weight: 80}),
      ({model}) => set(model, {doseUnit: 'mg'}),
    ], (nodes) => [validations(nodes), meta(nodes.model, 'weight')]);
    expectDeepEqual(snapshots, [
      [{weight: errors('Body weight is required for a dose in mg/kg')}, {hidden: false}],
      [{dose: warnings('Total dose 2400 mg is above the usual maximum of 2000 mg')}, {hidden: false}],
      [{}, {hidden: true}],
    ]);
  });

  test('A func source calls a platform function', async () => {
    const allometric = {
      id: 'allometric',
      type: 'rule' as const,
      from: 'weight:model/weight',
      to: 'cl:model/clearance',
      sources: {scaled: 'func("LibTests:AllometricClearance", weight: weight)'},
      when: 'not(missing(weight))',
      effects: ['set(cl, scaled, restriction: "restricted")'],
    };
    const pconf: any = await getProcessedConfig(pk([allometric]));
    expectDeepEqual(pconf.links.map((link: any) => link.id), ['allometric::data']);
    const controller = {getFirst: () => undefined} as any;
    const {scaled} = await resolveSources(controller, {scaled: compileSource(allometric.sources.scaled)},
      {$all: {}, weight: 70}) as any;
    expect(scaled, 5);
    const {scaled: none} = await resolveSources(controller, {scaled: compileSource(allometric.sources.scaled)},
      {$all: {}, weight: null}) as any;
    expect(none === null, true);
  });

  test('A validators source reports the function verdicts', async () => {
    const advice = {
      id: 'toleranceAdvice',
      type: 'rule' as const,
      from: ['tol:model/tolerance', 'method:model/method'],
      to: 't:model/tolerance',
      sources: {advice: 'validators(tol, names: ["LibTests:CheckTolerance"])'},
      when: 'eq(method, "LSODA")',
      effects: ['verdicts(t, advice)'],
    };
    const pconf: any = await getProcessedConfig(pk([advice]));
    expectDeepEqual(pconf.links.map((link: any) => link.id), ['toleranceAdvice::validator']);
    const controllerFor = (value: any) => ({getFirst: () => value}) as any;
    const source = {advice: compileSource(advice.sources.advice)};
    expectDeepEqual(await resolveSources(controllerFor(0.01), source, {$all: {}}), {advice: [
      {message: 'A tolerance above 0.001 can miss the absorption peak', isError: true, isHelper: false},
    ]});
    expectDeepEqual(await resolveSources(controllerFor(0.0001), source, {$all: {}}), {advice: []});
  });

  test('A query source binds its arguments', async () => {
    const pconf: any = await getProcessedConfig(pk([{
      id: 'compoundPresets',
      type: 'rule',
      from: 'key:model/compound',
      to: ['k:model/compound', '_(template):model/inputs(LibTests:PkModel, compound|$nonscalar|$linked)'],
      sources: {presets: 'query("Pkg:Compounds", `select compound, clearance, volume from compounds ' +
        'where compound = @key`, key: key)'},
      effects: ['assign(row(presets, "compound", key), restriction: "restricted", when: bool(row(presets, "compound", key)))'],
    }]));
    expectDeepEqual(pconf.links[0].params.sources.presets.query.args, {key: {var: 'key'}});
  });

  test('Validate an uploaded profile', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'profileQuality',
      type: 'rule',
      from: 'profile:model/profile',
      to: 'p:model/profile',
      when: 'not(missing(profile))',
      effects: [
        'error(p, map(columnsMissing(profile, [["subject", "string"], ["time", "numerical"]]), cat("Missing column ", $it)))',
        'warning(p, map(filter(column(profile, "time"), lt($it, 0)), cat("Negative sampling time ", $it)))',
        'warning(p, "At least three samples are needed for a fit", when: lt(len(profile), 3))',
      ],
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {profile: profile('id,t,parent\nS1,1,6.1\nS1,2,5.3\nS1,4,3.4')}),
      ({model}) => set(model, {profile: profile('subject,time,parent\nS1,-0.5,4.2\nS1,1,6.1\nS1,2,5.3')}),
      ({model}) => set(model, {profile: profile('subject,time,parent\nS1,1,6.1')}),
      ({model}) => set(model, {profile: profile()}),
    ], validations);
    expectDeepEqual(snapshots, [
      {profile: errors('Missing column subject (string)', 'Missing column time (numerical)')},
      {profile: warnings('Negative sampling time -0.5')},
      {profile: warnings('At least three samples are needed for a fit')},
      {},
    ]);
  });

  const COMPOUNDS = 'compound,clearance,volume\nCMP-0001,5.2,40\nCMP-0002,1.3,12';

  test('The compound table loads from the file share', async () => {
    const controller = {getFirst: () => undefined} as any;
    const {presets} = await resolveSources(controller,
      {presets: compileSource('file("System:AppData/LibTests/compounds.csv")')}) as any;
    expectDeepEqual(presets.col('compound').toList(), ['CMP-0001', 'CMP-0002']);
    expectDeepEqual(presets.columns.names(), DG.DataFrame.fromCsv(COMPOUNDS).columns.names());
  });

  test('Compound presets fill the model inputs', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'compoundPresets',
      type: 'rule',
      runOnInit: true,
      from: 'key:model/compound',
      to: ['k:model/compound', '_(template):model/inputs(LibTests:PkModel, compound|$nonscalar|$linked)'],
      sources: {presets: 'table("compound,clearance,volume\\nCMP-0001,5.2,40\\nCMP-0002,1.3,12")'},
      effects: [
        'warning(k, "Not in the compound table", when: and(not(missing(key)), not(in(key, column(presets, "compound")))))',
        'assign(row(presets, "compound", key), restriction: "restricted", when: in(key, column(presets, "compound")))',
      ],
    }]));
    const read = (nodes: Nodes) => [
      ['clearance', 'volume', 'dose'].map((io) => nodes.model.getStateStore().getState(io)),
      ['clearance', 'volume', 'dose'].map((io) => restriction(nodes.model, io)?.type),
      validations(nodes),
    ];
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {dose: 100, compound: 'CMP-0002'}),
      ({model}) => set(model, {compound: 'CMP-9999'}),
    ], read);
    expectDeepEqual(snapshots, [
      [[1.3, 12, 100], ['restricted', 'restricted', undefined], {}],
      [[1.3, 12, 100], [undefined, undefined, undefined], {compound: warnings('Not in the compound table')}],
    ]);
  });

  test('A default computed from the upstream step', async () => {
    const withCode = await getProcessedConfig(pk([{
      id: 'durationDefault',
      type: 'rule',
      runOnInit: true,
      from: 'profile:data/profile',
      to: 'd:model/duration',
      sources: {last: {js: {args: ['profile'], fn: (df?: DG.DataFrame) => df?.col('time')?.stats.max}}},
      when: 'not(missing(last))',
      effects: ['set(d, last, restriction: "restricted")'],
    }]));
    const withFormula = await getProcessedConfig(pk([{
      id: 'durationDefault',
      type: 'rule',
      runOnInit: true,
      from: 'profile:data/profile',
      to: 'd:model/duration',
      when: 'not(missing(profile))',
      effects: ['set(d, reduce(column(profile, "time"), max(current, accumulator), 0), restriction: "restricted")'],
    }]));
    for (const pconf of [withCode, withFormula]) {
      const snapshots = scenario(pconf, [
        () => {},
        ({data}) => set(data, {profile: profile()}),
      ], ({model}) => [model.getStateStore().getState('duration'), restriction(model, 'duration')?.type]);
      expectDeepEqual(snapshots, [[null, undefined], [8, 'restricted']]);
    }
  });

  test('Route-dependent absorption', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'absorption',
      type: 'rule',
      from: ['route:model/route', 'ka:model/ka'],
      to: 'k:model/ka',
      effects: [
        'show(k, when: eq(route, "oral"))',
        'error(k, "Absorption rate is required for oral dosing", when: and(eq(route, "oral"), missing(ka)))',
        'clear(k, when: and(eq(route, "iv"), not(missing(ka))))',
      ],
    }]));
    const snapshots = scenario(pconf, [
      ({model}) => set(model, {route: 'oral'}),
      ({model}) => set(model, {ka: 1.2}),
      ({model}) => set(model, {route: 'iv'}),
    ], (nodes) => [validations(nodes), meta(nodes.model, 'ka'), nodes.model.getStateStore().getState('ka')]);
    expectDeepEqual(snapshots, [
      [{ka: errors('Absorption rate is required for oral dosing')}, {hidden: false}, null],
      [{}, {hidden: false}, 1.2],
      [{}, {hidden: true}, null],
    ]);
  });

  test('One rule per model step of a dynamic workflow', async () => {
    const pconf = await getProcessedConfig({
      id: 'scenarios',
      type: 'dynamic',
      stepTypes: [
        {id: 'data', nqName: 'LibTests:LoadProfile'},
        {id: 'model', nqName: 'LibTests:PkModel'},
      ],
      initialSteps: ['data', 'model', 'model', 'data', 'model'],
      links: [{
        id: 'durationDefault',
        type: 'rule',
        runOnInit: true,
        base: 'base:expand(model)',
        from: 'profile:before(@base, data)/profile',
        to: 'd:same(@base)/duration',
        sources: {last: {js: {args: ['profile'], fn: (df?: DG.DataFrame) => df?.col('time')?.stats.max}}},
        when: 'not(missing(last))',
        effects: ['set(d, last, restriction: "restricted")'],
      }],
    });
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      cold('-a').subscribe(() => {
        set(nodeAt(tree, 0), {profile: profile()});
        set(nodeAt(tree, 3), {profile: profile('subject,time,parent\nS2,1,3\nS2,24,0.2')});
      });
      cold('--a').subscribe(() => values.push([1, 2, 4].map((i) => nodeAt(tree, i).getStateStore().getState('duration'))));
    });
    expectDeepEqual(values, [[8, 8, 24]]);
  });

  test('The object form is the same rule', async () => {
    const formula = await getProcessedConfig(pk([{
      id: 'oralOnly',
      type: 'rule',
      from: 'route:model/route',
      to: 'k:model/ka',
      when: 'eq(route, "oral")',
      effects: ['show(k)'],
    }]));
    const object = await getProcessedConfig(pk([{
      id: 'oralOnly',
      type: 'rule',
      from: 'route:model/route',
      to: 'k:model/ka',
      when: {'==': [{var: 'route'}, 'oral']},
      effects: [{effect: 'show', targets: 'k'}],
    }]));
    const params = (pconf: any) => pconf.links.map((link: any) => [link.id, link.params]);
    expectDeepEqual(params(formula), params(object));
  });
});

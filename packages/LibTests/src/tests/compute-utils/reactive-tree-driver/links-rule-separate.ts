import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {compileEffect} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/rule-formula';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode, PipelineNodeBase} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FuncCallInstancesBridge} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {DriverLogger} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineLinkConfigurationInput} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrows, expectThrowsAsync} from '../../../test-utils';

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

const scenarios = (initialSteps: string[], links: Links) => getProcessedConfig({
  id: 'scenarios',
  type: 'dynamic',
  stepTypes: [
    {id: 'data', nqName: 'LibTests:LoadProfile'},
    {id: 'model', nqName: 'LibTests:PkModel'},
  ],
  initialSteps,
  links,
} as PipelineConfiguration);

const nodeAt = (tree: StateTree, ...idx: number[]) =>
  tree.nodeTree.getNode(idx.map((i) => ({idx: i}))).getItem() as FuncCallNode;
const bridge = (node: FuncCallNode) => node.getStateStore() as FuncCallInstancesBridge;
const set = (node: FuncCallNode, values: Record<string, any>) => {
  for (const [io, value] of Object.entries(values))
    node.getStateStore().editState(io, value);
};
const get = (node: FuncCallNode, io: string) => node.getStateStore().getState(io);
const meta = (node: FuncCallNode, io: string) => bridge(node).meta[io].value;

// extra links or fields fail too
const expectExactLinks = (pconf: any, expected: any) =>
  expectDeepEqual(linksOf(pconf), expected, {forbidAdditionalProps: true});
const linksOf = (pconf: any) => Object.fromEntries(pconf.links.map((link: any) =>
  [link.id, {from: link.from.map((io: any) => io.name), to: link.to.map((io: any) => io.name)}]));

const solver = (setEffect: string): Links => [{
  id: 'solver',
  type: 'rule',
  from: ['method:model/method', 'unit:model/doseUnit'],
  to: ['t:model/tolerance', 'k:model/ka'],
  effects: [setEffect, 'hide(k, when: eq(unit, "mg"))'],
}];
const TOLERANCE = 'if(eq(method, "LSODA"), 0.000001, 0.0001)';

category('ComputeUtils: Driver rule separate effects', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  function run(pconf: any, body: (tree: StateTree, flush: () => void) => void, options: Record<string, any> = {}) {
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, ...options});
      tree.init().subscribe();
      flush();
      body(tree, flush);
    });
  }

  // the tolerance after each edit of the solver rule's inputs
  function solverEdits(pconf: any, options: Record<string, any> = {}) {
    const values: any[] = [];
    run(pconf, (tree, flush) => {
      const model = nodeAt(tree, 1);
      for (const step of [{method: 'RK45'}, {tolerance: 0.01}, {doseUnit: 'mg/kg'}, {method: 'LSODA'}]) {
        set(model, step);
        flush();
        values.push(get(model, 'tolerance'));
      }
    }, options);
    return values;
  }

  test('Without separate an input only a hide reads rewrites the set target', async () => {
    const values = solverEdits(await getProcessedConfig(pk(solver(`set(t, ${TOLERANCE})`))));
    expectDeepEqual(values, [0.0001, 0.01, 0.0001, 0.000001]);
  });

  test('A separate set keeps an edit when an input it does not read changes', async () => {
    const values = solverEdits(await getProcessedConfig(pk(solver(`set(t, ${TOLERANCE}, separate: true)`))));
    expectDeepEqual(values, [0.0001, 0.01, 0.01, 0.000001]);
  });

  test('A separate set keeps an edit with batched links', async () => {
    const values = solverEdits(await getProcessedConfig(pk(solver(`set(t, ${TOLERANCE}, separate: true)`))),
      {batchLinks: true});
    expectDeepEqual(values, [0.0001, 0.01, 0.01, 0.000001]);
  });

  test('The family link of a rule with a separate effect still applies its effects', async () => {
    const pconf = await getProcessedConfig(pk(solver(`set(t, ${TOLERANCE}, separate: true)`)));
    const hidden: any[] = [];
    run(pconf, (tree, flush) => {
      const model = nodeAt(tree, 1);
      for (const doseUnit of ['mg', 'mg/kg']) {
        set(model, {doseUnit});
        flush();
        hidden.push(meta(model, 'ka')?.hidden);
      }
    });
    expectDeepEqual(hidden, [true, false]);
  });

  test('A separate restricted set stays inconsistent when an input it does not read changes', async () => {
    const pconf =
      await getProcessedConfig(pk(solver(`set(t, ${TOLERANCE}, restriction: "restricted", separate: true)`)));
    const infos: any[] = [];
    run(pconf, (tree, flush) => {
      const model = nodeAt(tree, 1);
      for (const values of [{method: 'RK45'}, {tolerance: 0.01}, {doseUnit: 'mg/kg'}, {method: 'LSODA'}]) {
        set(model, values);
        flush();
        infos.push([get(model, 'tolerance'), model.consistencyInfo$.value['tolerance']?.inconsistent]);
      }
    });
    expectDeepEqual(infos, [[0.0001, false], [0.01, true], [0.01, true], [0.000001, false]]);
  });

  test('Additional triggers reset a target on inputs the effect does not read', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'subjectReset',
      type: 'rule',
      from: ['p:data/profile', 'c:model/compound', 'r:model/route', 'u:model/doseUnit', 'w:model/weight'],
      to: 's:model/subject',
      effects: ['clear(s, when: not(missing(p)), additionalTriggers: [c, r, u])'],
    }]));
    const subjects: any[] = [];
    run(pconf, (tree, flush) => {
      const [data, model] = [nodeAt(tree, 0), nodeAt(tree, 1)];
      const step = (node: FuncCallNode, values: Record<string, any>) => {
        set(node, values);
        flush();
        subjects.push(get(model, 'subject'));
      };
      step(data, {profile: DG.DataFrame.fromCsv('subject,time\nS1,1')});
      step(model, {subject: 'S1'});
      step(model, {weight: 70});
      step(model, {compound: 'C1'});
      step(model, {subject: 'S2'});
      step(model, {route: 'oral'});
      step(model, {subject: 'S3'});
      step(model, {doseUnit: 'mg/kg'});
      step(data, {profile: null});
      step(model, {subject: 'S4'});
      step(model, {compound: 'C2'});
    });
    expectDeepEqual(subjects, [null, 'S1', 'S1', null, 'S2', null, 'S3', null, null, 'S4', 'S4']);
  });

  test('A separate effect inherits the rule condition and its inputs trigger it', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'oralKa',
      type: 'rule',
      from: ['route:model/route', 'method:model/method'],
      to: 'k:model/ka',
      when: 'eq(route, "oral")',
      effects: ['set(k, 1.5, separate: true)'],
    }]));
    expectExactLinks(pconf, {'oralKa::data:1': {from: ['route'], to: ['k']}});
    const values: any[] = [];
    run(pconf, (tree, flush) => {
      const model = nodeAt(tree, 1);
      for (const step of [{method: 'RK45'}, {route: 'iv'}, {route: 'oral'}, {ka: 2}, {method: 'LSODA'}]) {
        set(model, step);
        flush();
        values.push(get(model, 'ka'));
      }
    });
    expectDeepEqual(values, [undefined, undefined, 1.5, 2, 2]);
  });

  test('A separate set without inputs writes its default once with runOnInit', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'defaults',
      type: 'rule',
      runOnInit: true,
      from: 'method:model/method',
      to: 't:model/tolerance',
      effects: ['set(t, 0.0001, separate: true)'],
    }]));
    expectExactLinks(pconf, {'defaults::data:1': {from: [], to: ['t']}});
    const values: any[] = [];
    run(pconf, (tree, flush) => {
      const model = nodeAt(tree, 1);
      values.push(get(model, 'tolerance'));
      for (const step of [{tolerance: 0.01}, {method: 'RK45'}]) {
        set(model, step);
        flush();
        values.push(get(model, 'tolerance'));
      }
    });
    expectDeepEqual(values, [0.0001, 0.01, 0.01]);
  });

  test('A source is resolved only by the links of the effects using it', async () => {
    let calls = 0;
    const rule = (separate: boolean) => getProcessedConfig(pk([{
      id: 'scaled',
      type: 'rule',
      from: ['w:model/weight', 'method:model/method'],
      to: ['cl:model/clearance', 't:model/tolerance'],
      sources: {scaled: {js: {args: ['w'], fn: (w: number) => {
        calls++;
        return w == null ? null : w * 2;
      }}}},
      effects: [`set(cl, scaled, separate: ${separate})`, `set(t, ${TOLERANCE}, separate: ${separate})`],
    }]));
    const counts = async (separate: boolean) => {
      const pconf = await rule(separate);
      const seen: any[] = [];
      run(pconf, (tree, flush) => {
        const model = nodeAt(tree, 1);
        calls = 0;
        for (const step of [{method: 'RK45'}, {weight: 70}]) {
          set(model, step);
          flush();
          seen.push([calls, get(model, 'clearance')]);
        }
      });
      return seen;
    };
    expectDeepEqual(await counts(false), [[1, null], [2, 140]], {prefix: 'family'});
    expectDeepEqual(await counts(true), [[0, undefined], [1, 140]], {prefix: 'separate'});
    const pconf = await rule(true);
    expectDeepEqual(pconf.links!.map((link) => Object.keys(link.params?.sources ?? {})), [['scaled'], []]);
  });

  test('Meta of a separate effect and of the family link merge on one target', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'kaInput',
      type: 'rule',
      from: ['route:model/route', 'unit:model/doseUnit'],
      to: 'k:model/ka',
      effects: ['items(k, [0.5, 1])', 'hide(k, when: eq(route, "iv"), separate: true)'],
    }]));
    const metas: any[] = [];
    run(pconf, (tree, flush) => {
      const model = nodeAt(tree, 1);
      for (const step of [{doseUnit: 'mg'}, {route: 'iv'}, {route: 'oral'}]) {
        set(model, step);
        flush();
        metas.push([meta(model, 'ka')?.items, meta(model, 'ka')?.hidden]);
      }
    });
    expectDeepEqual(metas, [[[0.5, 1], undefined], [[0.5, 1], true], [[0.5, 1], false]]);
  });

  test('Separate data effects of one rule may write the same input', async () => {
    const links = (rules: Links) => getProcessedConfig(pk(rules));
    const duplicates = (logger: DriverLogger) =>
      logger.errors.filter((error) => error.message.includes('Duplicate data link target')).length;
    const pconf = await links([{
      id: 'byMethod',
      type: 'rule',
      from: ['method:model/method', 'unit:model/doseUnit'],
      to: 't:model/tolerance',
      effects: [
        'set(t, 0.0001, restriction: "restricted", when: eq(method, "RK45"))',
        'set(t, 0.000001, restriction: "restricted", when: eq(method, "LSODA"), separate: true)',
      ],
    }]);
    const values: any[] = [];
    const logger = new DriverLogger();
    run(pconf, (tree, flush) => {
      const model = nodeAt(tree, 1);
      for (const method of ['RK45', 'LSODA']) {
        set(model, {method});
        flush();
        values.push(get(model, 'tolerance'));
      }
    }, {logger});
    expectDeepEqual(values, [0.0001, 0.000001]);
    expectDeepEqual(duplicates(logger), 0, {prefix: 'one rule'});
    const two = await links(['a', 'b'].map((id) => ({
      id, type: 'rule' as const, from: 'method:model/method', to: 't:model/tolerance', effects: ['set(t, 1)'],
    })));
    const twoLogger = new DriverLogger();
    run(two, () => {}, {logger: twoLogger});
    expectDeepEqual(duplicates(twoLogger) > 0, true, {prefix: 'two rules'});
  });

  test('A separate effect activates on its own inputs', async () => {
    const pconf = await scenarios(['model'], [{
      id: 'partial',
      type: 'rule',
      from: ['p:first(data)/profile', 'm:first(model)/method'],
      to: ['t:first(model)/tolerance', 'k:first(model)/ka'],
      effects: ['hide(k, when: missing(p))', 'set(t, if(eq(m, "LSODA"), 1, 2), separate: true)'],
    }]);
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      tree.init().subscribe();
      flush();
      set(nodeAt(tree, 0), {method: 'LSODA'});
      flush();
      expectDeepEqual([get(nodeAt(tree, 0), 'tolerance'), meta(nodeAt(tree, 0), 'ka')?.hidden], [1, undefined],
        {prefix: 'without a data step'});
      tree.addSubTree(tree.nodeTree.root.getItem().uuid, 'data', 0).subscribe();
      flush();
      expectDeepEqual(meta(nodeAt(tree, 1), 'ka')?.hidden, true, {prefix: 'with a data step'});
    });
  });

  test('A separate node message reading only the step count runs on rewire', async () => {
    const pconf = await scenarios([], [{
      id: 'stepsRequired',
      type: 'rule',
      from: 'm:first(model)/method',
      to: 'scenarios(node)',
      effects: ['error(scenarios, "Add a step", when: eq($steps, 0), separate: true)'],
    }]);
    expectDeepEqual(pconf.links!.map((link) => [link.id, link.from.length, link.params?.readsSteps]),
      [['stepsRequired::pipelineValidator:1', 0, true]]);
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      tree.init().subscribe();
      flush();
      const root = tree.nodeTree.root.getItem() as PipelineNodeBase;
      const messages = () => Object.values(root.pipelineValidations$.value).flatMap((r) => r?.errors ?? []);
      expectDeepEqual(messages(), [{description: 'Add a step'}], {prefix: 'empty'});
      tree.addSubTree(root.uuid, 'data', 0).subscribe();
      flush();
      expectDeepEqual(messages(), [], {prefix: 'with a step'});
    });
  });
});

category('ComputeUtils: Driver rule separate effects expansion', async () => {
  test('A separate effect gets its own link with the inputs it reads', async () => {
    const pconf = await getProcessedConfig(pk(solver(`set(t, ${TOLERANCE}, separate: true)`)));
    expectExactLinks(pconf, {
      'solver::meta': {from: ['method', 'unit'], to: ['k']},
      'solver::data:1': {from: ['method'], to: ['t']},
    });
  });

  test('A separate link reads the rule condition and its own condition', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'r',
      type: 'rule',
      from: ['m:model/method', 'w:model/weight', 'r:model/route', 'u:model/doseUnit'],
      to: 't:model/tolerance',
      when: 'not(missing(w))',
      effects: ['set(t, 1, when: eq(r, "oral"), separate: true)'],
    }]));
    expectExactLinks(pconf, {'r::data:1': {from: ['w', 'r'], to: ['t']}});
  });

  test('Additional triggers imply separate and keep the rule input order', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'r',
      type: 'rule',
      from: ['c:model/compound', 'r:model/route', 'u:model/doseUnit', 'p:data/profile'],
      to: 's:model/subject',
      effects: [{effect: 'clear', targets: 's', when: {'!': [{missing: ['p']}]}, additionalTriggers: ['u', 'c']}],
    }]));
    expectExactLinks(pconf, {'r::data:1': {from: ['c', 'u', 'p'], to: ['s']}});
  });

  test('A rule whose effects are all separate has no family links', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'r',
      type: 'rule',
      from: ['m:model/method', 'u:model/doseUnit'],
      to: ['t:model/tolerance', 'k:model/ka'],
      runOnInit: true,
      debounce: 0,
      effects: [
        'hide(k, when: eq(u, "mg"), separate: true)',
        'set(t, if(eq(m, "LSODA"), 1, 2), separate: true)',
        'error(t, "Too small", when: lt(m, 0), separate: true)',
      ],
    }]));
    const links = pconf.links! as any[];
    expectDeepEqual(links.map((link) => link.id), ['r::meta:1', 'r::data:2', 'r::validator:3']);
    expectDeepEqual(links.map((link) => link.type), ['meta', 'data', 'validator']);
    expectDeepEqual([links[1].runOnInit, links[2].debounce], [true, 0]);
  });

  test('A template input query is narrowed per alias', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'r',
      type: 'rule',
      from: '_(template):model/inputs(LibTests:PkModel, profile|time)',
      to: 'tol:model/tolerance',
      effects: ['set(tol, if(eq(method, "LSODA"), 1, 2), separate: true)', 'hide(tol, when: eq(route, "iv"))'],
    }]));
    const links = linksOf(pconf);
    expectDeepEqual(links['r::data:1'].from, ['method']);
    expectDeepEqual(links['r::meta'].from.length, 12);
  });

  test('Sources bring their inputs and the call into a separate link', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'r',
      type: 'rule',
      from: ['c:model/compound', 'r:model/route', 'u:model/doseUnit', 'm:model/method', 'tol:model/tolerance'],
      to: ['k:model/ka', 't:model/tolerance'],
      sources: {opts: 'choices(c)', advice: 'validators(tol, names: ["LibTests:CheckTolerance"])'},
      effects: [
        'items(k, opts.items, additionalTriggers: [r, u])',
        'verdicts(t, advice, separate: true)',
        'hide(t, when: eq(m, "x"))',
      ],
    }]));
    const links = pconf.links! as any[];
    const byId = Object.fromEntries(links.map((link) => [link.id, link]));
    expectDeepEqual(linksOf(pconf)['r::meta:1'].from, ['c', 'r', 'u', '$call']);
    expectDeepEqual(Object.keys(byId['r::meta:1'].params.sources), ['opts']);
    expectDeepEqual(linksOf(pconf)['r::validator:2'].from, ['tol']);
    expectDeepEqual(Object.keys(byId['r::validator:2'].params.sources), ['advice']);
    expectDeepEqual(linksOf(pconf)['r::meta'].from, ['c', 'r', 'u', 'm', 'tol', '$call']);
    expectDeepEqual(Object.keys(byId['r::meta'].params.sources), ['opts', 'advice']);
  });

  test('Data frame mutations and step count reads follow the separate link', async () => {
    const pconf = await getProcessedConfig(pk([{
      id: 'r',
      type: 'rule',
      from: ['p:data/profile', 'm:model/method'],
      to: ['t:model/tolerance', 'k:model/ka'],
      dataFrameMutations: ['p', 'm'],
      effects: ['set(t, if(eq(m, "LSODA"), 1, 2), separate: true)', 'hide(k, when: eq($steps, 0), separate: true)'],
    }]));
    const links = pconf.links! as any[];
    expectDeepEqual(links.map((link) => [link.id, link.dataFrameMutations, link.params.readsSteps]), [
      ['r::data:1', ['m'], undefined],
      ['r::meta:2', [], true],
    ]);
  });

  test('Reject invalid separate effects', async () => {
    const rule = (effect: any, extra: Record<string, any> = {}) => getProcessedConfig(pk([{
      id: 'r',
      type: 'rule',
      from: ['m:model/method', 'u:model/doseUnit'],
      to: 't:model/tolerance',
      sources: {src: {js: {args: ['m'], fn: (m: any) => m}}},
      effects: [effect],
      ...extra,
    }]));
    await expectThrowsAsync(() => rule({effect: 'clear', targets: 't', additionalTriggers: ['zz']}),
      /additionalTriggers of clear names unknown input alias zz/);
    await expectThrowsAsync(() => rule({effect: 'clear', targets: 't', additionalTriggers: 'src'}),
      /names unknown input alias src/);
    await expectThrowsAsync(() => rule({effect: 'clear', targets: 't', separate: false, additionalTriggers: 'm'}),
      /additionalTriggers of effect clear needs a separate link/);
    await expectThrowsAsync(() => rule({effect: 'clear', targets: 't', separate: 'yes'}),
      /separate of effect clear must be true or false/);
    await expectThrowsAsync(() => rule('hide(t, separate: true)'), /separate effect hide reads no input/);
    await expectThrowsAsync(() => rule('set(t, 1, separate: true)'),
      /separate effect set reads no input, add additionalTriggers or runOnInit/);
    await rule('set(t, 1, separate: true)', {runOnInit: true});
  });

  test('Formulas take separate and additional triggers on every effect', async () => {
    expectDeepEqual(compileEffect('set(t, 1, separate: true, additionalTriggers: [a, b])'),
      {effect: 'set', targets: 't', value: 1, separate: true, additionalTriggers: ['a', 'b']},
      {forbidAdditionalProps: true});
    expectDeepEqual(compileEffect('clear(t, additionalTriggers: a)'),
      {effect: 'clear', targets: 't', additionalTriggers: 'a'}, {forbidAdditionalProps: true});
    expectDeepEqual(compileEffect('meta(t, units: "K", separate: true)'),
      {effect: 'meta', targets: 't', meta: {units: 'K'}, separate: true}, {forbidAdditionalProps: true});
    expectThrows(() => compileEffect('set(t, 1, separate: 1)'), /separate must be true or false/);
    expectThrows(() => compileEffect('clear(t, additionalTriggers: "a")'),
      /additionalTriggers are an input alias or a list of them/);
  });
});

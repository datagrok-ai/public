import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallInstancesBridge} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineLinkConfigurationInput} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import * as DG from 'datagrok-api/dg';
import {createTestScheduler} from '../../../test-utils';

// the advanced examples of help/compute/workflows/link-query-language-advanced.mdx, on LibTests mocks:
// load = TestDF1 (out res), solver and analysis = TestAnnotatedInputs (in a, b, c, v, mode, df; out res),
// reset and metrics = TestAdd2

const DEFAULT_STEPS = ['load', 'solver', 'analysis', 'analysis', 'summary'];

const screening = (
  links: PipelineLinkConfigurationInput<string | string[]>[], initialSteps = DEFAULT_STEPS,
): PipelineConfiguration => ({
  id: 'screening',
  type: 'dynamic',
  stepTypes: [
    {id: 'load', nqName: 'LibTests:TestDF1'},
    {id: 'solver', nqName: 'LibTests:TestAnnotatedInputs', tags: ['report']},
    {id: 'analysis', nqName: 'LibTests:TestAnnotatedInputs'},
    {id: 'reset', nqName: 'LibTests:TestAdd2'},
    {id: 'summary', type: 'static', steps: [{id: 'metrics', nqName: 'LibTests:TestAdd2', tags: ['report']}]},
  ],
  initialSteps,
  links,
});

const table = (...names: string[]) =>
  DG.DataFrame.fromColumns(names.map((name) => DG.Column.fromList('double', name, [1, 2, 3])));

const bridgeAt = (tree: StateTree, ...idx: number[]) =>
  tree.nodeTree.getNode(idx.map((i) => ({idx: i}))).getItem().getStateStore() as FuncCallInstancesBridge;

category('ComputeUtils: Driver docs cases', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Chain each analysis to the next one', async () => {
    const pconf = await getProcessedConfig(screening([{
      id: 'chain',
      base: 'base:expand(analysis)',
      from: 'prev:same(@base)/res',
      to: 'next:after+(@base, analysis)/c',
    }]));
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const [first, second] = [bridgeAt(tree, 2), bridgeAt(tree, 3)];
      cold('-a').subscribe(() => {
        first.setState('res', 4);
        second.setState('res', 5);
      });
      cold('--a').subscribe(() => values.push([first.getState('c'), second.getState('c')]));
    });
    expectDeepEqual(values, [[null, 4]]);
  });

  test('Feed each analysis from the load of its section', async () => {
    const pconf = await getProcessedConfig(screening([{
      id: 'feed',
      base: 'base:expand(analysis)',
      from: 'table:before(@base, load, reset)/res',
      to: 'input:same(@base)/df',
    }], ['load', 'analysis', 'reset', 'analysis']));
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      cold('-a').subscribe(() => bridgeAt(tree, 0).setState('res', table('x')));
      cold('--a').subscribe(() => values.push([1, 3].map((i) => bridgeAt(tree, i).getState('df')?.columns.names())));
    });
    expectDeepEqual(values, [[['x'], undefined]]);
  });

  test('Collect every score before the summary', async () => {
    const pconf = await getProcessedConfig(screening([{
      id: 'collect',
      base: 'base:expand(summary)',
      from: 'scores:before*(@base, analysis)/res',
      to: 'target:same(@base)/metrics/a',
      handler: ({controller}) => controller.setAll('target', controller.getAll('scores')),
    }]));
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      cold('-a').subscribe(() => {
        bridgeAt(tree, 2).setState('res', 1);
        bridgeAt(tree, 3).setState('res', 2);
      });
      cold('--a').subscribe(() => values.push(bridgeAt(tree, 4, 0).getState('a')));
    });
    expectDeepEqual(values, [[1, 2]]);
  });

  test('Tags match across nesting', async () => {
    const pconf = await getProcessedConfig(screening([{
      id: 'reportsReady',
      type: 'pipelineValidator',
      debounce: 0,
      from: 'results:#all(report)/res',
      to: 'self',
      handler: ({controller}) => {
        const ready = (controller.getAll('results') ?? []).filter((result: any) => result != null).length;
        controller.setValidation('self', ready === 2 ? undefined : {errors: [`${ready} of 2 reports ready`]});
      },
    }]));
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const snap = () => values.push((tree.toState({skipFuncCalls: true}) as any).structureCheckResults?.errors);
      cold('-a').subscribe(() => bridgeAt(tree, 1).setState('res', 1));
      cold('--a').subscribe(snap);
      cold('---a').subscribe(() => bridgeAt(tree, 4, 0).setState('res', 2));
      cold('----a').subscribe(snap);
    });
    expectDeepEqual(values, [['1 of 2 reports ready'], undefined]);
  });

  test('Wildcard io selectors expand at config time', async () => {
    const pconf: any = await getProcessedConfig(screening([{
      id: 'loadToSolver',
      from: 'in_(template):load/outputs(LibTests:TestDF1)',
      to: 'out_(template):solver/inputs(LibTests:TestAnnotatedInputs, a|b|c|v|code|mode|col|mol)',
    }]));
    const names = (ios: any[]) => ios.map((io) => io.name);
    expectDeepEqual(names(pconf.links[0].from), ['in_res']);
    expectDeepEqual(names(pconf.links[0].to), ['out_df']);
  });

  test('Computed default from an upstream step', async () => {
    const pconf = await getProcessedConfig(screening([{
      id: 'rowsDefault',
      type: 'rule',
      runOnInit: true,
      base: 'base:expand(analysis)',
      from: 'table:before(@base, load)/res',
      to: 't:same(@base)/v',
      sources: {rows: {js: {args: ['table'], fn: (df?: DG.DataFrame) => df?.rowCount}}},
      when: {'!': {missing: ['table']}},
      effects: [{effect: 'set', targets: 't', value: {var: 'rows'}, restriction: 'restricted'}],
    }]));
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const analyses = [bridgeAt(tree, 2), bridgeAt(tree, 3)];
      const snap = () => values.push(analyses.map((b) => [b.getState('v'), b.inputRestrictions$.value['v']]));
      cold('-a').subscribe(snap);
      cold('--a').subscribe(() => bridgeAt(tree, 0).setState('res', table('x')));
      cold('---a').subscribe(snap);
    });
    expectDeepEqual(values, [
      [[null, undefined], [null, undefined]],
      [[3, {assignedValue: 3, type: 'restricted'}], [3, {assignedValue: 3, type: 'restricted'}]],
    ]);
  });

  test('Dynamic option list from an upstream table', async () => {
    const pconf = await getProcessedConfig(screening([{
      id: 'columnItems',
      type: 'rule',
      base: 'base:expand(analysis)',
      from: ['table:before(@base, load)/res', 'column:same(@base)/mode'],
      to: 'c:same(@base)/mode',
      sources: {columns: {js: {args: ['table'], fn: (df?: DG.DataFrame) => df ? df.columns.names() : []}}},
      effects: [
        {effect: 'items', targets: 'c', items: {var: 'columns'}},
        {
          effect: 'clear', targets: 'c',
          when: {and: [{'!': {missing: ['column']}}, {'!': {in: [{var: 'column'}, {var: 'columns'}]}}]},
        },
      ],
    }]));
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const load = bridgeAt(tree, 0);
      const analyses = [bridgeAt(tree, 2), bridgeAt(tree, 3)];
      const snap = () => values.push(analyses.map((b) => [b.getState('mode'), b.meta.mode.value?.items]));
      cold('-a').subscribe(() => load.setState('res', table('x', 's')));
      cold('--a').subscribe(() => analyses[0].setState('mode', 'x'));
      cold('---a').subscribe(snap);
      cold('----a').subscribe(() => load.setState('res', table('y')));
      cold('-----a').subscribe(snap);
    });
    expectDeepEqual(values, [
      [['x', ['x', 's']], [null, ['x', 's']]],
      [[null, ['y']], [null, ['y']]],
    ]);
  });

  test('Lookup table fills sibling inputs', async () => {
    const presets = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'mode', ['fast', 'exact']),
      DG.Column.fromList('double', 'a', [1, 10]),
      DG.Column.fromList('double', 'b', [2, 20]),
      DG.Column.fromList('double', 'c', [3, 30]),
    ]);
    const pconf = await getProcessedConfig(screening([{
      id: 'preset',
      type: 'rule',
      runOnInit: true,
      from: 'key:solver/mode',
      to: ['k:solver/mode', '_(template):solver/inputs(LibTests:TestAnnotatedInputs, mode|$nonscalar|$linked)'],
      sources: {presets: {js: {args: [], fn: () => presets}}},
      effects: [
        {effect: 'items', targets: 'k', items: {column: [{var: 'presets'}, 'mode']}},
        {
          effect: 'assign', values: {row: [{var: 'presets'}, 'mode', {var: 'key'}]}, restriction: 'restricted',
          when: {in: [{var: 'key'}, {column: [{var: 'presets'}, 'mode']}]},
        },
      ],
    }]));
    const values: any[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const solver = bridgeAt(tree, 1);
      const snap = () => values.push([
        ['a', 'b', 'c'].map((io) => solver.getState(io)),
        ['a', 'b', 'c'].map((io) => solver.inputRestrictions$.value[io]?.type),
        solver.meta.mode.value?.items,
      ]);
      cold('-a').subscribe(snap);
      cold('--a').subscribe(() => solver.setState('mode', 'exact'));
      cold('---a').subscribe(snap);
      cold('----a').subscribe(() => solver.setState('mode', 'custom'));
      cold('-----a').subscribe(snap);
    });
    expectDeepEqual(values, [
      [[null, null, null], [undefined, undefined, undefined], undefined],
      [[10, 20, 30], ['restricted', 'restricted', 'restricted'], ['fast', 'exact']],
      [[10, 20, 30], [undefined, undefined, undefined], ['fast', 'exact']],
    ]);
  });
});

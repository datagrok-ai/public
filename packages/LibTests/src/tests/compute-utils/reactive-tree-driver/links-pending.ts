import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {
  HeadlessResult, StepRunState, WorkflowRunState, runWorkflowConfig,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import {TestScheduler} from 'rxjs/testing';
import {of} from 'rxjs';
import {delay, distinctUntilChanged, tap} from 'rxjs/operators';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';

type CheckOptions = {debounce?: number, sequential?: boolean, asyncMs?: number, throws?: boolean};

// step2/a gets an error above 2
const makeConfig = (opts: CheckOptions): PipelineConfiguration => ({
  id: 'pipeline',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {id: 'step2', nqName: 'LibTests:TestMul2'},
  ],
  links: [
    {id: 'data', from: 'in:step1/res', to: 'out:step2/a'},
    {
      id: 'check', type: 'validator', from: 'in:step2/a', to: 'out:step2/a',
      ...(opts.debounce != null ? {debounce: opts.debounce} : {}),
      ...(opts.sequential ? {sequential: true} : {}),
      handler({controller}) {
        const set = () => {
          const a = controller.getFirst<number>('in');
          if (opts.throws && a != null)
            throw new Error('check failed');
          controller.setValidation('out', (a ?? 0) > 2 ? {errors: [{description: `a is ${a}`}]} : undefined);
        };
        if (opts.asyncMs == null)
          return set();
        return of(undefined).pipe(delay(opts.asyncMs), tap(set));
      },
    },
  ],
});

category('ComputeUtils: Driver pending links', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  async function runHeadless(conf: PipelineConfiguration) {
    const pconf = await getProcessedConfig(conf);
    const calls: string[] = [];
    let result: HeadlessResult | undefined;
    let doneAt = -1;
    testScheduler.run(() => {
      runWorkflowConfig({
        config: pconf, mockMode: true,
        instanceConfig: {steps: [{id: 'step1', initialValues: {a: 1, b: 2}}, {id: 'step2', initialValues: {b: 10}}]},
        mockRun: (id, state) => {
          calls.push(`${id} a=${state.a}`);
          return {res: id === 'step1' ? state.a + state.b : state.a * state.b};
        },
      }).subscribe((res) => {
        result = res;
        doneAt = testScheduler.now();
      });
    });
    const step2 = (result!.state as WorkflowRunState).steps[1] as StepRunState;
    return {calls, step2, doneAt};
  }

  async function expectBlocked(opts: CheckOptions, doneAt: number) {
    const res = await runHeadless(makeConfig(opts));
    expectDeepEqual(res.calls, ['step1 a=1']);
    expectDeepEqual(res.step2.validations.a?.errors, [{description: 'a is 3'}]);
    expectDeepEqual(res.doneAt, doneAt, {prefix: 'done at'});
  }

  test('A debounced validator blocks its step without waiting for the debounce', async () =>
    expectBlocked({debounce: 100}, 0));

  test('An async validator blocks its step', async () =>
    expectBlocked({asyncMs: 50}, 100));

  test('An async sequential validator blocks its step', async () =>
    expectBlocked({asyncMs: 50, sequential: true}, 100));

  test('An async debounced validator blocks its step without waiting for the debounce', async () =>
    expectBlocked({asyncMs: 50, debounce: 100}, 100));

  test('An async data link sets the input before its step runs', async () => {
    const res = await runHeadless({
      id: 'pipeline',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'data', from: 'in:step1/res', to: 'out:step2/a',
        handler({controller}) {
          return of(undefined).pipe(delay(50), tap(() => controller.setAll('out', controller.getFirst('in'))));
        },
      }],
    });
    expectDeepEqual(res.calls, ['step1 a=1', 'step2 a=3']);
  });

  async function editTimeline(opts: CheckOptions, edits: [number, number][]) {
    const pconf = await getProcessedConfig(makeConfig(opts));
    const locked: [number, boolean][] = [];
    let tree: StateTree;
    testScheduler.run(({cold}) => {
      tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      tree.init().subscribe();
      cold('5ms a').subscribe(() => {
        tree.treeMutationsLocked$.pipe(distinctUntilChanged())
          .subscribe((x) => locked.push([testScheduler.now(), x]));
      });
      for (const [time, a] of edits) {
        cold(`${time}ms a`).subscribe(() =>
          (tree.nodeTree.getNode([{idx: 1}]).getItem() as FuncCallNode).getStateStore().setState('a', a));
      }
    });
    const step2 = tree!.nodeTree.getNode([{idx: 1}]).getItem() as FuncCallNode;
    return {locked, validations: step2.validationInfo$.value};
  }

  test('An edit locks the tree until its debounced validator runs', async () => {
    const res = await editTimeline({debounce: 100}, [[10, 5]]);
    expectDeepEqual(res.locked, [[5, false], [10, true], [110, false]]);
    expectDeepEqual(res.validations, {a: {errors: [{description: 'a is 5'}], warnings: [], notifications: []}});
  });

  test('An edit locks the tree until its async validator finishes', async () => {
    const res = await editTimeline({asyncMs: 50}, [[10, 5]]);
    expectDeepEqual(res.locked, [[5, false], [10, true], [60, false]]);
  });

  test('A superseded async validator run keeps the tree locked until the latest one finishes', async () => {
    const res = await editTimeline({asyncMs: 50}, [[10, 5], [30, 7]]);
    expectDeepEqual(res.locked, [[5, false], [10, true], [80, false]]);
    expectDeepEqual(res.validations, {a: {errors: [{description: 'a is 7'}], warnings: [], notifications: []}});
  });

  test('A throwing debounced validator unlocks the tree', async () => {
    const res = await editTimeline({debounce: 100, throws: true}, [[10, 5]]);
    expectDeepEqual(res.locked, [[5, false], [10, true], [110, false]]);
  });

  test('A mutation waits for a pending debounced validator', async () => {
    const pconf = await getProcessedConfig(makeConfig({debounce: 100}));
    const events: string[] = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      tree.init().subscribe();
      const step2 = () => tree.nodeTree.getNode([{idx: 1}]).getItem() as FuncCallNode;
      cold('10ms a').subscribe(() => step2().getStateStore().setState('a', 5));
      cold('20ms a').subscribe(() => {
        step2().validationInfo$.subscribe((v) => {
          if (Object.keys(v).length)
            events.push(`validated ${testScheduler.now()}`);
        });
        tree.runMutateTree().subscribe(() => events.push(`mutated ${testScheduler.now()}`));
      });
    });
    expectDeepEqual(events, ['validated 110', 'mutated 110']);
  });
});

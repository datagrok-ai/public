import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';
import {of} from 'rxjs';
import {delay, map, switchMap, tap} from 'rxjs/operators';
import {FuncCallNode, StaticPipelineNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';


category('ComputeUtils: Driver hooks running', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  const config1: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'static',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestAdd2',
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestMul2',
      },
    ],
    states: [{
      id: 'meta1',
    }],
    onInit: {
      id: 'link1',
      from: 'in1:step1/b',
      to: 'out1:meta1',
      handler({controller}) {
        return of(undefined).pipe(
          delay(250),
          tap(() => controller.setAll('out1', 10)),
        );
      },
    },
  };

  const config2: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'static',
    steps: [
      {
        ...config1,
      },
    ],
  };

  const config3: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'parallel',
    stepTypes: [{
      ...config1,
    }],
  };

  const config4: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'static',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestAdd2',
      },
    ],
    onReturn: {
      id: 'link1',
      from: 'in1:step1/res',
      type: 'return',
      to: [],
      handler({controller}) {
        const res = controller.getFirst('in1');
        controller.returnResult(res);
      },
    },
  };

  test('Run onInit', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const rnode = tree.nodeTree.root;
      expectObservable(rnode.getItem().getStateStore().getStateChanges('meta1')).toBe('a 249ms b', {a: undefined, b: 10});
    });
  });

  test('Run nested onInit', async () => {
    const pconf = await getProcessedConfig(config2);

    testScheduler.run((helpers) => {
      const {expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const item = tree.nodeTree.getItem([{idx: 0}]) as StaticPipelineNode;
      expectObservable(item.getStateStore().getStateChanges('meta1')).toBe('a 249ms b', {a: undefined, b: 10});
    });
  });

  test('Run nested onInit for dynamic items', async () => {
    const pconf = await getProcessedConfig(config3);

    testScheduler.run((helpers) => {
      const {expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const rnode = tree.nodeTree.root;
      tree.addSubTree(rnode.getItem().uuid, 'pipeline1', 0).subscribe();
      const item$ = tree.makeStateRequests$.pipe(
        map(() => tree.nodeTree.getItem([{idx: 0}])),
        switchMap((item) => item.getStateStore().getStateChanges('meta1')),
      );
      expectObservable(item$).toBe('250ms b', {b: 10});
    });
  });

  test('Run onInit with states shorthand', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      states: ['meta1'],
      onInit: {
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:meta1',
        handler({controller}) {
          return of(undefined).pipe(
            delay(250),
            tap(() => controller.setAll('out1', 10)),
          );
        },
      },
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const rnode = tree.nodeTree.root;
      expectObservable(rnode.getItem().getStateStore().getStateChanges('meta1')).toBe('a 249ms b', {a: undefined, b: 10});
    });
  });

  const initChainConfig = (runOnInit: boolean): PipelineConfiguration => ({
    id: 'pipeline1',
    type: 'static',
    steps: [
      {id: 'step1', nqName: 'LibTests:TestAdd2'},
      {id: 'step2', nqName: 'LibTests:TestMul2'},
    ],
    states: ['s'],
    onInit: {
      id: 'init',
      from: [],
      to: 'out:s',
      handler({controller}) {
        controller.setAll('out', 5);
      },
    },
    links: [
      {id: 'l1', from: 'in1:s', to: 'out1:step1/a', runOnInit},
      {id: 'l2', from: 'in2:step1/a', to: 'out2:step2/a', runOnInit},
    ],
  });

  const runOnInitConfig: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'static',
    steps: [
      {id: 'step1', nqName: 'LibTests:TestAdd2'},
      {id: 'step2', nqName: 'LibTests:TestMul2'},
    ],
    links: [
      {
        id: 'l0', from: [], to: 'out0:step1/a', runOnInit: true,
        handler({controller}) {
          controller.setAll('out0', 7);
        },
      },
      {id: 'l2', from: 'in2:step1/a', to: 'out2:step2/a'},
    ],
  };

  test('onInit writes reach plain data links', async () => {
    const pconf = await getProcessedConfig(initChainConfig(false));
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const root = tree.nodeTree.root.getItem().getStateStore();
      const step1 = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const step2 = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore();
      const snap = () => snapshots.push([root.getState('s'), step1.getState('a'), step2.getState('a')]);
      cold('-a').subscribe(snap);
      cold('--a').subscribe(() => root.setState('s', 9));
      cold('---a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [[5, 5, 5], [9, 9, 9]]);
  });

  test('onInit writes reach runOnInit links and chain through them', async () => {
    const pconf = await getProcessedConfig(initChainConfig(true));
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const root = tree.nodeTree.root.getItem().getStateStore();
      const step1 = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const step2 = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore();
      cold('-a').subscribe(() => snapshots.push([root.getState('s'), step1.getState('a'), step2.getState('a')]));
    });
    expectDeepEqual(snapshots, [[5, 5, 5]]);
  });

  test('runOnInit link writes reach plain data links', async () => {
    const pconf = await getProcessedConfig(runOnInitConfig);
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const step1 = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const step2 = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore();
      const snap = () => snapshots.push([step1.getState('a'), step2.getState('a')]);
      cold('-a').subscribe(snap);
      cold('--a').subscribe(() => step1.setState('a', 8));
      cold('---a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [[7, 7], [8, 8]]);
  });

  // the same pipeline as the root, as a static nested step, and as a dynamic item added after init
  const placements = ['root', 'nested', 'added'] as const;
  type Placement = typeof placements[number];

  const placed = (inner: PipelineConfiguration, placement: Placement): PipelineConfiguration =>
    placement === 'root' ? inner :
      placement === 'nested' ? {id: 'outer', type: 'static', steps: [inner]} :
        {id: 'outer', type: 'parallel', stepTypes: [inner]};

  async function valuesAfterInit(inner: PipelineConfiguration, read: (stores: any[]) => any[]) {
    const results: Record<string, any[]> = {};
    for (const placement of placements) {
      const pconf = await getProcessedConfig(placed(inner, placement));
      testScheduler.run((helpers) => {
        const {cold} = helpers;
        const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
        tree.init().subscribe();
        if (placement === 'added')
          cold('-a').subscribe(() => tree.addSubTree(tree.nodeTree.root.getItem().uuid, inner.id!, 0).subscribe());
        cold('100ms a').subscribe(() => {
          const base = placement === 'root' ? [] : [{idx: 0}];
          const stores = [base, [...base, {idx: 0}], [...base, {idx: 1}]]
            .map((path) => tree.nodeTree.getNode(path).getItem().getStateStore());
          results[placement] = read(stores);
        });
      });
    }
    return results;
  }

  const sameEverywhere = (expected: any[]) => Object.fromEntries(placements.map((placement) => [placement, expected]));

  test('onInit writes reach plain data links at root, nested and added', async () => {
    const results = await valuesAfterInit(initChainConfig(false),
      ([pipeline, step1, step2]) => [pipeline.getState('s'), step1.getState('a'), step2.getState('a')]);
    expectDeepEqual(results, sameEverywhere([5, 5, 5]));
  });

  test('runOnInit link writes reach plain data links at root, nested and added', async () => {
    const results = await valuesAfterInit(runOnInitConfig,
      ([, step1, step2]) => [step1.getState('a'), step2.getState('a')]);
    expectDeepEqual(results, sameEverywhere([7, 7]));
  });

  test('onInit writes chain through runOnInit links at root, nested and added', async () => {
    const results = await valuesAfterInit(initChainConfig(true),
      ([pipeline, step1, step2]) => [pipeline.getState('s'), step1.getState('a'), step2.getState('a')]);
    expectDeepEqual(results, sameEverywhere([5, 5, 5]));
  });

  test('A plain link from the workflow fills an added item', async () => {
    const pconf = await getProcessedConfig({
      id: 'outer',
      type: 'parallel',
      stepTypes: [{
        id: 'pipeline1',
        type: 'static',
        steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
      }],
      states: ['src'],
      links: [{id: 'feed', base: 'base:expand(pipeline1)', from: 'in:src', to: 'out:same(@base)/step1/a'}],
    });
    let value: any;
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      cold('-a').subscribe(() => tree.nodeTree.root.getItem().getStateStore().setState('src', 3));
      cold('--a').subscribe(() => tree.addSubTree(tree.nodeTree.root.getItem().uuid, 'pipeline1', 0).subscribe());
      cold('100ms a').subscribe(() => {
        value = tree.nodeTree.getNode([{idx: 0}, {idx: 0}]).getItem().getStateStore().getState('a');
      });
    });
    expectDeepEqual(value, 3);
  });

  test('Run onReturn hook', async () => {
    const pconf = await getProcessedConfig(config4);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.root.getChild({idx: 0}).getItem() as FuncCallNode;
      node.getStateStore().run({res: 10}, 10).subscribe();
      cold('20ms a').subscribe(() => {
        tree.returnResult().subscribe();
      });
      expectObservable(tree.result$).toBe('20ms a', {a: 10});
    });
  });
});

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

  test('onInit writes do not reach plain data links', async () => {
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
    expectDeepEqual(snapshots, [[5, undefined, undefined], [9, 9, 9]]);
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

  test('runOnInit link writes do not reach plain data links', async () => {
    const pconf = await getProcessedConfig({
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
    });
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
    expectDeepEqual(snapshots, [[7, undefined], [8, 8]]);
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

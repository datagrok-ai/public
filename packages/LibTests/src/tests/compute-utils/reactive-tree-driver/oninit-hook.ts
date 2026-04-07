import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {of} from 'rxjs';
import {delay, map, switchMap, tap} from 'rxjs/operators';
import {FuncCallNode, StaticPipelineNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';


category('ComputeUtils: Driver hooks running', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
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

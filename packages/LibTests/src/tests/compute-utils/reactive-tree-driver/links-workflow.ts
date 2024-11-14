import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {of} from 'rxjs';
import {concatMap, delay, map, switchMap} from 'rxjs/operators';
import {makeValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {RunHelpers} from 'rxjs/internal/testing/TestScheduler';
import {NodeAddress} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/BaseTree';


category('ComputeUtils: Driver workflow test', async () => {
  let testScheduler: TestScheduler;

  const config1: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'static',
    steps: [
      {
        id: 'stepAdd',
        nqName: 'LibTests:TestAdd2',
      },
      {
        id: 'stepSub',
        nqName: 'LibTests:TestSub2',
      },
      {
        id: 'pipeline3',
        type: 'sequential',
        stepTypes: [
          {
            id: 'stepAdd',
            nqName: 'LibTests:TestAdd2',
          },
          {
            id: 'stepSub',
            nqName: 'LibTests:TestSub2',
          },
          {
            id: 'stepMul',
            nqName: 'LibTests:TestMul2',
          },
          {
            id: 'stepDiv',
            nqName: 'LibTests:TestDiv2',
          },
        ],
        initialSteps: [
          {
            id: 'stepAdd',
          },
          {
            id: 'stepSub',
          },
          {
            id: 'stepMul',
          },
          {
            id: 'stepDiv',
          },
        ],
        links: [
          {
            id: 'link31',
            base: 'base:expand(stepSub)',
            from: 'in1:same(@base,stepSub)/res',
            to: 'out1:after+(@base,stepMul)/a',
          },
          {
            id: 'link32',
            type: 'validator',
            base: 'base:expand(stepSub)',
            from: 'in1:same(@base,stepSub)/a',
            to: 'out1:same(@base,stepSub)/a',
            handler({controller}) {
              controller.setValidation('out1', makeValidationResult({warnings: ['warning from link 32']}));
            },
          },
          {
            id: 'link33',
            type: 'meta',
            base: 'base:expand(stepSub)',
            from: 'in1:same(@base,stepSub)/a',
            to: 'out1:same(@base,stepSub)/a',
            handler({controller}) {
              controller.setViewMeta('out1', {payload: 'meta from link33'});
            },
          },
          {
            id: 'link34',
            type: 'validator',
            base: 'base:expand(stepDiv)',
            from: 'in1:same(@base,stepDiv)/a',
            to: 'out1:same(@base,stepDiv)/a',
            handler({controller}) {
              controller.setValidation('out1', makeValidationResult({warnings: ['warning from link 34']}));
            },
          },
          {
            id: 'link35',
            type: 'meta',
            base: 'base:expand(stepDiv)',
            from: 'in1:same(@base,stepDiv)/a',
            to: 'out1:same(@base,stepDiv)/a',
            handler({controller}) {
              controller.setViewMeta('out1', {payload: 'meta from link35'});
            },
          },
        ],
      },
      {
        id: 'stepMul',
        nqName: 'LibTests:TestMul2',
      },
      {
        id: 'pipeline5',
        type: 'parallel',
        stepTypes: [
          {
            id: 'stepAdd',
            nqName: 'LibTests:TestAdd2',
          },
          {
            id: 'stepSub',
            nqName: 'LibTests:TestSub2',
          },
          {
            id: 'stepMul',
            nqName: 'LibTests:TestMul2',
          },
          {
            id: 'stepDiv',
            nqName: 'LibTests:TestDiv2',
          },
        ],
        initialSteps: [
          {
            id: 'stepAdd',
          },
          {
            id: 'stepSub',
          },
          {
            id: 'stepMul',
          },
          {
            id: 'stepDiv',
          },
        ],
      },
      {
        id: 'stepDiv',
        nqName: 'LibTests:TestDiv2',
      },
    ],
    links: [
      {
        id: 'link1',
        from: 'in1:stepAdd/res',
        to: 'out1:stepSub/a',
      }, {
        id: 'link2',
        from: 'in1:stepSub/res',
        to: 'out1:pipeline3/all(stepSub|stepDiv)/a',
      }, {
        id: 'link3',
        from: 'in1:pipeline3/last(stepAdd|stepSub|stepMul|stepDiv)/res',
        to: 'out1:stepMul/a',
      }, {
        id: 'link4',
        from: 'in1:stepMul/res',
        to: 'out1:pipeline5/all(stepAdd|stepSub|stepMul|stepDiv)/a',
      }, {
        id: 'link5',
        from: 'in1:pipeline5/all(stepAdd|stepSub|stepMul|stepDiv)/res',
        to: 'out1:stepDiv/a',
        handler({controller}) {
          const inputs = controller.getAll<number>('in1');
          const res = inputs?.reduce((acc, num) => acc + num, 0);
          controller.setAll('out1', res);
        },
      }, {
        id: 'link6',
        type: 'meta',
        from: 'in1:stepAdd/res',
        to: 'out1:stepDiv/b',
        handler({controller}) {
          controller.setViewMeta('out1', {payload: 'meta from link6'});
        },
      }, {
        id: 'link7',
        type: 'validator',
        from: 'in1:stepAdd/res',
        to: 'out1:stepDiv/b',
        handler({controller}) {
          controller.setValidation('out1', makeValidationResult({warnings: ['warning from link 7']}));
        },
      }, {
        id: 'link8',
        type: 'meta',
        from: 'in1:stepAdd/res',
        to: 'out1:stepSub/b',
        handler({controller}) {
          controller.setViewMeta('out1', {payload: 'meta from link8'});
        },
      }, {
        id: 'link9',
        type: 'validator',
        from: 'in1:stepAdd/res',
        to: 'out1:stepSub/b',
        handler({controller}) {
          controller.setValidation('out1', makeValidationResult({warnings: ['warning from link 9']}));
        },
      }, {
        id: 'link10',
        type: 'meta',
        from: 'in1:stepSub/res',
        to: 'out1:pipeline3/all(stepSub|stepDiv)/b',
        handler({controller}) {
          controller.setViewMeta('out1', {payload: 'meta from link10'});
        },
      }, {
        id: 'link11',
        type: 'meta',
        from: 'in1:pipeline3/last(stepAdd|stepSub|stepMul|stepDiv)/res',
        to: 'out1:stepMul/a',
        handler({controller}) {
          controller.setViewMeta('out1', {payload: 'meta from link11'});
        },
      }, {
        id: 'link12',
        from: 'in1:stepAdd/res',
        to: 'out1:stepDiv/b',
      },
    ],
  };

  function mockWorkflow(tree: StateTree) {
    return of(null).pipe(
      delay(10),
      // step1
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 0}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('a', 2);
        item.getStateStore().setState('b', 2);
        return item.instancesWrapper.run({res: 4}, 10);
      }),
      delay(10),
      // step2
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 1}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('b', 1);
        return item.instancesWrapper.run({res: 3}, 10);
      }),
      delay(10),
      // step3-1
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 2}, {idx: 0}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('a', 2);
        item.getStateStore().setState('b', 2);
        return item.instancesWrapper.run({res: 4}, 10);
      }),
      delay(10),
      // step3-2
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 2}, {idx: 1}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('b', 1);
        return item.instancesWrapper.run({res: 2}, 10);
      }),
      delay(10),
      // step3-3
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 2}, {idx: 2}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('a', 2);
        item.getStateStore().setState('b', 3);
        return item.instancesWrapper.run({res: 6}, 10);
      }),
      delay(10),
      // step3-4
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 2}, {idx: 3}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('b', 6);
        return item.instancesWrapper.run({res: 0.5}, 10);
      }),
      delay(10),
      // step4
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 3}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('b', 10);
        return item.instancesWrapper.run({res: 5}, 10);
      }),
      delay(10),
      // step5-1
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 4}, {idx: 0}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('b', 2);
        return item.instancesWrapper.run({res: 7}, 10);
      }),
      delay(10),
      // step5-2
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 4}, {idx: 1}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('b', 1);
        return item.instancesWrapper.run({res: 4}, 10);
      }),
      delay(10),
      // step5-3
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 4}, {idx: 2}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('b', 2);
        return item.instancesWrapper.run({res: 10}, 10);
      }),
      delay(10),
      // step5-4
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 4}, {idx: 3}]);
        const item = inNode.getItem() as FuncCallNode;
        item.getStateStore().setState('b', 2);
        return item.instancesWrapper.run({res: 2.5}, 10);
      }),
      delay(10),
      // step6
      concatMap(() => {
        const inNode = tree.nodeTree.getNode([{idx: 5}]);
        const item = inNode.getItem() as FuncCallNode;
        return item.instancesWrapper.run({res: 5.875}, 10);
      }),
    );
  }

  interface ExpectNode {
    address: NodeAddress,
    io: Record<string, [string, Record<string, any>, string?]>,
  }

  function expectNodes(tree: StateTree, helpers: RunHelpers, payload: ExpectNode[], sub?: string) {
    for (const nodeData of payload) {
      const node = tree.nodeTree.getNode(nodeData.address);
      for (const [ioName, spec] of Object.entries(nodeData.io)) {
        const stateChanges$ = node.getItem().getStateStore().getStateChanges(ioName);
        helpers.expectObservable(stateChanges$, sub).toBe(spec[0], spec[1]);
      }
    }
  }

  function expectMeta(tree: StateTree, helpers: RunHelpers, payload: ExpectNode[], sub?: string) {
    for (const nodeData of payload) {
      const node = tree.nodeTree.getNode(nodeData.address);
      for (const [ioName, spec] of Object.entries(nodeData.io)) {
        const stateChanges$ = (node.getItem() as FuncCallNode).metaInfo$.pipe(
          switchMap((meta) => {
            return meta[ioName] ?? of(undefined);
          }),
        );
        helpers.expectObservable(stateChanges$, sub).toBe(spec[0], spec[1]);
      }
    }
  }

  function expectValidations(tree: StateTree, helpers: RunHelpers, payload: ExpectNode[], sub?: string) {
    for (const nodeData of payload) {
      const node = tree.nodeTree.getNode(nodeData.address);
      for (const [ioName, spec] of Object.entries(nodeData.io)) {
        const stateChanges$ = (node.getItem() as FuncCallNode).validationInfo$.pipe(
          map((val) => {
            return val[ioName];
          }),
        );
        helpers.expectObservable(stateChanges$, sub).toBe(spec[0], spec[1]);
      }
    }
  }

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Run static workflow', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      mockWorkflow(tree).subscribe();
      expectNodes(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a 9ms b', {a: undefined, b: 2}],
            b: ['a 9ms b', {a: undefined, b: 2}],
            res: ['a 19ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 19ms b', {a: undefined, b: 4}],
            b: ['a 29ms b', {a: undefined, b: 1}],
            res: ['a 39ms b', {a: undefined, b: 3}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a 49ms b', {a: undefined, b: 2}],
            b: ['a 49ms b', {a: undefined, b: 2}],
            res: ['a 59ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 69ms b', {a: undefined, b: 1}],
            res: ['a 79ms b', {a: undefined, b: 2}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a 79ms b 9ms c', {a: undefined, b: 2, c: 2}],
            b: ['a 89ms b', {a: undefined, b: 3}],
            res: ['a 99ms b', {a: undefined, b: 6}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 109ms b', {a: undefined, b: 6}],
            res: ['a 119ms b', {a: undefined, b: 0.5}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: 0.5}],
            b: ['a 129ms b', {a: undefined, b: 10}],
            res: ['a 139ms b', {a: undefined, b: 5}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 149ms b', {a: undefined, b: 2}],
            res: ['a 159ms b', {a: undefined, b: 7}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 169ms b', {a: undefined, b: 1}],
            res: ['a 179ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 189ms b', {a: undefined, b: 2}],
            res: ['a 199ms b', {a: undefined, b: 10}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 209ms b', {a: undefined, b: 2}],
            res: ['a 219ms b', {a: undefined, b: 2.5}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 159ms b 19ms c 19ms d 19ms e', {a: undefined, b: NaN, c: NaN, d: NaN, e: 23.5}],
            b: ['a 19ms b', {a: undefined, b: 4}],
            res: ['a 239ms b', {a: undefined, b: 5.875}],
          },
        },
      ]);
      expectMeta(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link8',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link33',
            }}],
            b: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link35',
            }}],
            b: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: {
              'payload': 'meta from link11',
            }}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link6',
            }}],
            res: ['a', {a: undefined}],
          },
        },
      ]);
      expectValidations(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {a: undefined, b: {
              'errors': [],
              'warnings': [
                {
                  'description': 'warning from link 9',
                },
              ],
              'notifications': [],
            }}],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 289ms b', {a: undefined, b: {
              'errors': [],
              'warnings': [
                {
                  'description': 'warning from link 32',
                },
              ],
              'notifications': [],
            }}],
            b: ['a 289ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 289ms b', {a: undefined, b: {
              'errors': [],
              'warnings': [
                {
                  'description': 'warning from link 34',
                },
              ],
              'notifications': [],
            }}],
            b: ['a 289ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {a: undefined, b: {
              'errors': [],
              'warnings': [
                {
                  'description': 'warning from link 7',
                },
              ],
              'notifications': [],
            }}],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
      ]);
    });
  });

  test('Simulate adding node 1', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      mockWorkflow(tree).subscribe();
      helpers.cold('500ms a').subscribe(() => {
        tree.runMutateTree({mutationRootPath: [], addIdx: 0}).subscribe();
      });
      expectNodes(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a 9ms b', {a: undefined, b: 2}],
            b: ['a 9ms b', {a: undefined, b: 2}],
            res: ['a 19ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 19ms b 479ms c', {a: undefined, b: 4, c: 4}],
            b: ['a 29ms b', {a: undefined, b: 1}],
            res: ['a 39ms b', {a: undefined, b: 3}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a 49ms b', {a: undefined, b: 2}],
            b: ['a 49ms b', {a: undefined, b: 2}],
            res: ['a 59ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 69ms b', {a: undefined, b: 1}],
            res: ['a 79ms b', {a: undefined, b: 2}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a 79ms b 9ms c', {a: undefined, b: 2, c: 2}],
            b: ['a 89ms b', {a: undefined, b: 3}],
            res: ['a 99ms b', {a: undefined, b: 6}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 109ms b', {a: undefined, b: 6}],
            res: ['a 119ms b', {a: undefined, b: 0.5}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: 0.5}],
            b: ['a 129ms b', {a: undefined, b: 10}],
            res: ['a 139ms b', {a: undefined, b: 5}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 149ms b', {a: undefined, b: 2}],
            res: ['a 159ms b', {a: undefined, b: 7}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 169ms b', {a: undefined, b: 1}],
            res: ['a 179ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 189ms b', {a: undefined, b: 2}],
            res: ['a 199ms b', {a: undefined, b: 10}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 209ms b', {a: undefined, b: 2}],
            res: ['a 219ms b', {a: undefined, b: 2.5}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 159ms b 19ms c 19ms d 19ms e', {a: undefined, b: NaN, c: NaN, d: NaN, e: 23.5}],
            b: ['a 19ms b 479ms c', {a: undefined, b: 4, c: 4}],
            res: ['a 239ms b', {a: undefined, b: 5.875}],
          },
        },
      ]);
      expectMeta(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b 479ms b', {a: undefined, b: {
              'payload': 'meta from link8',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link33',
            }}],
            b: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link35',
            }}],
            b: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b 379ms b', {a: undefined, b: {
              'payload': 'meta from link11',
            }}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b 479ms b', {a: undefined, b: {
              'payload': 'meta from link6',
            }}],
            res: ['a', {a: undefined}],
          },
        },
      ]);
      expectValidations(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b 229ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 9',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 289ms b 209ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 32',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b 209ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b 209ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 289ms b 209ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 34',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b 209ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b 209ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b 229ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 7',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
          },
        },
      ]);
    });
  });

  test('Simulate adding node 2', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      mockWorkflow(tree).subscribe();
      helpers.cold('500ms a').subscribe(() => {
        tree.runMutateTree({mutationRootPath: [], addIdx: 1}).subscribe();
      });
      expectNodes(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a 9ms b', {a: undefined, b: 2}],
            b: ['a 9ms b', {a: undefined, b: 2}],
            res: ['a 19ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 19ms b 479ms c', {a: undefined, b: 4, c: 4}],
            b: ['a 29ms b', {a: undefined, b: 1}],
            res: ['a 39ms b', {a: undefined, b: 3}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a 49ms b', {a: undefined, b: 2}],
            b: ['a 49ms b', {a: undefined, b: 2}],
            res: ['a 59ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b 459ms b', {a: undefined, b: 3}],
            b: ['a 69ms b', {a: undefined, b: 1}],
            res: ['a 79ms b', {a: undefined, b: 2}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a 79ms b 9ms c', {a: undefined, b: 2, c: 2}],
            b: ['a 89ms b', {a: undefined, b: 3}],
            res: ['a 99ms b', {a: undefined, b: 6}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b 459ms b', {a: undefined, b: 3}],
            b: ['a 109ms b', {a: undefined, b: 6}],
            res: ['a 119ms b', {a: undefined, b: 0.5}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: 0.5}],
            b: ['a 129ms b', {a: undefined, b: 10}],
            res: ['a 139ms b', {a: undefined, b: 5}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 149ms b', {a: undefined, b: 2}],
            res: ['a 159ms b', {a: undefined, b: 7}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 169ms b', {a: undefined, b: 1}],
            res: ['a 179ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 189ms b', {a: undefined, b: 2}],
            res: ['a 199ms b', {a: undefined, b: 10}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 209ms b', {a: undefined, b: 2}],
            res: ['a 219ms b', {a: undefined, b: 2.5}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 159ms b 19ms c 19ms d 19ms e', {a: undefined, b: NaN, c: NaN, d: NaN, e: 23.5}],
            b: ['a 19ms b 479ms c', {a: undefined, b: 4, c: 4}],
            res: ['a 239ms b', {a: undefined, b: 5.875}],
          },
        },
      ]);
      expectMeta(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b 479ms b', {a: undefined, b: {
              'payload': 'meta from link8',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b 459ms (bb)', {a: undefined, b: {
              'payload': 'meta from link33',
            }}],
            b: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b 459ms (bb)', {a: undefined, b: {
              'payload': 'meta from link35',
            }}],
            b: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b 379ms b', {a: undefined, b: {
              'payload': 'meta from link11',
            }}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b 479ms b', {a: undefined, b: {
              'payload': 'meta from link6',
            }}],
            res: ['a', {a: undefined}],
          },
        },
      ]);
      expectValidations(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b 229ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 9',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 289ms b 209ms b 249ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 32',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 289ms b 209ms b 249ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 34',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b 229ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 7',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
          },
        },
      ]);
    });
  });

  test('Simulate adding pipeline 3', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      mockWorkflow(tree).subscribe();
      helpers.cold('500ms a').subscribe(() => {
        tree.runMutateTree({mutationRootPath: [], addIdx: 2}).subscribe();
      });
      // expectNodes(tree, helpers, [
      //   {
      //     address: [{idx: 0}],
      //     io: {
      //       a: ['a 9ms b', {a: undefined, b: 2}],
      //       b: ['a 9ms b', {a: undefined, b: 2}],
      //       res: ['a 19ms b', {a: undefined, b: 4}],
      //     },
      //   },
      //   {
      //     address: [{idx: 1}],
      //     io: {
      //       a: ['a 19ms b', {a: undefined, b: 4}],
      //       b: ['a 29ms b', {a: undefined, b: 1}],
      //       res: ['a 39ms b', {a: undefined, b: 3}],
      //     },
      //   },
      //   {
      //     address: [{idx: 2}, {idx: 0}],
      //     io: {
      //       a: ['a 49ms b', {a: undefined, b: 2}],
      //       b: ['a 49ms b', {a: undefined, b: 2}],
      //       res: ['a 59ms b', {a: undefined, b: 4}],
      //     },
      //   },
      //   {
      //     address: [{idx: 2}, {idx: 1}],
      //     io: {
      //       a: ['a 39ms b 459ms b', {a: undefined, b: 3}],
      //       b: ['a 69ms b', {a: undefined, b: 1}],
      //       res: ['a 79ms b', {a: undefined, b: 2}],
      //     },
      //   },
      //   {
      //     address: [{idx: 2}, {idx: 2}],
      //     io: {
      //       a: ['a 79ms b 9ms c', {a: undefined, b: 2, c: 2}],
      //       b: ['a 89ms b', {a: undefined, b: 3}],
      //       res: ['a 99ms b', {a: undefined, b: 6}],
      //     },
      //   },
      //   {
      //     address: [{idx: 2}, {idx: 3}],
      //     io: {
      //       a: ['a 39ms b 459ms b', {a: undefined, b: 3}],
      //       b: ['a 109ms b', {a: undefined, b: 6}],
      //       res: ['a 119ms b', {a: undefined, b: 0.5}],
      //     },
      //   },
      //   {
      //     address: [{idx: 3}],
      //     io: {
      //       a: ['a 119ms b 379ms b', {a: undefined, b: 0.5}],
      //       b: ['a 129ms b', {a: undefined, b: 10}],
      //       res: ['a 139ms b', {a: undefined, b: 5}],
      //     },
      //   },
      //   {
      //     address: [{idx: 4}, {idx: 0}],
      //     io: {
      //       a: ['a 139ms b', {a: undefined, b: 5}],
      //       b: ['a 149ms b', {a: undefined, b: 2}],
      //       res: ['a 159ms b', {a: undefined, b: 7}],
      //     },
      //   },
      //   {
      //     address: [{idx: 4}, {idx: 1}],
      //     io: {
      //       a: ['a 139ms b', {a: undefined, b: 5}],
      //       b: ['a 169ms b', {a: undefined, b: 1}],
      //       res: ['a 179ms b', {a: undefined, b: 4}],
      //     },
      //   },
      //   {
      //     address: [{idx: 4}, {idx: 2}],
      //     io: {
      //       a: ['a 139ms b', {a: undefined, b: 5}],
      //       b: ['a 189ms b', {a: undefined, b: 2}],
      //       res: ['a 199ms b', {a: undefined, b: 10}],
      //     },
      //   },
      //   {
      //     address: [{idx: 4}, {idx: 3}],
      //     io: {
      //       a: ['a 139ms b', {a: undefined, b: 5}],
      //       b: ['a 209ms b', {a: undefined, b: 2}],
      //       res: ['a 219ms b', {a: undefined, b: 2.5}],
      //     },
      //   },
      //   {
      //     address: [{idx: 5}],
      //     io: {
      //       a: ['a 159ms b 19ms c 19ms d 19ms e', {a: undefined, b: NaN, c: NaN, d: NaN, e: 23.5}],
      //       b: ['a 19ms b 479ms c', {a: undefined, b: 4, c: 4}],
      //       res: ['a 239ms b', {a: undefined, b: 5.875}],
      //     },
      //   },
      // ]);
      // expectMeta(tree, helpers, [
      //   {
      //     address: [{idx: 0}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a', {a: undefined}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 1}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a 19ms b', {a: undefined, b: {
      //         'payload': 'meta from link8',
      //       }}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 2}, {idx: 0}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a', {a: undefined}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 2}, {idx: 1}],
      //     io: {
      //       a: ['a 39ms b 459ms (bb)', {a: undefined, b: {
      //         'payload': 'meta from link33',
      //       }}],
      //       b: ['a 39ms b 459ms b', {a: undefined, b: {
      //         'payload': 'meta from link10',
      //       }}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 2}, {idx: 2}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a', {a: undefined}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 2}, {idx: 3}],
      //     io: {
      //       a: ['a 39ms b 459ms (bb)', {a: undefined, b: {
      //         'payload': 'meta from link35',
      //       }}],
      //       b: ['a 39ms b 459ms b', {a: undefined, b: {
      //         'payload': 'meta from link10',
      //       }}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 3}],
      //     io: {
      //       a: ['a 119ms b 379ms b', {a: undefined, b: {
      //         'payload': 'meta from link11',
      //       }}],
      //       b: ['a', {a: undefined}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 4}, {idx: 0}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a', {a: undefined}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 4}, {idx: 1}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a', {a: undefined}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 4}, {idx: 2}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a', {a: undefined}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 4}, {idx: 3}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a', {a: undefined}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      //   {
      //     address: [{idx: 5}],
      //     io: {
      //       a: ['a', {a: undefined}],
      //       b: ['a 19ms b 479ms b', {a: undefined, b: {
      //         'payload': 'meta from link6',
      //       }}],
      //       res: ['a', {a: undefined}],
      //     },
      //   },
      // ]);
      expectValidations(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 9',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 289ms b 209ms b 249ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 32',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 289ms b 209ms b 249ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 34',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b 229ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 7',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b 229ms b', {a: undefined, b: undefined}],
          },
        },
      ]);
    });
  });

  test('Simulate adding node 3-3', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      mockWorkflow(tree).subscribe();
      helpers.cold('500ms a').subscribe(() => {
        tree.runMutateTree({mutationRootPath: [{idx: 2}], addIdx: 2}).subscribe();
      });
      expectNodes(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a 9ms b', {a: undefined, b: 2}],
            b: ['a 9ms b', {a: undefined, b: 2}],
            res: ['a 19ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 19ms b', {a: undefined, b: 4}],
            b: ['a 29ms b', {a: undefined, b: 1}],
            res: ['a 39ms b', {a: undefined, b: 3}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a 49ms b', {a: undefined, b: 2}],
            b: ['a 49ms b', {a: undefined, b: 2}],
            res: ['a 59ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 69ms b', {a: undefined, b: 1}],
            res: ['a 79ms b', {a: undefined, b: 2}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a 79ms b 9ms c 409ms c', {a: undefined, b: 2, c: 2}],
            b: ['a 89ms b', {a: undefined, b: 3}],
            res: ['a 99ms b', {a: undefined, b: 6}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b 459ms b', {a: undefined, b: 3}],
            b: ['a 109ms b', {a: undefined, b: 6}],
            res: ['a 119ms b', {a: undefined, b: 0.5}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: 0.5}],
            b: ['a 129ms b', {a: undefined, b: 10}],
            res: ['a 139ms b', {a: undefined, b: 5}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 149ms b', {a: undefined, b: 2}],
            res: ['a 159ms b', {a: undefined, b: 7}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 169ms b', {a: undefined, b: 1}],
            res: ['a 179ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 189ms b', {a: undefined, b: 2}],
            res: ['a 199ms b', {a: undefined, b: 10}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 209ms b', {a: undefined, b: 2}],
            res: ['a 219ms b', {a: undefined, b: 2.5}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 159ms b 19ms c 19ms d 19ms e', {a: undefined, b: NaN, c: NaN, d: NaN, e: 23.5}],
            b: ['a 19ms b', {a: undefined, b: 4, c: 4}],
            res: ['a 239ms b', {a: undefined, b: 5.875}],
          },
        },
      ]);
      expectMeta(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link8',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link33',
            }}],
            b: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b 459ms (bb)', {a: undefined, b: {
              'payload': 'meta from link35',
            }}],
            b: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b 379ms b', {a: undefined, b: {
              'payload': 'meta from link11',
            }}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link6',
            }}],
            res: ['a', {a: undefined}],
          },
        },
      ]);
      expectValidations(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 9',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 289ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 32',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 289ms b 209ms b 249ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 34',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 7',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
      ]);
    });
  });

  test('Simulate removing at node 3-3', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      mockWorkflow(tree).subscribe();
      helpers.cold('500ms a').subscribe(() => {
        tree.runMutateTree({mutationRootPath: [{idx: 2}], removeIdx: 2}).subscribe();
      });
      expectNodes(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a 9ms b', {a: undefined, b: 2}],
            b: ['a 9ms b', {a: undefined, b: 2}],
            res: ['a 19ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 19ms b', {a: undefined, b: 4}],
            b: ['a 29ms b', {a: undefined, b: 1}],
            res: ['a 39ms b', {a: undefined, b: 3}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a 49ms b', {a: undefined, b: 2}],
            b: ['a 49ms b', {a: undefined, b: 2}],
            res: ['a 59ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 69ms b', {a: undefined, b: 1}],
            res: ['a 79ms b', {a: undefined, b: 2}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a 79ms b 9ms c 409ms c', {a: undefined, b: 2, c: 2}],
            b: ['a 89ms b', {a: undefined, b: 3}],
            res: ['a 99ms b', {a: undefined, b: 6}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b 459ms b', {a: undefined, b: 3}],
            b: ['a 109ms b', {a: undefined, b: 6}],
            res: ['a 119ms b', {a: undefined, b: 0.5}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b 379ms b', {a: undefined, b: 0.5}],
            b: ['a 129ms b', {a: undefined, b: 10}],
            res: ['a 139ms b', {a: undefined, b: 5}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 149ms b', {a: undefined, b: 2}],
            res: ['a 159ms b', {a: undefined, b: 7}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 169ms b', {a: undefined, b: 1}],
            res: ['a 179ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 189ms b', {a: undefined, b: 2}],
            res: ['a 199ms b', {a: undefined, b: 10}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 209ms b', {a: undefined, b: 2}],
            res: ['a 219ms b', {a: undefined, b: 2.5}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 159ms b 19ms c 19ms d 19ms e', {a: undefined, b: NaN, c: NaN, d: NaN, e: 23.5}],
            b: ['a 19ms b', {a: undefined, b: 4, c: 4}],
            res: ['a 239ms b', {a: undefined, b: 5.875}],
          },
        },
      ]);
      expectMeta(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link8',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link33',
            }}],
            b: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b 459ms (bb)', {a: undefined, b: {
              'payload': 'meta from link35',
            }}],
            b: ['a 39ms b 459ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b 379ms b', {a: undefined, b: {
              'payload': 'meta from link11',
            }}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link6',
            }}],
            res: ['a', {a: undefined}],
          },
        },
      ]);
      expectValidations(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 9',
                  },
                ],
                'notifications': [],
              },
              c: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 9',
                  },
                  {
                    'description': 'warning from link 9',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 289ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 32',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 289ms b 209ms b 249ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 34',
                  },
                ],
                'notifications': [],
              },
              c: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 34',
                  },
                  {
                    'description': 'warning from link 34',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b 209ms b 249ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 7',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
      ]);
    });
  });

  test('Simulate adding node 5-3', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      mockWorkflow(tree).subscribe();
      helpers.cold('500ms a').subscribe(() => {
        tree.runMutateTree({mutationRootPath: [{idx: 4}], addIdx: 2}).subscribe();
      });
      expectNodes(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a 9ms b', {a: undefined, b: 2}],
            b: ['a 9ms b', {a: undefined, b: 2}],
            res: ['a 19ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 19ms b', {a: undefined, b: 4}],
            b: ['a 29ms b', {a: undefined, b: 1}],
            res: ['a 39ms b', {a: undefined, b: 3}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a 49ms b', {a: undefined, b: 2}],
            b: ['a 49ms b', {a: undefined, b: 2}],
            res: ['a 59ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 69ms b', {a: undefined, b: 1}],
            res: ['a 79ms b', {a: undefined, b: 2}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a 79ms b 9ms c', {a: undefined, b: 2, c: 2}],
            b: ['a 89ms b', {a: undefined, b: 3}],
            res: ['a 99ms b', {a: undefined, b: 6}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 109ms b', {a: undefined, b: 6}],
            res: ['a 119ms b', {a: undefined, b: 0.5}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: 0.5}],
            b: ['a 129ms b', {a: undefined, b: 10}],
            res: ['a 139ms b', {a: undefined, b: 5}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 149ms b', {a: undefined, b: 2}],
            res: ['a 159ms b', {a: undefined, b: 7}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 169ms b', {a: undefined, b: 1}],
            res: ['a 179ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a 139ms b 359ms b', {a: undefined, b: 5}],
            b: ['a 189ms b', {a: undefined, b: 2}],
            res: ['a 199ms b', {a: undefined, b: 10}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a 139ms b 359ms b', {a: undefined, b: 5}],
            b: ['a 209ms b', {a: undefined, b: 2}],
            res: ['a 219ms b', {a: undefined, b: 2.5}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 159ms b 19ms c 19ms d 19ms e 279ms e', {a: undefined, b: NaN, c: NaN, d: NaN, e: 23.5}],
            b: ['a 19ms b', {a: undefined, b: 4, c: 4}],
            res: ['a 239ms b', {a: undefined, b: 5.875}],
          },
        },
      ]);
      expectMeta(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link8',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link33',
            }}],
            b: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link35',
            }}],
            b: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: {
              'payload': 'meta from link11',
            }}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link6',
            }}],
            res: ['a', {a: undefined}],
          },
        },
      ]);
      expectValidations(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 9',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 289ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 32',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 289ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 34',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 7',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
      ]);
    });
  });

  test('Simulate removing node 5-3', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      mockWorkflow(tree).subscribe();
      helpers.cold('500ms a').subscribe(() => {
        tree.runMutateTree({mutationRootPath: [{idx: 4}], removeIdx: 2}).subscribe();
      });
      expectNodes(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a 9ms b', {a: undefined, b: 2}],
            b: ['a 9ms b', {a: undefined, b: 2}],
            res: ['a 19ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 19ms b', {a: undefined, b: 4}],
            b: ['a 29ms b', {a: undefined, b: 1}],
            res: ['a 39ms b', {a: undefined, b: 3}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a 49ms b', {a: undefined, b: 2}],
            b: ['a 49ms b', {a: undefined, b: 2}],
            res: ['a 59ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 69ms b', {a: undefined, b: 1}],
            res: ['a 79ms b', {a: undefined, b: 2}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a 79ms b 9ms c', {a: undefined, b: 2, c: 2}],
            b: ['a 89ms b', {a: undefined, b: 3}],
            res: ['a 99ms b', {a: undefined, b: 6}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b', {a: undefined, b: 3}],
            b: ['a 109ms b', {a: undefined, b: 6}],
            res: ['a 119ms b', {a: undefined, b: 0.5}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: 0.5}],
            b: ['a 129ms b', {a: undefined, b: 10}],
            res: ['a 139ms b', {a: undefined, b: 5}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 149ms b', {a: undefined, b: 2}],
            res: ['a 159ms b', {a: undefined, b: 7}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a 139ms b', {a: undefined, b: 5}],
            b: ['a 169ms b', {a: undefined, b: 1}],
            res: ['a 179ms b', {a: undefined, b: 4}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a 139ms b 359ms b', {a: undefined, b: 5}],
            b: ['a 189ms b', {a: undefined, b: 2}],
            res: ['a 199ms b', {a: undefined, b: 10}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a 139ms b 359ms b', {a: undefined, b: 5}],
            b: ['a 209ms b', {a: undefined, b: 2}],
            res: ['a 219ms b', {a: undefined, b: 2.5}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 159ms b 19ms c 19ms d 19ms e 279ms e', {a: undefined, b: NaN, c: NaN, d: NaN, e: 23.5}],
            b: ['a 19ms b', {a: undefined, b: 4, c: 4}],
            res: ['a 239ms b', {a: undefined, b: 5.875}],
          },
        },
      ]);
      expectMeta(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link8',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link33',
            }}],
            b: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link35',
            }}],
            b: ['a 39ms b', {a: undefined, b: {
              'payload': 'meta from link10',
            }}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a 119ms b', {a: undefined, b: {
              'payload': 'meta from link11',
            }}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a 19ms b', {a: undefined, b: {
              'payload': 'meta from link6',
            }}],
            res: ['a', {a: undefined}],
          },
        },
      ]);
      expectValidations(tree, helpers, [
        {
          address: [{idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 1}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 9',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 1}],
          io: {
            a: ['a 289ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 32',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 2}, {idx: 3}],
          io: {
            a: ['a 289ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 34',
                  },
                ],
                'notifications': [],
              },
            }],
            b: ['a 289ms b', {a: undefined, b: undefined}],
            res: ['a 289ms b', {a: undefined, b: undefined}],
          },
        },
        {
          address: [{idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 0}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 1}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 2}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 4}, {idx: 3}],
          io: {
            a: ['a', {a: undefined}],
            b: ['a', {a: undefined}],
            res: ['a', {a: undefined}],
          },
        },
        {
          address: [{idx: 5}],
          io: {
            a: ['a 269ms b', {a: undefined, b: undefined}],
            b: ['a 269ms b', {
              a: undefined,
              b: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 7',
                  },
                ],
                'notifications': [],
              },
              c: {
                'errors': [],
                'warnings': [
                  {
                    'description': 'warning from link 7',
                  },
                  {
                    'description': 'warning from link 7',
                  },
                ],
                'notifications': [],
              },
            }],
            res: ['a 269ms b', {a: undefined, b: undefined}],
          },
        },
      ]);
    });
  });

});

import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {of} from 'rxjs';
import {concatMap, delay} from 'rxjs/operators';
import {makeValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {expectTreeConsistency, expectTreeData, expectTreeMeta, expectTreeValidations, getTreeStates, runRXTreeSnapshotTest} from '../../../test-utils';

category('ComputeUtils: Driver workflow test', async () => {
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
            defaultRestrictions: {
              'out1': 'restricted',
            },
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
        defaultRestrictions: {
          'out1': 'restricted',
        },
      }, {
        id: 'link2',
        from: 'in1:stepSub/res',
        to: 'out1:pipeline3/all(stepSub|stepDiv)/a',
        defaultRestrictions: {
          'out1': 'restricted',
        },
      }, {
        id: 'link3',
        from: 'in1:pipeline3/last(stepAdd|stepSub|stepMul|stepDiv)/res',
        to: 'out1:stepMul/a',
        defaultRestrictions: {
          'out1': 'restricted',
        },

      }, {
        id: 'link4',
        from: 'in1:stepMul/res',
        to: 'out1:pipeline5/all(stepAdd|stepSub|stepMul|stepDiv)/a',
        defaultRestrictions: {
          'out1': 'restricted',
        },
      }, {
        id: 'link5',
        from: 'in1:pipeline5/all(stepAdd|stepSub|stepMul|stepDiv)/res',
        to: 'out1:stepDiv/a',
        handler({controller}) {
          const inputs = controller.getAll<number>('in1');
          const res = inputs?.reduce((acc, num) => acc + num, 0);
          controller.setAll('out1', res == null || isNaN(res) ? null : res, 'restricted');
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
        defaultRestrictions: {
          'out1': 'restricted',
        },
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

  test('Static workflow data', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Static workflow data', (expectObservable) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeData(tree$, expectObservable);
    });
  });

  test('Static workflow meta', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Static workflow meta', (expectObservable) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeMeta(tree$, expectObservable);
    });
  });

  test('Static workflow validations', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Static workflow validations', (expectObservable) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeValidations(tree$, expectObservable);
    });
  });

  test('Static workflow consistency', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Static workflow consistency', (expectObservable) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeConsistency(tree$, expectObservable);
    });
  });

  test('Remove node 1 data', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Remove node 1 data', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeData(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 1}]);
        tree.removeSubtree(target.uuid).subscribe();
      });
    });
  });

  test('Remove node 1 meta', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Remove node 1 meta', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeMeta(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 1}]);
        tree.removeSubtree(target.uuid).subscribe();
      });
    });
  });

  test('Remove node 1 validations', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Remove node 1 validations', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeValidations(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 1}]);
        tree.removeSubtree(target.uuid).subscribe();
      });
    });
  });

  test('Remove node 1 consistency', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Remove node 1 consistency', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeConsistency(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 1}]);
        tree.removeSubtree(target.uuid).subscribe();
      });
    });
  });

  test('Remove node 4 data', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Remove node 4 data', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeData(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 3}]);
        tree.removeSubtree(target.uuid).subscribe();
      });
    });
  });

  test('Remove node 4 meta', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Remove node 4 meta', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeMeta(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 3}]);
        tree.removeSubtree(target.uuid).subscribe();
      });
    });
  });

  test('Remove node 4 validations', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Remove node 4 validations', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeValidations(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 3}]);
        tree.removeSubtree(target.uuid).subscribe();
      });
    });
  });

  test('Remove node 4 consistency', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Remove node 4 consistency', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeConsistency(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 3}]);
        tree.removeSubtree(target.uuid).subscribe();
      });
    });
  });

  test('Add node 3-5 data', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Add node 3-5 data', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeData(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 2}]);
        tree.addSubTree(target.uuid, 'stepDiv', 4).subscribe();
      });
    });
  });

  test('Add node 3-5 meta', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Add node 3-5 meta', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeMeta(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 2}]);
        tree.addSubTree(target.uuid, 'stepDiv', 4).subscribe();
      });
    });
  });

  test('Add node 3-5 validations', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Add node 3-5 validations', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeValidations(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 2}]);
        tree.addSubTree(target.uuid, 'stepDiv', 4).subscribe();
      });
    });
  });

  test('Add node 3-5 consistency', async () => {
    const pconf = await getProcessedConfig(config1);
    await runRXTreeSnapshotTest('Add node 3-5 consistency', (expectObservable, cold) => {
      const [tree, tree$] = getTreeStates(pconf);
      mockWorkflow(tree).subscribe();
      expectTreeConsistency(tree$, expectObservable);
      cold('500ms a').subscribe(() => {
        const target = tree.nodeTree.getItem([{idx: 2}]);
        tree.addSubTree(target.uuid, 'stepDiv', 4).subscribe();
      });
    });
  });
});

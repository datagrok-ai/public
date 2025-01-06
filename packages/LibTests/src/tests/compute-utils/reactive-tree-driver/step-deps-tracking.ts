import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';


category('ComputeUtils: Driver steps dependencies tracking', async () => {
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
        nqName: 'LibTests:TestSub2',
      },
      {
        id: 'step3',
        nqName: 'LibTests:TestMul2',
      },
      {
        id: 'step4',
        nqName: 'LibTests:TestDiv2',
      },
    ],
    links: [
      {
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
        defaultRestrictions: {
          out1: 'restricted',
        },
      },
      {
        id: 'link2',
        from: 'in1:step2/res',
        to: 'out1:step3/a',
        defaultRestrictions: {
          out1: 'restricted',
        },
      },
      {
        id: 'link3',
        from: 'in1:step3/res',
        to: 'out1:step4/a',
        defaultRestrictions: {
          out1: 'restricted',
        },
      },
    ],
  };

  test('Track direct dependencies', async () => {
    const pconf = await getProcessedConfig(config1);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node2 = tree.nodeTree.getNode([{idx: 1}]);
      const node3 = tree.nodeTree.getNode([{idx: 2}]);
      const node4 = tree.nodeTree.getNode([{idx: 3}]);
      const state1$ = tree.getFuncCallStates()[node1.getItem().uuid];
      const state2$ = tree.getFuncCallStates()[node2.getItem().uuid];
      const state3$ = tree.getFuncCallStates()[node3.getItem().uuid];
      const state4$ = tree.getFuncCallStates()[node4.getItem().uuid];

      cold('-a').subscribe(() => {
        tree.runStep(node1.getItem().uuid, {res: 1}).subscribe();
      });
      cold('---a').subscribe(() => {
        tree.runStep(node2.getItem().uuid, {res: 2}).subscribe();
      });
      cold('-----a').subscribe(() => {
        tree.runStep(node3.getItem().uuid, {res: 3}).subscribe();
      });

      expectObservable(state1$).toBe('a(bcdef)', {
        a: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        b: {
          'isRunning': true,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        c: {
          'isRunning': true,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        d: {
          'isRunning': true,
          'isRunnable': false,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
        e: {
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
        f: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
      });
      expectObservable(state2$).toBe('ax-(bcdef)', {
        a: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node1.getItem().uuid],
        },
        x: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        b: {
          'isRunning': true,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        c: {
          'isRunning': true,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        d: {
          'isRunning': true,
          'isRunnable': false,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
        e: {
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
        f: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
      });
      expectObservable(state3$).toBe('ax-y-(bcdef)', {
        a: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node2.getItem().uuid],
        },
        x: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node2.getItem().uuid],
        },
        y: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        b: {
          'isRunning': true,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        c: {
          'isRunning': true,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        d: {
          'isRunning': true,
          'isRunnable': false,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
        e: {
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
        f: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
      });
      expectObservable(state4$).toBe('ax-y-z', {
        a: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        x: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        y: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        z: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
      });
    });
  });

  test('Track transitive dependencies', async () => {
    const pconf = await getProcessedConfig(config1);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node2 = tree.nodeTree.getNode([{idx: 1}]);
      const node3 = tree.nodeTree.getNode([{idx: 2}]);
      const node4 = tree.nodeTree.getNode([{idx: 3}]);
      const state4$ = tree.getFuncCallStates()[node4.getItem().uuid];

      cold('-a').subscribe(() => {
        tree.runStep(node1.getItem().uuid, {res: 1}).subscribe();
      });
      cold('---a').subscribe(() => {
        tree.runStep(node2.getItem().uuid, {res: 2}).subscribe();
      });
      cold('-----a').subscribe(() => {
        tree.runStep(node3.getItem().uuid, {res: 3}).subscribe();
      });
      cold('-------a').subscribe(() => {
        tree.runStep(node1.getItem().uuid, {res: 4}).subscribe();
      });
      expectObservable(state4$).toBe('ax-y-z', {
        a: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        x: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        y: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        z: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
      });
      expectObservable(((node4.getItem()) as FuncCallNode).consistencyInfo$).toBe('a 4ms b', {
        a: {},
        b: {
          a: {
            'restriction': 'restricted',
            'inconsistent': false,
            'assignedValue': 3,
          },
        },
      });
    });
  });

  test('Update restriction on runned steps', async () => {
    const pconf = await getProcessedConfig(config1);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node2 = tree.nodeTree.getNode([{idx: 1}]);
      const node3 = tree.nodeTree.getNode([{idx: 2}]);
      const node4 = tree.nodeTree.getNode([{idx: 3}]);
      const state4$ = tree.getFuncCallStates()[node4.getItem().uuid];

      cold('-a').subscribe(() => {
        tree.runStep(node1.getItem().uuid, {res: 1}).subscribe();
      });
      cold('---a').subscribe(() => {
        tree.runStep(node2.getItem().uuid, {res: 2}).subscribe();
      });
      cold('-----a').subscribe(() => {
        tree.runStep(node3.getItem().uuid, {res: 3}).subscribe();
      });
      cold('-------a').subscribe(() => {
        tree.runStep(node1.getItem().uuid, {res: 4}).subscribe();
      });
      cold('---------a').subscribe(() => {
        tree.runStep(node4.getItem().uuid, {res: 5}).subscribe();
      });
      cold('-----------a').subscribe(() => {
        tree.runStep(node3.getItem().uuid, {res: 10}).subscribe();
      });

      expectObservable(state4$).toBe('ax-y-z---(bcdef)', {
        a: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        x: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        y: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },
        z: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        b:
         {
           'isRunning': true,
           'isRunnable': true,
           'isOutputOutdated': true,
           'pendingDependencies': [],
         },
        c: {
          'isRunning': true,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        d: {
          'isRunning': true,
          'isRunnable': false,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
        e: {
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': false,
          'pendingDependencies': [],
        },
        f:
          {
            'isRunning': false,
            'isRunnable': true,
            'isOutputOutdated': false,
            'pendingDependencies': [],
          },
      });
      expectObservable(((node4.getItem()) as FuncCallNode).consistencyInfo$).toBe('a 4ms b---b-c', {
        a: {},
        b: {
          a: {
            'restriction': 'restricted',
            'inconsistent': false,
            'assignedValue': 3,
          },
        },
        c: {
          a: {
            'restriction': 'restricted',
            'inconsistent': true,
            'assignedValue': 10,
          },
        },
      });
    });
  });
});

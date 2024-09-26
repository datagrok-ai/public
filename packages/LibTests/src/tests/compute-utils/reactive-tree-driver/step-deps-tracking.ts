import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {of, Subject} from 'rxjs';
import {delay, mapTo, switchMap} from 'rxjs/operators';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {makeValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
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
      },
      {
        id: 'link2',
        from: 'in1:step2/res',
        to: 'out1:step3/a',
      },
      {
        id: 'link3',
        from: 'in1:step3/res',
        to: 'out1:step4/a',
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
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        e: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
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
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        e: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
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
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        e: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
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
      expectObservable(state4$).toBe('ax-y-z-b', {
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
        b: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [node3.getItem().uuid],
        },

      });
    });
  });

});

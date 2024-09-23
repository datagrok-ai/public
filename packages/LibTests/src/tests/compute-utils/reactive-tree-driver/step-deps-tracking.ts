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
      console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Track dependencies', async () => {
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
      links: [{
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step2/a',
      }],
    };

    const pconf = await getProcessedConfig(config1);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        tree.runStep(inNode.getItem().uuid).subscribe();
      });
      expectObservable((outNode.getItem() as FuncCallNode).pendingDependencies$).toBe('a b', {
        a: [inNode.getItem().uuid],
        b: [],
      });
    });
  }, {skipReason: 'TODO'});

  test('Track and update chain dependencies', async () => {
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
        {
          id: 'step2',
          nqName: 'LibTests:TestSub2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
      }, {
        id: 'link2',
        from: 'in1:step2/res',
        to: 'out1:step3/res',
      }],
    };

    const pconf = await getProcessedConfig(config1);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const midNode = tree.nodeTree.getNode([{idx: 1}]);
      const outNode = tree.nodeTree.getNode([{idx: 2}]);
      cold('-a').subscribe(() => {
        tree.runStep(inNode.getItem().uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        tree.runStep(midNode.getItem().uuid).subscribe();
      });
      expectObservable((outNode.getItem() as FuncCallNode).pendingDependencies$).toBe('a bc', {
        a: [inNode.getItem().uuid, midNode.getItem().uuid],
        b: [midNode.getItem().uuid],
        c: [],
      });
    });
  }, {skipReason: 'TODO'});

});

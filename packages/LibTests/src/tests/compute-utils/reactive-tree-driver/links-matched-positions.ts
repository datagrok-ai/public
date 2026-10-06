import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {MatchedNodeInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/RuntimeControllers';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';

category('ComputeUtils: Driver links matched positions', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Report positions for data link ios', async () => {
    let in1Pos: MatchedNodeInfo[] | undefined;
    let in2Pos: MatchedNodeInfo[] | undefined;
    let out1Pos: MatchedNodeInfo[] | undefined;
    let basePos: MatchedNodeInfo | undefined | 'unset' = 'unset';
    const config: PipelineConfiguration = {
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
          id: 'step3',
          nqName: 'LibTests:TestSub2',
        },
      ],
      links: [{
        id: 'link1',
        from: ['in1:step1/res', 'in2:step2/res'],
        to: 'out1:step3/a',
        handler({controller}) {
          in1Pos = controller.getMatchedPositions('in1');
          in2Pos = controller.getMatchedPositions('in2');
          out1Pos = controller.getMatchedPositions('out1');
          basePos = controller.getBasePosition();
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('res', 1);
      });
      cold('--z').subscribe(() => {
        expectDeepEqual(in1Pos, [{path: [{id: 'step1', idx: 0}], position: 0, ioName: 'res'}], {prefix: 'in1'});
        expectDeepEqual(in2Pos, [{path: [{id: 'step2', idx: 1}], position: 1, ioName: 'res'}], {prefix: 'in2'});
        expectDeepEqual(out1Pos, [{path: [{id: 'step3', idx: 2}], position: 2, ioName: 'a'}], {prefix: 'out1'});
        expectDeepEqual(basePos, undefined, {prefix: 'base'});
      });
    });
  });

  test('Report positions for multiple matched nodes', async () => {
    let inPos: MatchedNodeInfo[] | undefined;
    let inVals: number[] | undefined;
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'pipelinePar',
          type: 'parallel',
          stepTypes: [
            {
              id: 'stepAdd',
              nqName: 'LibTests:TestAdd2',
            },
            {
              id: 'stepMul',
              nqName: 'LibTests:TestMul2',
            },
          ],
          initialSteps: [
            {id: 'stepAdd'},
            {id: 'stepMul'},
            {id: 'stepAdd'},
          ],
        },
        {
          id: 'step3',
          nqName: 'LibTests:TestSub2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:pipelinePar/all(stepAdd)/res',
        to: 'out1:step3/a',
        handler({controller}) {
          inPos = controller.getMatchedPositions('in1');
          inVals = controller.getAll('in1');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const add1 = tree.nodeTree.getNode([{idx: 0}, {idx: 0}]);
      const add2 = tree.nodeTree.getNode([{idx: 0}, {idx: 2}]);
      cold('-a').subscribe(() => {
        add1.getItem().getStateStore().setState('res', 10);
      });
      cold('--a').subscribe(() => {
        add2.getItem().getStateStore().setState('res', 30);
      });
      cold('---z').subscribe(() => {
        expectDeepEqual(inPos, [
          {path: [{id: 'pipelinePar', idx: 0}, {id: 'stepAdd', idx: 0}], position: 0, ioName: 'res'},
          {path: [{id: 'pipelinePar', idx: 0}, {id: 'stepAdd', idx: 2}], position: 2, ioName: 'res'},
        ], {prefix: 'in1 positions'});
        expectDeepEqual(inVals, [10, 30], {prefix: 'in1 values'});
      });
    });
  });

  test('Report base position for base instantiated links', async () => {
    const basePositions: MatchedNodeInfo[] = [];
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'sequential',
      stepTypes: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      initialSteps: [
        {id: 'step1'},
        {id: 'step2'},
        {id: 'step2'},
      ],
      links: [{
        id: 'link1',
        base: 'base:expand(step2)',
        from: 'from:same(@base, step2)/a',
        to: 'to:same(@base, step2)/b',
        handler({controller}) {
          const basePos = controller.getBasePosition();
          if (basePos)
            basePositions.push(basePos);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const mul1 = tree.nodeTree.getNode([{idx: 1}]);
      const mul2 = tree.nodeTree.getNode([{idx: 2}]);
      cold('-a').subscribe(() => {
        mul1.getItem().getStateStore().setState('a', 1);
      });
      cold('--a').subscribe(() => {
        mul2.getItem().getStateStore().setState('a', 2);
      });
      cold('---z').subscribe(() => {
        const sorted = [...basePositions].sort((p1, p2) => p1.position - p2.position);
        expectDeepEqual(sorted[0], {path: [{id: 'step2', idx: 1}], position: 1, ioName: undefined},
          {prefix: 'first base'});
        expectDeepEqual(sorted[sorted.length - 1], {path: [{id: 'step2', idx: 2}], position: 2, ioName: undefined},
          {prefix: 'last base'});
      });
    });
  });

  test('Report positions in selector links', async () => {
    let outPos: MatchedNodeInfo[] | undefined;
    const config: PipelineConfiguration = {
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
        id: 'selector',
        type: 'selector',
        from: 'in:step2/a',
        to: 'out1:title',
        handler({controller}) {
          outPos = controller.getMatchedPositions('out1');
          const [inPos] = controller.getMatchedPositions('in');
          controller.setDescriptionItem('out1', `Step ${inPos.position + 1}`);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 1}]);
      const pipeline = tree.nodeTree.root;
      cold('-a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 1);
      });
      cold('--z').subscribe(() => {
        expectDeepEqual(outPos, [{path: [], position: -1, ioName: 'title'}], {prefix: 'out1'});
      });
      expectObservable(pipeline.getItem().nodeDescription.getStateChanges('title')).toBe('ab',
        {a: undefined, b: 'Step 2'});
    });
  });

  test('Throw on unknown io names', async () => {
    let error: Error | undefined;
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/a',
        to: 'out1:step1/b',
        handler({controller}) {
          try {
            controller.getMatchedPositions('nope');
          } catch (e: any) {
            error = e;
          }
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 1);
      });
      cold('--z').subscribe(() => {
        expectDeepEqual(!!error?.message?.includes('unknown io nope'), true, {prefix: 'error message'});
      });
    });
  });
});

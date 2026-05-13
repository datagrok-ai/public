import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';
import {of, Subject} from 'rxjs';
import {delay, filter, mapTo, take} from 'rxjs/operators';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {DriverLogger} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';


// ============================================================
// Multi-hop data propagation
// ============================================================

category('ComputeUtils: Driver multi-hop data propagation', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('4-step chain propagates data end-to-end via default handlers', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestSub2'},
        {id: 'step4', nqName: 'LibTests:TestDiv2'},
      ],
      links: [
        {id: 'link1', from: 'in1:step1/b', to: 'out1:step2/a'},
        {id: 'link2', from: 'in1:step2/a', to: 'out1:step3/a'},
        {id: 'link3', from: 'in1:step3/a', to: 'out1:step4/a'},
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({expectObservable, cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node4 = tree.nodeTree.getNode([{idx: 3}]);

      cold('-a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 42);
      });

      expectObservable(node4.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 42});
    });
  });

  test('4-step chain with custom handlers transforms data at each hop', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestSub2'},
        {id: 'step4', nqName: 'LibTests:TestDiv2'},
      ],
      links: [
        {
          id: 'link1', from: 'in1:step1/b', to: 'out1:step2/a',
          handler({controller}) {
            controller.setAll('out1', controller.getFirst('in1')! * 2);
          },
        },
        {
          id: 'link2', from: 'in1:step2/a', to: 'out1:step3/a',
          handler({controller}) {
            controller.setAll('out1', controller.getFirst('in1')! + 10);
          },
        },
        {
          id: 'link3', from: 'in1:step3/a', to: 'out1:step4/a',
          handler({controller}) {
            controller.setAll('out1', controller.getFirst('in1')! - 1);
          },
        },
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({expectObservable, cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node4 = tree.nodeTree.getNode([{idx: 3}]);

      // 5 * 2 = 10 -> 10 + 10 = 20 -> 20 - 1 = 19
      cold('-a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 5);
      });

      expectObservable(node4.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 19});
    });
  });

  test('Re-trigger source in 4-step chain replaces final value', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestSub2'},
        {id: 'step4', nqName: 'LibTests:TestDiv2'},
      ],
      links: [
        {id: 'link1', from: 'in1:step1/b', to: 'out1:step2/a'},
        {id: 'link2', from: 'in1:step2/a', to: 'out1:step3/a'},
        {id: 'link3', from: 'in1:step3/a', to: 'out1:step4/a'},
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({expectObservable, cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node4 = tree.nodeTree.getNode([{idx: 3}]);

      cold('-a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 1);
      });
      cold('--a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 99);
      });

      expectObservable(node4.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a b c', {a: undefined, b: 1, c: 99});
    });
  });

  test('4-step chain with restrictions propagates restriction metadata', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestSub2'},
        {id: 'step4', nqName: 'LibTests:TestDiv2'},
      ],
      links: [
        {id: 'link1', from: 'in1:step1/b', to: 'out1:step2/a', defaultRestrictions: {out1: 'restricted'}},
        {id: 'link2', from: 'in1:step2/a', to: 'out1:step3/a', defaultRestrictions: {out1: 'restricted'}},
        {id: 'link3', from: 'in1:step3/a', to: 'out1:step4/a', defaultRestrictions: {out1: 'restricted'}},
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({expectObservable, cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node4 = tree.nodeTree.getNode([{idx: 3}]);

      cold('-a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 7);
      });

      expectObservable((node4.getItem().getStateStore() as FuncCallInstancesBridge).inputRestrictions$)
        .toBe('a b', {
          a: {},
          b: {a: {type: 'restricted', assignedValue: 7}},
        });
    });
  });
});

// ============================================================
// Chain error recovery
// ============================================================

category('ComputeUtils: Driver chain error recovery', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Error in middle of 3-step chain does not block recovery', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestSub2'},
      ],
      links: [
        {
          id: 'link1', from: 'in1:step1/b', to: 'out1:step2/a',
          handler({controller}) {
            const v = controller.getFirst('in1');
            if (v === 'bad') throw new Error('bad value');
            controller.setAll('out1', v);
          },
        },
        {id: 'link2', from: 'in1:step2/a', to: 'out1:step3/a'},
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({expectObservable, cold}) => {
      const logger = new DriverLogger();
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, logger});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node3 = tree.nodeTree.getNode([{idx: 2}]);

      // First: bad value — link1 throws, step3 never sees it
      cold('-a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 'bad');
      });
      // Second: good value — link1 recovers, propagates through link2 to step3
      cold('--a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 'good');
      });

      expectObservable(node3.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a-b', {a: undefined, b: 'good'});
    });
  });

  test('Error in async handler mid-chain: next good value completes full chain', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestSub2'},
      ],
      links: [
        {
          id: 'link1', from: 'in1:step1/b', to: 'out1:step2/a',
          handler({controller}) {
            const v = controller.getFirst('in1');
            if (v === 'bad') throw new Error('async fail');
            controller.setAll('out1', v);
            return of(null).pipe(delay(50), mapTo(void(0)));
          },
        },
        {id: 'link2', from: 'in1:step2/a', to: 'out1:step3/a'},
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({expectObservable, cold}) => {
      const logger = new DriverLogger();
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, logger});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node3 = tree.nodeTree.getNode([{idx: 2}]);

      cold('-a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 'bad');
      });
      cold('--a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 'recovered');
      });

      // The error is synchronous (throw before the async return), so the second
      // trigger's async handler completes at cold('--a') + 50ms delay
      expectObservable(node3.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a 51ms b', {a: undefined, b: 'recovered'});
    });
  });
});

// ============================================================
// Config processing: action steps
// ============================================================

category('ComputeUtils: Driver config processing: action steps', async () => {
  before(async () => {});

  test('Process config with action step', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'actions', type: 'action', friendlyName: 'Actions'} as any,
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config as any);
    // Action step should be converted to static pipeline
    const actionStep = (pconf as any).steps[1];
    expectDeepEqual(actionStep.type, 'static', {prefix: 'Action step converted to static'});
    expectDeepEqual(actionStep.isActionStep, true, {prefix: 'isActionStep flag set'});
    expectDeepEqual(actionStep.steps.length, 0, {prefix: 'No child steps'});
    expectDeepEqual(actionStep.disableHistory, true, {prefix: 'History disabled'});
    expectDeepEqual(actionStep.friendlyName, 'Actions', {prefix: 'friendlyName preserved'});
  });

  test('Action step creates valid tree node', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'actions', type: 'action'} as any,
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
    };
    const pconf = await getProcessedConfig(config as any);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    // Tree should have 3 children: step1, action, step2
    const root = tree.nodeTree.root;
    expectDeepEqual(root.getChildren().length, 3, {prefix: 'Root has 3 children'});
    // Action step node should be a StaticPipelineNode
    const actionNode = tree.nodeTree.getNode([{idx: 1}]);
    expectDeepEqual(actionNode.getItem().nodeType, 'static', {prefix: 'Action node type'});
    expectDeepEqual(actionNode.getChildren().length, 0, {prefix: 'Action node has no children'});
  });

  test('Action step config carries isActionStep through to node', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'actions', type: 'action', friendlyName: 'My Actions'} as any,
      ],
    };
    const pconf = await getProcessedConfig(config as any);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    // Verify via the node's config (avoids toState which needs FuncCall)
    const actionNode = tree.nodeTree.getNode([{idx: 1}]);
    const nodeConfig = actionNode.getItem().config as any;
    expectDeepEqual(nodeConfig.isActionStep, true, {prefix: 'isActionStep on node config'});
    expectDeepEqual(nodeConfig.disableHistory, true, {prefix: 'disableHistory on node config'});
    expectDeepEqual(nodeConfig.type, 'static', {prefix: 'type is static on node config'});
  });

  test('Links skip action step (no IO to match)', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'actions', type: 'action'} as any,
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config as any);
    const testScheduler = createTestScheduler();

    testScheduler.run(({expectObservable, cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node1 = tree.nodeTree.getNode([{idx: 0}]);
      const node3 = tree.nodeTree.getNode([{idx: 2}]);

      cold('-a').subscribe(() => {
        node1.getItem().getStateStore().setState('b', 42);
      });

      expectObservable(node3.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 42});
    });
  });
});

// ============================================================
// Link retention across multiple mutations
// ============================================================

category('ComputeUtils: Driver link retention: multiple mutations', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Links retained across add-then-add mutations', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'sequential',
      stepTypes: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      initialSteps: [
        {id: 'step1'},
      ],
      links: [
        {
          id: 'link1',
          base: 'base:expand(step1)',
          from: 'from:same(@base, step1)/res',
          to: 'to:after+(@base, step2)/a',
        },
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();

      const initialLinks = [...tree.linksState.links.values()];
      expectDeepEqual(initialLinks.length, 0, {prefix: 'Initial: no expanded links (no step2 yet)'});

      const updateFinished$ = new Subject<true>();

      // First mutation: add step2
      cold('-a').subscribe(() => {
        tree.addSubTree(tree.nodeTree.root.getItem().uuid, 'step2', 1).subscribe();
        tree.globalROLocked$.pipe(filter((x) => !x), take(1)).subscribe(() => updateFinished$.next(true));
      });

      let firstMutationLinkUuid: string;
      updateFinished$.pipe(take(1)).subscribe(() => {
        const links = [...tree.linksState.links.values()];
        expectDeepEqual(links.length, 1, {prefix: 'After first add: 1 link'});
        firstMutationLinkUuid = links[0].uuid;
      });

      // Second mutation: add another step1
      cold('---a').subscribe(() => {
        tree.addSubTree(tree.nodeTree.root.getItem().uuid, 'step1', 0).subscribe();
        tree.globalROLocked$.pipe(filter((x) => !x), take(1)).subscribe(() => updateFinished$.next(true));
      });

      updateFinished$.pipe(take(1)).subscribe(() => {
        const links = [...tree.linksState.links.values()];
        expectDeepEqual(links.length >= 1, true, {prefix: 'After second add: at least 1 link'});
        // Original link should be retained (same UUID)
        const retained = links.find((l) => l.uuid === firstMutationLinkUuid);
        expectDeepEqual(retained != null, true, {prefix: 'Original link retained'});
      });
    });
  });
});

// ============================================================
// Cross-pipeline data propagation
// ============================================================

category('ComputeUtils: Driver cross-pipeline data propagation', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Data propagates from step inside nested pipeline to outer step', async () => {
    const config: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [
        {
          id: 'inner',
          type: 'static',
          steps: [
            {id: 'step1', nqName: 'LibTests:TestAdd2'},
          ],
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'crossLink',
        from: 'in1:inner/step1/b',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({expectObservable, cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const innerStep = tree.nodeTree.getNode([{idx: 0}, {idx: 0}]);
      const outerStep = tree.nodeTree.getNode([{idx: 1}]);

      cold('-a').subscribe(() => {
        innerStep.getItem().getStateStore().setState('b', 77);
      });

      expectObservable(outerStep.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 77});
    });
  });

  test('Data propagates from outer step into nested pipeline step', async () => {
    const config: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'inner',
          type: 'static',
          steps: [
            {id: 'step2', nqName: 'LibTests:TestMul2'},
          ],
        },
      ],
      links: [{
        id: 'crossLink',
        from: 'in1:step1/b',
        to: 'out1:inner/step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run(({expectObservable, cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const outerStep = tree.nodeTree.getNode([{idx: 0}]);
      const innerStep = tree.nodeTree.getNode([{idx: 1}, {idx: 0}]);

      cold('-a').subscribe(() => {
        outerStep.getItem().getStateStore().setState('b', 33);
      });

      expectObservable(innerStep.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 33});
    });
  });
});


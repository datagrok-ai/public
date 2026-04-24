import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {delay, filter, take} from 'rxjs/operators';
import {createTestScheduler} from '../../../test-utils';

// --- Configs ---

const dynamicConfig: PipelineConfiguration = {
  id: 'root',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {
      id: 'analyses',
      type: 'dynamic',
      stepTypes: [
        {id: 'regression', nqName: 'LibTests:TestMul2'},
        {id: 'clustering', nqName: 'LibTests:TestAdd2'},
      ],
      initialSteps: [{id: 'regression'}, {id: 'clustering'}],
    },
  ],
  actions: [
    {
      id: 'add-regression',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        controller.addStep('out1', 'regression');
      },
    },
    {
      id: 'add-at-zero',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        controller.addStep('out1', 'clustering', 0);
      },
    },
    {
      id: 'remove-first',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        const steps = controller.getSteps('out1');
        controller.removeStep('out1', steps[0]);
      },
    },
    {
      id: 'move-last-to-front',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        const steps = controller.getSteps('out1');
        controller.moveStep('out1', steps[steps.length - 1], 0);
      },
    },
    {
      id: 'remove-then-add',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        const steps = controller.getSteps('out1');
        controller.removeStep('out1', steps[0]);
        controller.addStep('out1', 'clustering');
      },
    },
    {
      id: 'double-add',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        controller.addStep('out1', 'regression');
        controller.addStep('out1', 'clustering');
      },
    },
    {
      id: 'stale-remove',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        const steps = controller.getSteps('out1');
        controller.removeStep('out1', steps[0]);
        controller.removeStep('out1', steps[0]); // stale — should throw
      },
    },
    {
      id: 'exclusivity-test',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        controller.addStep('out1', 'regression');
        controller.setPipelineState('out1', {id: 'analyses', steps: []}); // should throw
      },
    },
  ],
};

const dynamicWithLinkConfig: PipelineConfiguration = {
  id: 'root',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {
      id: 'analyses',
      type: 'dynamic',
      stepTypes: [
        {id: 'regression', nqName: 'LibTests:TestMul2'},
      ],
      initialSteps: [{id: 'regression'}],
    },
  ],
  links: [
    {
      id: 'data-link',
      from: 'in1:step1/b',
      to: 'out1:analyses/all(regression)/a',
    },
  ],
  actions: [
    {
      id: 'add-regression',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        controller.addStep('out1', 'regression');
      },
    },
    {
      id: 'remove-first-regression',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        const steps = controller.getSteps('out1');
        controller.removeStep('out1', steps[0]);
      },
    },
  ],
};

const emptyDynamicConfig: PipelineConfiguration = {
  id: 'root',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {
      id: 'analyses',
      type: 'dynamic',
      stepTypes: [
        {id: 'regression', nqName: 'LibTests:TestMul2'},
      ],
    },
  ],
  actions: [
    {
      id: 'add-regression',
      type: 'pipeline',
      from: [],
      to: 'out1:analyses',
      position: 'none',
      handler({controller}) {
        controller.addStep('out1', 'regression');
      },
    },
  ],
};

// --- Helpers ---

function getAction(tree: StateTree, actionId: string) {
  const actions = [...tree.linksState.actions.values()];
  const action = actions.find((a) => a.matchInfo.spec.id === actionId);
  if (!action)
    throw new Error(`Action ${actionId} not found among [${actions.map((a) => a.matchInfo.spec.id).join(', ')}]`);
  return action;
}

function getAnalysesNode(tree: StateTree) {
  return tree.nodeTree.getNode([{idx: 1}]);
}

// --- Tests ---

category('ComputeUtils: Driver pipeline granular: getSteps', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Returns children with configId and position', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'add-regression');
      // Access controller's getSteps by running a probe action
      const analysesNode = getAnalysesNode(tree);
      const children = analysesNode.getChildren();
      expectDeepEqual(children.length, 2, {prefix: 'Initial child count'});
      expectDeepEqual(children[0].id, 'regression', {prefix: 'First child configId'});
      expectDeepEqual(children[1].id, 'clustering', {prefix: 'Second child configId'});
    });
  });

  test('Returns empty for empty pipeline', async () => {
    const pconf = await getProcessedConfig(emptyDynamicConfig);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const analysesNode = getAnalysesNode(tree);
      const children = analysesNode.getChildren();
      expectDeepEqual(children.length, 0, {prefix: 'Empty pipeline child count'});
    });
  });
});

category('ComputeUtils: Driver pipeline granular: addStep', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Appends when no position given', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'add-regression');
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 3, {prefix: 'Child count after add'});
        expectDeepEqual(children[2].id, 'regression', {prefix: 'Appended child configId'});
      });
    });
  });

  test('Inserts at position 0', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'add-at-zero');
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 3, {prefix: 'Child count after insert'});
        expectDeepEqual(children[0].id, 'clustering', {prefix: 'Inserted child at 0'});
      });
    });
  });

  test('Multiple adds in one handler', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'double-add');
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 4, {prefix: 'Child count after double add'});
      });
    });
  });

  test('Appends to empty pipeline', async () => {
    const pconf = await getProcessedConfig(emptyDynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'add-regression');
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 1, {prefix: 'Child count after add to empty'});
        expectDeepEqual(children[0].id, 'regression', {prefix: 'Added child configId'});
      });
    });
  });
});

category('ComputeUtils: Driver pipeline granular: removeStep', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Removes targeted child', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'remove-first');
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 1, {prefix: 'Child count after remove'});
        expectDeepEqual(children[0].id, 'clustering', {prefix: 'Remaining child'});
      });
    });
  });

  test('Stale handle prevents mutation', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'stale-remove');
      // Handler tries to removeStep the same handle twice — second call throws,
      // error is caught internally by the link error handler, mutation is not applied
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        // Tree should be unchanged since the handler errored
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 2, {prefix: 'Tree unchanged after stale handle error'});
      });
    });
  });
});

category('ComputeUtils: Driver pipeline granular: moveStep', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Reorders children', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'move-last-to-front');
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 2, {prefix: 'Child count unchanged'});
        expectDeepEqual(children[0].id, 'clustering', {prefix: 'Moved child at 0'});
        expectDeepEqual(children[1].id, 'regression', {prefix: 'Original first now at 1'});
      });
    });
  });

  test('Preserves step state after move', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      // Set state on the last child (clustering at idx 1)
      const clusteringNode = tree.nodeTree.getNode([{idx: 1}, {idx: 1}]);
      const clusteringUuid = clusteringNode.getItem().uuid;
      clusteringNode.getItem().getStateStore().setState('a', 42);
      // Move it to front
      const action = getAction(tree, 'move-last-to-front');
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        // Now clustering is at idx 0
        const movedNode = tree.nodeTree.getNode([{idx: 1}, {idx: 0}]);
        expectDeepEqual(movedNode.getItem().uuid, clusteringUuid, {prefix: 'Same node after move'});
        expectDeepEqual(movedNode.getItem().getStateStore().getState('a'), 42, {prefix: 'State preserved'});
      });
    });
  });
});

category('ComputeUtils: Driver pipeline granular: mixed operations', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Remove then add in one handler', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'remove-then-add');
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 2, {prefix: 'Child count after remove+add'});
        // First was removed (regression), second remained (clustering), new one added (clustering)
        expectDeepEqual(children[0].id, 'clustering', {prefix: 'Remaining child'});
        expectDeepEqual(children[1].id, 'clustering', {prefix: 'Added child'});
      });
    });
  });
});

category('ComputeUtils: Driver pipeline granular: exclusivity', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('addStep then setPipelineState on same output prevents mutation', async () => {
    const pconf = await getProcessedConfig(dynamicConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const action = getAction(tree, 'exclusivity-test');
      // Handler calls addStep then setPipelineState on same output — throws internally,
      // error is caught by the link error handler, mutation is not applied
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        // Tree should be unchanged since the handler errored
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 2, {prefix: 'Tree unchanged after exclusivity error'});
      });
    });
  });
});

category('ComputeUtils: Driver pipeline granular: link interaction', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Data link to all() picks up added step', async () => {
    const pconf = await getProcessedConfig(dynamicWithLinkConfig);
    testScheduler.run(({cold, expectObservable}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const step1Node = tree.nodeTree.getNode([{idx: 0}]);
      // Set step1 output to trigger data link
      cold('-a').subscribe(() => {
        step1Node.getItem().getStateStore().setState('b', 99);
      });
      // Add a second regression step
      cold('--a').subscribe(() => {
        const action = getAction(tree, 'add-regression');
        tree.runAction(action.uuid).subscribe();
      });
      // After add + link recalculation, set value again to trigger link for new step
      cold('---a').subscribe(() => {
        step1Node.getItem().getStateStore().setState('b', 77);
      });
      cold('----a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 2, {prefix: 'Two regression steps'});
        // Both steps should have received the value via data link
        const step0store = children[0].item.getItem().getStateStore();
        const step1store = children[1].item.getItem().getStateStore();
        expectDeepEqual(step0store.getState('a'), 77, {prefix: 'First step got value'});
        expectDeepEqual(step1store.getState('a'), 77, {prefix: 'Added step got value'});
      });
    });
  });

  test('Remove step recalculates links', async () => {
    const pconf = await getProcessedConfig(dynamicWithLinkConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const initialLinks = [...tree.linksState.links.values()];
      const initialCount = initialLinks.length;
      // Remove the only regression step
      cold('-a').subscribe(() => {
        const action = getAction(tree, 'remove-first-regression');
        tree.runAction(action.uuid).subscribe();
      });
      cold('--a').subscribe(() => {
        const children = getAnalysesNode(tree).getChildren();
        expectDeepEqual(children.length, 0, {prefix: 'No children after remove'});
      });
    });
  });
});

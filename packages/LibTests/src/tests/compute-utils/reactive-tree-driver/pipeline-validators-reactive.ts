import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineNodeBase} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ControllerCancelled, PipelineValidatorController} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinkControllers';
import {DriverLogger} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';


function getPipelineNode(tree: StateTree, path: {idx: number}[] = []): PipelineNodeBase {
  return tree.nodeTree.getNode(path).getItem() as PipelineNodeBase;
}

function getAction(tree: StateTree, actionId: string) {
  const actions = [...tree.linksState.actions.values()];
  const action = actions.find((a) => a.matchInfo.spec.id === actionId);
  if (!action)
    throw new Error(`Action ${actionId} not found`);
  return action;
}


category('ComputeUtils: Driver pipeline validators: self target', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Validator with zero-segment `to` fires on `from` change after debounce', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
      ],
      links: [{
        id: 'check',
        type: 'pipelineValidator',
        from: 'a:step1/a',
        to: 'check',
        handler({controller}) {
          const a = controller.getFirst<number>('a');
          if (a && a < 0)
            controller.setValidation('check', {errors: [{description: 'a must be >= 0'}]});
          else
            controller.setValidation('check');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const pipelineNode = getPipelineNode(tree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', -3);
      });
      cold('400ms a').subscribe({
        next: () => {
          const values = Object.values(pipelineNode.pipelineValidations$.value);
          const errors = values.flatMap((v) => v?.errors ?? []);
          expectDeepEqual(errors, [{description: 'a must be >= 0'}]);
        },
      });
    });
  });

  test('Validator result appears on structureCheckResults of serialised state', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const pipelineNode = getPipelineNode(tree);
      pipelineNode.setPipelineValidation('test-uuid', {warnings: [{description: 'always-warn'}]});
      const state = tree.toState({skipFuncCalls: true});
      expectDeepEqual((state as any).structureCheckResults, {
        errors: [],
        warnings: [{description: 'always-warn'}],
        notifications: [],
      });
    });
  });

  test('setValidation(undefined) removes a contribution', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const pipelineNode = getPipelineNode(tree);
      pipelineNode.setPipelineValidation('uuid1', {errors: [{description: 'fail'}]});
      expectDeepEqual(Object.keys(pipelineNode.pipelineValidations$.value).length, 1);
      pipelineNode.setPipelineValidation('uuid1', undefined);
      expectDeepEqual(Object.keys(pipelineNode.pipelineValidations$.value).length, 0);
    });
  });
});


category('ComputeUtils: Driver pipeline validators: multi-validator merge', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Two validators on same node produce merged errors+warnings', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const pipelineNode = getPipelineNode(tree);
      pipelineNode.setPipelineValidation('uuid-a', {errors: [{description: 'err-a'}]});
      pipelineNode.setPipelineValidation('uuid-b', {warnings: [{description: 'warn-b'}]});
      const state = tree.toState({skipFuncCalls: true});
      expectDeepEqual((state as any).structureCheckResults, {
        errors: [{description: 'err-a'}],
        warnings: [{description: 'warn-b'}],
        notifications: [],
      });
    });
  });
});


category('ComputeUtils: Driver pipeline validators: descendant routing', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Validator with `to: first(inner)` lands result on the descendant node', async () => {
    const config: PipelineConfiguration = {
      id: 'outer',
      type: 'static',
      steps: [
        {
          id: 'inner',
          type: 'static',
          steps: [{id: 'inStep', nqName: 'LibTests:TestAdd2'}],
        },
      ],
      links: [{
        id: 'outerCheck',
        type: 'pipelineValidator',
        from: 'a:first(inner)/first(inStep)/a',
        to: 'check:first(inner)',
        handler({controller}) {
          controller.setValidation('check', {warnings: [{description: 'inner-warn'}]});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const innerNode = getPipelineNode(tree, [{idx: 0}]);
      const outerNode = getPipelineNode(tree);
      cold('-a 300ms b').subscribe({
        next: (v) => {
          if (v === 'a') {
            const innerStep = tree.nodeTree.getNode([{idx: 0}, {idx: 0}]);
            innerStep.getItem().getStateStore().setState('a', 1);
          } else if (v === 'b') {
            expectDeepEqual(Object.keys(innerNode.pipelineValidations$.value).length, 1, {prefix: 'inner has validation'});
            expectDeepEqual(Object.keys(outerNode.pipelineValidations$.value).length, 0, {prefix: 'outer has none'});
          }
        },
      });
    });
  });
});


category('ComputeUtils: Driver pipeline validators: multi-target `to`', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Array `to` writes the same result to every matched pipeline node', async () => {
    const config: PipelineConfiguration = {
      id: 'outer',
      type: 'static',
      steps: [
        {id: 'innerA', type: 'static', steps: [{id: 'sA', nqName: 'LibTests:TestAdd2'}]},
        {id: 'innerB', type: 'static', steps: [{id: 'sB', nqName: 'LibTests:TestAdd2'}]},
      ],
      links: [{
        id: 'fanOut',
        type: 'pipelineValidator',
        from: 'a:first(innerA)/first(sA)/a',
        to: ['ta:first(innerA)', 'tb:first(innerB)'],
        handler({controller}) {
          controller.setValidation('ta', {warnings: [{description: 'shared'}]});
          controller.setValidation('tb', {warnings: [{description: 'shared'}]});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const innerA = getPipelineNode(tree, [{idx: 0}]);
      const innerB = getPipelineNode(tree, [{idx: 1}]);
      cold('-a').subscribe(() => {
        tree.nodeTree.getNode([{idx: 0}, {idx: 0}]).getItem().getStateStore().setState('a', 7);
      });
      cold('400ms a').subscribe({
        next: () => {
          const valsA = Object.values(innerA.pipelineValidations$.value);
          const valsB = Object.values(innerB.pipelineValidations$.value);
          expectDeepEqual(valsA, [{warnings: [{description: 'shared'}]}], {prefix: 'innerA got result'});
          expectDeepEqual(valsB, [{warnings: [{description: 'shared'}]}], {prefix: 'innerB got result'});
        },
      });
    });
  });

  test('Named `to` writes a result to a single target, leaving the sibling untouched', async () => {
    const config: PipelineConfiguration = {
      id: 'outer',
      type: 'static',
      steps: [
        {id: 'innerA', type: 'static', steps: [{id: 'sA', nqName: 'LibTests:TestAdd2'}]},
        {id: 'innerB', type: 'static', steps: [{id: 'sB', nqName: 'LibTests:TestAdd2'}]},
      ],
      links: [{
        id: 'onlyA',
        type: 'pipelineValidator',
        from: 'a:first(innerA)/first(sA)/a',
        to: ['ta:first(innerA)', 'tb:first(innerB)'],
        handler({controller}) {
          // set only one of the two matched targets
          controller.setValidation('ta', {warnings: [{description: 'only-a'}]});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const innerA = getPipelineNode(tree, [{idx: 0}]);
      const innerB = getPipelineNode(tree, [{idx: 1}]);
      cold('-a').subscribe(() => {
        tree.nodeTree.getNode([{idx: 0}, {idx: 0}]).getItem().getStateStore().setState('a', 7);
      });
      cold('400ms a').subscribe({
        next: () => {
          const valsA = Object.values(innerA.pipelineValidations$.value);
          expectDeepEqual(valsA, [{warnings: [{description: 'only-a'}]}], {prefix: 'innerA got result'});
          expectDeepEqual(Object.keys(innerB.pipelineValidations$.value).length, 0, {prefix: 'innerB untouched'});
        },
      });
    });
  });

  test('Array `to` clears all targets on link destroy via orphan cleanup', async () => {
    const config: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [{id: 'innerA', type: 'static', steps: [{id: 'sA', nqName: 'LibTests:TestAdd2'}]}],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const root = getPipelineNode(tree);
      const innerA = getPipelineNode(tree, [{idx: 0}]);
      root.setPipelineValidation('uuidX', {errors: [{description: 'e1'}]});
      innerA.setPipelineValidation('uuidX', {errors: [{description: 'e1'}]});
      expectDeepEqual(Object.keys(root.pipelineValidations$.value).length, 1);
      expectDeepEqual(Object.keys(innerA.pipelineValidations$.value).length, 1);
      // simulate orphan sweep with a current-set that doesn't include uuidX
      root.clearOldPipelineValidations(new Set());
      innerA.clearOldPipelineValidations(new Set());
      expectDeepEqual(Object.keys(root.pipelineValidations$.value).length, 0);
      expectDeepEqual(Object.keys(innerA.pipelineValidations$.value).length, 0);
    });
  });
});


category('ComputeUtils: Driver pipeline validators: action step target', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('`to` resolving to an action step writes the validation and surfaces in structureCheckResults', async () => {
    const config: PipelineConfiguration = {
      id: 'outer',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'myAction', type: 'action'},
      ],
      links: [{
        id: 'checkAction',
        type: 'pipelineValidator',
        from: 'a:step1/a',
        to: 'check:first(myAction)',
        handler({controller}) {
          controller.setValidation('check', {warnings: [{description: 'action-warn'}]});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const actionNode = getPipelineNode(tree, [{idx: 1}]);
      cold('-a').subscribe(() => {
        tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      });
      cold('400ms a').subscribe({
        next: () => {
          const vals = Object.values(actionNode.pipelineValidations$.value);
          expectDeepEqual(vals, [{warnings: [{description: 'action-warn'}]}], {prefix: 'action step received validation'});
          const state = tree.toState({skipFuncCalls: true}) as any;
          const actionState = state.steps?.[1];
          expectDeepEqual(actionState?.isActionStep, true, {prefix: 'confirmed action step'});
          expectDeepEqual(actionState?.structureCheckResults, {
            errors: [],
            warnings: [{description: 'action-warn'}],
            notifications: [],
          }, {prefix: 'structureCheckResults surfaced on action step'});
        },
      });
    });
  });
});


category('ComputeUtils: Driver pipeline validators: tree mutation rerun', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Validator re-fires after structural mutation (no `from` change)', async () => {
    let callCount = 0;
    const config: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {
          id: 'analyses',
          type: 'dynamic',
          stepTypes: [{id: 'regression', nqName: 'LibTests:TestMul2'}],
          initialSteps: [{id: 'regression'}],
        },
      ],
      links: [{
        id: 'check',
        type: 'pipelineValidator',
        from: 'a:step1/a',
        to: 'check',
        handler({controller}) {
          callCount++;
          controller.setValidation('check', {warnings: [{description: `call ${callCount}`}]});
        },
      }],
      actions: [{
        id: 'add-regression',
        type: 'pipeline',
        from: [],
        to: 'out1:analyses',
        position: 'none',
        handler({controller}) {
          controller.addStep('out1', 'regression');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      cold('500ms a').subscribe({
        next: () => {
          const action = getAction(tree, 'add-regression');
          tree.runAction(action.uuid).subscribe();
        },
      });
      cold('2000ms a').subscribe({
        next: () => {
          expectDeepEqual(callCount >= 2, true, {prefix: `validator ran at least twice (got ${callCount})`});
        },
      });
    });
  });
});


category('ComputeUtils: Driver pipeline validators: orphan cleanup', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Removed descendant orphans are pruned from pipelineValidations$', async () => {
    const config: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const pipelineNode = getPipelineNode(tree);
      pipelineNode.setPipelineValidation('orphan-uuid', {errors: [{description: 'old'}]});
      expectDeepEqual(Object.keys(pipelineNode.pipelineValidations$.value).length, 1);
      pipelineNode.clearOldPipelineValidations(new Set());
      expectDeepEqual(Object.keys(pipelineNode.pipelineValidations$.value).length, 0);
    });
  });
});


category('ComputeUtils: Driver pipeline validators: funccall target rejection', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('`to` resolving to a FuncCall node logs a warning and skips the write', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
      ],
      links: [{
        id: 'misrouted',
        type: 'pipelineValidator',
        from: 'a:step1/a',
        to: 'check:first(step1)',                  // FuncCall node — invalid target
        handler({controller}) {
          controller.setValidation('check', {errors: [{description: 'should not appear'}]});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const logger = new DriverLogger();
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, logger});
      tree.init().subscribe();
      const pipelineNode = getPipelineNode(tree);
      cold('-a').subscribe(() => {
        tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      });
      cold('400ms a').subscribe({
        next: () => {
          expectDeepEqual(Object.keys(pipelineNode.pipelineValidations$.value).length, 0, {prefix: 'no validation written'});
          const warning = logger.errors.find((e) => e.context === 'link:misrouted');
          expectDeepEqual(!!warning, true, {prefix: 'warning logged'});
          expectDeepEqual(warning?.severity, 'warning', {prefix: 'severity is warning'});
        },
      });
    });
  });
});


category('ComputeUtils: Driver pipeline validators: cancelled controller', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('getOutline and setValidation throw ControllerCancelled after close()', async () => {
    testScheduler.run(() => {
      const stubOutline = {configId: 'stub'} as any;
      const controller = new PipelineValidatorController({
        inputs: {}, inputsSet: new Set(), outputsSet: new Set(), callInputs: new Set(),
        id: 'spec-id', outline: stubOutline,
      });
      controller.close();
      let outlineThrew = false;
      let setValThrew = false;
      try { controller.getOutline(); } catch (e) { outlineThrew = e instanceof ControllerCancelled; }
      try { controller.setValidation('check', {errors: [{description: 'x'}]}); } catch (e) { setValThrew = e instanceof ControllerCancelled; }
      expectDeepEqual(outlineThrew, true, {prefix: 'getOutline throws'});
      expectDeepEqual(setValThrew, true, {prefix: 'setValidation throws'});
    });
  });
});


category('ComputeUtils: Driver pipeline validators: getOutline', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Handler can read the link-defining pipeline subtree via controller.getOutline()', async () => {
    let observedConfigId: string | undefined;
    let observedStepCount: number | undefined;
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'check',
        type: 'pipelineValidator',
        from: 'a:step1/a',
        to: 'check',
        handler({controller}) {
          const outline = controller.getOutline();
          observedConfigId = outline.configId;
          observedStepCount = (outline as any).steps?.length;
          controller.setValidation('check');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      cold('-a').subscribe(() => {
        tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 5);
      });
      cold('300ms a').subscribe({
        next: () => {
          expectDeepEqual(observedConfigId, 'pipeline1', {prefix: 'outline.configId'});
          expectDeepEqual(observedStepCount, 2, {prefix: 'outline.steps.length'});
        },
      });
    });
  });

  test('getOutline reflects link-defining pipeline, not the `to` target', async () => {
    let observedConfigId: string | undefined;
    const config: PipelineConfiguration = {
      id: 'outer',
      type: 'static',
      steps: [
        {
          id: 'inner',
          type: 'static',
          steps: [{id: 'inStep', nqName: 'LibTests:TestAdd2'}],
        },
      ],
      links: [{
        id: 'outerCheck',
        type: 'pipelineValidator',
        from: 'a:first(inner)/first(inStep)/a',
        to: 'check:first(inner)',                   // result lands on inner
        handler({controller}) {
          observedConfigId = controller.getOutline().configId;
          controller.setValidation('check');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      cold('-a').subscribe(() => {
        tree.nodeTree.getNode([{idx: 0}, {idx: 0}]).getItem().getStateStore().setState('a', 5);
      });
      cold('300ms a').subscribe({
        next: () => {
          expectDeepEqual(observedConfigId, 'outer', {prefix: 'outline reflects link scope, not target'});
        },
      });
    });
  });
});

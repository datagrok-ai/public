import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {of, Subject} from 'rxjs';
import {delay, filter, mapTo, take} from 'rxjs/operators';
import {snapshotCompare, createTestScheduler} from '../../../test-utils';


category('ComputeUtils: Driver links reactivity: actions', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Get and run pipeline validation actions', async () => {
    const s = new Subject<string>();
    const config3: PipelineConfiguration = {
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
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        actions: 'actions',
        type: 'validator',
        handler({controller}) {
          const v = controller.getFirst('in1');
          if (v === 1) {
            const action = controller.getValidationAction('actions', 'action1');
            s.next(action);
          }
          return;
        },
      }],
      actions: [{
        id: 'action1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        position: 'none',
        handler({controller}) {
          controller.setAll('out1', 10);
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config3);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      s.pipe(
        delay(1),
      ).subscribe((uuid) => {
        tree.runAction(uuid);
      });
      expectObservable(inNode.getItem().getStateStore().getStateChanges('a')).toBe('ab 250ms c', {a: undefined, b: 1, c: 10});
    });
  });

  test('Get and run pipeline validation actions with params', async () => {
    const s = new Subject<string>();
    const config3: PipelineConfiguration = {
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
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        actions: 'actions',
        type: 'validator',
        handler({controller}) {
          const v = controller.getFirst('in1');
          if (v === 1) {
            const action = controller.getValidationAction('actions', 'action1');
            s.next(action);
          }
          return;
        },
      }],
      actions: [{
        id: 'action1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        position: 'none',
        handler({controller}) {
          const val = controller.getAdditionalParam('val');
          controller.setAll('out1', val);
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config3);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      s.pipe(
        delay(1),
      ).subscribe((uuid) => {
        tree.runAction(uuid, {val: 100});
      });
      expectObservable(inNode.getItem().getStateStore().getStateChanges('a')).toBe('ab 250ms c', {a: undefined, b: 1, c: 100});
    });
  });

  test('Get and run funcall validation actions', async () => {
    const s = new Subject<string>();
    const config3: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          actions: [{
            id: 'action1',
            from: 'in1:a',
            to: 'out1:a',
            position: 'none',
            handler({controller}) {
              controller.setAll('out1', 10);
              return;
            },
          }],
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        actions: 'actions:step1',
        type: 'validator',
        handler({controller}) {
          const v = controller.getFirst('in1');
          if (v === 1) {
            const action = controller.getValidationAction('actions', 'action1');
            s.next(action);
          }
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config3);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      s.pipe(
        delay(1),
      ).subscribe((uuid) => {
        tree.runAction(uuid);
      });
      expectObservable(inNode.getItem().getStateStore().getStateChanges('a')).toBe('ab 250ms c', {a: undefined, b: 1, c: 10});
    });
  });

  test('Run pipeline mutation actions', async () => {
    const config4: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'nestedPipeline',
          type: 'parallel',
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
            {
              id: 'step1',
            },
          ],
        },
        {
          id: 'stepr',
          nqName: 'LibTests:TestSub2',
        },
      ],
      actions: [{
        id: 'action1',
        from: [],
        position: 'none',
        to: 'out1:nestedPipeline',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1', {
            id: 'nestedPipeline',
            steps: [
              {
                id: 'step2',
                initialValues: {
                  a: 5,
                },
              },
              {
                id: 'step1',
                initialValues: {
                  a: 10,
                },
              },
            ],
          });
        },
      }],
      links: [{
        id: 'link1',
        from: 'in1:nestedPipeline/all(step1|step2)/a',
        to: 'out1:stepr/a',
        handler({controller}) {
          const v = controller.getAll<number>('in1')!;
          const r = v.reduce((acc, val) => acc + val, 0);
          controller.setAll('out1', r);
        },
      }],
    };

    const pconf = await getProcessedConfig(config4);
    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const out = tree.nodeTree.getNode([{idx: 1}]);
      const action = [...tree.linksState.actions.values()][0];
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      expectObservable(out.getItem().getStateStore().getStateChanges('a')).toBe('ab', {a: undefined, b: 15});
    });
  });

  test('Run pipeline mutation actions on root pipeline', async () => {
    const config5: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'parallel',
      stepTypes: [
        {
          id: 'nestedPipeline1',
          type: 'parallel',
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
            {
              id: 'step1',
            },
          ],
        },
        {
          id: 'step',
          nqName: 'LibTests:TestSub2',
        },
      ],
      actions: [{
        id: 'action1',
        from: [],
        position: 'none',
        to: 'out1',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1', {
            id: 'pipeline1',
            steps: [{
              id: 'nestedPipeline1',
              steps: [
                {
                  id: 'step2',
                  initialValues: {
                    a: 5,
                  },
                },
                {
                  id: 'step1',
                  initialValues: {
                    a: 10,
                  },
                },
              ],
            }],
          });
        },
      }],
    };

    const pconf = await getProcessedConfig(config5);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    await snapshotCompare(state, 'Run pipeline mutation actions on root pipeline');
  });

  test('Run pipeline mutation actions on pipeline', async () => {
    const config5: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'nestedPipeline1',
          type: 'parallel',
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
            {
              id: 'step1',
            },
          ],
        },
        {
          id: 'step',
          nqName: 'LibTests:TestSub2',
        },
      ],
      actions: [{
        id: 'action1',
        from: [],
        position: 'none',
        to: 'out1:nestedPipeline1',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1',
            {
              id: 'nestedPipeline1',
              steps: [
                {
                  id: 'step2',
                  initialValues: {
                    a: 5,
                  },
                },
                {
                  id: 'step1',
                  initialValues: {
                    a: 10,
                  },
                },
              ],
            },
          );
        },
      }],
    };

    const pconf = await getProcessedConfig(config5);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    await snapshotCompare(state, 'Run pipeline mutation actions on pipeline');
  });

  test('Select state descriptions', async () => {
    const config5: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      links: [{
        id: 'selector',
        type: 'selector',
        from: 'in:step1/a',
        to: ['out1:title', 'out2:description', 'out3:tags'],
        handler({controller}) {
          const val = controller.getFirst('in');
          controller.setDescriptionItem('out1', `Title ${val}`);
          controller.setDescriptionItem('out2', `Description ${val}`);
          controller.setDescriptionItem('out3', [`tag ${val}`]);
        },
      }],
    };
    const pconf = await getProcessedConfig(config5);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]);
      const link = [...tree.linksState.links.values()][0];
      const pipeline = tree.nodeTree.root;
      cold('-a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 1);
      });
      cold('--a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 2);
      });
      expectObservable(pipeline.getItem().nodeDescription.getStateChanges('title')).toBe('abc',
        {a: undefined, b: 'Title 1', c: 'Title 2'});
      expectObservable(pipeline.getItem().nodeDescription.getStateChanges('description')).toBe('abc',
        {a: undefined, b: 'Description 1', c: 'Description 2'});
      expectObservable(pipeline.getItem().nodeDescription.getStateChanges('tags')).toBe('abc',
        {a: undefined, b: {[link.uuid]: ['tag 1']}, c: {[link.uuid]: ['tag 2']}});
    });
  });

  test('FuncCall pipeline actions', async () => {
    const config1: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          initialValues: {
            a: 1,
            b: 2,
          },
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
          initialValues: {
            a: 3,
            b: 4,
          },
        },
      ],
      actions: [{
        id: 'action1',
        type: 'funccall',
        from: ['a:step1/a', 'b:step1/b', 'fc(call):step2'],
        to: 'out(call):step2',
        position: 'none',
        async handler({controller}) {
          const a = controller.getFirst('a');
          const b = controller.getFirst('b');
          const fc = controller.getFirst('fc');
          expectDeepEqual(a, 1);
          expectDeepEqual(b, 2);
          expectDeepEqual(fc instanceof DG.FuncCall, true);
          const func: DG.Func = await grok.functions.eval('LibTests:simpleInputs');
          const nfc = func.prepare({a: 33, b: 44});
          controller.setFuncCall('out', nfc);
        },
      }],
    };

    const pconf = await getProcessedConfig(config1);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    await tree.treeMutationsLocked$.pipe(filter((x) => !x), take(1)).toPromise();
    const node = tree.nodeTree.getNode([{idx: 1}]);
    const a = node.getItem().getStateStore().getState('a');
    const b = node.getItem().getStateStore().getState('b');
    expectDeepEqual(a, 33);
    expectDeepEqual(b, 44);
  });

  test('FuncCall step actions', async () => {
    const config1: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          initialValues: {
            a: 1,
            b: 2,
          },
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
          initialValues: {
            a: 3,
            b: 4,
          },
          actions: [{
            id: 'action1',
            type: 'funccall',
            from: 'fc(call)',
            to: 'out(call)',
            position: 'none',
            async handler({controller}) {
              const fc = controller.getFirst('fc');
              expectDeepEqual(fc instanceof DG.FuncCall, true);
              const func: DG.Func = await grok.functions.eval('LibTests:simpleInputs');
              const nfc = func.prepare({a: 33, b: 44});
              controller.setFuncCall('out', nfc);
            },
          }],
        },
      ],
    };

    const pconf = await getProcessedConfig(config1);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    await tree.treeMutationsLocked$.pipe(filter((x) => !x), take(1)).toPromise();
    const node = tree.nodeTree.getNode([{idx: 1}]);
    const a = node.getItem().getStateStore().getState('a');
    const b = node.getItem().getStateStore().getState('b');
    expectDeepEqual(a, 33);
    expectDeepEqual(b, 44);
  });

  test('Links priority ordering per node', async () => {
    const config1: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestMultiarg5',
        },
      ],
      links: [
        {
          id: 'd-e',
          from: 'in:step1/d',
          to: 'out:step1/e',
          nodePriority: 1,
          handler({controller}) {
            const in1 = controller.getFirst('in');
            controller.setAll('out', in1);
            return of(void(0)).pipe(delay(100));
          },
        },
        {
          id: 'a-b',
          from: 'in:step1/a',
          to: 'out:step1/b',
          nodePriority: 4,
          handler({controller}) {
            const in1 = controller.getFirst('in');
            controller.setAll('out', in1);
            return of(void(0)).pipe(delay(100));
          },
        },
        {
          id: 'c-d',
          from: 'in:step1/c',
          to: 'out:step1/d',
          nodePriority: 2,
          handler({controller}) {
            const in1 = controller.getFirst('in');
            controller.setAll('out', in1);
            return of(void(0)).pipe(delay(100));
          },
        },
        {
          id: 'b-c',
          from: 'in:step1/b',
          to: 'out:step1/c',
          nodePriority: 3,
          handler({controller}) {
            const in1 = controller.getFirst('in');
            controller.setAll('out', in1);
            return of(void(0)).pipe(delay(100));
          },
        },
      ],
    };
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 1);
      });
      expectObservable(node.getItem().getStateStore().getStateChanges('a')).toBe('ab', {a: undefined, b: 1});
      expectObservable(node.getItem().getStateStore().getStateChanges('b')).toBe('a 100ms b', {a: undefined, b: 1});
      expectObservable(node.getItem().getStateStore().getStateChanges('c')).toBe('a 200ms b', {a: undefined, b: 1});
      expectObservable(node.getItem().getStateStore().getStateChanges('d')).toBe('a 300ms b', {a: undefined, b: 1});
      expectObservable(node.getItem().getStateStore().getStateChanges('e')).toBe('a 400ms b', {a: undefined, b: 1});
    });
  });

  test('Links optional match test', async () => {
    const inputs$ = new Subject<Set<string>>();
    const outputs$ = new Subject<Set<string>>();
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
            {
              id: 'stepAdd',
            },
          ],
        },
        {
          id: 'step3',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: ['in1:pipelinePar/stepAdd/a', 'in2(optional):pipelinePar/stepMul/a'],
        to: ['out1:step3/a', 'out2(optional):step2/a'],
        handler({controller}) {
          inputs$.next(controller.getMatchedInputs());
          outputs$.next(controller.getMatchedOutputs());
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}, {idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 10);
      });
      expectObservable(inputs$, '^ 1000ms !').toBe('-a', {a: new Set(['in1'])});
      expectObservable(outputs$, '^ 1000ms !').toBe('-a', {a: new Set(['out1'])});
    });
  });

  test('Pipeline action with visibleOn routes to child step', async () => {
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
      actions: [{
        id: 'action1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
        position: 'buttons',
        visibleOn: 'step2',
        handler({controller}) {
          controller.setAll('out1', controller.getFirst('in1'));
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();

      // The action should be visible on step2, not on the pipeline root
      const step2Node = tree.nodeTree.getNode([{idx: 1}]);
      const step2Actions = tree.linksState.getNodeActionsData(step2Node.getItem().uuid);
      expectDeepEqual(step2Actions?.length, 1);
      expectDeepEqual(step2Actions?.[0].id, 'action1');
      expectDeepEqual(step2Actions?.[0].position, 'buttons');

      // Pipeline root should have no actions
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(rootActions, undefined);

      // The action should still work — it resolves from/to at pipeline level
      const action = [...tree.linksState.actions.values()][0];
      const step1Node = tree.nodeTree.getNode([{idx: 0}]);

      // Set step1 output, then run the action to copy it to step2
      cold('-a').subscribe(() => {
        step1Node.getItem().getStateStore().editState('res', 42);
      });
      cold('--a').subscribe(() => {
        tree.runAction(action.uuid);
      });
      expectObservable(step2Node.getItem().getStateStore().getStateChanges('a')).toBe('a-b', {a: undefined, b: 42});
    });
  });

  test('Pipeline action with visibleOn routes to deeply nested step', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'nestedPipeline',
          type: 'static',
          steps: [
            {
              id: 'innerStep1',
              nqName: 'LibTests:TestAdd2',
            },
            {
              id: 'innerStep2',
              nqName: 'LibTests:TestMul2',
            },
          ],
        },
      ],
      actions: [{
        id: 'action1',
        from: 'in1:nestedPipeline/innerStep1/res',
        to: 'out1:nestedPipeline/innerStep2/a',
        position: 'buttons',
        visibleOn: 'innerStep2',
        handler({controller}) {
          controller.setAll('out1', controller.getFirst('in1'));
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();

      // innerStep2 is at path [nestedPipeline, innerStep2]
      const innerStep2Node = tree.nodeTree.getNode([{idx: 0}, {idx: 1}]);
      const innerStep2Actions = tree.linksState.getNodeActionsData(innerStep2Node.getItem().uuid);
      expectDeepEqual(innerStep2Actions?.length, 1);
      expectDeepEqual(innerStep2Actions?.[0].id, 'action1');

      // Root and nestedPipeline should have no actions
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(rootActions, undefined);
      const nestedActions = tree.linksState.getNodeActionsData(tree.nodeTree.getNode([{idx: 0}]).getItem().uuid);
      expectDeepEqual(nestedActions, undefined);

      // Action should still work across the nested boundary
      const innerStep1Node = tree.nodeTree.getNode([{idx: 0}, {idx: 0}]);
      cold('-a').subscribe(() => {
        innerStep1Node.getItem().getStateStore().editState('res', 77);
      });
      cold('--a').subscribe(() => {
        const action = [...tree.linksState.actions.values()][0];
        tree.runAction(action.uuid);
      });
      expectObservable(innerStep2Node.getItem().getStateStore().getStateChanges('a')).toBe('a-b', {a: undefined, b: 77});
    });
  });

  test('FuncCall action with visibleOn routes to child step', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          initialValues: {
            a: 1,
            b: 2,
          },
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
          initialValues: {
            a: 3,
            b: 4,
          },
        },
      ],
      actions: [{
        id: 'action1',
        type: 'funccall',
        from: ['a:step1/a', 'b:step1/b', 'fc(call):step2'],
        to: 'out(call):step2',
        position: 'buttons',
        visibleOn: 'step1',
        async handler({controller}) {
          const a = controller.getFirst('a');
          const b = controller.getFirst('b');
          const fc = controller.getFirst('fc');
          expectDeepEqual(fc instanceof DG.FuncCall, true);
          const func: DG.Func = await grok.functions.eval('LibTests:simpleInputs');
          const nfc = func.prepare({a: a + 10, b: b + 10});
          controller.setFuncCall('out', nfc);
        },
      }],
    };

    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();

    // Action should be visible on step1, not on pipeline root
    const step1Node = tree.nodeTree.getNode([{idx: 0}]);
    const step1Actions = tree.linksState.getNodeActionsData(step1Node.getItem().uuid);
    expectDeepEqual(step1Actions?.length, 1);
    expectDeepEqual(step1Actions?.[0].id, 'action1');
    expectDeepEqual(step1Actions?.[0].position, 'buttons');

    const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
    expectDeepEqual(rootActions, undefined);

    // Action should still work — replaces step2's FuncCall
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    await tree.treeMutationsLocked$.pipe(filter((x) => !x), take(1)).toPromise();
    const step2Node = tree.nodeTree.getNode([{idx: 1}]);
    const a = step2Node.getItem().getStateStore().getState('a');
    const b = step2Node.getItem().getStateStore().getState('b');
    expectDeepEqual(a, 11);
    expectDeepEqual(b, 12);
  });

  test('Pipeline mutation action with visibleOn routes to child step', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'nestedPipeline1',
          type: 'parallel',
          stepTypes: [
            {
              id: 'step1',
              nqName: 'LibTests:TestAdd2',
            },
          ],
          initialSteps: [
            {
              id: 'step1',
            },
          ],
        },
        {
          id: 'step',
          nqName: 'LibTests:TestSub2',
        },
      ],
      actions: [{
        id: 'action1',
        from: [],
        position: 'buttons',
        to: 'out1:nestedPipeline1',
        type: 'pipeline',
        visibleOn: 'step',
        handler({controller}) {
          controller.setPipelineState('out1',
            {
              id: 'nestedPipeline1',
              steps: [
                {
                  id: 'step1',
                  initialValues: {a: 5},
                },
                {
                  id: 'step1',
                  initialValues: {a: 10},
                },
              ],
            },
          );
        },
      }],
    };

    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();

    // Action should be visible on 'step', not on pipeline root
    const stepNode = tree.nodeTree.getNode([{idx: 1}]);
    const stepActions = tree.linksState.getNodeActionsData(stepNode.getItem().uuid);
    expectDeepEqual(stepActions?.length, 1);
    expectDeepEqual(stepActions?.[0].id, 'action1');
    expectDeepEqual(stepActions?.[0].position, 'buttons');

    const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
    expectDeepEqual(rootActions, undefined);

    // Action should still work — mutates nestedPipeline1
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    const nestedNode = tree.nodeTree.getNode([{idx: 0}]);
    const children = nestedNode.getChildren();
    expectDeepEqual(children.length, 2);
    const child1Val = tree.nodeTree.getNode([{idx: 0}, {idx: 0}]).getItem().getStateStore().getState('a');
    const child2Val = tree.nodeTree.getNode([{idx: 0}, {idx: 1}]).getItem().getStateStore().getState('a');
    expectDeepEqual(child1Val, 5);
    expectDeepEqual(child2Val, 10);
  });

  test('Pipeline action without visibleOn stays on pipeline node', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      actions: [{
        id: 'action1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        position: 'buttons',
        handler({controller}) {
          controller.setAll('out1', 99);
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();

      // Without visibleOn, action stays on pipeline root
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(rootActions?.length, 1);
      expectDeepEqual(rootActions?.[0].id, 'action1');

      // step1 should have no actions
      const step1Node = tree.nodeTree.getNode([{idx: 0}]);
      const step1Actions = tree.linksState.getNodeActionsData(step1Node.getItem().uuid);
      expectDeepEqual(step1Actions, undefined);
    });
  });
});

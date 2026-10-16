import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {of, Subject} from 'rxjs';
import {delay, filter, map, mapTo, take} from 'rxjs/operators';
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

  test('setPipelineState accepts string-shorthand steps', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'nestedPipeline',
          type: 'parallel',
          stepTypes: [
            {id: 'stepAdd', nqName: 'LibTests:TestAdd2'},
            {id: 'stepMul', nqName: 'LibTests:TestMul2'},
          ],
          initialSteps: ['stepAdd'],
        },
        {id: 'stepr', nqName: 'LibTests:TestSub2'},
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
            steps: ['stepMul', {id: 'stepAdd', initialValues: {a: 4, b: 6}}],
          });
        },
      }],
      links: [{
        id: 'link1',
        from: 'in1:nestedPipeline/all(stepAdd|stepMul)/a',
        to: 'out1:stepr/a',
        handler({controller}) {
          const v = controller.getAll<number>('in1')!;
          const r = v.reduce((acc, val) => acc + (val ?? 0), 0);
          controller.setAll('out1', r);
        },
      }],
    };

    const pconf = await getProcessedConfig(config);
    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const out = tree.nodeTree.getNode([{idx: 1}]);
      const action = [...tree.linksState.actions.values()][0];
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      // After mutation: nestedPipeline has [stepMul (no a), stepAdd (a=4)] — sum of 'a' = 4
      expectObservable(out.getItem().getStateStore().getStateChanges('a')).toBe('ab', {a: undefined, b: 4});
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

  test('hasCall (call,optional) branches handler output: guard present', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestDiv2'},
      ],
      links: [{
        id: 'link1',
        from: ['in1:step1/b', 'guard(call,optional):step3'],
        to: 'out1:step2/a',
        handler({controller}) {
          const in1 = controller.getFirst<number>('in1')!;
          const v = controller.hasCall('guard') ? in1 * 100 : in1;
          controller.setAll('out1', v);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 3);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 300});
    });
  });

  test('hasCall (call,optional) branches handler output: guard absent', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: ['in1:step1/b', 'guard(call,optional):missingStep'],
        to: 'out1:step2/a',
        handler({controller}) {
          const in1 = controller.getFirst<number>('in1')!;
          const v = controller.hasCall('guard') ? in1 * 100 : in1;
          controller.setAll('out1', v);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 3);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 3});
    });
  });

  test('hasCall throws on unknown call input name (typo guard)', async () => {
    const errors$ = new Subject<unknown>();
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestDiv2'},
      ],
      links: [{
        id: 'link1',
        from: ['in1:step1/b', 'guard(call,optional):step3'],
        to: 'out1:step2/a',
        handler({controller}) {
          try {
            controller.hasCall('typo');
          } catch (e) {
            errors$.next(e);
          }
          controller.setAll('out1', controller.getFirst('in1'));
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
      });
      expectObservable(errors$.pipe(map((e: any) => e instanceof Error && /unknown call input typo/.test(e.message))), '^ 500ms !').toBe('-a', {a: true});
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

  test('Pipeline action with visibleOn pointing to unknown id is hidden', async () => {
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
        visibleOn: 'doesNotExist',
        handler({controller}) {
          controller.setAll('out1', 99);
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();

      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(rootActions, undefined);

      const step1Node = tree.nodeTree.getNode([{idx: 0}]);
      const step1Actions = tree.linksState.getNodeActionsData(step1Node.getItem().uuid);
      expectDeepEqual(step1Actions, undefined);
    });
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

  test('Action with showWhen matching is visible', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
      ],
      actions: [{
        id: 'action1',
        position: 'none',
        from: [],
        to: 'out1:step1/a',
        showWhen: 'guard: step1',
        handler({controller}) {
          controller.setAll('out1', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(rootActions?.length, 1);
      expectDeepEqual(rootActions?.[0].id, 'action1');
      expectDeepEqual(tree.linksState.isActionVisible(rootActions![0].uuid), true);
    });
  });

  test('Action with unmatched showWhen is hidden but still present', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
      ],
      actions: [{
        id: 'action1',
        position: 'none',
        from: [],
        to: 'out1:step1/a',
        showWhen: 'guard: missingStep',
        handler({controller}) {
          controller.setAll('out1', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      // Action is still in nodesActions (RTD does not filter), but hidden
      expectDeepEqual(rootActions?.length, 1);
      expectDeepEqual(rootActions?.[0].id, 'action1');
      expectDeepEqual(tree.linksState.isActionVisible(rootActions![0].uuid), false);
      // Action is still runnable (no runtime gating)
      const actionUuid = rootActions![0].uuid;
      expectDeepEqual(tree.linksState.actions.has(actionUuid), true);
    });
  });

  const driverMode = {batchLinks: true, annotationChecks: true};

  const hiddenByHideWhen = (mode: {batchLinks?: boolean, annotationChecks?: boolean}) => async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
      ],
      actions: [{
        id: 'action1',
        position: 'none',
        from: [],
        to: 'out1:step1/a',
        hideWhen: 'guard: step1',
        handler({controller}) {
          controller.setAll('out1', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, ...mode});
      tree.init().subscribe();
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(rootActions?.length, 1);
      expectDeepEqual(tree.linksState.isActionVisible(rootActions![0].uuid), false);
    });
  };

  test('Action with matching hideWhen is hidden', hiddenByHideWhen({}));

  test('Action with matching hideWhen is hidden in Driver mode', hiddenByHideWhen(driverMode));

  test('hideWhen wins over showWhen when both match', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
      ],
      actions: [{
        id: 'action1',
        position: 'none',
        from: [],
        to: 'out1:step1/a',
        showWhen: 'guard: step1',
        hideWhen: 'guard: step1',
        handler({controller}) {
          controller.setAll('out1', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(tree.linksState.isActionVisible(rootActions![0].uuid), false);
    });
  });

  test('showWhen with (optional) flag does not gate', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
      ],
      actions: [{
        id: 'action1',
        position: 'none',
        from: [],
        to: 'out1:step1/a',
        showWhen: 'guard(optional): missingStep',
        handler({controller}) {
          controller.setAll('out1', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(tree.linksState.isActionVisible(rootActions![0].uuid), true);
    });
  });

  test('Action with no visibility condition is visible', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
      ],
      actions: [{
        id: 'action1',
        position: 'none',
        from: [],
        to: 'out1:step1/a',
        handler({controller}) {
          controller.setAll('out1', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const rootActions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid);
      expectDeepEqual(tree.linksState.isActionVisible(rootActions![0].uuid), true);
    });
  });

  test('Visibility re-evaluates after tree mutation', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'nestedPipeline',
          type: 'parallel',
          stepTypes: [
            {id: 'step1', nqName: 'LibTests:TestAdd2'},
            {id: 'step2', nqName: 'LibTests:TestMul2'},
          ],
          initialSteps: [{id: 'step1'}],
        },
      ],
      actions: [{
        id: 'addStep2',
        position: 'none',
        from: [],
        to: 'out1:nestedPipeline',
        type: 'pipeline',
        showWhen: 'guard: nestedPipeline/first(step1)',
        hideWhen: 'guard: nestedPipeline/first(step2)',
        handler({controller}) {
          controller.setPipelineState('out1', {
            id: 'nestedPipeline',
            steps: [{id: 'step1'}, {id: 'step2'}],
          });
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const rootUuid = tree.nodeTree.root.getItem().uuid;
    const before = tree.linksState.getNodeActionsData(rootUuid);
    expectDeepEqual(before?.length, 1);
    expectDeepEqual(tree.linksState.isActionVisible(before![0].uuid), true);

    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();

    const after = tree.linksState.getNodeActionsData(rootUuid);
    expectDeepEqual(after?.length, 1);
    // step2 now exists → hideWhen matches → hidden
    expectDeepEqual(tree.linksState.isActionVisible(after![0].uuid), false);
  });

  test('Validator gets a hidden action through getValidationAction', async () => {
    const observed = new Subject<string | undefined>();
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
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
            observed.next(controller.getValidationAction('actions', 'action1'));
          }
        },
      }],
      actions: [{
        id: 'action1',
        position: 'none',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        showWhen: 'guard: missingStep',
        handler({controller}) {
          controller.setAll('out1', 99);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      expectObservable(observed.pipe(map((uuid) => !!uuid))).toBe('- 250ms a', {a: true});
    });
  });

  test('Validator getValidationAction returns undefined for unknown name/actionId', async () => {
    const observed = new Subject<{known: boolean, unknownName: boolean, unknownId: boolean}>();
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
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
            observed.next({
              known: !!controller.getValidationAction('actions', 'action1'),
              unknownName: !!controller.getValidationAction('nope', 'action1'),
              unknownId: !!controller.getValidationAction('actions', 'nope'),
            });
          }
        },
      }],
      actions: [{
        id: 'action1',
        position: 'none',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        handler({controller}) {
          controller.setAll('out1', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      expectObservable(observed).toBe('- 250ms a', {a: {known: true, unknownName: false, unknownId: false}});
    });
  });

  const daysConfig = (
    initialDays: number, links: any[], daysActions: any[], rootActions: any[] = [],
  ): PipelineConfiguration => ({
    id: 'pipeline1',
    type: 'static',
    steps: [
      {id: 'cooling', nqName: 'LibTests:TestAdd2'},
      {
        id: 'days',
        type: 'parallel',
        stepTypes: [{id: 'day', nqName: 'LibTests:TestAdd2'}],
        initialSteps: Array.from({length: initialDays}, () => ({id: 'day'})),
        actions: daysActions,
      },
    ],
    links,
    actions: rootActions,
  });

  const recordingValidator = (actions: string, last: {uuid?: string, runs: number}) => ({
    id: 'check',
    type: 'validator' as const,
    debounce: 0,
    from: 'in:cooling/a',
    to: 'out:cooling/a',
    actions,
    handler({controller}: any) {
      last.uuid = controller.getValidationAction('act', 'action1');
      last.runs++;
    },
  });

  const daysAction = {
    id: 'action1',
    position: 'none' as const,
    from: 'in:all(day)/a',
    to: 'out:all(day)/b',
    handler({controller}: any) {
      controller.setAll('out', 1);
    },
  };

  test('Validation action stays runnable after its action is recreated', async () => {
    const last = {runs: 0} as {uuid?: string, runs: number};
    const config = daysConfig(1, [recordingValidator('act:days', last)], [daysAction]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const cooling = tree.nodeTree.getNode([{idx: 0}]);
      cooling.getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.actions.has(last.uuid!), true, {prefix: 'before mutation'});

      const daysUuid = tree.nodeTree.getNode([{idx: 1}]).getItem().uuid;
      tree.addSubTree(daysUuid, 'day', 1).subscribe();
      flush();
      expectDeepEqual(tree.linksState.actions.has(last.uuid!), true, {prefix: 'after mutation'});

      cooling.getItem().getStateStore().setState('a', 2);
      flush();
      expectDeepEqual(tree.linksState.actions.has(last.uuid!), true, {prefix: 'after input change'});
    });
  });

  test('Validation action is found once its action first matches', async () => {
    const last = {runs: 0} as {uuid?: string, runs: number};
    const config = daysConfig(0, [recordingValidator('act:days', last)], [daysAction]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const cooling = tree.nodeTree.getNode([{idx: 0}]);
      cooling.getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(last.uuid === undefined, true, {prefix: 'no action before mutation'});

      const daysUuid = tree.nodeTree.getNode([{idx: 1}]).getItem().uuid;
      tree.addSubTree(daysUuid, 'day', 0).subscribe();
      flush();
      expectDeepEqual(tree.linksState.actions.has(last.uuid!), true, {prefix: 'after mutation'});

      cooling.getItem().getStateStore().setState('a', 2);
      flush();
      expectDeepEqual(tree.linksState.actions.has(last.uuid!), true, {prefix: 'after input change'});
    });
  });

  const rootAction = (handler?: () => void, extra: Record<string, any> = {}) => ({
    id: 'action1',
    position: 'none' as const,
    from: [],
    to: 'out:cooling/b',
    handler: handler ?? (() => {}),
    ...extra,
  });

  const hidingLink = (
    id: string, from: string, actions: string, visible: (v: any) => boolean | undefined, runs?: {count: number},
  ) => ({
    id,
    type: 'nodemeta' as const,
    from,
    to: [],
    actions,
    handler({controller}: any) {
      if (runs)
        runs.count++;
      const value = visible(controller.getFirst('in'));
      if (value !== undefined)
        controller.setActionVisible('act', 'action1', value);
    },
  });

  const actionUuid = (tree: StateTree, nodePath: {idx: number}[] = []) =>
    tree.linksState.getNodeActionsData(tree.nodeTree.getNode(nodePath).getItem().uuid)![0].uuid;

  const hidesFromInputs = (mode: {batchLinks?: boolean, annotationChecks?: boolean}) => async () => {
    const config = daysConfig(0, [hidingLink('vis', 'in:cooling/a', 'act', (v) => v !== 1)], [], [rootAction()]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, ...mode});
      tree.init().subscribe();
      flush();
      const cooling = tree.nodeTree.getNode([{idx: 0}]);
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), true, {prefix: 'initial'});

      cooling.getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'hidden'});
      expectDeepEqual(tree.linksState.hiddenActions$.value.has(actionUuid(tree)), true, {prefix: 'published'});

      cooling.getItem().getStateStore().setState('a', 2);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), true, {prefix: 'shown again'});
    });
  };

  test('Node meta link hides an action from input values', hidesFromInputs({}));

  test('Node meta link hides an action from input values in Driver mode', hidesFromInputs(driverMode));

  test('Node meta link withdraws its hiding when it stops setting it', async () => {
    const config = daysConfig(0, [hidingLink('vis', 'in:cooling/a', 'act', (v) => v === 1 ? false : undefined)], [],
      [rootAction()]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const cooling = tree.nodeTree.getNode([{idx: 0}]);
      cooling.getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'hidden'});

      cooling.getItem().getStateStore().setState('a', 2);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), true, {prefix: 'withdrawn'});
    });
  });

  test('Later node meta link can show an action hidden by hideWhen', async () => {
    const runs = {count: 0};
    const config = daysConfig(0, [hidingLink('vis', 'in:cooling/a', 'act', (v) => v === 1 ? true : undefined, runs)],
      [], [rootAction(undefined, {hideWhen: 'guard: cooling'})]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'hidden by hideWhen'});

      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(runs.count > 0, true, {prefix: 'user link ran'});
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), true, {prefix: 'shown by the later link'});
    });
  });

  const laterLinkWins = async (first: boolean, second: boolean) => {
    const config = daysConfig(0, [
      hidingLink('first', 'in:cooling/a', 'act', () => first),
      hidingLink('second', 'in:cooling/a', 'act', () => second),
    ], [], [rootAction()]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();
      return tree.linksState.isActionVisible(actionUuid(tree));
    });
  };

  test('Later node meta link overrides an earlier one', async () => {
    expectDeepEqual(await laterLinkWins(true, false), false, {prefix: 'shown then hidden'});
    expectDeepEqual(await laterLinkWins(false, true), true, {prefix: 'hidden then shown'});
  });

  test('Node meta link addresses an action shown elsewhere via visibleOn', async () => {
    const config = daysConfig(0, [hidingLink('vis', 'in:cooling/a', 'act:cooling', () => false)], [],
      [rootAction(undefined, {visibleOn: 'cooling'})]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree, [{idx: 0}])), false);
    });
  });

  test('Hidden action stays hidden after its action is recreated', async () => {
    const config = daysConfig(1, [hidingLink('vis', 'in:cooling/a', 'act:days', () => false)], [daysAction]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();
      const before = actionUuid(tree, [{idx: 1}]);
      expectDeepEqual(tree.linksState.isActionVisible(before), false, {prefix: 'before mutation'});

      const daysUuid = tree.nodeTree.getNode([{idx: 1}]).getItem().uuid;
      tree.addSubTree(daysUuid, 'day', 1).subscribe();
      flush();
      const after = actionUuid(tree, [{idx: 1}]);
      expectDeepEqual(after !== before, true, {prefix: 'action recreated'});
      expectDeepEqual(tree.linksState.isActionVisible(after), false, {prefix: 'after mutation'});
    });
  });

  test('Removed node meta link stops hiding its action', async () => {
    const config = daysConfig(1, [hidingLink('vis', 'in:days/first(day)/a', 'act', () => false)], [], [rootAction()]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const day = tree.nodeTree.getNode([{idx: 1}, {idx: 0}]);
      day.getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'hidden'});

      tree.removeSubtree(day.getItem().uuid).subscribe();
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), true, {prefix: 'after link removal'});
    });
  });

  test('Node meta hiding stays when hideWhen stops matching', async () => {
    const config = daysConfig(1, [hidingLink('vis', 'in:cooling/a', 'act', (v) => v === 1 ? false : undefined)], [],
      [rootAction(undefined, {hideWhen: 'guard: days/first(day)'})]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'hidden by hideWhen'});
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();

      tree.removeSubtree(tree.nodeTree.getNode([{idx: 1}, {idx: 0}]).getItem().uuid).subscribe();
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'still hidden by the link'});
    });
  });

  test('hideWhen that starts matching overrides an earlier node meta link', async () => {
    const config = daysConfig(0, [hidingLink('vis', 'in:cooling/a', 'act', (v) => v === 1 ? true : undefined)], [],
      [rootAction(undefined, {hideWhen: 'guard: days/first(day)'})]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), true, {prefix: 'shown by the link'});

      tree.addSubTree(tree.nodeTree.getNode([{idx: 1}]).getItem().uuid, 'day', 0).subscribe();
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'hidden by hideWhen'});
    });
  });

  test('Node meta hiding survives an unrelated tree change', async () => {
    const config = daysConfig(0, [hidingLink('vis', 'in:cooling/a', 'act', (v) => v === 1 ? false : undefined)], [],
      [rootAction()]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'hidden'});

      tree.addSubTree(tree.nodeTree.getNode([{idx: 1}]).getItem().uuid, 'day', 0).subscribe();
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree)), false, {prefix: 'after adding a step'});
    });
  });

  test('hideWhen hides an action shown elsewhere via visibleOn', async () => {
    const config = daysConfig(0, [], [], [rootAction(undefined, {visibleOn: 'cooling', hideWhen: 'guard: cooling'})]);
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(actionUuid(tree, [{idx: 0}])), false);
    });
  });

  test('hideWhen hides every copy of a base-expanded action', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'parallel',
      stepTypes: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
      initialSteps: [{id: 'step1'}, {id: 'step1'}, {id: 'step1'}],
      actions: [{
        id: 'act',
        base: 'base:expand(step1)',
        position: 'none',
        from: 'in:same(@base)/a',
        to: 'out:same(@base)/b',
        hideWhen: 'guard: first(step1)',
        handler() {},
      }],
    };
    const pconf = await getProcessedConfig(config);
    return testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const actions = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid) ?? [];
      expectDeepEqual(actions.length, 3, {prefix: 'copies'});
      expectDeepEqual(actions.map((action) => tree.linksState.isActionVisible(action.uuid)), [false, false, false]);
    });
  });
});

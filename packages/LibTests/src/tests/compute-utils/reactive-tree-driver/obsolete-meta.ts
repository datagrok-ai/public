import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {switchMap} from 'rxjs/operators';
import {makeValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {FuncCallNode, PipelineNodeBase} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';


category('ComputeUtils: Driver obsolete meta cleanup', async () => {
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
        nqName: 'LibTests:TestDiv2',
      },
    ],
    links: [{
      id: 'link1',
      from: 'in1:step1/res',
      to: 'out1:step3/a',
      handler({controller}) {
        controller.setAll('out1', 10, 'info');
      },
    }, {
      id: 'link2',
      from: 'in1:step1/res',
      to: 'out1:step3/a',
      type: 'validator',
      handler({controller}) {
        controller.setValidation('out1', makeValidationResult({warnings: ['test warning']}));
      },
    }, {
      id: 'link3',
      from: 'in1:step1/res',
      to: 'out1:step3/a',
      type: 'meta',
      handler({controller}) {
        controller.setViewMeta('out1', {key: 'val'});
      },
    }],
  };

  let testScheduler: TestScheduler;
  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Remove obsolete consistency', async () => {
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const outNode = tree.nodeTree.getNode([{idx: 2}]);
      const [link1] = tree.linksState.links.values();

      cold('-a').subscribe(() => {
        outNode.getItem().getStateStore().setState('b', 1, 'restricted');
        link1.trigger();
      });
      cold('5ms a').subscribe(() => {
        tree.runMutateTree().subscribe();
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('ab', {a: undefined, b: 10});
      expectObservable(outNode.getItem().getStateStore().getStateChanges('b'), '^ 1000ms !').toBe('ab', {a: undefined, b: 1});
      expectObservable(tree.getConsistency()[outNode.getItem().uuid], '^ 1000ms !').toBe('a (bc)d', {
        a: {},
        b: {
          'b': {
            'restriction': 'restricted',
            'inconsistent': false,
            'assignedValue': 1,
          },
        },
        c: {
          'b': {
            'restriction': 'restricted',
            'inconsistent': false,
            'assignedValue': 1,
          },
          'a': {
            'restriction': 'info',
            'inconsistent': false,
            'assignedValue': 10,
          },
        },
        d: {
          'a': {
            'restriction': 'info',
            'inconsistent': false,
            'assignedValue': 10,
          },
        },
      });
    });
  });

  test('Remove obsolete meta', async () => {
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const outNode = tree.nodeTree.getNode([{idx: 2}]);
      const [,, link3] = tree.linksState.links.values();

      cold('-a').subscribe(() => {
        (outNode.getItem() as FuncCallNode).instancesWrapper.setMeta('b', {key: 'val2'});
        link3.trigger();
      });
      cold('5ms a').subscribe(() => {
        tree.runMutateTree().subscribe();
      });
      const aMeta$ = (outNode.getItem() as FuncCallNode).metaInfo$.pipe(
        switchMap((x) => x.a),
      );
      const bMeta$ = (outNode.getItem() as FuncCallNode).metaInfo$.pipe(
        switchMap((x) => x.b),
      );
      const outA$ = tree.getMeta()[outNode.getItem().uuid].pipe(
        switchMap((x) => x.a),
      );
      const outB$ = tree.getMeta()[outNode.getItem().uuid].pipe(
        switchMap((x) => x.b),
      );
      expectObservable(aMeta$, '^ 1000ms !').toBe('ab', {a: undefined, b: {key: 'val'}});
      expectObservable(bMeta$, '^ 1000ms !').toBe('ab---c', {a: undefined, b: {key: 'val2'}, c: undefined});
      expectObservable(outA$, '^ 1000ms !').toBe('ab', {
        a: undefined,
        b: {key: 'val'},
      });
      expectObservable(outB$, '^ 1000ms !').toBe('ab---c', {
        a: undefined,
        b: {key: 'val2'},
        c: undefined,
      });
    });
  });

  test('Remove obsolete validation', async () => {
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.linksState.forceInitialMetaRun = true;
      tree.init().subscribe();
      const outNode = tree.nodeTree.getNode([{idx: 2}]);

      cold('-a').subscribe(() => {
        (outNode.getItem() as FuncCallNode).instancesWrapper.setValidation('a', 'asdf', makeValidationResult({warnings: ['test warning2']}));
      });
      cold('10ms a').subscribe(() => {
        tree.runMutateTree().subscribe();
      });
      expectObservable(tree.getValidations()[outNode.getItem().uuid], '^ 1000ms !').toBe('ab 8ms (ba)', {
        a: {
          'a': {
            'errors': [],
            'warnings': [
              {
                'description': 'test warning',
              },
            ],
            'notifications': [],
          },
        },
        b: {
          'a': {
            'errors': [],
            'warnings': [
              {
                'description': 'test warning',
              },
              {
                'description': 'test warning2',
              },
            ],
            'notifications': [],
          },
        },
      });
    });
  });

  test('Remove obsolete tags', async () => {
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.linksState.forceInitialMetaRun = true;
      tree.init().subscribe();
      const outNode = tree.nodeTree.root;
      cold('-a').subscribe(() => {
        (outNode.getItem() as PipelineNodeBase).nodeDescription.setState('tags', {'1234': ['tag1', 'tag2']});
      });
      cold('10ms a').subscribe(() => {
        tree.runMutateTree().subscribe();
      });
      expectObservable((outNode.getItem() as PipelineNodeBase).nodeDescription.getStateChanges('tags')).toBe('ab 8ms c', {
        a: undefined,
        b: {
          '1234': [
            'tag1',
            'tag2',
          ],
        },
        c: {},
      });
    });
  });
});

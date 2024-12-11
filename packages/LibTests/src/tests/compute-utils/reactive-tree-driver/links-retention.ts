import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {Subject} from 'rxjs';
import { filter, take} from 'rxjs/operators';

category('ComputeUtils: Driver links retention', async () => {
  let testScheduler: TestScheduler;
  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Should retain unchanged links after tree mutations', async () => {
    const config1: PipelineConfiguration = {
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
        {
          id: 'step3',
          nqName: 'LibTests:TestSub2',
        },
      ],
      initialSteps: [
        { id: 'step1' },
        { id: 'step2' }
      ],
      links: [
        {
          id: 'link1',
          base: 'base:expand(step1)',
          from: 'from:same(@base, step1)/res',
          to: 'to:after+(@base, step2)/a',
        },
        {
          id: 'link1',
          base: 'base:expand(step2)',
          from: 'from:same(@base, step2)/res',
          to: 'to:after+(@base, step3)/a',
        }
      ]
    };
    const pconf = await getProcessedConfig(config1);
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const links = [...tree.linksState.links.values()];
      expectDeepEqual(links.length, 1, { prefix: 'Initial links count' });
      const updateFinished$ = new Subject<true>();
      cold('-a').subscribe(() => {
        tree.addSubTree(tree.nodeTree.root.getItem().uuid, 'step3', 2).subscribe();
        tree.globalROLocked$.pipe(filter(x => !x), take(1)).subscribe(() => updateFinished$.next(true));
      });
      updateFinished$.pipe(take(1)).subscribe(() => {
        const nlinks = [...tree.linksState.links.values()];
        expectDeepEqual(nlinks.length, 2, { prefix: 'Final links count' });
        expectDeepEqual(links[0].uuid, nlinks[0].uuid, { prefix: 'Retained link' });
      });
    });
  });

  test('Should retain unchanged links with actions after tree mutations', async () => {
    const config1: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'sequential',
      stepTypes: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          actions: [{
            id: 'action1',
            from: 'from:a',
            to: 'to:b',
            position: 'none',
            handler({controller}) {
              controller.setAll('to', 10);
            },
          }]
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
          actions: [{
            id: 'action2',
            from: 'from:a',
            to: 'to:b',
            position: 'none',
            handler({controller}) {
              controller.setAll('to', 10);
            },
          }]
        },
        {
          id: 'step3',
          nqName: 'LibTests:TestSub2',
        },
      ],
      initialSteps: [
        { id: 'step1' },
        { id: 'step2' },
        { id: 'step3' },
      ],
      links: [
        {
          id: 'link1',
          type: 'validator',
          base: 'base:expand(step1)',
          from: 'from:same(@base, step1)/res',
          to: 'to:after+(@base, step2)/a',
          actions: 'actions:same(@base, step1)',
          handler({controller}) {
            console.log('link1');
            const act = controller.getValidationAction('actions', 'action1');
            expectDeepEqual(typeof act, 'string', { prefix: 'Link1 should get action1' })
            controller.setValidation('to');
          }
        },
        {
          id: 'link1',
          type: 'validator',
          base: 'base:expand(step2)',
          from: 'from:same(@base, step2)/res',
          to: 'to:after+(@base, step3)/a',
          actions: 'actions:same(@base, step2)',
          handler({controller}) {
            console.log('link2');
            const act = controller.getValidationAction('actions', 'action2');
            expectDeepEqual(typeof act, 'string', { prefix: 'Link1 should get action2' });
            controller.setValidation('to');
          }
        }
      ]
    };
    const pconf = await getProcessedConfig(config1);
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const links = [...tree.linksState.links.values()];
      const actions = [...tree.linksState.actions.values()];
      expectDeepEqual(links.length, 2, { prefix: 'Initial links count' });
      expectDeepEqual(actions.length, 2, { prefix: 'Initial actions count' });

      const updateFinished$ = new Subject<true>();
      cold('-a').subscribe(() => {
        const node = tree.nodeTree.getNode([{idx: 2}]);
        tree.removeSubtree(node.getItem().uuid).subscribe();
        tree.globalROLocked$.pipe(filter(x => !x), take(1)).subscribe(() => updateFinished$.next(true));
      });
      updateFinished$.pipe(take(1)).subscribe(() => {
        const nlinks = [...tree.linksState.links.values()];
        const nActions = [...tree.linksState.actions.values()];
        expectDeepEqual(nlinks.length, 1, { prefix: 'Final links count' });
        expectDeepEqual(nActions.length, 2, { prefix: 'Final actions count' });
        expectDeepEqual(links[0].uuid, nlinks[0].uuid, { prefix: 'Retained link' });
        expectDeepEqual(nActions[0].uuid, nActions[0].uuid, { prefix: 'Retained action1' });
        expectDeepEqual(nActions[1].uuid, nActions[1].uuid, { prefix: 'Retained action2' });
      });
    });
  });

});

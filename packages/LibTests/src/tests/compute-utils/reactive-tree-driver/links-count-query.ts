import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {parseLinkIO} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/LinkSpec';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {IControllerBase} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/RuntimeControllers';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrows} from '../../../test-utils';

const shifts = (links: any[]): PipelineConfiguration => ({
  id: 'pipeline1',
  type: 'static',
  steps: [{
    id: 'shifts',
    type: 'parallel',
    stepTypes: [{id: 'day', nqName: 'LibTests:TestAdd2'}, {id: 'night', nqName: 'LibTests:TestMul2'}],
    initialSteps: [{id: 'day'}, {id: 'night'}, {id: 'day'}],
    links,
  }],
});

category('ComputeUtils: Driver links count query', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('count parses as the last segment of an input query', async () => {
    expectDeepEqual(parseLinkIO('n:count()', 'input')[0].segments,
      [{type: 'selector', selector: 'count', ids: [], stopIds: [], ref: undefined}]);
    expectDeepEqual(parseLinkIO('n:shifts/count(day|night)', 'input')[0].segments[1],
      {type: 'selector', selector: 'count', ids: ['day', 'night'], stopIds: [], ref: undefined});
    expectThrows(() => parseLinkIO('n:count()', 'output'), /count selector is allowed only in input queries/);
    expectThrows(() => parseLinkIO('n:count(day)/a', 'input'), /count selector must be the last segment/);
    expectThrows(() => parseLinkIO('n(call):count()', 'input'), /\(call\) flag is not allowed with count/);
    expectThrows(() => parseLinkIO('n(template):count(day)', 'input'), /\(template\) flag is not allowed with count/);
    expectThrows(() => parseLinkIO('n:count(@ref, day)', 'input'), /non-ref selector with ref/);
  });

  test('A count input is the number of matching child steps', async () => {
    let counts: any;
    const pconf = await getProcessedConfig(shifts([{
      id: 'counter',
      from: ['all:count()', 'days:count(day)', 'none:count(evening)'],
      to: 'out:first(night)/a',
      handler({controller}: {controller: IControllerBase}) {
        counts = ['all', 'days', 'none'].map((name) => controller.getAll(name));
      },
    }]));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      tree.linksState.links.forEach((link) => link.trigger());
      flush();
      expectDeepEqual(counts, [[3], [2], [0]]);
    });
  });

  test('A count gives one value per matched parent', async () => {
    let counts: any;
    const pconf = await getProcessedConfig({
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'a', type: 'parallel', stepTypes: [{id: 'day', nqName: 'LibTests:TestAdd2'}], initialSteps: [{id: 'day'}]},
        {id: 'b', type: 'parallel', stepTypes: [{id: 'day', nqName: 'LibTests:TestAdd2'}],
          initialSteps: [{id: 'day'}, {id: 'day'}]},
      ],
      links: [{
        id: 'counter', from: 'n:all(a|b)/count(day)', to: 'out:a/first(day)/a',
        handler({controller}: {controller: IControllerBase}) {
          counts = controller.getAll('n');
        },
      }],
    });
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      tree.linksState.links.forEach((link) => link.trigger());
      flush();
      expectDeepEqual(counts, [1, 2]);
    });
  });

  test('A count link runs again only when the number changes', async () => {
    const seen: number[] = [];
    const pconf = await getProcessedConfig(shifts([{
      id: 'counter', type: 'nodemeta', from: 'n:count(day)', to: 'title:title',
      handler({controller}: {controller: IControllerBase}) {
        seen.push(controller.getFirst<number>('n')!);
      },
    }]));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const shiftsNode = tree.nodeTree.getNode([{idx: 0}]);
      tree.addSubTree(shiftsNode.getItem().uuid, 'day', 0).subscribe();
      flush();
      expectDeepEqual(seen, [3], {prefix: 'after adding a day'});
      tree.moveSubtree(shiftsNode.getChildren()[0].item.getItem().uuid, 3).subscribe();
      flush();
      tree.addSubTree(shiftsNode.getItem().uuid, 'night', 0).subscribe();
      flush();
      expectDeepEqual(seen, [3], {prefix: 'after a move and a night'});
      tree.removeSubtree(shiftsNode.getChildren()[1].item.getItem().uuid).subscribe();
      flush();
      expectDeepEqual(seen, [3, 2], {prefix: 'after removing a day'});
    });
  });
});

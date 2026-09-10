import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {BaseTree, NodeAddress} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/BaseTree';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';

function addr(...idxs: number[]): NodeAddress {
  return idxs.map((idx) => ({id: '', idx}));
}

let matchedInputs: string[] = [];

const orderedIOConfig: PipelineConfiguration = {
  id: 'root',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {id: 'step2', nqName: 'LibTests:TestMul2'},
  ],
  links: [{
    id: 'order-link',
    from: ['inB:step2/a', 'inA:step1/a'],
    to: 'out1:step2/b',
    handler({controller}) {
      matchedInputs = [...(controller as any).getMatchedInputs()];
      controller.setAll('out1', 42);
    },
  }],
};

category('ComputeUtils: Driver link io order', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('compareAddresses orders a prefix before its descendants', async () => {
    expectDeepEqual(BaseTree.compareAddresses(addr(0), addr(0, 1)) < 0, true, {prefix: 'prefix before descendant'});
    expectDeepEqual(BaseTree.compareAddresses(addr(0, 1), addr(0)) > 0, true, {prefix: 'descendant after prefix'});
  });

  test('compareAddresses orders siblings and equal addresses', async () => {
    expectDeepEqual(BaseTree.compareAddresses(addr(0, 1), addr(0, 2)) < 0, true, {prefix: 'sibling order'});
    expectDeepEqual(BaseTree.compareAddresses(addr(2), addr(1, 5)) > 0, true, {prefix: 'first level wins'});
    expectDeepEqual(BaseTree.compareAddresses(addr(1, 2), addr(1, 2)), 0, {prefix: 'equal addresses'});
    expectDeepEqual(BaseTree.compareAddresses(addr(), addr(0)) < 0, true, {prefix: 'root before children'});
  });

  test('Matched inputs iterate in tree order', async () => {
    const pconf = await getProcessedConfig(orderedIOConfig);
    matchedInputs = [];
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      cold('---a').subscribe(() => {
        expectDeepEqual(matchedInputs, ['inA', 'inB'], {prefix: 'Matched inputs order'});
      });
    });
  });
});

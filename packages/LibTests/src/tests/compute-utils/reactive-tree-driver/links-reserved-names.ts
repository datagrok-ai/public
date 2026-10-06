import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrowsAsync} from '../../../test-utils';

const RESERVED = 'LibTests:TestReservedNames';

category('ComputeUtils: Driver reserved names', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('An input named call does not break evaluated choices', async () => {
    const pconf: any =
      await getProcessedConfig({id: 'pipeline1', type: 'static', steps: [{id: 'step', nqName: RESERVED}]});
    const choices = pconf.steps[0].links.find((link: any) => link.id === '::city:choices::meta');
    const names = choices.from.map((item: any) => item.name);
    expectDeepEqual(['call', '$call'].map((name) => names.includes(name)), [true, true]);
  });

  test('A rule reads an alias named all', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'r', type: 'rule', from: 'all:step1/a', to: 't:step2/a',
        effects: [{effect: 'set', targets: 't', value: {'+': [{var: 'all'}, {var: '$all.all.0'}]}}],
      }],
    };
    const pconf = await getProcessedConfig(config);
    let value: any;
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const step1 = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const step2 = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore();
      cold('-a').subscribe(() => step1.setState('a', 3));
      cold('--a').subscribe(() => value = step2.getState('a'));
    });
    expectDeepEqual(value, 6);
  });

  test('Expression checks see inputs named like driver aliases', async () => {
    const pconf = await getProcessedConfig({id: 'pipeline1', type: 'static', steps: [{id: 'step', nqName: RESERVED}]});
    let inputs: string[] = [];
    testScheduler.run(() => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, defaultValidators: true});
      tree.init().subscribe();
      const link = [...tree.linksState.links.values()].find((item) => item.matchInfo.spec.id === '::x:validator')!;
      inputs = Object.keys(link.matchInfo.inputs);
    });
    expectDeepEqual(['call', 'all', 'table', 'target', 'literals'].map((name) => inputs.includes(name)),
      [true, true, true, true, true]);
  });

  test('Rules and checks reject aliases starting with $', async () => {
    const steps = [{id: 'step1', nqName: 'LibTests:TestAdd2'}, {id: 'step2', nqName: 'LibTests:TestMul2'}];
    await expectThrowsAsync(() => getProcessedConfig({id: 'p', type: 'static', steps, links: [{
      id: 'r', type: 'rule', from: '$x:step1/a', to: 't:step2/a', effects: [{effect: 'set', targets: 't', value: 1}],
    }]}), /reserved for the driver/);
    await expectThrowsAsync(() => getProcessedConfig({id: 'p', type: 'static', steps, links: [{
      id: 'c', type: 'check', io: 'step2/a', check: {validator: 'k > 1'}, vars: {$k: 'step1/a'},
    }]} as any), /reserved for the driver/);
  });
});

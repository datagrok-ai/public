import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

category('ComputeUtils: Driver links inbound outgoing', async () => {
  test('Root changes', async () => {
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
      links: [{
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const link = [...tree.linksState.links.values()][0];
    expectDeepEqual(tree.linksState.isAffected([], link), true, {prefix: 'isAffected'});
    expectDeepEqual(tree.linksState.isInbound([], link), false, {prefix: 'isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link), false, {prefix: 'isOutgoing'});
  }, {skipReason: 'TODO'});

  test('Root changes with an offset', async () => {
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
      links: [{
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const link = [...tree.linksState.links.values()][0];
    expectDeepEqual(tree.linksState.isAffected([], link, 1), true, {prefix: 'isAffected'});
    expectDeepEqual(tree.linksState.isInbound([], link, 1), true, {prefix: 'isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link, 1), false, {prefix: 'isOutgoing'});
  }, {skipReason: 'TODO'});

  test('Nested changes', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'pipeline1',
          type: 'static',
          steps: [
            {
              id: 'step1',
              nqName: 'LibTests:TestAdd2',
            },
          ],
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:pipeline1/step1/res',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const link = [...tree.linksState.links.values()][0];
    expectDeepEqual(tree.linksState.isAffected([{idx: 0}], link), true, {prefix: 'isAffected'});
    expectDeepEqual(tree.linksState.isInbound([{idx: 0}], link), false, {prefix: 'isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([{idx: 0}], link), true, {prefix: 'isOutgoing'});
  }, {skipReason: 'TODO'});

  test('Nested changes with an offset', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'pipeline1',
          type: 'static',
          steps: [
            {
              id: 'step1',
              nqName: 'LibTests:TestAdd2',
            },
            {
              id: 'step2',
              nqName: 'LibTests:TestAdd2',
            },
          ],
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:pipeline1/step2/res',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const link = [...tree.linksState.links.values()][0];
    expectDeepEqual(tree.linksState.isAffected([{idx: 0}], link, 1), true, {prefix: 'isAffected'});
    expectDeepEqual(tree.linksState.isInbound([{idx: 0}], link, 1), false, {prefix: 'isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([{idx: 0}], link, 1), true, {prefix: 'isOutgoing'});
  }, {skipReason: 'TODO'});

  test('Nested changes with multiple io', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'pipeline1',
          type: 'static',
          steps: [
            {
              id: 'step1',
              nqName: 'LibTests:TestAdd2',
            },
            {
              id: 'step2',
              nqName: 'LibTests:TestAdd2',
            },
          ],
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
        {
          id: 'step3',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: ['in1:pipeline1/step2/res', 'in2:pipeline1/step2/res'],
        to: ['out1:step2/a', 'out2:step3/a'],
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const link = [...tree.linksState.links.values()][0];
    expectDeepEqual(tree.linksState.isAffected([{idx: 0}], link, 1), true, {prefix: 'isAffected'});
    expectDeepEqual(tree.linksState.isInbound([{idx: 0}], link, 1), false, {prefix: 'isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([{idx: 0}], link, 1), true, {prefix: 'isOutgoing'});
  }, {skipReason: 'TODO'});

});

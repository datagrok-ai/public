import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

category('ComputeUtils: Driver links inbound outgoing', async () => {
  const config1: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'static',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestAdd2',
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestSub2',
      },
      {
        id: 'pipeline3',
        type: 'static',
        steps: [
          {
            id: 'step1',
            nqName: 'LibTests:TestAdd2',
          },
          {
            id: 'step2',
            nqName: 'LibTests:TestSub2',
          },
          {
            id: 'step3',
            nqName: 'LibTests:TestMul2',
          },
          {
            id: 'step4',
            nqName: 'LibTests:TestDiv2',
          },
        ],
        links: [{
          id: 'link5',
          from: 'in1:step2/res',
          to: 'out1:step3/a',
        }],
      },
      {
        id: 'step4',
        nqName: 'LibTests:TestMul2',
      },
      {
        id: 'step5',
        nqName: 'LibTests:TestDiv2',
      },
    ],
    links: [{
      id: 'link1',
      from: 'in1:step2/res',
      to: 'out1:pipeline3/step1/a',
    }, {
      id: 'link2',
      from: 'in1:pipeline3/step4/res',
      to: 'out1:step4/a',
    }, {
      id: 'link3',
      from: 'in1:step1/res',
      to: 'out1:step5/a',
    }, {
      id: 'link4',
      from: 'in1:step4/res',
      to: 'out1:step5/b',
    }],
  };

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
        {
          id: 'step3',
          nqName: 'LibTests:TestDiv2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step3/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const [link] = [...tree.linksState.links.values()];
    expectDeepEqual(tree.linksState.isInbound([], link), false, {prefix: 'isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link), false, {prefix: 'isOutgoing'});
  });

  test('Root add item', async () => {
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
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const [link] = [...tree.linksState.links.values()];
    expectDeepEqual(tree.linksState.isInbound([], link, 1), true, {prefix: 'isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link, 1), false, {prefix: 'isOutgoing'});
  });

  test('Root remove item', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
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
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const [link] = [...tree.linksState.links.values()];
    expectDeepEqual(tree.linksState.isInbound([], link, undefined, 1), true, {prefix: 'isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link, undefined), false, {prefix: 'isOutgoing'});
  });

  test('Root move item', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'step3',
          nqName: 'LibTests:TestMul2',
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestDiv2',
        },
        {
          id: 'step4',
          nqName: 'LibTests:TestSub2',
        },
      ],
      links: [
        {
          id: 'link1',
          from: 'in1:step1/res',
          to: 'out1:step3/a',
        },
        {
          id: 'link2',
          from: 'in1:step3/res',
          to: 'out1:step4/a',
        },
      ],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const [link1, link2] = [...tree.linksState.links.values()];
    expectDeepEqual(tree.linksState.isInbound([], link1, 1, 2), true, {prefix: 'link1 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link1, 1), false, {prefix: 'link1 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link2, 1, 2), false, {prefix: 'link2 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link2, 1), true, {prefix: 'link2 isOutgoing'});
  });


  test('Nested add nested item', async () => {
    const pconf = await getProcessedConfig(config1);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const [link1, link2, link3, link4, link5] = [...tree.linksState.links.values()];
    expectDeepEqual(tree.linksState.isInbound([{idx: 2}], link1, 1), false, {prefix: 'link1 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([{idx: 2}], link1, 1), false, {prefix: 'link1 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([{idx: 2}], link2, 1), false, {prefix: 'link2 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([{idx: 2}], link2, 1), false, {prefix: 'link2 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([{idx: 2}], link3, 1), false, {prefix: 'link3 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([{idx: 2}], link3, 1), false, {prefix: 'link3 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([{idx: 2}], link4, 1), false, {prefix: 'link4 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([{idx: 2}], link4, 1), false, {prefix: 'link4 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([{idx: 2}], link5, 1), false, {prefix: 'link5 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([{idx: 2}], link5, 1), true, {prefix: 'link5 isOutgoing'});
  });

  test('Nested remove item', async () => {
    const pconf = await getProcessedConfig(config1);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const [link1, link2, link3, link4, link5] = [...tree.linksState.links.values()];
    expectDeepEqual(tree.linksState.isInbound([], link1, undefined, 2), true, {prefix: 'link1 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link1, undefined), false, {prefix: 'link1 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link2, undefined, 2), false, {prefix: 'link2 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link2, undefined), false, {prefix: 'link2 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link3, undefined, 2), true, {prefix: 'link3 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link3, undefined), false, {prefix: 'link3 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link4, undefined, 2), false, {prefix: 'link4 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link4, undefined), false, {prefix: 'link4 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link5, undefined, 2), false, {prefix: 'link5 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link5, undefined), false, {prefix: 'link5 isOutgoing'});
  });

  test('Nested move item', async () => {
    const pconf = await getProcessedConfig(config1);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    tree.init().subscribe();
    const [link1, link2, link3, link4, link5] = [...tree.linksState.links.values()];
    expectDeepEqual(tree.linksState.isInbound([], link1, 2, 1), true, {prefix: 'link1 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link1, 2), false, {prefix: 'link1 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link2, 2, 1), false, {prefix: 'link2 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link2, 2), true, {prefix: 'link2 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link3, 2, 1), true, {prefix: 'link3 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link3, 2), false, {prefix: 'link3 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link4, 2, 1), false, {prefix: 'link4 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link4, 2), false, {prefix: 'link4 isOutgoing'});
    expectDeepEqual(tree.linksState.isInbound([], link5, 2, 1), false, {prefix: 'link5 isInbound'});
    expectDeepEqual(tree.linksState.isOutgoing([], link5, 2), false, {prefix: 'link5 isOutgoing'});
  });
});

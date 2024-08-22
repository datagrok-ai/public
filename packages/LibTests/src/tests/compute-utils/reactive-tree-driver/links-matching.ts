import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import { PipelineConfiguration } from '@datagrok-libraries/compute-utils';
import { getProcessedConfig, PipelineConfigurationStaticProcessed } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import { StateTree } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import { matchLink } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/link-matching';
import { PipelineLinkConfiguration } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import { LinkSpecString } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import { snapshotCompare } from '../../../test-utils';

category('ComputeUtils: Driver links matching', async () => {
  test('Links match static pipeline', async () => {
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
        from: 'in1:step1/res',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig({config: pconf});
    await tree.initAll().toPromise();
    const links = (tree.getRoot().getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfo = matchLink(tree, [], links[0]);
    await snapshotCompare(matchInfo, 'Links match static pipeline');
  });

  test('Links match nested static pipeline', async () => {
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
            }
          ]
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
    const tree = StateTree.fromConfig({config: pconf});
    await tree.initAll().toPromise();
    const links = (tree.getRoot().getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfo = matchLink(tree, [], links[0]);
    await snapshotCompare(matchInfo, 'Links match nested static pipeline');
  });

  test('Links with multiple io match nested static pipeline', async () => {
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
            }
          ]
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: ['in1:pipeline1/step1/res', 'in2:pipeline1/step1/a'],
        to: ['out1:step2/a', 'out2:step2/b'],
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig({config: pconf});
    await tree.initAll().toPromise();
    const links = (tree.getRoot().getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfo = matchLink(tree, [], links[0]);
    await snapshotCompare(matchInfo, 'Links with multiple io match nested static pipeline');
  });

  test('Links matching all parallel pipeline items', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
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
              id: 'stepAdd'
            },
            {
              id: 'stepMul'
            },
            {
              id: 'stepMul'
            },
            {
              id: 'stepAdd'
            },
          ]
        },
        {
          id: 'step3',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: ['out1:pipelinePar/all(stepAdd)/a', 'out2:pipelinePar/all(stepMul)/a'],
      }, {
        id: 'link2',
        from: ['in1:pipelinePar/all(stepAdd)/res', 'in2:pipelinePar/all(stepMul)/res'],
        to: 'out1:step3/a'
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig({config: pconf});
    await tree.initAll().toPromise();
    const links = (tree.getRoot().getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map(link => matchLink(tree, [], link));
    await snapshotCompare(matchInfos, 'Links matching all parallel pipeline items');
  });

});

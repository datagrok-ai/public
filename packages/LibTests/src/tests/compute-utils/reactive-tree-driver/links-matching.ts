import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig, PipelineConfigurationStaticProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {parseLinkIO} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/LinkSpec';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {MatchInfo, matchLink} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/link-matching';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {snapshotCompare} from '../../../test-utils';

function cleanMatchInfos(data: (MatchInfo[] | undefined)[]) {
  for (const infos of data) {
    if (infos) {
      infos.forEach((info: any) => {
        delete info.basePathUUID;
        delete info.inputsUUID;
        delete info.outputsUUID;
      });
    }
  }
}

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
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfo = matchLink(tree, [], links[0]);
    cleanMatchInfos([matchInfo]);
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
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfo = matchLink(tree, [], links[0]);
    cleanMatchInfos([matchInfo]);
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
        from: ['in1:pipeline1/step1/res', 'in2:pipeline1/step1/a'],
        to: ['out1:step2/a', 'out2:step2/b'],
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfo = matchLink(tree, [], links[0]);
    cleanMatchInfos([matchInfo]);
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
              id: 'stepAdd',
            },
            {
              id: 'stepMul',
            },
            {
              id: 'stepMul',
            },
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
        from: 'in1:step1/res',
        to: ['out1:pipelinePar/all(stepAdd)/a', 'out2:pipelinePar/all(stepMul)/a'],
      }, {
        id: 'link2',
        from: ['in1:pipelinePar/all(stepAdd)/res', 'in2:pipelinePar/all(stepMul)/res'],
        to: 'out1:step3/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links matching all parallel pipeline items');
  });

  test('Links matching multiple ids', async () => {
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
              id: 'stepAdd',
            },
            {
              id: 'stepMul',
            },
            {
              id: 'stepMul',
            },
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
        from: 'in1:step1/res',
        to: ['out1:pipelinePar/all(stepAdd|stepMul)/a'],
      }, {
        id: 'link2',
        from: ['in1:pipelinePar/all(stepAdd|stepMul)/res'],
        to: 'out1:step3/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links matching multiple ids');
  });

  test('Links templates matching', async () => {
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
              id: 'stepAdd',
            },
            {
              id: 'stepMul',
            },
            {
              id: 'stepMul',
            },
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
        from: 'in1:step1/res',
        to: ['out1(template):pipelinePar/all(stepAdd|stepMul)/a|b'],
      }, {
        id: 'link2',
        from: ['in1(template):pipelinePar/all(stepAdd|stepMul)/a|b'],
        to: 'out1:step3/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links templates matching');
  });

  test('Links expaning base path', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'sequential',
      stepTypes: [
        {
          id: 'stepAdd',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'stepMul',
          nqName: 'LibTests:TestMul2',
        },
        {
          id: 'stepSub',
          nqName: 'LibTests:TestSub2',
        },
        {
          id: 'stepDiv',
          nqName: 'LibTests:TestDiv2',
        },
      ],
      initialSteps: [
        {id: 'stepAdd'},
        {id: 'stepMul'},
        {id: 'stepSub'},
        {id: 'stepDiv'},
        {id: 'stepAdd'},
        {id: 'stepMul'},
        {id: 'stepSub'},
        {id: 'stepDiv'},
      ],
      links: [{
        id: 'link1',
        from: [],
        to: [],
        base: 'base:expand(stepMul)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links expaning base path');
  });

  test('Links referencing adjacent to base segment', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'sequential',
      stepTypes: [
        {
          id: 'stepAdd',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'stepMul',
          nqName: 'LibTests:TestMul2',
        },
        {
          id: 'stepSub',
          nqName: 'LibTests:TestSub2',
        },
        {
          id: 'stepDiv',
          nqName: 'LibTests:TestDiv2',
        },
      ],
      initialSteps: [
        {id: 'stepAdd'},
        {id: 'stepMul'},
        {id: 'stepSub'},
        {id: 'stepDiv'},
        {id: 'stepAdd'},
        {id: 'stepMul'},
        {id: 'stepSub'},
        {id: 'stepDiv'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:before+(@base,stepAdd)/res',
        to: 'out1:after+(@base,stepSub)/a',
        base: 'base:expand(stepMul)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links referencing adjacent to base segment');
  });

  test('Links referencing remote from base path', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'sequential',
      stepTypes: [
        {
          id: 'stepAdd',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'stepMul',
          nqName: 'LibTests:TestMul2',
        },
        {
          id: 'stepSub',
          nqName: 'LibTests:TestSub2',
        },
        {
          id: 'stepDiv',
          nqName: 'LibTests:TestDiv2',
        },
      ],
      initialSteps: [
        {id: 'stepAdd'},
        {id: 'stepSub'},
        {id: 'stepMul'},
        {id: 'stepSub'},
        {id: 'stepDiv'},
        {id: 'stepAdd'},
        {id: 'stepSub'},
        {id: 'stepMul'},
        {id: 'stepSub'},
        {id: 'stepDiv'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:before(@base,stepAdd)/res',
        to: 'out1:after(@base,stepDiv)/a',
        base: 'base:expand(stepMul)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links referencing remote from base path');
  });

  test('Links referencing all match limits', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'sequential',
      stepTypes: [
        {
          id: 'stepAdd',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'stepMul',
          nqName: 'LibTests:TestMul2',
        },
        {
          id: 'stepSub',
          nqName: 'LibTests:TestSub2',
        },
        {
          id: 'stepDiv',
          nqName: 'LibTests:TestDiv2',
        },
      ],
      initialSteps: [
        {id: 'stepAdd'},
        {id: 'stepSub'},

        {id: 'stepAdd'},
        {id: 'stepAdd'},
        {id: 'stepMul'},
        {id: 'stepDiv'},
        {id: 'stepDiv'},

        {id: 'stepSub'},

        {id: 'stepAdd'},
        {id: 'stepAdd'},
        {id: 'stepMul'},
        {id: 'stepDiv'},
        {id: 'stepDiv'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:before*(@base,stepAdd,stepSub)/res',
        to: 'out1:after*(@base,stepDiv,stepSub)/a',
        base: 'base:expand(stepMul)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links referencing all match limits');
  });

  test('Links deep base matching', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [{
        id: 'stepAdd',
        nqName: 'LibTests:TestAdd2',
      }, {
        id: 'pipeline1',
        type: 'sequential',
        stepTypes: [
          {
            id: 'stepAdd',
            nqName: 'LibTests:TestAdd2',
          },
          {
            id: 'stepMul',
            nqName: 'LibTests:TestMul2',
          },
          {
            id: 'stepSub',
            nqName: 'LibTests:TestSub2',
          },
          {
            id: 'stepDiv',
            nqName: 'LibTests:TestDiv2',
          },
          {
            id: 'pipeline1',
            type: 'sequential',
            stepTypes: [
              {
                id: 'stepAdd',
                nqName: 'LibTests:TestAdd2',
              },
              {
                id: 'stepMul',
                nqName: 'LibTests:TestMul2',
              },
              {
                id: 'stepSub',
                nqName: 'LibTests:TestSub2',
              },
              {
                id: 'stepDiv',
                nqName: 'LibTests:TestDiv2',
              },
            ],
            initialSteps: [
              {id: 'stepAdd'},
              {id: 'stepSub'},

              {id: 'stepAdd'},
              {id: 'stepAdd'},
              {id: 'stepMul'},
              {id: 'stepDiv'},
              {id: 'stepDiv'},

              {id: 'stepSub'},

              {id: 'stepAdd'},
              {id: 'stepAdd'},
              {id: 'stepMul'},
              {id: 'stepDiv'},
              {id: 'stepDiv'},
            ],
          },
        ],
        initialSteps: [
          {id: 'stepAdd'},
          {id: 'stepMul'},
          {id: 'pipeline1'},
          {id: 'stepSub'},
          {id: 'stepDiv'},
          {id: 'stepAdd'},
          {id: 'stepMul'},
          {id: 'pipeline1'},
          {id: 'stepSub'},
          {id: 'stepDiv'},
        ],
      }],
      links: [{
        id: 'link1',
        from: [],
        to: [],
        base: 'base:pipeline1/expand(pipeline1)/expand(stepMul)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links deep base matching');
  });

  test('Links expanding deep base path', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [{
        id: 'stepAdd',
        nqName: 'LibTests:TestAdd2',
      }, {
        id: 'pipeline1',
        type: 'sequential',
        stepTypes: [
          {
            id: 'stepAdd',
            nqName: 'LibTests:TestAdd2',
          },
          {
            id: 'stepMul',
            nqName: 'LibTests:TestMul2',
          },
          {
            id: 'stepSub',
            nqName: 'LibTests:TestSub2',
          },
          {
            id: 'stepDiv',
            nqName: 'LibTests:TestDiv2',
          },
          {
            id: 'pipeline1',
            type: 'sequential',
            stepTypes: [
              {
                id: 'stepAdd',
                nqName: 'LibTests:TestAdd2',
              },
              {
                id: 'stepMul',
                nqName: 'LibTests:TestMul2',
              },
              {
                id: 'stepSub',
                nqName: 'LibTests:TestSub2',
              },
              {
                id: 'stepDiv',
                nqName: 'LibTests:TestDiv2',
              },
            ],
            initialSteps: [
              {id: 'stepAdd'},
              {id: 'stepSub'},

              {id: 'stepAdd'},
              {id: 'stepAdd'},
              {id: 'stepMul'},
              {id: 'stepDiv'},
              {id: 'stepDiv'},

              {id: 'stepSub'},

              {id: 'stepAdd'},
              {id: 'stepAdd'},
              {id: 'stepMul'},
              {id: 'stepDiv'},
              {id: 'stepDiv'},
            ],
          },
        ],
        initialSteps: [
          {id: 'stepAdd'},
          {id: 'stepMul'},
          {id: 'pipeline1'},
          {id: 'stepSub'},
          {id: 'stepDiv'},
          {id: 'stepAdd'},
          {id: 'stepMul'},
          {id: 'pipeline1'},
          {id: 'stepSub'},
          {id: 'stepDiv'},
        ],
      }],
      links: [{
        id: 'link1',
        from: 'in1:pipeline1/same(@base,pipeline1)/before*(@base,stepAdd,stepSub)/res',
        to: 'out1:pipeline1/same(@base,pipeline1)/after*(@base,stepDiv,stepSub)/a',
        base: 'base:pipeline1/expand(pipeline1)/expand(stepMul)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links expanding deep base path');
  });

  test('same(@base) without ids matches identically to same(@base, ids)', async () => {
    const buildConfig = (fromSpec: string, toSpec: string): PipelineConfiguration => ({
      id: 'pipeline1',
      type: 'static',
      steps: [{
        id: 'pipeline1',
        type: 'sequential',
        stepTypes: [
          {id: 'stepAdd', nqName: 'LibTests:TestAdd2'},
          {id: 'stepMul', nqName: 'LibTests:TestMul2'},
        ],
        initialSteps: [
          {id: 'stepAdd'},
          {id: 'stepMul'},
          {id: 'stepAdd'},
          {id: 'stepMul'},
        ],
      }],
      links: [{
        id: 'link1',
        from: fromSpec,
        to: toSpec,
        base: 'base:pipeline1/expand(stepAdd)',
      }],
    });
    const runMatch = async (fromSpec: string, toSpec: string) => {
      const pconf = await getProcessedConfig(buildConfig(fromSpec, toSpec));
      const tree = StateTree.fromPipelineConfig({config: pconf});
      await tree.init().toPromise();
      const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
      const info = links.map((link) => matchLink(tree, [], link));
      cleanMatchInfos(info);
      // Strip spec — it parses to different ASTs (with vs without ids) and isn't part of the matching outcome.
      return info.map((linkInfos) => linkInfos?.map(({spec: _spec, ...rest}: any) => rest));
    };
    const withIds = await runMatch(
      'in1:same(@base,pipeline1)/same(@base,stepAdd)/res',
      'out1:same(@base,pipeline1)/after+(@base,stepMul)/a',
    );
    const withoutIds = await runMatch(
      'in1:same(@base)/same(@base)/res',
      'out1:same(@base)/after+(@base,stepMul)/a',
    );
    expectDeepEqual(withoutIds, withIds);
  });

  test('Links not matching', async () => {
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
        to: 'out1:step2/a',
        not: 'not:step3',
      }],
    };

    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links not matching');
  });

  test('Links optional matching', async () => {
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
        from: ['in1:pipelinePar/stepAdd/res', 'in2(optional):pipelinePar/stepMul/res'],
        to: ['out1:step3/a'],
        handler({controller}) {
          controller.setAll('out', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links optional matching');
  });

  test('Links single tag matching', async () => {
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
              tags: ['step1'],
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
        from: ['in1:#first(step1)/a'],
        to: ['out1:step3/a'],
        handler({controller}) {
          controller.setAll('out', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links single tag matching');
  });

  test('Links single tag with prefix matching', async () => {
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
              tags: ['step1'],
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
        from: ['in1:pipelinePar/#first(step1)/a'],
        to: ['out1:step3/a'],
        handler({controller}) {
          controller.setAll('out', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links single tag with prefix matching');
  });

  test('Links multiple tags matching', async () => {
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
              tags: ['step1'],
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
            {
              id: 'stepMul',
            },
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
        from: ['in1:#all(step1)/a'],
        to: ['out1:step3/a'],
        handler({controller}) {
          controller.setAll('out', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links multiple tags matching');
  });

  test('Links relative tags matching', async () => {
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
              tags: ['taggedStep'],
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
            {
              id: 'stepMul',
            },
            {
              id: 'stepAdd',
            },
            {
              id: 'stepMul',
            },
            {
              id: 'stepMul',
            },
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
        base: 'base:#expand(taggedStep)',
        from: 'in1:#same(@base,taggedStep)/res',
        to: 'out1:#after(@base,taggedStep)/a',
        handler({controller}) {
          controller.setAll('out', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links relative tags matching');
  });

  test('Links nested relative tags matching', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'pipelinePar',
          type: 'sequential',
          stepTypes: [
            {
              id: 'stepAdd',
              nqName: 'LibTests:TestAdd2',
              tags: ['taggedStep'],
            },
            {
              id: 'pipelineNested1',
              type: 'sequential',
              stepTypes: [
                {
                  id: 'stepAdd',
                  nqName: 'LibTests:TestAdd2',
                  tags: ['taggedStep'],
                },
                {
                  id: 'stepMul',
                  nqName: 'LibTests:TestMul2',
                },
              ],
              initialSteps: [
                {
                  id: 'stepMul',
                },
                {
                  id: 'stepAdd',
                },
              ],
            },
            {
              id: 'pipelineNested2',
              type: 'sequential',
              stepTypes: [
                {
                  id: 'stepAdd',
                  nqName: 'LibTests:TestAdd2',
                  tags: ['taggedStep'],
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
                {
                  id: 'stepMul',
                },
                {
                  id: 'stepAdd',
                },
              ],
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
            {
              id: 'stepMul',
            },
            {
              id: 'pipelineNested1',
            },
            {
              id: 'stepAdd',
            },
            {
              id: 'stepMul',
            },
            {
              id: 'pipelineNested2',
            },
            {
              id: 'stepMul',
            },
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
        base: 'base:#expand(taggedStep)',
        from: 'in1:#same(@base,taggedStep)/res',
        to: 'out1:#after(@base,taggedStep)/a',
        handler({controller}) {
          controller.setAll('out', 1);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links nested relative tags matching');
  });

  test('Links selector and tags matching', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'pipeline2',
          type: 'static',
          tags: ['taggedPipeline'],
          steps: [
            {
              id: 'step1',
              nqName: 'LibTests:TestMul2',
              tags: ['taggedStep'],
            },
            {
              id: 'pipeline3',
              type: 'static',
              tags: ['taggedPipeline'],
              steps: [
                {
                  id: 'step1',
                  nqName: 'LibTests:TestMul2',
                  tags: ['taggedStep'],
                }
              ]
            },
            {
              id: 'step2',
              nqName: 'LibTests:TestMul2',
              tags: ['taggedStep'],
            },
          ]
        },
        {
          id: 'pipeline4',
          type: 'static',
          tags: ['taggedPipeline'],
          steps: [
            {
              id: 'step1',
              nqName: 'LibTests:TestMul2',
              tags: ['taggedStep'],
            },
            {
              id: 'pipeline5',
              type: 'static',
              tags: ['taggedPipeline'],
              steps: [
                {
                  id: 'step1',
                  nqName: 'LibTests:TestMul2',
                  tags: ['taggedStep'],
                }
              ]
            },
            {
              id: 'step2',
              nqName: 'LibTests:TestMul2',
              tags: ['taggedStep'],
            },
          ]
        },
        {
          id: 'step3',
          nqName: 'LibTests:TestMul2',
          tags: ['taggedStep'],

        },
      ],
      links: [{
        id: 'link1',
        base: 'base:expand(pipeline2|pipeline4)/expand(pipeline3|pipeline5)/#expand(taggedStep)',
        from: 'in1:#same(@base,taggedStep)/res',
        to: 'out1:#after(@base,taggedStep)/a',
        type: 'meta',
        handler() {}
      }, {
        id: 'link2',
        base: 'base:#expand(taggedPipeline)/expand(pipeline3|pipeline5)/#expand(taggedStep)',
        from: 'in1:#same(@base,taggedStep)/res',
        to: 'out1:#after(@base,taggedStep)/a',
        type: 'meta',
        handler() {}
      }]
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const matchInfos = links.map((link) => matchLink(tree, [], link));
    cleanMatchInfos(matchInfos);
    await snapshotCompare(matchInfos, 'Links selector and tags matching');
  });

  test('not(call):step3 snapshot equals not:step3', async () => {
    const make = (notSpec: string): PipelineConfiguration => ({
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
        {id: 'step3', nqName: 'LibTests:TestDiv2'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
        not: notSpec,
      }],
    });
    const treeA = StateTree.fromPipelineConfig({config: await getProcessedConfig(make('guard:step3'))});
    await treeA.init().toPromise();
    const treeB = StateTree.fromPipelineConfig({config: await getProcessedConfig(make('guard(call):step3'))});
    await treeB.init().toPromise();
    const linksA = (treeA.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const linksB = (treeB.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const a = linksA.map((link) => matchLink(treeA, [], link));
    const b = linksB.map((link) => matchLink(treeB, [], link));
    cleanMatchInfos(a);
    cleanMatchInfos(b);
    expectDeepEqual(a, b);
  });

  test('Reject (call) without (optional) in data link', async () => {
    let threw = false;
    try {
      await getProcessedConfig({
        id: 'pipeline1',
        type: 'static',
        steps: [
          {id: 'step1', nqName: 'LibTests:TestAdd2'},
          {id: 'step2', nqName: 'LibTests:TestMul2'},
        ],
        links: [{
          id: 'link1',
          from: 'fc(call):step1',
          to: 'out1:step2/a',
        }],
      });
    } catch {
      threw = true;
    }
    expectDeepEqual(threw, true);
  });

  test('Reject (call) on `to` in data link', async () => {
    let threw = false;
    try {
      await getProcessedConfig({
        id: 'pipeline1',
        type: 'static',
        steps: [
          {id: 'step1', nqName: 'LibTests:TestAdd2'},
          {id: 'step2', nqName: 'LibTests:TestMul2'},
        ],
        links: [{
          id: 'link1',
          from: 'in1:step1/res',
          to: 'out1(call):step2',
        }],
      });
    } catch {
      threw = true;
    }
    expectDeepEqual(threw, true);
  });
});

category('ComputeUtils: Driver link path shorthands', async () => {
  test('Self-ref forms in, in:, in:., in:./ parse to identical structure', async () => {
    const forms = ['in', 'in:', 'in:.', 'in:./'];
    const parsed = forms.map((f) => parseLinkIO(f, 'input'));
    for (let i = 1; i < parsed.length; i++)
      expectDeepEqual(parsed[i], parsed[0]);
    expectDeepEqual(parsed[0], [{name: 'in', segments: [], flags: []}]);
  });

  test('Self-ref equivalence holds across all ioTypes', async () => {
    const forms = ['in', 'in:', 'in:.', 'in:./'];
    for (const ioType of ['input', 'output', 'base', 'actions'] as const) {
      const parsed = forms.map((f) => parseLinkIO(f, ioType));
      for (let i = 1; i < parsed.length; i++)
        expectDeepEqual(parsed[i], parsed[0]);
    }
  });

  test('./path prefix parses identically to bare path', async () => {
    expectDeepEqual(parseLinkIO('in1:./step1/res', 'input'), parseLinkIO('in1:step1/res', 'input'));
    expectDeepEqual(parseLinkIO('out1:./step2/a', 'output'), parseLinkIO('out1:step2/a', 'output'));
    expectDeepEqual(
      parseLinkIO('in1:./pipeline1/step1/res', 'input'),
      parseLinkIO('in1:pipeline1/step1/res', 'input'),
    );
  });

  test('./ prefix works with flags and selectors', async () => {
    expectDeepEqual(
      parseLinkIO('in1(optional):./pipelinePar/stepMul/res', 'input'),
      parseLinkIO('in1(optional):pipelinePar/stepMul/res', 'input'),
    );
    expectDeepEqual(
      parseLinkIO('out1:./pipelinePar/all(stepAdd)/a', 'output'),
      parseLinkIO('out1:pipelinePar/all(stepAdd)/a', 'output'),
    );
  });

  test('End-to-end: ./path link matches identically to bare path link', async () => {
    const buildConfig = (fromSpec: string): PipelineConfiguration => ({
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{id: 'link1', from: fromSpec, to: 'out1:step2/a'}],
    });
    const runMatch = async (spec: string) => {
      const pconf = await getProcessedConfig(buildConfig(spec));
      const tree = StateTree.fromPipelineConfig({config: pconf});
      await tree.init().toPromise();
      const links = (tree.nodeTree.root.getItem().config as PipelineConfigurationStaticProcessed).links!;
      const info = matchLink(tree, [], links[0]);
      cleanMatchInfos([info]);
      return info;
    };
    expectDeepEqual(await runMatch('in1:./step1/res'), await runMatch('in1:step1/res'));
  });

  test('End-to-end: self-ref forms produce identical pipelineValidator links', async () => {
    const buildConfig = (toSpec: string): PipelineConfiguration => ({
      id: 'pipeline1',
      type: 'static',
      steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
      links: [{
        id: 'check',
        type: 'pipelineValidator',
        from: 'a:step1/a',
        to: toSpec,
        handler() {},
      }],
    });
    const forms = ['check', 'check:', 'check:.', 'check:./'];
    const processed = await Promise.all(forms.map((f) => getProcessedConfig(buildConfig(f))));
    for (let i = 1; i < processed.length; i++) {
      const a = (processed[i] as PipelineConfigurationStaticProcessed).links![0];
      const b = (processed[0] as PipelineConfigurationStaticProcessed).links![0];
      expectDeepEqual(a.to, b.to);
    }
  });

  const expectThrow = (spec: string, ioType: 'input' | 'output' = 'input') => {
    try {
      parseLinkIO(spec, ioType);
    } catch (e) {
      return;
    }
    throw new Error(`Expected parseLinkIO(${JSON.stringify(spec)}) to throw`);
  };

  test('Reject leading slash in:/path', async () => {
    expectThrow('in1:/step1/res');
  });

  test('Reject double dot in:..', async () => {
    expectThrow('in1:..');
  });

  test('Reject mid-path dot foo/./bar', async () => {
    expectThrow('in1:step1/./res');
  });

  test('Reject trailing dot foo/.', async () => {
    expectThrow('in1:step1/.');
  });

  test('Reject dot-prefix identifier in:.foo', async () => {
    expectThrow('in1:.foo');
  });

  test('Reject in:./. and in:./..', async () => {
    expectThrow('in1:./.');
    expectThrow('in1:./..');
  });

  test('Context-level constraint: self-ref forms still rejected where bare in is rejected', async () => {
    // In ioType='input'/'output' contexts, a path that ends with a non-`first` selector throws.
    // Bare self-ref skips that check today. The new alias forms must skip it identically (not bypass any *added* checks).
    // We verify equivalence at parse time across all relevant ioTypes — if a downstream check rejects bare 'in',
    // identical AST guarantees the new forms are rejected the same way.
    const baseForms = ['in', 'in:', 'in:.', 'in:./'];
    for (const ioType of ['input', 'output', 'base', 'actions'] as const) {
      const parsed = baseForms.map((f) => parseLinkIO(f, ioType));
      // every form yields one parsed entry with name='in', empty segments, empty flags
      for (const p of parsed) {
        expectDeepEqual(p.length, 1);
        expectDeepEqual(p[0].name, 'in');
        expectDeepEqual(p[0].segments, []);
      }
    }
  });

  test('same(@base) parses with empty ids', async () => {
    const [parsed] = parseLinkIO('in1:same(@base)/step1/res', 'input');
    expectDeepEqual(parsed.segments[0], {type: 'selector', selector: 'same', ids: [], stopIds: [], ref: 'base'});
  });

  test('#same(@base) parses with empty tags', async () => {
    const [parsed] = parseLinkIO('in1:#same(@base)/res', 'input');
    expectDeepEqual(parsed.segments[0], {type: 'tag', selector: 'same', tags: [], ref: 'base'});
  });

  test('same(@base, x|y) still parses with ids', async () => {
    const [parsed] = parseLinkIO('in1:same(@base,x|y)/step1/res', 'input');
    expectDeepEqual(parsed.segments[0], {type: 'selector', selector: 'same', ids: ['x', 'y'], stopIds: [], ref: 'base'});
  });

  test('Reject before(@base) without ids', async () => {
    expectThrow('in1:before(@base)/res');
  });

  test('Reject after(@base) without ids', async () => {
    expectThrow('in1:after(@base)/res');
  });

  test('Reject before*(@base) without ids', async () => {
    expectThrow('in1:before*(@base)/res');
  });

  test('Reject after+(@base) without ids', async () => {
    expectThrow('in1:after+(@base)/res');
  });

  test('Reject #before(@base) without tags', async () => {
    expectThrow('in1:#before(@base)/res');
  });

  test('Reject (call,template) combination', async () => {
    expectThrow('fc(call,template):step1|step2');
  });

  test('Reject (call) in base query', async () => {
    try {
      parseLinkIO('base(call):expand(step1)', 'base');
    } catch {
      return;
    }
    throw new Error('Expected parseLinkIO to throw on (call) in base');
  });
});

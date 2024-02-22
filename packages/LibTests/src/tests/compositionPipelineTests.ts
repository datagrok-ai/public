import * as DG from 'datagrok-api/dg';
import {serialize, deserialize, applyTransformations} from '@datagrok-libraries/utils/src/json-serialization';
import {category, test, before, delay} from '@datagrok-libraries/utils/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {CompositionPipeline, PipelineConfiguration} from '@datagrok-libraries/compute-utils/composition-pipeline/src/composition-pipeline';

function pickPipelineConfData(obj: any) {
  const res = {};
  Object.assign(res, ...['config', 'ioInfo', 'nodes', 'links', 'hooks', 'steps'].map((key) => ({[key]: obj[key]})));
  return res;
}

category('CompositionPipeline single config', async () => {
  before(async () => {});

  test('Simple config', async () => {
    const config: PipelineConfiguration = {
      id: 'testPipeline',
      nqName: 'LibTests:MockWrapper',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:AddMock',
        },
        {
          id: 'step2',
          nqName: 'LibTests:MulMock',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'a'],
        to: ['step2', 'a'],
      }]
    };
    const pipeline = new CompositionPipeline(config);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('Complex Config', async () => {
    const config: PipelineConfiguration = {
      id: 'testPipeline',
      nqName: 'LibTests:MockWrapper',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:AddMock',
          states: [{id: 'state11'}, {id: 'state12'}],
          actions: [{
            id: 'action1',
            from: [['step1', 'state11'], ['step1', 'state12']],
            to: ['step2', 'state21'],
            handler: 'LibTests:MockHandler1',
            position: 'buttons',
          }],
        },
        {
          id: 'step2',
          nqName: 'LibTests:MulMock',
          states: [{id: 'state21'}, {id: 'state22'}],
          popups: [{
            id: 'popup1',
            nqName: 'LibTests:TestFn1',
            position: 'menu',
          }],
        },
      ],
      hooks: {
        beforeFuncCallReady: [
          {
            id: 'hook1',
            handler: 'LibTests:MockHook1',
            from: ['step1', 'a'],
            to: ['step1', 'b'],
          },
        ],
      },
      states: [{id: 'state1'}, {id: 'state2'}],
      links: [
        {
          id: 'link1',
          from: ['step2', 'a'],
          to: ['step2', 'popup1', 'a'],
        },
        {
          id: 'link2',
          from: ['step2', 'b'],
          to: ['step2', 'popup1', 'b'],
        },
      ],
    };
    const pipeline = new CompositionPipeline(config);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });
});

category('CompositionPipeline composition config', async () => {
  test('2 simple configs', async () => {
    const config1: PipelineConfiguration = {
      id: 'testPipeline1',
      nqName: 'LibTests:MockWrapper',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:AddMock',
        },
        {
          id: 'step2',
          nqName: 'LibTests:MulMock',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'a'],
        to: ['step2', 'a'],
      }]
    };
    const config2: PipelineConfiguration = {
      id: 'testPipeline2',
      nqName: 'LibTests:MockWrapper',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestFn1',
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestFn2',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'a'],
        to: ['step2', 'a'],
      }]

    };
    const composedConfig = CompositionPipeline.compose([config1], config2);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });
});

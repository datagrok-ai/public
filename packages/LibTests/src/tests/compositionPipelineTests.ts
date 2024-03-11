import * as DG from 'datagrok-api/dg';
import {serialize, deserialize, applyTransformations} from '@datagrok-libraries/utils/src/json-serialization';
import {category, test, before, delay} from '@datagrok-libraries/utils/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import { CompositionPipeline, PipelineCompositionConfiguration, PipelineConfiguration} from '@datagrok-libraries/compute-utils/composition-pipeline/src/composition-pipeline';

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
            friendlyName: 'action1',
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
            friendlyName: 'popup1',
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
  const sconfig1: PipelineConfiguration = {
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
  const sconfig2: PipelineConfiguration = {
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

  const sconfig3: PipelineConfiguration = {
    id: 'testPipeline3',
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

  const conf1: PipelineConfiguration = {
    id: 'testPipeline1',
    nqName: 'LibTests:MockWrapper',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestFn1',
        states: [{id: 'state11'}, {id: 'state12'}],
        actions: [{
          id: 'action1',
          friendlyName: 'action1',
          from: [['step1', 'a'], ['step1', 'b']],
          to: ['step2', 'b'],
          handler: 'LibTests:MockHandler1',
          position: 'buttons',
        }],
        popups: [{
          id: 'popup1',
          friendlyName: 'popup1',
          nqName: 'LibTests:MockPopup1',
          position: 'menu',
          states: [{id: 'state111'}, {id: 'state112'}],
          actions: [{
            id: 'action11',
            friendlyName: 'action11',
            from: [['step1', 'popup1', 'a'], ['step1', 'popup1', 'b']],
            to: ['step2', 'b'],
            handler: 'LibTests:MockHandler1',
            position: 'buttons',
          }],
        }],
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestFn2',
      },
    ],
    states: [{id: 'state1'}, {id: 'state2'}],
    links: [{
      id: 'link1',
      from: ['step1', 'a'],
      to: ['step2', 'a'],
    }],
    hooks: {
      beforeFuncCallReady: [
        {
          id: 'hook1',
          handler: 'LibTests:MockHook1',
          from: ['step2', 'a'],
          to: ['step2', 'b'],
        },
      ],
    },
  };

  const conf2: PipelineConfiguration = {
    id: 'testPipeline2',
    nqName: 'LibTests:MockWrapper',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestFn3',
        states: [{id: 'state11'}, {id: 'state12'}],
        actions: [{
          id: 'action1',
          friendlyName: 'action1',
          from: [['step1', 'a'], ['step1', 'b']],
          to: ['step2', 'b'],
          handler: 'LibTests:MockHandler1',
          position: 'buttons',
        }],
        popups: [{
          id: 'popup1',
          friendlyName: 'popup1',
          nqName: 'LibTests:MockPopup1',
          position: 'menu',
          states: [{id: 'state111'}, {id: 'state112'}],
          actions: [{
            id: 'action11',
            friendlyName: 'action11',
            from: [['step1', 'popup1', 'a'], ['step1', 'popup1', 'b']],
            to: ['step2', 'b'],
            handler: 'LibTests:MockHandler1',
            position: 'buttons',
          }],
        }],
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestFn4',
      },
    ],
    states: [{id: 'state1'}, {id: 'state2'}],
    links: [{
      id: 'link1',
      from: ['step1', 'a'],
      to: ['step2', 'a'],
    }],
    hooks: {
      beforeFuncCallReady: [
        {
          id: 'hook1',
          handler: 'LibTests:MockHook1',
          from: ['step2', 'a'],
          to: ['step2', 'b'],
        },
      ],
    },
  };

  const conf3: PipelineConfiguration = {
    id: 'testPipeline3',
    nqName: 'LibTests:MockWrapper',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestFn5',
        states: [{id: 'state11'}, {id: 'state12'}],
        actions: [{
          id: 'action1',
          friendlyName: 'action1',
          from: [['step1', 'a'], ['step1', 'b']],
          to: ['step2', 'b'],
          handler: 'LibTests:MockHandler1',
          position: 'buttons',
        }],
        popups: [{
          id: 'popup1',
          friendlyName: 'popup1',
          nqName: 'LibTests:MockPopup1',
          position: 'menu',
          states: [{id: 'state111'}, {id: 'state112'}],
          actions: [{
            id: 'action11',
            friendlyName: 'action11',
            from: [['step1', 'popup1', 'a'], ['step1', 'popup1', 'b']],
            to: ['step2', 'b'],
            handler: 'LibTests:MockHandler1',
            position: 'buttons',
          }],
        }],
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestFn6',
      },
    ],
    states: [{id: 'state1'}, {id: 'state2'}],
    links: [{
      id: 'link1',
      from: ['step1', 'a'],
      to: ['step2', 'a'],
    }],
    hooks: {
      beforeFuncCallReady: [
        {
          id: 'hook1',
          handler: 'LibTests:MockHook1',
          from: ['step2', 'a'],
          to: ['step2', 'b'],
        },
      ],
    },
  };

  test('2 configs simple', async () => {
    const composedConfig = CompositionPipeline.compose(sconfig2, [sconfig1]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs simple configs', async () => {
    const composedConfig = CompositionPipeline.compose(sconfig3, [sconfig1, sconfig2]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs simple nested composition',  async () => {
    const composedConfig = CompositionPipeline.compose(sconfig3, [CompositionPipeline.compose(sconfig2, [sconfig1])]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('2 configs', async () => {
    const composedConfig = CompositionPipeline.compose(conf2, [conf1]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs', async () => {
    const composedConfig = CompositionPipeline.compose(conf3, [conf1, conf2]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs nested composition',  async () => {
    const composedConfig = CompositionPipeline.compose(conf3, [CompositionPipeline.compose(conf2, [conf1])]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs remove items', async () => {
    const cconf: PipelineCompositionConfiguration = {
      ...conf3,
      itemsToRemove: [
        ['testPipeline1', 'step1', 'popup1', 'action11'],
        ['testPipeline1', 'step2'],
        ['testPipeline2', 'step1', 'popup1'],
        ['testPipeline2', 'link1']
      ]
    };
    const composedConfig = CompositionPipeline.compose(cconf, [conf1, conf2]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs nested remove items', async () => {
    const cconf: PipelineCompositionConfiguration = {
      ...conf3,
      itemsToRemove: [
        ['testPipeline2', 'testPipeline1', 'step1', 'popup1', 'action11'],
        ['testPipeline2', 'testPipeline1', 'step2'],
        ['testPipeline2', 'step1', 'popup1'],
        ['testPipeline2', 'link1']
      ]
    };
    const composedConfig = CompositionPipeline.compose(cconf, [CompositionPipeline.compose(conf2, [conf1])]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs add items', async () => {
    const cconf: PipelineCompositionConfiguration = {
      ...conf3,
      actionsToAdd: [
        [{
          id: 'addedAction1',
          friendlyName: 'addedAction1',
          handler: 'LibTests:MockHandler2',
          position: 'buttons',
        }, ['step2']],
        [{
          id: 'addedAction2',
          friendlyName: 'addedAction2',
          handler: 'LibTests:MockHandler2',
          position: 'menu',
        }, ['testPipeline2', 'step2', 'addedPopup1']]
      ],
      popupsToAdd: [
        [{
          id: 'addedPopup1',
          friendlyName: 'addedPopup1',
          nqName: 'LibTests:MockPopup2',
          position: 'buttons',
        }, ['testPipeline2', 'step2']]
      ]
    };
    const composedConfig = CompositionPipeline.compose(cconf, [conf1, conf2]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs nested add items', async () => {
    const cconf3: PipelineCompositionConfiguration = {
      ...conf3,
      actionsToAdd: [
        [{
          id: 'addedAction1',
          friendlyName: 'addedAction1',
          handler: 'LibTests:MockHandler2',
          position: 'buttons',
        }, ['step2']],
        [{
          id: 'addedAction2',
          friendlyName: 'addedAction2',
          handler: 'LibTests:MockHandler2',
          position: 'menu',
        }, ['testPipeline2', 'testPipeline1', 'step2', 'addedPopup1']]
      ],
    };
    const cconf2: PipelineCompositionConfiguration = {
      ...conf2,
      popupsToAdd: [
        [{
          id: 'addedPopup1',
          friendlyName: 'addedPopup2',
          nqName: 'LibTests:MockPopup2',
          position: 'buttons',
        }, ['testPipeline1', 'step2']]
      ]
    };
    const composedConfig = CompositionPipeline.compose(cconf3, [CompositionPipeline.compose(cconf2, [conf1])]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('2 configs add steps', async () => {
    const cconf: PipelineCompositionConfiguration = {
      ...sconfig2,
      stepsToAdd: [[
        {
          id: 'addedStep1',
          nqName: 'LibTests:MockPopup2'
        },
        ['testPipeline1', 'step1'],
      ]]
    };
    const composedConfig = CompositionPipeline.compose(cconf, [sconfig1]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('2 configs add/remove steps', async () => {
    const cconf: PipelineCompositionConfiguration = {
      ...sconfig2,
      stepsToAdd: [[
        {
          id: 'addedStep1',
          nqName: 'LibTests:MockPopup2'
        },
        ['testPipeline1', 'step2'],
      ]],
      itemsToRemove: [[
        'testPipeline1', 'step2'
      ]]
    };
    const composedConfig = CompositionPipeline.compose(cconf, [sconfig1]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

  test('3 configs nested add steps',  async () => {
    const cconf3: PipelineCompositionConfiguration = {
      ...conf3,
      stepsToAdd: [[
        {
          id: 'addedStep1',
          nqName: 'LibTests:MockPopup2'
        },
        ['testPipeline2', 'step2'],
      ]],
    };
    const cconf2: PipelineCompositionConfiguration = {
      ...conf2,
      stepsToAdd: [[
        {
          id: 'addedStep2',
          nqName: 'LibTests:MockPopup2'
        },
        ['testPipeline1', 'step2'],
      ]],
    };
    const composedConfig = CompositionPipeline.compose(cconf3, [CompositionPipeline.compose(cconf2, [conf1])]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
  });

});

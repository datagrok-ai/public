import * as DG from 'datagrok-api/dg';
import {serialize, applyTransformations} from '@datagrok-libraries/utils/src/json-serialization';
import {category, test, before, delay} from '@datagrok-libraries/utils/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {CompositionPipeline, PipelineCompositionConfiguration, PipelineConfiguration} from '@datagrok-libraries/compute-utils/composition-pipeline/src/composition-pipeline';
import cloneDeepWith from 'lodash.clonedeepwith';
import {BehaviorSubject, Observable, Subject} from 'rxjs';
import simpleConfData from './snapshots/simple-conf-data.json';
import complexConfData from './snapshots/complex-conf-data.json';
import twoConfsSimple from './snapshots/2-confs-simple.json';
import threeConfsSimple from './snapshots/3-confs-simple.json';
import threeConfsSimpleNested from './snapshots/3-confs-simple-nested.json';
import twoConfs from './snapshots/2-confs.json';
import threeConfs from './snapshots/3-confs.json';
import threeConfsNested from './snapshots/3-confs-nested.json';
import threeConfsRemove from './snapshots/3-confs-remove-items.json';
import threeConfsRemoveNested from './snapshots/3-confs-remove-nested-items.json';
import threeConfsAdd from './snapshots/3-confs-add-items.json';
import threeConfsAddNested from './snapshots/3-confs-add-nested-items.json';
import twoConfsAddSteps from './snapshots/2-confs-add-steps.json';
import twoConfsAddRemoveSteps from './snapshots/2-confs-add-remove-steps.json';
import threeConfsAddNestedSteps from './snapshots/3-confs-nested-add-steps.json';

export function removeObservables<T>(config: T): T {
  return cloneDeepWith(config, (val) => {
    if (val instanceof Map) {
      const entries = removeObservables([...val.entries()]);
      return new Map(entries);
    }
    if ((val instanceof Subject) || (val instanceof BehaviorSubject) || (val instanceof Observable))
      return '$observable';
    if (val instanceof Function)
      return 'function';
  });
}

function pickPipelineConfData(obj: any) {
  const res = {};
  Object.assign(res, ...['config', 'ioInfo', 'nodes', 'links', 'hooks', 'steps'].map((key) => ({[key]: obj[key]})));
  // return serialize(removeObservables(res));
  return removeObservables(res);
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
    pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(simpleConfData));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(complexConfData));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(twoConfsSimple));
  });

  test('3 configs simple configs', async () => {
    const composedConfig = CompositionPipeline.compose(sconfig3, [sconfig1, sconfig2]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(threeConfsSimple));
  });

  test('3 configs simple nested composition',  async () => {
    const composedConfig = CompositionPipeline.compose(sconfig3, [CompositionPipeline.compose(sconfig2, [sconfig1])]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    console.log(actual);
    // expectDeepEqual(actual, applyTransformations(threeConfsSimpleNested));
  });

  test('2 configs', async () => {
    const composedConfig = CompositionPipeline.compose(conf2, [conf1]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(twoConfs));
  });

  test('3 configs', async () => {
    const composedConfig = CompositionPipeline.compose(conf3, [conf1, conf2]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(threeConfs));
  });

  test('3 configs nested composition',  async () => {
    const composedConfig = CompositionPipeline.compose(conf3, [CompositionPipeline.compose(conf2, [conf1])]);
    const pipeline = new CompositionPipeline(composedConfig);
    const _view = pipeline.makePipelineView();
    await pipeline.init();
    const actual = pickPipelineConfData(pipeline);
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(threeConfsNested));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(threeConfsRemove));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(threeConfsRemoveNested));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(threeConfsAdd));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(threeConfsAddNested));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(twoConfsAddSteps));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(twoConfsAddRemoveSteps));
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
    // console.log(actual);
    expectDeepEqual(actual, applyTransformations(threeConfsAddNestedSteps));
  });

});

category('CompositionPipeline reactivity', async () => {
  test('simple link',  async () => {
    const pipeline = new CompositionPipeline({
      id: 'testPipeline',
      nqName: 'LibTests:MockWrapper1',
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
        from: ['step1', 'res'],
        to: ['step2', 'a'],
      }]
    });
    const view = pipeline.makePipelineView();
    await pipeline.init();
    const rfv1 = view.getStepView('LibTests:AddMock');
    rfv1.setInput('a', 2);
    rfv1.setInput('b', 3);
    await delay(300);
    await rfv1.run();
    await delay(300);
    const rfv2 = view.getStepView('LibTests:MulMock');
    rfv2.setInput('b', 4);
    await delay(300);
    await rfv2.run();
    await delay(300);
    const a = rfv2.getParamValue('a');
    const res = rfv2.getParamValue('res');
    expectDeepEqual(a, 5);
    expectDeepEqual(res, 20);
  });

  test('composition links',  async () => {
    const conf1: PipelineCompositionConfiguration = {
      id: 'testPipeline1',
      nqName: 'LibTests:MockWrapper7',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:AddMock',
        },
      ],
      stepsToAdd: [[
        {
          id: 'step2',
          nqName: 'LibTests:MulMock',
        },
        ['testPipeline2', 'step1']
      ]],
      links: [{
        id: 'link1',
        from: ['step1', 'res'],
        to: ['testPipeline2', 'step1', 'a'],
      }, {
        id: 'link2',
        from: ['testPipeline2', 'step1', 'res'],
        to: ['testPipeline2', 'step2', 'a'],
      }]
    };
    const conf2 = {
      id: 'testPipeline2',
      nqName: 'LibTests:MockWrapper7',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestFn1',
        },
      ]
    }
    const pipeline = new CompositionPipeline(CompositionPipeline.compose(conf1, [conf2]));
    const view = pipeline.makePipelineView();
    await pipeline.init();
    const rfv1 = view.getStepView('LibTests:AddMock');
    rfv1.setInput('a', 1);
    rfv1.setInput('b', 2);
    await delay(300);
    await rfv1.run();
    await delay(300);
    const rfv2 = view.getStepView('LibTests:TestFn1');
    rfv2.setInput('b', 4);
    rfv2.setInput('c', 5);
    await delay(300);
    await rfv2.run();
    await delay(300);
    const rfv3 = view.getStepView('LibTests:MulMock');
    rfv3.setInput('b', 10);
    await delay(300);
    await rfv3.run();
    await delay(300);
    const a = rfv3.getParamValue('a');
    const res = rfv3.getParamValue('res');
    expectDeepEqual(a, 12);
    expectDeepEqual(res, 120);
  });

});

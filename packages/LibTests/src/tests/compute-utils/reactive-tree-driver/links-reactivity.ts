import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {of} from 'rxjs';
import {delay, mapTo} from 'rxjs/operators';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import { makeValidationResult } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';

category('ComputeUtils: Driver links reactivity', async () => {
  let testScheduler: TestScheduler;

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
        nqName: 'LibTests:TestMul2',
      },
    ],
    links: [{
      id: 'link1',
      from: 'in1:step1/b',
      to: 'out1:step2/a',
    }],
  };

  const config2: PipelineConfiguration = {
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
      from: 'in1:step1/a',
      to: 'out1:step1/a',
      isValidator: true,
      handler({controller}) {
        controller.setValidation('out1', makeValidationResult({warnings: ['test warn']}));
        return;
      },
    }],
  };

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      expectDeepEqual(actual, expected);
      // console.log(actual, expected);
    });
  });

  test('Run default handler on trigger', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 1});
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
        link.trigger();
      });
    });
  });

  test('Run enabled default handler', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 1});
      cold('-a').subscribe(() => {
        link.setActive();
        inNode.getItem().getStateStore().setState('b', 1);
      });
    });
  });

  test('Dont run disabled handlers', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a', {a: undefined});
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
      });
    });
  });

  test('Multiple changes buffering', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 2});
      cold('-a').subscribe(() => {
        link.setActive();
        inNode.getItem().getStateStore().setState('b', 1);
        inNode.getItem().getStateStore().setState('b', 2);
      });
    });
  });

  test('Run custom handlers', async () => {
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
        handler({controller}) {
          const in1 = controller.getFirst('in1');
          controller.setAll('out1', in1 + 1);
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 2});
      cold('-a').subscribe(() => {
        link.setActive();
        inNode.getItem().getStateStore().setState('b', 1);
      });
    });
  });

  test('Run custom async handlers', async () => {
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
        handler({controller}) {
          const in1 = controller.getFirst('in1');
          controller.setAll('out1', in1 + 1);
          return of(null).pipe(delay(100), mapTo(void(0)));
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a 100ms b', {a: undefined, b: 2});
      cold('-a').subscribe(() => {
        link.setActive();
        inNode.getItem().getStateStore().setState('b', 1);
      });
    });
  });

  test('Properly cancel async handlers', async () => {
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
        handler({controller}) {
          const in1 = controller.getFirst('in1');
          controller.setAll('out1', in1 + 1);
          return of(null).pipe(delay(100), mapTo(void(0)));
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a 149ms b', {a: undefined, b: 3});
      cold('-a').subscribe(() => {
        link.setActive();
        inNode.getItem().getStateStore().setState('b', 1);
      });
      cold('50ms a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 2);
      });
    });
  });


  test('Run validators with debounce', async () => {
    const pconf = await getProcessedConfig(config2);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      link.wire(tree);
      expectObservable(link.isRunning$, '^ 1000ms !').toBe('a 250ms (bc)', {a: false, b: true, c: false});
      cold('-a').subscribe(() => {
        link.setActive();
        inNode.getItem().getStateStore().setState('a', 1);
      });

      cold('252ms a').subscribe(() => {
        const validators = (inNode.getItem().getStateStore() as FuncCallInstancesBridge).validations$.value;
        const expected = {
          'a': {
            'warnings': [
              {
                'description': 'test warn',
              },
            ],
          },
        };
        expectDeepEqual(validators[link.uuid], expected);
      });
    });
  });

  test('Run validators on trigger', async () => {
    const pconf = await getProcessedConfig(config2);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      link.wire(tree);
      expectObservable(link.isRunning$, '^ 1000ms !').toBe('a (bc)', {a: false, b: true, c: false});
      cold('-a').subscribe(() => {
        link.trigger();
      });

      cold('--a').subscribe(() => {
        const validators = (inNode.getItem().getStateStore() as FuncCallInstancesBridge).validations$.value;
        const expected = {
          'a': {
            'warnings': [
              {
                'description': 'test warn',
              },
            ],
          },
        };
        expectDeepEqual(validators[link.uuid], expected);
      });
    });
  });

  test('Propagate restriction info in a default handler', async () => {
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
        defaultRestrictions: {
          out1: 'restricted',
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 1});
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
        link.trigger();
      });
      cold('--a').subscribe(() => {
        const restrictions = (outNode.getItem().getStateStore() as FuncCallInstancesBridge).inputRestrictions$.value;
        const expected = {
          'a': {
            'type': 'restricted',
            'assignedValue': 1,
          },
        };
        expectDeepEqual(restrictions, expected);
      });
    });
  });

  test('Propagate restriction info in a custom handler', async () => {
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
        handler({controller}) {
          const in1 = controller.getFirst('in1');
          controller.setAll('out1', in1 + 1, 'restricted');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 2});
      cold('-a').subscribe(() => {
        link.setActive();
        inNode.getItem().getStateStore().setState('b', 1);
      });
      cold('--a').subscribe(() => {
        const restrictions = (outNode.getItem().getStateStore() as FuncCallInstancesBridge).inputRestrictions$.value;
        const expected = {
          'a': {
            'type': 'restricted',
            'assignedValue': 2,
          },
        };
        expectDeepEqual(restrictions, expected);
      });
    });
  });

  test('Links propagate multiple validators results', async () => {
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
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        isValidator: true,
        handler({controller}) {
          controller.setValidation('out1', makeValidationResult({warnings: ['some warn']}));
          return;
        },
      }, {
        id: 'link2',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        isValidator: true,
        handler({controller}) {
          controller.setValidation('out1', makeValidationResult({warnings: ['another warn']}));
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link1, link2] = ls.createAutoLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      link1.wire(tree);
      link2.wire(tree);
      expectObservable(link1.isRunning$, '^ 1000ms !').toBe('a (bc)', {a: false, b: true, c: false});
      expectObservable(link2.isRunning$, '^ 1000ms !').toBe('a (bc)', {a: false, b: true, c: false});
      cold('-a').subscribe(() => {
        link1.trigger();
        link2.trigger();
      });

      cold('--a').subscribe(() => {
        const validators = (inNode.getItem().getStateStore() as FuncCallInstancesBridge).validations$.value;
        const expected = {
          [link1.uuid]: {
            a: {'warnings': [
              {
                'description': 'some warn',
              }],
            }},
          [link2.uuid]: {
            a: {'warnings': [
              {
                'description': 'another warn',
              }],
            },
          },
        };
        expectDeepEqual(validators, expected);
      });
    });
  });


});

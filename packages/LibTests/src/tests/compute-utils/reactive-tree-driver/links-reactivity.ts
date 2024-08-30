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

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      expectDeepEqual(actual, expected);
    });
  });

  test('Run default handler on trigger', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      const outNode = tree.getNode([{idx: 1}]);
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
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      const outNode = tree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 1});
      cold('-a').subscribe(() => {
        link.enable();
        inNode.getItem().getStateStore().setState('b', 1);
      });
    });
  });

  test('Dont run disabled handlers', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      const outNode = tree.getNode([{idx: 1}]);
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
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      const outNode = tree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 2});
      cold('-a').subscribe(() => {
        link.enable();
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
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      const outNode = tree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 2});
      cold('-a').subscribe(() => {
        link.enable();
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
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      const outNode = tree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a 100ms b', {a: undefined, b: 2});
      cold('-a').subscribe(() => {
        link.enable();
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
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      const outNode = tree.getNode([{idx: 1}]);
      link.wire(tree);
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a 149ms b', {a: undefined, b: 3});
      cold('-a').subscribe(() => {
        link.enable();
        inNode.getItem().getStateStore().setState('b', 1);
      });
      cold('50ms a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 2);
      });
    });
  });


});

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {createTestScheduler} from '../../../test-utils';


category('ComputeUtils: Driver links reactivity: meta', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Propagate meta info', async () => {
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
        type: 'meta',
        handler({controller}) {
          controller.setViewMeta('out1', {key: 'val'});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
      });
      expectObservable((outNode.getItem().getStateStore() as FuncCallInstancesBridge).meta.a
      ).toBe('ab', {
        a: undefined,
        b: {
          'key': 'val',
        },
      });
    });
  });

  test('Propagate multiple meta info', async () => {
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
        type: 'meta',
        handler({controller}) {
          controller.setViewMeta('out1', {key1: 'val1'});
        },
      }, {
        id: 'link2',
        from: 'in1:step1/b',
        to: 'out1:step2/a',
        type: 'meta',
        handler({controller}) {
          controller.setViewMeta('out1', {key2: 'val1'});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
      });
      expectObservable((outNode.getItem().getStateStore() as FuncCallInstancesBridge).meta.a
      ).toBe('a(bc)', {
        a: undefined,
        b: {
          'key1': 'val1',
        },
        c: {
          'key1': 'val1',
          'key2': 'val1',
        },
      });
    });
  });

  test('Run no inputs meta on init', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      links: [{
        id: 'link1',
        from: [],
        to: 'out1:step1/a',
        type: 'meta',
        handler({controller}) {
          controller.setViewMeta('out1', {key: 'val'});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);
    testScheduler.run((helpers) => {
      const {expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, defaultValidators: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]);
      expectObservable((node.getItem().getStateStore() as FuncCallInstancesBridge).meta.a
      ).toBe('a', {
        a: {
          'key': 'val',
        },
      });
    });
  });
});

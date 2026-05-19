import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig, PipelineConfigurationStaticProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {parseLinkIO} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/LinkSpec';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';

function expectThrows(fn: () => unknown, match?: RegExp) {
  let threw = false;
  let err: unknown = undefined;
  try {
    fn();
  } catch (e) {
    threw = true;
    err = e;
  }
  expectDeepEqual(threw, true);
  if (match) expectDeepEqual(match.test(String((err as Error)?.message ?? err)), true);
}

async function expectThrowsAsync(fn: () => Promise<unknown>, match?: RegExp) {
  let threw = false;
  let err: unknown = undefined;
  try {
    await fn();
  } catch (e) {
    threw = true;
    err = e;
  }
  expectDeepEqual(threw, true);
  if (match) expectDeepEqual(match.test(String((err as Error)?.message ?? err)), true);
}

category('ComputeUtils: Driver template expansion grammar', async () => {
  test('outputs(nqName) parses with ioExpand + nqName', async () => {
    const [parsed] = parseLinkIO('in_(template):stepA/outputs(LibTests:TestAdd2)', 'input');
    const last = parsed.segments[parsed.segments.length - 1] as any;
    expectDeepEqual(last.ioExpand, 'outputs');
    expectDeepEqual(last.nqName, 'LibTests:TestAdd2');
    expectDeepEqual(last.excludeIds, []);
  });

  test('inputs(nqName, excludes) parses both args', async () => {
    const [parsed] = parseLinkIO('out_(template):stepB/inputs(LibTests:TestMul2, a)', 'output');
    const last = parsed.segments[parsed.segments.length - 1] as any;
    expectDeepEqual(last.ioExpand, 'inputs');
    expectDeepEqual(last.nqName, 'LibTests:TestMul2');
    expectDeepEqual(last.excludeIds, ['a']);
  });

  test('Pipe-separated excludes parsed', async () => {
    const [parsed] = parseLinkIO('in_(template):stepA/outputs(LibTests:TestAdd2, a|b|c)', 'input');
    const last = parsed.segments[parsed.segments.length - 1] as any;
    expectDeepEqual(last.excludeIds, ['a', 'b', 'c']);
  });

  test('No direction validation: outputs() on to and inputs() on from parse', async () => {
    parseLinkIO('out_(template):stepB/outputs(LibTests:TestAdd2)', 'output');
    parseLinkIO('in_(template):stepA/inputs(LibTests:TestAdd2)', 'input');
  });

  test('Reject IO selector without nqName arg', async () => {
    expectThrows(() => parseLinkIO('in_(template):stepA/outputs()', 'input'));
  });

  test('Reject single-arg without colon (not a valid nqName)', async () => {
    expectThrows(() => parseLinkIO('in_(template):stepA/outputs(stepA)', 'input'));
    expectThrows(() => parseLinkIO('in_(template):stepA/outputs(a|b)', 'input'));
  });

  test('Reject IO selector without (template) flag', async () => {
    expectThrows(() => parseLinkIO('in_:stepA/outputs(LibTests:TestAdd2)', 'input'));
  });

  test('Reject IO selector outside last segment', async () => {
    expectThrows(() => parseLinkIO('in_(template):outputs(LibTests:TestAdd2)/res', 'input'));
  });

  test('Reject IO selector in not: / base: / actions: queries', async () => {
    // not: queries match nodes; they have no IO semantics. Deferred expansion
    // only runs for from/to, so a deferred selector in not: would silently
    // survive to the runtime — reject at parse time.
    expectThrows(() => parseLinkIO('in_(template):outputs(LibTests:TestAdd2)', 'not'));
    expectThrows(() => parseLinkIO('in_(template):outputs(LibTests:TestAdd2)', 'base'));
    expectThrows(() => parseLinkIO('in_(template):outputs(LibTests:TestAdd2)', 'actions'));
  });

  test('Reject (template, call) with IO selector', async () => {
    expectThrows(() => parseLinkIO('in_(template,call):stepA/outputs(LibTests:TestAdd2)', 'input'));
  });

  test('all(...) is rejected as deferred form (only inputs/outputs supported)', async () => {
    // all(...) is parsed as a regular node selector now; LibTests:TestAdd2 has a colon
    // which the node-selector grammar does not accept → parse error.
    expectThrows(() => parseLinkIO('in_(template):stepA/all(LibTests:TestAdd2)', 'input'));
  });

  test('Existing bare-list template still parses positionally', async () => {
    const parsed = parseLinkIO('in_(template):step1/a|b', 'input');
    expectDeepEqual(parsed.length, 2);
    expectDeepEqual(parsed[0].name, 'in_a');
    expectDeepEqual(parsed[1].name, 'in_b');
  });
});

category('ComputeUtils: Driver template expansion config-time', async () => {
  test('Expansion: outputs(nqName) replaces deferred entry with bare-list entries', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in_(template):step1/outputs(LibTests:TestAdd2)',
        to: 'out_(template):step2/inputs(LibTests:TestMul2)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const links = (pconf as PipelineConfigurationStaticProcessed).links!;
    // TestAdd2 outputs: ['res']; TestMul2 inputs: ['a', 'b']
    expectDeepEqual(links[0].from.map((io) => io.name), ['in_res']);
    expectDeepEqual(links[0].to.map((io) => io.name), ['out_a', 'out_b']);
    for (const io of [...links[0].from, ...links[0].to]) {
      const last = io.segments[io.segments.length - 1] as any;
      expectDeepEqual(last.ioExpand, undefined);
    }
  });

  test('Expansion respects exclusion', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMultiarg5'},
      ],
      links: [{
        id: 'link1',
        from: 'in_(template):step1/outputs(LibTests:TestAdd2)',
        to: 'out_(template):step2/inputs(LibTests:TestMultiarg5, c|d|e)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const links = (pconf as PipelineConfigurationStaticProcessed).links!;
    expectDeepEqual(links[0].to.map((io) => io.name), ['out_a', 'out_b']);
  });

  test('Error: unknown nqName', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
      links: [{
        id: 'link1',
        from: 'in_(template):step1/outputs(LibTests:DoesNotExist)',
        to: 'out_:step1/a',
      }],
    };
    await expectThrowsAsync(() => getProcessedConfig(config), /DoesNotExist/);
  });
});

category('ComputeUtils: Driver template expansion runtime', async () => {
  let testScheduler: TestScheduler;
  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Positional pairing transfers values when source/target IO orders match', async () => {
    // TestIONamesA outputs ['x', 'y']; TestIONamesB inputs ['x', 'y']. Same order;
    // default handler positional pairing → x→x.
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestIONamesA'},
        {id: 'step2', nqName: 'LibTests:TestIONamesB'},
      ],
      links: [{
        id: 'link1',
        from: 'in_(template):step1/outputs(LibTests:TestIONamesA)',
        to: 'out_(template):step2/inputs(LibTests:TestIONamesB)',
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
        inNode.getItem().getStateStore().setState('x', 11);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('x'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 11});
    });
  });

  test('Existing bare-list template still pairs by index (regression)', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in_:step1/b',
        to: 'out_:step2/a',
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
        inNode.getItem().getStateStore().setState('b', 5);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 5});
    });
  });
});

void DG;

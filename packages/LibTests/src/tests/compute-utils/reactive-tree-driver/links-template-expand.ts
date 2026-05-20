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

  test('Marker: templateName is the prefix string for prefixed templates; absent on bare entries', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestIONamesA'},
        {id: 'step2', nqName: 'LibTests:TestIONamesB'},
      ],
      links: [{
        id: 'link1',
        from: ['in_(template):step1/outputs(LibTests:TestIONamesA)', 'bare_:step1/x'],
        to: ['out_(template):step2/inputs(LibTests:TestIONamesB)', 'plain_:step2/x'],
      }],
    };
    const pconf = await getProcessedConfig(config);
    const link = (pconf as PipelineConfigurationStaticProcessed).links![0];
    for (const io of link.from)
      expectDeepEqual((io as any).templateName, io.name.startsWith('in_') ? 'in_' : undefined);
    for (const io of link.to)
      expectDeepEqual((io as any).templateName, io.name.startsWith('out_') ? 'out_' : undefined);
  });

  test('Name matching propagates regardless of target IO declaration order', async () => {
    // TestIONamesA outputs [x, y]; TestIONamesBReversed inputs [y, x].
    // Positional pairing would miswire (x→y); name matching pairs x→x.
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestIONamesA'},
        {id: 'step2', nqName: 'LibTests:TestIONamesBReversed'},
      ],
      links: [{
        id: 'link1',
        from: 'in_(template):step1/outputs(LibTests:TestIONamesA)',
        to: 'out_(template):step2/inputs(LibTests:TestIONamesBReversed)',
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

  test('Unmatched script IOs inside a template are dropped', async () => {
    // TestIONamesAExtra outputs [x, y, z]; TestIONamesB inputs [x, y]. z is dropped.
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestIONamesAExtra'},
        {id: 'step2', nqName: 'LibTests:TestIONamesB'},
      ],
      links: [{
        id: 'link1',
        from: 'in_(template):step1/outputs(LibTests:TestIONamesAExtra)',
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
        inNode.getItem().getStateStore().setState('y', 22);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('y'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 22});
    });
  });

  test('Interleaved bare and template slots', async () => {
    // bare slot pair (a→a) plus template slot pair (x,y matched by name).
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestIONamesA'},
        {id: 'step2', nqName: 'LibTests:TestIONamesBReversed'},
      ],
      links: [{
        id: 'link1',
        from: ['seed_:step1/seed', 'in_(template):step1/outputs(LibTests:TestIONamesA)'],
        to: ['drop_:step2/y', 'out_(template):step2/inputs(LibTests:TestIONamesBReversed)'],
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
        inNode.getItem().getStateStore().setState('x', 7);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('x'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 7});
    });
  });

  test('No-op when no templates: bare entries still pair positionally', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in_:step1/res',
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
        inNode.getItem().getStateStore().setState('res', 13);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !')
        .toBe('a b', {a: undefined, b: 13});
    });
  });

  test('_(template) on source with bare on target: templateName is numeric 0', async () => {
    // _(template) produces link-io names equal to script IO ids verbatim (x, y) —
    // no prefix. The slot identity is the numeric index 0, distinct from any
    // user-typed string template name (grammar identifiers can never be pure digits).
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestIONamesA'},
        {id: 'step2', nqName: 'LibTests:TestIONamesB'},
      ],
      links: [{
        id: 'link1',
        from: '_(template):step1/outputs(LibTests:TestIONamesA)',
        to: 'plain_:step2/x',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const link = (pconf as PipelineConfigurationStaticProcessed).links![0];
    expectDeepEqual(link.from.map((io) => io.name), ['x', 'y']);
    for (const io of link.from)
      expectDeepEqual((io as any).templateName, 0);
    expectDeepEqual((link.to[0] as any).templateName, undefined);
  });

  test('Two _(template) on the same side index to 0, 1', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestIONamesA'},
        {id: 'step2', nqName: 'LibTests:TestIONamesB'},
      ],
      links: [{
        id: 'link1',
        from: [
          '_(template):step1/outputs(LibTests:TestIONamesA)',
          '_(template):step1/inputs(LibTests:TestIONamesA, seed)',
        ],
        to: 'out_(template):step2/inputs(LibTests:TestIONamesB)',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const link = (pconf as PipelineConfigurationStaticProcessed).links![0];
    // First anonymous template expands TestIONamesA's outputs (x, y); second
    // expands its inputs minus `seed` (so nothing). The counter increments
    // per anonymous operator regardless of how many entries it emits.
    const firstTplNames = link.from.filter((io: any) => io.templateName === 0).map((io) => io.name);
    expectDeepEqual(firstTplNames, ['x', 'y']);
    expectDeepEqual(link.from.every((io: any) => io.templateName !== 1), true);
  });

  test('Repeated non-_ prefix on the same side still throws (manual disambiguation is the right answer)', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestIONamesA'},
        {id: 'step2', nqName: 'LibTests:TestIONamesB'},
      ],
      links: [{
        id: 'link1',
        from: [
          'dup_(template):step1/outputs(LibTests:TestIONamesA)',
          'dup_(template):step1/inputs(LibTests:TestIONamesA, seed)',
        ],
        to: 'out_(template):step2/inputs(LibTests:TestIONamesB)',
      }],
    };
    await expectThrowsAsync(() => getProcessedConfig(config), /distinct prefixes/);
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

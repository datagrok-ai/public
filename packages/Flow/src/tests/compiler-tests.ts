import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {topologicalSort} from '../compiler/topological-sort';
import {emitScript} from '../compiler/script-emitter';
import {validateGraph} from '../compiler/validator';
import {makeEditor, destroyEditor, addNode} from './test-utils';

const SETTINGS = {name: 'TestFlow', description: 'test', tags: ['funcflow']};

category('Flow: topological sort', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('orders source before target', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      const sorted = topologicalSort(e.flow);
      expect(sorted.length, 2);
      expect(sorted.indexOf(input.id) < sorted.indexOf(output.id), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('chained pass-through preserves order', async () => {
    const e = makeEditor();
    try {
      const c = await addNode(e.flow, 'Constants/String');
      const t1 = await addNode(e.flow, 'Utilities/ToString');
      const t2 = await addNode(e.flow, 'Utilities/ToString');
      await e.flow.addConnectionByKeys(c.id, 'value', t1.id, 'value');
      await e.flow.addConnectionByKeys(t1.id, 'text', t2.id, 'value');
      const sorted = topologicalSort(e.flow);
      expect(sorted.indexOf(c.id) < sorted.indexOf(t1.id), true);
      expect(sorted.indexOf(t1.id) < sorted.indexOf(t2.id), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('disjoint subgraphs execute top-first, fully', async () => {
    const e = makeEditor();
    try {
      // Bottom chain created FIRST — canvas position must win over creation order.
      const bottomIn = await addNode(e.flow, 'Constants/String', 0, 500);
      const bottomOut = await addNode(e.flow, 'Utilities/ToString', 300, 500);
      await e.flow.addConnectionByKeys(bottomIn.id, 'value', bottomOut.id, 'value');
      const topIn = await addNode(e.flow, 'Constants/String', 0, 50);
      const topOut = await addNode(e.flow, 'Utilities/ToString', 300, 50);
      await e.flow.addConnectionByKeys(topIn.id, 'value', topOut.id, 'value');

      const sorted = topologicalSort(e.flow);
      const at = (id: string): number => sorted.indexOf(id);
      // The whole top component drains before the bottom one starts —
      // lower disjoint paths may implicitly consume what upper ones produced.
      expect(Math.max(at(topIn.id), at(topOut.id)) < Math.min(at(bottomIn.id), at(bottomOut.id)), true,
        'top chain must fully precede the bottom chain');
    } finally {
      destroyEditor(e);
    }
  });

  test('within a component, ready nodes process top-to-bottom', async () => {
    const e = makeEditor();
    try {
      // Two independent sources merging into one sink; the lower source was
      // created first, but the upper one must come first in the order.
      const lower = await addNode(e.flow, 'Constants/String', 0, 400);
      const upper = await addNode(e.flow, 'Constants/String', 0, 100);
      const sink = await addNode(e.flow, 'Comparisons/Equals (==)', 300, 250);
      await e.flow.addConnectionByKeys(upper.id, 'value', sink.id, 'a');
      await e.flow.addConnectionByKeys(lower.id, 'value', sink.id, 'b');

      const sorted = topologicalSort(e.flow);
      expect(sorted.indexOf(upper.id) < sorted.indexOf(lower.id), true, 'upper source first');
      expect(sorted.indexOf(sink.id), sorted.length - 1, 'sink last');
    } finally {
      destroyEditor(e);
    }
  });

  test('detects cycles', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Utilities/ToString');
      const b = await addNode(e.flow, 'Utilities/ToString');
      await e.flow.addConnectionByKeys(a.id, 'text', b.id, 'value');
      await e.flow.addConnectionByKeys(b.id, 'text', a.id, 'value');
      let threw = false;
      try {
        topologicalSort(e.flow);
      } catch {
        threw = true;
      }
      expect(threw, true);
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: script emitter', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('emits input/output headers and pass-through body', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      input.properties['paramName'] = 'myTable';
      const output = await addNode(e.flow, 'Outputs/Table Output');
      output.properties['paramName'] = 'res';
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('//input: dataframe myTable'), true, 'input header');
      expect(script.includes('//output: dataframe res'), true, 'output header');
      expect(script.includes('res = myTable;'), true, 'output assignment');
      expect(script.includes('//language: javascript'), true, 'language header');
    } finally {
      destroyEditor(e);
    }
  });

  test('constant string is inlined into Value Output', async () => {
    const e = makeEditor();
    try {
      const c = await addNode(e.flow, 'Constants/String');
      c.properties['value'] = 'hello';
      // Constant nodes title themselves after their value — emission must
      // dispatch on the registered type, not the (user-editable) label.
      c.label = 'const: hello';
      const out = await addNode(e.flow, 'Outputs/Value Output');
      out.properties['paramName'] = 'greeting';
      await e.flow.addConnectionByKeys(c.id, 'value', out.id, 'value');

      const script = emitScript(e.flow, SETTINGS);
      // The constant is declared as its own variable and the output references it:
      //   let constHello = "hello";  …  greeting = constHello;
      expect(script.includes('"hello"'), true, 'string literal present');
      expect(/greeting\s*=\s*\w+;/.test(script), true, 'output assigned from upstream variable');
    } finally {
      destroyEditor(e);
    }
  });

  test('instrumented mode wraps steps with run events', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');

      const script = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'run-123'});
      expect(script.includes('run-123'), true, 'run id embedded');
      expect(script.includes('run-complete'), true, 'completion event emitted');
    } finally {
      destroyEditor(e);
    }
  });

  test('side-effect-only utilities (Log) emit no output summary when instrumented', async () => {
    // Log/Info/Warning declare no variable — the instrumented wrapper used to
    // emit `__ff_summarize(log)` anyway, so any flow with a Log node failed
    // at run time with "log is not defined" (ReferenceError).
    const e = makeEditor();
    try {
      const c = await addNode(e.flow, 'Constants/String');
      c.properties['value'] = 'hello';
      const log = await addNode(e.flow, 'Utilities/Log');
      await e.flow.addConnectionByKeys(c.id, 'value', log.id, 'value');
      const toStr = await addNode(e.flow, 'Utilities/ToString');
      await e.flow.addConnectionByKeys(c.id, 'value', toStr.id, 'value');

      const script = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'run-log'});
      expect(script.includes('console.log('), true, 'the Log body is emitted');
      expect(script.includes('__ff_summarize(log'), false, 'no summary of an undeclared variable');
      expect(script.includes(`__ff_stash(${JSON.stringify(log.id)}`), false,
        'no stash of an undeclared variable (the compiler gives Log no phantom output)');
      expect(script.includes(`__ff_emit('node-complete', '${log.id}');`), true,
        'the Log node completes with no outputs payload');
      // A value-producing utility on the same canvas still summarizes normally.
      // In instrumented mode the `let` is hoisted out of the try block, so the
      // body line is a bare assignment.
      const toStringVar = script.match(/(\w+) = \(.*\)\.toString\(\);/)?.[1];
      expect(!!toStringVar, true, 'the ToString step declares its variable');
      expect(script.includes(`__ff_summarize(${toStringVar})`), true, 'ToString keeps its output summary');
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: validator', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('empty graph yields a warning', async () => {
    const e = makeEditor();
    try {
      const results = validateGraph(e.flow);
      expect(results.length >= 1, true);
      expect(results[0].severity, 'warning');
    } finally {
      destroyEditor(e);
    }
  });

  test('column input without table input is an error', async () => {
    const e = makeEditor();
    try {
      await addNode(e.flow, 'Inputs/Column Input');
      const results = validateGraph(e.flow);
      expect(results.some((r) => r.severity === 'error'), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('valid input → output graph has no errors', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      const results = validateGraph(e.flow);
      expect(results.some((r) => r.severity === 'error'), false);
    } finally {
      destroyEditor(e);
    }
  });
});

import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, createNode} from '../rete/node-factory';
import {EXEC_IN_KEY, EXEC_OUT_KEY} from '../rete/scheme';
import {areTypesCompatible} from '../types/type-map';
import {topologicalSort} from '../compiler/topological-sort';
import {emitScript} from '../compiler/script-emitter';
import {validateGraph} from '../compiler/validator';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

const SETTINGS = {name: 'OrderFlow', description: '', tags: ['funcflow']};

category('Flow: order edges', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('order sockets connect only to each other (isolated from data + wildcards)', async () => {
    expect(areTypesCompatible('order', 'order'), true);
    expect(areTypesCompatible('order', 'dataframe'), false);
    expect(areTypesCompatible('dataframe', 'order'), false);
    // The dynamic/object wildcards must NOT swallow order ports.
    expect(areTypesCompatible('order', 'dynamic'), false);
    expect(areTypesCompatible('dynamic', 'order'), false);
    expect(areTypesCompatible('order', 'object'), false);
  });

  test('every created node gets an exec-in and exec-out port of type order', async () => {
    for (const type of ['Utilities/Info', 'Inputs/Table Input', 'Outputs/Value Output', 'Constants/String']) {
      const node = createNode(type)!;
      expect(EXEC_IN_KEY in node.inputs, true, `${type} has exec-in`);
      expect(EXEC_OUT_KEY in node.outputs, true, `${type} has exec-out`);
      expect((node.inputs[EXEC_IN_KEY] as {socket: {dgType: string}}).socket.dgType, 'order');
      expect((node.outputs[EXEC_OUT_KEY] as {socket: {dgType: string}}).socket.dgType, 'order');
    }
  });

  test('exec squares are hover-only until wired: data-wired flips with an order edge', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Utilities/Info', 0, 0);
      const b = await addNode(e.flow, 'Utilities/Info', 300, 0);
      const row = (id: string): HTMLElement | null =>
        e.container.querySelector(`.ff-node[data-node-id="${id}"] .ff-node-exec-row`);
      await until(() => row(a.id) != null && row(b.id) != null);
      // Unwired → the row is marked not-wired (CSS keeps it invisible except
      // on node hover / during an order drag).
      expect(row(a.id)!.dataset.wired, 'false', 'unwired row hidden by default');
      // The tooltip explains the port in plain words, not just "order" — and
      // the inner socket dot must NOT carry its own title (a nested title
      // would shadow the wrapper's explanation with the bare type name).
      const port = e.container.querySelector(`.ff-node[data-node-id="${a.id}"] .ff-exec-out`);
      expect((port?.getAttribute('title') ?? '').includes('Run order'), true, 'plain-language tooltip');
      expect(port?.querySelector('.ff-socket')?.hasAttribute('title'), false, 'order dot has no shadowing title');

      await e.flow.addConnectionByKeys(a.id, EXEC_OUT_KEY, b.id, EXEC_IN_KEY);
      await e.flow.updateNode(a.id);
      await e.flow.updateNode(b.id);
      const wired = await until(() =>
        row(a.id)?.dataset.wired === 'true' && row(b.id)?.dataset.wired === 'true');
      expect(wired, true, 'both ends stay visible once an order edge exists');
    } finally {
      destroyEditor(e);
    }
  });

  test('an order drag lights every other node\'s opposite exec square and dims the rest', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Utilities/Info', 0, 0);
      const b = await addNode(e.flow, 'Utilities/Info', 300, 0);
      const c = await addNode(e.flow, 'Utilities/Info', 600, 0);
      const nodeEl = (id: string): HTMLElement | null =>
        e.container.querySelector(`.ff-node[data-node-id="${id}"]`);
      await until(() => nodeEl(a.id) != null && nodeEl(b.id) != null && nodeEl(c.id) != null);

      (e.flow as unknown as {beginConnectHints(n: string, k: string, s: string): void})
        .beginConnectHints(a.id, EXEC_OUT_KEY, 'output');
      expect(e.container.classList.contains('ff-connecting'), true, 'connect mode on');
      expect(e.container.classList.contains('ff-connecting-order'), true, 'order drag marked (squares stay visible)');
      expect(nodeEl(a.id)!.classList.contains('ff-node-source'), true, 'source undimmed');
      for (const other of [b, c]) {
        expect(nodeEl(other.id)!.classList.contains('ff-node-compat'), true, 'other node lit as a target');
        expect(nodeEl(other.id)!.querySelector('[data-testid="ff-exec-in"]')!
          .classList.contains('ff-socket-compat'), true, 'its exec-in square glows');
      }

      (e.flow as unknown as {endConnectHints(): void}).endConnectHints();
      expect(e.container.classList.contains('ff-connecting-order'), false, 'order class cleared on drop');
      expect(e.container.querySelectorAll('.ff-socket-compat, .ff-node-compat').length, 0, 'hints cleared');
    } finally {
      destroyEditor(e);
    }
  });

  test('an order edge overrides vertical position in the topological sort', async () => {
    const e = makeEditor();
    try {
      // `first` is placed LOW on the canvas, `second` HIGH — without an edge the
      // sort would run the higher (second) node first. The order edge flips it.
      const first = await addNode(e.flow, 'Utilities/Info', 0, 500);
      const second = await addNode(e.flow, 'Utilities/Info', 0, 10);

      const before = topologicalSort(e.flow);
      expect(before.indexOf(second.id) < before.indexOf(first.id), true, 'higher node first when unordered');

      await e.flow.addConnectionByKeys(first.id, EXEC_OUT_KEY, second.id, EXEC_IN_KEY);
      const after = topologicalSort(e.flow);
      expect(after.indexOf(first.id) < after.indexOf(second.id), true, 'order edge forces first → second');
    } finally {
      destroyEditor(e);
    }
  });

  test('order edges sequence the emitted script but add no data / variables', async () => {
    const e = makeEditor();
    try {
      // Two Log nodes, distinguishable by label, with `b` positioned ABOVE `a`.
      const a = await addNode(e.flow, 'Utilities/Log', 0, 400);
      a.properties['label'] = 'AAA';
      const b = await addNode(e.flow, 'Utilities/Log', 0, 20);
      b.properties['label'] = 'BBB';
      await e.flow.addConnectionByKeys(a.id, EXEC_OUT_KEY, b.id, EXEC_IN_KEY);

      const script = emitScript(e.flow, SETTINGS);
      // The ordering edge is invisible in the output — no exec ports leak in.
      expect(script.includes('__exec'), false, 'no exec port keys in generated code');
      // ...but the run order is enforced: AAA logged before BBB.
      const iA = script.indexOf(`'AAA:'`);
      const iB = script.indexOf(`'BBB:'`);
      expect(iA >= 0 && iB >= 0, true, 'both Log steps present');
      expect(iA < iB, true, 'order edge sequenced AAA before BBB');
    } finally {
      destroyEditor(e);
    }
  });

  test('an order edge does not feed data into the target', async () => {
    const e = makeEditor();
    try {
      // A constant ordered before a Value Output via exec ports only — the
      // output has no data input, so it must emit nothing for its value.
      const c = await addNode(e.flow, 'Constants/String', 0, 0);
      c.properties['value'] = 'hello';
      const out = await addNode(e.flow, 'Outputs/Value Output', 300, 0);
      out.properties['paramName'] = 'result';
      await e.flow.addConnectionByKeys(c.id, EXEC_OUT_KEY, out.id, EXEC_IN_KEY);

      const script = emitScript(e.flow, SETTINGS);
      // The order edge must NOT be mistaken for the output's value.
      expect(/result\s*=/.test(script), false, 'order edge is not wired as the output value');
      expect(script.includes('__exec'), false);
    } finally {
      destroyEditor(e);
    }
  });

  test('cycle via order edges is detected by validation', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Utilities/Info', 0, 0);
      const b = await addNode(e.flow, 'Utilities/Info', 0, 100);
      await e.flow.addConnectionByKeys(a.id, EXEC_OUT_KEY, b.id, EXEC_IN_KEY);
      await e.flow.addConnectionByKeys(b.id, EXEC_OUT_KEY, a.id, EXEC_IN_KEY);
      const errors = validateGraph(e.flow).filter((r) => r.severity === 'error');
      expect(errors.some((r) => /cycle/i.test(r.message)), true, 'order-edge cycle reported');
    } finally {
      destroyEditor(e);
    }
  });

  test('order edges round-trip through serialization', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Utilities/Info', 0, 0);
      const b = await addNode(e.flow, 'Utilities/Info', 0, 100);
      await e.flow.addConnectionByKeys(a.id, EXEC_OUT_KEY, b.id, EXEC_IN_KEY);

      const doc = serializeFlow(e.flow, {scriptName: 'S', scriptDescription: '', tags: []});
      expect(doc.connections.length, 1);
      expect(doc.connections[0].sourceOutput, EXEC_OUT_KEY);
      expect(doc.connections[0].targetInput, EXEC_IN_KEY);

      await deserializeFlow(doc, e.flow);
      expect(e.flow.getConnectionCount(), 1, 'order edge restored');
      const conn = e.flow.getConnections()[0];
      expect(String(conn.sourceOutput), EXEC_OUT_KEY);
      expect(String(conn.targetInput), EXEC_IN_KEY);
    } finally {
      destroyEditor(e);
    }
  });
});

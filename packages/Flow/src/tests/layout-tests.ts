import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, createNode} from '../rete/node-factory';
import {computeLayers, estimateNodeWidth, estimateNodeHeight} from '../rete/graph-layout';
import {buildCreationScriptGraph, applyGraphToEditor} from '../import/creation-script-importer';
import {makeEditor, destroyEditor, addNode, nodesByFunc} from './test-utils';
import {FlowNode} from '../rete/scheme';

function noOverlap(nodes: FlowNode[]): boolean {
  for (let i = 0; i < nodes.length; i++) {
    for (let j = i + 1; j < nodes.length; j++) {
      const a = nodes[i]; const b = nodes[j];
      const separated =
        a.pos.x + estimateNodeWidth(a) <= b.pos.x || b.pos.x + estimateNodeWidth(b) <= a.pos.x ||
        a.pos.y + estimateNodeHeight(a) <= b.pos.y || b.pos.y + estimateNodeHeight(b) <= a.pos.y;
      if (!separated) return false;
    }
  }
  return true;
}

category('Flow: layout', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('computeLayers: chain increments, diamond takes the longest path', async () => {
    const a = createNode('Constants/String')!;
    const b = createNode('Utilities/ToString')!;
    const c = createNode('Utilities/ToString')!;
    const d = createNode('Outputs/Value Output')!;
    // a → b → c → d, plus a short-circuit a → d: d must land at the longest path (3).
    const edges = [
      {source: a, target: b}, {source: b, target: c},
      {source: c, target: d}, {source: a, target: d},
    ];
    const layers = computeLayers([a, b, c, d], edges);
    expect(layers.get(a), 0);
    expect(layers.get(b), 1);
    expect(layers.get(c), 2);
    expect(layers.get(d), 3, 'sink takes the longest path, not the short-circuit');
  });

  test('computeLayers: a node with no edges is a source (layer 0)', async () => {
    const a = createNode('Constants/String')!;
    const b = createNode('Constants/Int')!;
    const layers = computeLayers([a, b], []);
    expect(layers.get(a), 0);
    expect(layers.get(b), 0);
  });

  test('autoLayout arranges a scrambled graph left-to-right with no overlap', async () => {
    const e = makeEditor();
    try {
      // Positions are deliberately scrambled — autoLayout must fix them.
      const c = await addNode(e.flow, 'Constants/String', 900, 700);
      const t1 = await addNode(e.flow, 'Utilities/ToString', 10, 10);
      const t2 = await addNode(e.flow, 'Utilities/ToString', 400, 999);
      await e.flow.addConnectionByKeys(c.id, 'value', t1.id, 'value');
      await e.flow.addConnectionByKeys(t1.id, 'text', t2.id, 'value');

      await e.flow.autoLayout();

      // Every edge points strictly rightward.
      for (const conn of e.flow.getConnections()) {
        const s = e.flow.getNodeById(conn.source)!;
        const t = e.flow.getNodeById(conn.target)!;
        expect(s.pos.x < t.pos.x, true, `${s.label} must be left of ${t.label}`);
      }
      expect(noOverlap(e.flow.getNodes()), true, 'no node boxes overlap');
    } finally {
      destroyEditor(e);
    }
  });

  test('autoLayout reproduces producer-above-consumer banding in the editor', async () => {
    const e = makeEditor();
    try {
      // Join (defined first) reads "Second" via a Select Table; the Second
      // producer path is defined later — autoLayout must still place the
      // producer band above the consumer band.
      const g = buildCreationScriptGraph([
        'Joined = JoinTables("First", "Second", ["Id"], ["Id"], ["Id"], ["Id"])',
        'Second = OpenFile("s.csv")',
        'AddNewColumn(Second, "1", "x")',
      ].join('\n'));
      await applyGraphToEditor(g, e.flow);

      // Scramble positions so the test exercises autoLayout, not the import layout.
      for (const n of e.flow.getNodes()) n.pos = {x: 0, y: 0};

      await e.flow.autoLayout();

      const join = nodesByFunc(g, 'JoinTables')[0];
      const openSecond = nodesByFunc(g, 'OpenFile')[0];
      const setSecond = nodesByFunc(g, 'SetVar').find((n) => n.inputValues['variableName'] === 'Second')!;
      expect(openSecond.pos.y < join.pos.y, true, 'producer OpenFile above the join');
      expect(setSecond.pos.y < join.pos.y, true, 'producer SetVar above the join');
      expect(noOverlap(e.flow.getNodes()), true, 'no node boxes overlap');
    } finally {
      destroyEditor(e);
    }
  });
});

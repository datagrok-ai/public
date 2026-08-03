import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {makeEditor, destroyEditor, addNode} from './test-utils';

const SETTINGS = {scriptName: 'RoundTrip', scriptDescription: 'desc', tags: ['a', 'b']};

category('Flow: serializer', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('serialize captures nodes, connections and settings', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input', 10, 20);
      input.properties['paramName'] = 'src';
      const output = await addNode(e.flow, 'Outputs/Table Output', 340, 20);
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');

      const doc = serializeFlow(e.flow, SETTINGS);
      expect(doc.version, '2.0');
      expect(doc.nodes.length, 2);
      expect(doc.connections.length, 1);
      expect(doc.metadata?.settings?.scriptName, 'RoundTrip');
      const inputDoc = doc.nodes.find((n) => n.typeName === 'Inputs/Table Input');
      expect(inputDoc?.properties['paramName'], 'src');
      expect(inputDoc?.pos.x, 10);
    } finally {
      destroyEditor(e);
    }
  });

  test('round-trip preserves topology and properties', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input', 10, 20);
      input.properties['paramName'] = 'src';
      const output = await addNode(e.flow, 'Outputs/Table Output', 340, 20);
      output.properties['paramName'] = 'finalResult';
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');

      const doc = serializeFlow(e.flow, SETTINGS);
      await deserializeFlow(doc, e.flow);

      expect(e.flow.getNodeCount(), 2);
      expect(e.flow.getConnectionCount(), 1);
      const nodes = e.flow.getNodes();
      const reInput = nodes.find((n) => n.dgTypeName === 'Inputs/Table Input');
      const reOutput = nodes.find((n) => n.dgTypeName === 'Outputs/Table Output');
      expect(reInput?.properties['paramName'], 'src');
      expect(reOutput?.properties['paramName'], 'finalResult');
      // The single connection still runs input → output.
      const conn = e.flow.getConnections()[0];
      expect(conn.source, reInput!.id);
      expect(conn.target, reOutput!.id);
    } finally {
      destroyEditor(e);
    }
  });

  test('unknown node types are skipped on deserialize', async () => {
    const e = makeEditor();
    try {
      const doc = serializeFlow(e.flow, SETTINGS);
      doc.nodes.push({
        id: 'ghost', typeName: 'Bogus/Type', label: 'Ghost',
        description: '', pos: {x: 0, y: 0}, properties: {}, inputValues: {},
      });
      await deserializeFlow(doc, e.flow);
      expect(e.flow.getNodeCount(), 0); // ghost skipped, nothing else added
    } finally {
      destroyEditor(e);
    }
  });
});

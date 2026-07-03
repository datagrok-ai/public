/** Manual viewer nodes: a `Viewers/<type>` node creates a DG viewer from a
 *  wired table via `table.plot.fromType(type, {})` + `setOptions(look)`. Covers
 *  node construction, the curated option specs, and the emitted code (clean +
 *  instrumented). */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import {registerBuiltinNodes, createNode, VIEWER_NODE_TYPES} from '../rete/node-factory';
import {emitScript} from '../compiler/script-emitter';
import {makeEditor, destroyEditor, addNode} from './test-utils';

const SETTINGS = {name: 'T', description: '', tags: []};

category('Flow: viewers', () => {
  before(async () => {
    registerBuiltinNodes();
  });

  test('core viewer node types are registered (Scatter Plot present)', async () => {
    const labels = VIEWER_NODE_TYPES.map((v) => v.label);
    expect(labels.includes('Scatter Plot'), true, 'Scatter Plot viewer registered');
    expect(VIEWER_NODE_TYPES.some((v) => v.nodeTypeName === 'Viewers/Scatter Plot'), true);
  });

  test('a viewer node has a table input, a viewer output, and its type/specs', async () => {
    const node = createNode('Viewers/Scatter Plot')!;
    expect(!!node, true, 'created');
    expect(node.dgNodeType, 'utility');
    expect(node.dgOutputType, 'viewer');
    expect(node.properties['viewerType'], 'Scatter plot');
    expect('table' in node.inputs, true, 'has table input');
    expect('viewer' in node.outputs, true, 'has viewer output');
    const specs = node.properties['viewerOptionSpecs'] as Array<{key: string}>;
    expect(specs.some((s) => s.key === 'xColumnName'), true, 'exposes X');
    expect(node.requiredInputs.includes('table'), true, 'table is required');
    // Column options are connectable input sockets, not panel-only.
    for (const key of ['xColumnName', 'yColumnName', 'colorColumnName', 'sizeColumnName'])
      expect(key in node.inputs, true, `${key} is a connectable input`);
  });

  test('a wired column option feeds setOptions via its .name', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const selCol = await addNode(e.flow, 'Utilities/Select Column');
      selCol.properties['columnName'] = 'age';
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot');
      await e.flow.addConnectionByKeys(input.id, 'table', selCol.id, 'table');
      await e.flow.addConnectionByKeys(input.id, 'table', viewer.id, 'table');
      await e.flow.addConnectionByKeys(selCol.id, 'column', viewer.id, 'xColumnName');

      const clean = emitScript(e.flow, SETTINGS);
      expect(/"xColumnName": \(.+\)\.name/.test(clean), true, 'wired column option uses its .name');
    } finally {
      destroyEditor(e);
    }
  });

  test('emits plot.fromType + setOptions from the wired table and stored look', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot');
      (viewer.properties['viewerLook'] as Record<string, unknown>) = {xColumnName: 'age', title: 'My Plot'};
      await e.flow.addConnectionByKeys(input.id, 'table', viewer.id, 'table');

      const clean = emitScript(e.flow, SETTINGS);
      expect(/\.plot\.fromType\("Scatter plot", \{\}\)/.test(clean), true, 'creates the viewer via fromType');
      expect(clean.includes('.setOptions({'), true, 'applies stored options');
      expect(/"xColumnName":\s*"age"/.test(clean), true, 'options include the exposed column');
      expect(clean.includes('await '), true, 'fromType is awaited');

      const inst = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'r1'});
      expect(inst.includes('__ff_summarize'), true, 'instrumented summarizes the viewer');
      expect(/__ff_summarize\(\w+, 'viewer'\)/.test(inst), true, 'declared as a viewer output');
    } finally {
      destroyEditor(e);
    }
  });

  test('a viewer node with no wired table emits nothing (can\'t plot)', async () => {
    const e = makeEditor();
    try {
      await addNode(e.flow, 'Viewers/Scatter Plot');
      const clean = emitScript(e.flow, SETTINGS);
      expect(clean.includes('plot.fromType'), false, 'no viewer code without a table');
    } finally {
      destroyEditor(e);
    }
  });
});

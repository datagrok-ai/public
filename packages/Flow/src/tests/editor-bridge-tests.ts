/** Regression tests for the per-node editor bridge (`FlowNode.editorBridge`).
 *
 *  The bridge used to be a page-level global (`window.__ff_editor`) that every
 *  `FlowEditor` construction rebound to itself and every `destroy()` deleted.
 *  With a second editor on the page (a file preview, the creation-script
 *  dialog, a detached compile editor) the first editor's collapse carets went
 *  dead and its collapsed nodes dropped their socket DOM on re-render, leaving
 *  connections dangling at stale positions. The bridge now lives on each node,
 *  stamped by the owning editor, so any number of editors coexist. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes} from '../rete/node-factory';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {makeEditor, destroyEditor, addNode} from './test-utils';

const SETTINGS = {scriptName: 'BridgeFlow', scriptDescription: '', tags: ['funcflow']};

category('Flow: editor bridge', () => {
  before(async () => {
    registerBuiltinNodes();
  });

  test('every added node is stamped with a working bridge', async () => {
    const e = makeEditor();
    try {
      const info = await addNode(e.flow, 'Utilities/Info');
      expect(info.editorBridge != null, true, 'bridge stamped on add');

      expect(info.collapsed, false);
      info.editorBridge!.toggleCollapsed(info.id);
      expect(info.collapsed, true, 'bridge toggles collapse on the owning editor');

      const input = await addNode(e.flow, 'Inputs/String Input');
      expect(info.editorBridge!.isSocketConnected(info.id, 'input', 'message'), false);
      await e.flow.addConnectionByKeys(input.id, 'value', info.id, 'message');
      expect(info.editorBridge!.isSocketConnected(info.id, 'input', 'message'), true,
        'bridge sees the owning editor\'s connections');
    } finally {
      destroyEditor(e);
    }
  });

  test('two live editors keep independent bridges', async () => {
    const a = makeEditor();
    const b = makeEditor(); // second construction must not hijack a's bridge
    try {
      const nodeA = await addNode(a.flow, 'Utilities/Info');
      const nodeB = await addNode(b.flow, 'Utilities/Info');
      expect(nodeA.editorBridge === nodeB.editorBridge, false, 'each editor stamps its own bridge');

      nodeA.editorBridge!.toggleCollapsed(nodeA.id);
      expect(nodeA.collapsed, true);
      expect(nodeB.collapsed, false, 'toggling in a never leaks into b');

      const inputA = await addNode(a.flow, 'Inputs/String Input');
      await a.flow.addConnectionByKeys(inputA.id, 'value', nodeA.id, 'message');
      expect(nodeA.editorBridge!.isSocketConnected(nodeA.id, 'input', 'message'), true);
      expect(nodeB.editorBridge!.isSocketConnected(nodeB.id, 'input', 'message'), false);
    } finally {
      destroyEditor(b);
      destroyEditor(a);
    }
  });

  test('destroying another editor leaves the first bridge intact', async () => {
    const a = makeEditor();
    try {
      const input = await addNode(a.flow, 'Inputs/String Input');
      const info = await addNode(a.flow, 'Utilities/Info');
      await a.flow.addConnectionByKeys(input.id, 'value', info.id, 'message');

      // The old global bridge was deleted by ANY editor's destroy() — the
      // detached compile editor did exactly this in normal use.
      const detached = makeEditor();
      destroyEditor(detached);

      info.editorBridge!.toggleCollapsed(info.id);
      expect(info.collapsed, true, 'caret still works after another editor is destroyed');
      expect(info.editorBridge!.isSocketConnected(info.id, 'input', 'message'), true,
        'collapsed-socket rendering still sees connections');
    } finally {
      destroyEditor(a);
    }
  });

  test('deserialized nodes are stamped by the target editor', async () => {
    const source = makeEditor();
    const target = makeEditor();
    try {
      const input = await addNode(source.flow, 'Inputs/String Input');
      const info = await addNode(source.flow, 'Utilities/Info');
      await source.flow.addConnectionByKeys(input.id, 'value', info.id, 'message');
      const doc = serializeFlow(source.flow, SETTINGS);
      expect(JSON.stringify(doc).includes('editorBridge'), false, 'bridge never serializes');

      await deserializeFlow(doc, target.flow);
      for (const node of target.flow.getNodes()) {
        expect(node.editorBridge != null, true, `deserialized "${node.label}" has a bridge`);
        node.editorBridge!.toggleCollapsed(node.id);
        expect(node.collapsed, true, `bridge of "${node.label}" targets the new editor`);
      }
    } finally {
      destroyEditor(target);
      destroyEditor(source);
    }
  });
});

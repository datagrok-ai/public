/** Duplicate / copy / paste: node sets copy with their internal connections. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes} from '../rete/node-factory';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

function key(init: KeyboardEventInit): void {
  window.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, cancelable: true, ...init}));
}

category('Flow: clipboard', () => {
  before(async () => {
    registerBuiltinNodes();
  });

  test('the node context menu offers Copy; copyNodes feeds paste without a selection', async () => {
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, 'Inputs/String Input', 0, 0);
      await until(() => !!e.container.querySelector(`.ff-node[data-node-id="${node.id}"]`));

      const nodeEl = e.container.querySelector<HTMLElement>(`.ff-node[data-node-id="${node.id}"]`)!;
      nodeEl.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: 10, clientY: 10}));
      const labels = (): string[] => Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
        .map((el) => el.textContent?.trim() ?? '');
      expect(await until(() => labels().includes('Copy')), true, `no Copy item; saw: ${labels().join(', ')}`);
      for (const el of Array.from(document.querySelectorAll('.d4-menu-popup, .d4-menu-dropdown')))
        el.remove();

      // What the menu item runs: copy THIS node (no selection needed), then paste.
      expect(e.flow.copyNodes([node.id]), 1, 'one node copied');
      const pasted = await e.flow.pasteClipboard();
      expect(pasted.length, 1, 'paste materialized the copy');
      expect(e.flow.getNodeCount(), 2);
      expect(e.flow.copyNodes([]), 0, 'nothing to copy leaves the clipboard untouched');
      expect((await e.flow.pasteClipboard()).length, 1, 'previous clipboard still pastes');
    } finally {
      destroyEditor(e);
    }
  });

  test('duplicateNodes copies nodes + internal connections and selects the copies', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const output = await addNode(e.flow, 'Outputs/Table Output', 300, 0);
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');

      const copies = await e.flow.duplicateNodes([input.id, output.id]);
      expect(copies.length, 2, 'both nodes duplicated');
      expect(e.flow.getNodeCount(), 4);
      expect(e.flow.getConnectionCount(), 2, 'the connection between the copies was duplicated');

      const copyIds = new Set(copies.map((n) => n.id));
      const copiedConn = e.flow.getConnections()
        .find((c) => copyIds.has(c.source) && copyIds.has(c.target));
      expect(!!copiedConn, true, 'the duplicated connection links the copies, not the originals');

      const selected = new Set(e.flow.getSelectedNodeIds());
      expect(copies.every((n) => selected.has(n.id)), true, 'the copies are selected');
      expect(selected.has(input.id) || selected.has(output.id), false, 'the originals are not');

      const outCopy = copies.find((n) => n.dgNodeType === 'output')!;
      expect(outCopy.properties.paramName !== output.properties.paramName, true,
        'the copied output got a unique paramName (no duplicate-variable error)');
      const inCopy = copies.find((n) => n.dgNodeType === 'input')!;
      expect(inCopy.pos.x, input.pos.x + 30, 'copies are offset from the originals');
    } finally {
      destroyEditor(e);
    }
  });

  test('right-clicked Duplicate of a single node still works via duplicateNodes', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/String Input', 0, 0);
      a.label = 'My Input';
      a.description = 'annotated';
      const [copy] = await e.flow.duplicateNodes([a.id]);
      expect(copy.label, 'My Input', 'label copied');
      expect(copy.description, 'annotated', 'description copied');
      expect(e.flow.getNodeCount(), 2);
    } finally {
      destroyEditor(e);
    }
  });

  test('copySelection + pasteClipboard round-trip; repeat pastes fan out', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const output = await addNode(e.flow, 'Outputs/Table Output', 300, 0);
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      await e.flow.selectNode(input.id);
      await e.flow.selectNode(output.id, true);

      expect(e.flow.copySelection(), 2, 'two nodes copied');
      const first = await e.flow.pasteClipboard();
      expect(first.length, 2, 'paste materializes both');
      expect(e.flow.getConnectionCount(), 2, 'paste re-creates the internal connection');
      const firstIn = first.find((n) => n.dgNodeType === 'input')!;
      expect(firstIn.pos.x, 30, 'first paste offsets by 30');

      const second = await e.flow.pasteClipboard();
      const secondIn = second.find((n) => n.dgNodeType === 'input')!;
      expect(secondIn.pos.x, 60, 'second paste fans out further');
      expect(e.flow.getNodeCount(), 6);

      const selected = new Set(e.flow.getSelectedNodeIds());
      expect(second.every((n) => selected.has(n.id)), true, 'the latest paste is the selection');
      expect(selected.size, 2, 'only the latest paste is selected');
    } finally {
      destroyEditor(e);
    }
  });

  test('Ctrl+C / Ctrl+V keyboard path copies the selection', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/String Input', 0, 0);
      await until(() => !!e.container.querySelector(`.ff-node[data-node-id="${a.id}"]`));
      await e.flow.selectNode(a.id);

      key({key: 'c', ctrlKey: true});
      key({key: 'v', ctrlKey: true});
      expect(await until(() => e.flow.getNodeCount() === 2), true, 'Ctrl+C, Ctrl+V duplicated the node');

      key({key: 'v', ctrlKey: true});
      expect(await until(() => e.flow.getNodeCount() === 3), true, 'the clipboard pastes repeatedly');
    } finally {
      destroyEditor(e);
    }
  });

  test('an empty selection leaves the clipboard untouched', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/String Input', 0, 0);
      await e.flow.selectNode(a.id);
      expect(e.flow.copySelection(), 1);
      await e.flow.unselectAllNodes();
      expect(e.flow.copySelection(), 0, 'nothing selected → nothing copied');
      const pasted = await e.flow.pasteClipboard();
      expect(pasted.length, 1, 'the previous clipboard survives an empty copy');
    } finally {
      destroyEditor(e);
    }
  });

  test('right-click on a selected node keeps the multi-selection', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/String Input', 0, 0);
      const b = await addNode(e.flow, 'Inputs/String Input', 250, 0);
      await until(() => !!e.container.querySelector(`.ff-node[data-node-id="${b.id}"]`));
      await e.flow.selectNode(a.id, true);
      await e.flow.selectNode(b.id, true);

      const el = e.container.querySelector(`.ff-node[data-node-id="${a.id}"]`) as HTMLElement;
      const r = el.getBoundingClientRect();
      const init = {bubbles: true, cancelable: true, button: 2,
        clientX: r.left + r.width / 2, clientY: r.top + 8};
      el.dispatchEvent(new PointerEvent('pointerdown', init));
      await new Promise((res) => setTimeout(res, 60));
      el.dispatchEvent(new PointerEvent('pointerup', init));
      await new Promise((res) => setTimeout(res, 60));
      expect(e.flow.getSelectedNodeIds().length, 2,
        'right-click must not clear the selection Duplicate acts on');
    } finally {
      destroyEditor(e);
    }
  });

  test('Escape clears the selection', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/String Input', 0, 0);
      const b = await addNode(e.flow, 'Inputs/String Input', 250, 0);
      await e.flow.selectNode(a.id, true);
      await e.flow.selectNode(b.id, true);
      expect(e.flow.getSelectedNodeIds().length, 2);

      key({key: 'Escape'});
      expect(await until(() => e.flow.getSelectedNodeIds().length === 0), true,
        'Escape deselects everything');
    } finally {
      destroyEditor(e);
    }
  });
});

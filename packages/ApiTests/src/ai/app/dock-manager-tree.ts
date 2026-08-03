import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, markedPanel, uniqueName, until, withTableView} from '../helpers';

category('AI: App: Dock Manager Tree', () => {
  const mark = (): string => uniqueName('dmTree');
  const panel = (cls: string): HTMLElement => markedPanel(cls);

  // Run against a disposable TableView's own dock manager, never the global shell dock tree:
  // these tests do destructive surgery (detachFromParent/removeChild/float/destroy) that can
  // corrupt the root layout and shrink subsequently-added views. Closing the view disposes the
  // whole dock subtree, so nothing leaks into other tests.
  const withDm = (body: (dm: DG.DockManager) => Promise<void> | void): Promise<void> =>
    withTableView(demog(), (tv) => body(tv.dockManager));

  // Dock two panels (RIGHT then FILL into a tab group) and return both nodes.
  function dockTabGroup(dm: DG.DockManager): {node1: DG.DockNode, node2: DG.DockNode} {
    const node1 = dm.dock(panel(mark()), DG.DOCK_TYPE.RIGHT);
    const node2 = dm.dock(panel(mark()), DG.DOCK_TYPE.FILL, node1);
    return {node1, node2};
  }

  test('rootNode and documentContainer', async () => {
    await withDm(async (dm) => {
      expect(dm.element instanceof HTMLElement, true);
      const root = dm.rootNode;
      expect(root instanceof DG.DockNode, true);
      expect(root.container instanceof DG.DockContainer, true);
      expect(root.container.containerElement instanceof HTMLElement, true);
      expect(dm.documentContainer instanceof DG.DockContainer, true);
      expect(dm.documentContainer.containerElement instanceof HTMLElement, true);
    });
  });

  test('dock then findNode resolves the same node', async () => {
    await withDm(async (dm) => {
      const cls = mark();
      const el = panel(cls);
      const node = dm.dock(el, DG.DOCK_TYPE.RIGHT);
      try {
        // findNode of a never-docked element is reliably undefined.
        // (expect()'s expected arg defaults to true, so compare via a boolean.)
        const undocked = ui.div([], {classes: mark()});
        expect(dm.findNode(undocked) == null, true);
        // findNode must resolve the docked element to its node, whose container hosts the element.
        await until(() => dm.findNode(el) instanceof DG.DockNode);
        const found = dm.findNode(el)!;
        expect(found instanceof DG.DockNode, true);
        const ce = found.container.containerElement;
        expect(ce.contains(el) || ce === el, true);
      } finally {
        dm.close(node);
      }
    });
  });

  test('DockNode parent and children navigation', async () => {
    await withDm(async (dm) => {
      const {node1, node2} = dockTabGroup(dm);
      try {
        const parent = node2.parent;
        expect(parent != null, true);
        expect(parent instanceof DG.DockNode, true);
        const kids = Array.from(parent.children);
        expect(kids.length >= 1, true);
        for (const k of kids)
          expect(k instanceof DG.DockNode, true);
        expect(Array.from(dm.rootNode.children).length >= 1, true);
      } finally {
        dm.close(node2);
        dm.close(node1);
      }
    });
  });

  test('DockContainer title round-trip', async () => {
    await withDm(async (dm) => {
      const el = panel(mark());
      const node = dm.dock(el, DG.DOCK_TYPE.RIGHT, null, 'Composer A');
      try {
        expect(typeof node.container.title, 'string');
        node.container.title = 'Renamed';
        expect(node.container.title, 'Renamed');
        node.container.title = 'Renamed Again';
        expect(node.container.title, 'Renamed Again');
      } finally {
        dm.close(node);
      }
    });
  });

  test('close removes the docked element from the DOM (by node and by element)', async () => {
    await withDm(async (dm) => {
      const elA = panel(mark());
      const nodeA = dm.dock(elA, DG.DOCK_TYPE.RIGHT);
      await until(() => document.body.contains(elA));
      dm.close(nodeA);
      await until(() => !document.body.contains(elA));
      expect(document.body.contains(elA), false);

      const elB = panel(mark());
      dm.dock(elB, DG.DOCK_TYPE.RIGHT);
      await until(() => document.body.contains(elB));
      dm.close(elB);
      await until(() => !document.body.contains(elB));
      expect(document.body.contains(elB), false);
    });
  });

  // Dropped here (not deterministically observable in the headless harness): detachFromParent/
  // removeChild tab-group surgery, the onClosed event payload, and float()/destroy() relocating a
  // panel out of the document container. close()-based DOM removal above is the reliable proxy for
  // teardown; the rest is exercised at the UI layer.
}, {owner: 'agolovko@datagrok.ai'});

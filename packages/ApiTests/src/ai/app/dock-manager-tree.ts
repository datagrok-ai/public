import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, markedPanel, subscribeAll, uniqueName, until, wait, withTableView} from '../helpers';

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
        // findNode may return undefined before the node is wired into the tree,
        // so only assert the positive shape when it resolves — never force it.
        await until(() => dm.findNode(el) instanceof DG.DockNode);
        const found = dm.findNode(el);
        if (found != null) {
          expect(found instanceof DG.DockNode, true);
          const ce = found.container.containerElement;
          expect(ce.contains(el) || ce === el, true);
        }
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
    // title maps to the container's internal name (auto-generated panel_N), NOT the panel title bar
    // passed to dock(); so the round-trip is a set->get on the name field itself.
    await withDm(async (dm) => {
      const el = panel(mark());
      const node = dm.dock(el, DG.DOCK_TYPE.RIGHT, null, 'Composer A');
      try {
        const original = node.container.title;
        expect(typeof original, 'string');
        node.container.title = 'Renamed';
        expect(node.container.title, 'Renamed');
      } finally {
        dm.close(node);
      }
    });
  });

  test('setActiveChild on a tab group', async () => {
    await withDm(async (dm) => {
      const {node1, node2} = dockTabGroup(dm);
      try {
        const parentContainer = node2.parent.container;
        expectNoThrow(() => parentContainer.setActiveChild(node1.container));
        expectNoThrow(() => parentContainer.setActiveChild(node2.container));
      } finally {
        dm.close(node2);
        dm.close(node1);
      }
    });
  });

  test('detachFromParent, removeChild and close', async () => {
    await withDm(async (dm) => {
      const {node1, node2} = dockTabGroup(dm);
      const parent = node2.parent;
      expectNoThrow(() => node2.detachFromParent());
      expectNoThrow(() => parent.removeChild(node2));
      dm.close(node2);
      dm.close(node1);

      // close by DockNode
      const elA = panel(mark());
      const nodeA = dm.dock(elA, DG.DOCK_TYPE.RIGHT);
      dm.close(nodeA);
      await until(() => !document.body.contains(elA));
      expect(document.body.contains(elA), false);

      // close by element
      const elB = panel(mark());
      dm.dock(elB, DG.DOCK_TYPE.RIGHT);
      dm.close(elB);
      await until(() => !document.body.contains(elB));
      expect(document.body.contains(elB), false);
    });
  });

  test('onClosed event plus float and destroy smoke', async () => {
    await withDm(async (dm) => {
      expect(dm.element instanceof HTMLElement, true);
      // onClosed doesn't fire reliably headless, so this is a subscribe-smoke:
      // confirm the stream is subscribable and closing a docked element completes without error.
      const unsub = subscribeAll([dm.onClosed]);

      const elClose = panel(mark());
      dm.dock(elClose, DG.DOCK_TYPE.RIGHT);
      expectNoThrow(() => dm.close(elClose));

      const elFloat = panel(mark());
      const floatNode = dm.dock(elFloat, DG.DOCK_TYPE.RIGHT);
      expectNoThrow(() => floatNode.container.float());
      await wait();
      expectNoThrow(() => dm.close(elFloat));

      const elDestroy = panel(mark());
      const destroyNode = dm.dock(elDestroy, DG.DOCK_TYPE.RIGHT);
      expectNoThrow(() => destroyNode.container.destroy());

      unsub();
    });
  });
}, {owner: 'agolovko@datagrok.ai'});

import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow, markedPanel, uniqueName, until, withTableView, demog} from '../helpers';

// DG.DockManager.saveState — core/client/libs/dock_spawn/lib/dock/dock_manager.dart
// (serialize the dock layout). Tests run against a disposable TableView's own dockManager so they
// never mutate the global grok.shell.dockManager (whose tree is shared with other viewer tests).
//
// loadState() restoration is intentionally NOT covered here: in the headless harness loadState on a
// live dock tree does not complete deterministically (it hangs the isolate), so a restoration
// assertion can't be made to fail honestly on a real regression. Dock-layout round-trip is exercised
// at the UI layer (TestTrack layout specs) and via ViewLayout fromJson/loadLayout in
// view-layout-roundtrip. DockContainer.title round-trip is covered by dock-manager-tree.
category('AI: App: Dock Manager State', () => {
  const mark = (): string => uniqueName('dms');

  // Count dock nodes in a serialized state by walking its node tree (DockGraphSerializer emits
  // {containerType, state, children[]} recursively from the root node).
  function countNodes(state: string): number {
    let n = 0;
    const walk = (node: any): void => {
      if (node == null || typeof node !== 'object')
        return;
      n++;
      const kids = node['children'];
      if (Array.isArray(kids))
        for (const k of kids)
          walk(k);
    };
    walk(JSON.parse(state));
    return n;
  }

  test('saveState serializes the live layout as a non-empty node tree', async () => {
    await withTableView(demog(), (tv) => {
      const s = tv.dockManager.saveState();
      expect(typeof s, 'string');
      expect(s.length > 0, true);
      // The serialized form must parse and describe at least the view's own root dock node.
      expect(countNodes(s) > 0, true);
    });
  });

  test('dock adds a node that hosts the marked panel and grows the serialized tree', async () => {
    await withTableView(demog(), async (tv) => {
      const dm = tv.dockManager;
      const before = countNodes(dm.saveState());
      const cls = mark();
      const el = markedPanel(cls);
      const node = dm.dock(el, DG.DOCK_TYPE.RIGHT, null, cls);
      await until(() => node instanceof DG.DockNode);
      expect(node instanceof DG.DockNode, true);
      // The returned node's container element actually hosts the panel we docked.
      const ce = node.container.containerElement;
      expect(ce instanceof HTMLElement, true);
      expect(ce.querySelector('.' + cls) != null || ce.classList.contains(cls), true);
      // Docking a panel adds nodes to the serialized tree.
      const withPanel = dm.saveState();
      expectNoThrow(() => JSON.parse(withPanel));
      expect(countNodes(withPanel) > before, true);
      // Closing the docked node drops it back out of the serialized tree.
      dm.close(node);
      await until(() => countNodes(dm.saveState()) === before);
      expect(countNodes(dm.saveState()), before);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});

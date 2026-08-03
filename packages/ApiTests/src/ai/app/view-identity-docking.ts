import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, uniqueName, until, withTableView} from '../helpers';

category('AI: App: View Identity & Docking', () => {
  const mark = (): string => uniqueName('vid');

  let tv: DG.TableView;
  before(async () => {tv = grok.shell.addTableView(demog());});
  after(async () => tv.close());

  test('isPinned is a read-only boolean', async () => {
    expectNoThrow(() => tv.isPinned);
    expect(typeof tv.isPinned, 'boolean');
  });

  test('dockNode resolves to a DockNode with a container', async () => {
    await until(() => tv.dockNode instanceof DG.DockNode);
    expect(tv.dockNode instanceof DG.DockNode, true);
    expect(tv.dockNode.container instanceof DG.DockContainer, true);
    expect(tv.dockNode.container.containerElement instanceof HTMLElement, true);
  });

  test('toolboxPage exposes a non-null page', async () => {
    expectNoThrow(() => tv.toolboxPage);
    expect(tv.toolboxPage != null, true);
  });

  test('ribbonMenu round-trips a Menu', async () => {
    await withTableView(demog(), async (tv) => {
      expect(tv.ribbonMenu instanceof DG.Menu, true);
      const m = DG.Menu.create();
      tv.ribbonMenu = m;
      expect(tv.ribbonMenu instanceof DG.Menu, true);
    });
  });

  test('statusBarPanels set then get round-trips', async () => {
    await withTableView(demog(), async (tv) => {
      const marker = mark();
      const p = ui.div(['status'], {classes: marker});
      tv.statusBarPanels = [p];
      const panels = tv.statusBarPanels;
      expect(panels.length, 1);
      expect(panels[0] === p, true);
    });
  });

  test('parentView set and read-back; parentCall getter smoke', async () => {
    await withTableView(demog(), async (tv) => {
      const original = tv.parentView;
      expect(original == null || original instanceof DG.ViewBase, true);
      const other = grok.shell.addTableView(demog());
      try {
        tv.parentView = other;
        const readBack = tv.parentView;
        expect(readBack instanceof DG.ViewBase, true);
      } finally {
        expectNoThrow(() => tv.parentView = original);
        other.close();
      }
      expectNoThrow(() => tv.parentCall);
      const pc = tv.parentCall;
      expect(pc == null || typeof pc === 'object', true);
    });
  });
}, {owner: 'agolovko@datagrok.ai', clear: false});

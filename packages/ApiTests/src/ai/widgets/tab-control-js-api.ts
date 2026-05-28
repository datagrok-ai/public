// DG.TabControl — core/client/d4/lib/src/widgets/tab_control/tab_control.dart (scenario: tab-control-js-api)
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectFiresWithin, subscribeAll, expectNoThrow} from '../helpers';

category('AI: Widgets: TabControl JS API', () => {
  test('create + addPane + currentPane getter reflects first activated pane', async () => {
    const tc = DG.TabControl.create();
    expect(tc instanceof DG.TabControl, true);
    const alpha = tc.addPane('alpha', () => ui.divText('alpha-content'));
    expect(alpha instanceof DG.TabPane, true);
    expect(alpha.name, 'alpha');
    expect(tc.panes.length, 1);
    // First addPane auto-activates the pane (currentPage ??= sheet).
    expect(tc.currentPane != null, true);
    expect(tc.currentPane.name, 'alpha');
  });

  test('currentPane setter switches pane and fires onTabChanged', async () => {
    const tc = DG.TabControl.create();
    tc.addPane('alpha', () => ui.divText('a'));
    const beta = tc.addPane('beta', () => ui.divText('b'));
    // currentPane is still alpha (??= does not overwrite); switching to beta is a real change.
    expect(tc.currentPane.name, 'alpha');
    await expectFiresWithin(tc.onTabChanged, () => {tc.currentPane = beta;});
    expect(tc.currentPane.name, 'beta');
  });

  test('onBeforeTabChanged fires before a real switch', async () => {
    const tc = DG.TabControl.create();
    const alpha = tc.addPane('alpha', () => ui.divText('a'));
    tc.addPane('beta', () => ui.divText('b'));
    // Active pane is alpha; switch back to a different pane to trigger the event.
    tc.currentPane = tc.getPane('beta');
    await expectFiresWithin(tc.onBeforeTabChanged, () => {tc.currentPane = alpha;});
    expect(tc.currentPane.name, 'alpha');
  });

  test('onTabAdded fires when a pane is added', async () => {
    const tc = DG.TabControl.create();
    tc.addPane('alpha', () => ui.divText('a'));
    await expectFiresWithin(tc.onTabAdded, () => {tc.addPane('gamma', () => ui.div());});
    expect(tc.panes.length, 2);
  });

  test('onTabRemoved is a healthy Observable', async () => {
    const tc = DG.TabControl.create();
    tc.addPane('alpha', () => ui.divText('a'));
    // Not cleanly triggerable from JS (no removePane wrapper; clear() does not fire it).
    subscribeAll([tc.onTabRemoved])();
  });

  test('TabPane parent + content', async () => {
    const tc = DG.TabControl.create();
    const pane = tc.addPane('alpha', () => ui.divText('alpha-body'));
    expect(pane.parent instanceof DG.TabControl, true);
    expect(pane.content instanceof HTMLElement, true);
    expect(pane.content.textContent != null && pane.content.textContent.indexOf('alpha-body') >= 0, true);
  });

  test('boundary: empty control + vertical factory', async () => {
    const empty = DG.TabControl.create();
    expect(empty.panes.length, 0);
    expectNoThrow(() => {const cur = empty.currentPane; expect(cur == null, true);});
    const vertical = DG.TabControl.create(true);
    expect(vertical instanceof DG.TabControl, true);
    const pane = vertical.addPane('v', () => ui.divText('v-body'));
    expect(pane != null, true);
    expect(vertical.currentPane != null, true);
    expect(vertical.currentPane.name, 'v');
  });
}, {owner: 'agolovko@datagrok.ai'});

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
    expect(tc.currentPane != null, true);
    expect(tc.currentPane.name, 'alpha');
  });

  test('currentPane setter switches pane and fires onTabChanged', async () => {
    const tc = DG.TabControl.create();
    tc.addPane('alpha', () => ui.divText('a'));
    const beta = tc.addPane('beta', () => ui.divText('b'));
    expect(tc.currentPane.name, 'alpha');
    await expectFiresWithin(tc.onTabChanged, () => {tc.currentPane = beta;});
    expect(tc.currentPane.name, 'beta');
  });

  test('onBeforeTabChanged fires before a real switch', async () => {
    const tc = DG.TabControl.create();
    const alpha = tc.addPane('alpha', () => ui.divText('a'));
    tc.addPane('beta', () => ui.divText('b'));
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
    const vertical = DG.TabControl.create({vertical: true});
    expect(vertical instanceof DG.TabControl, true);
    const pane = vertical.addPane('v', () => ui.divText('v-body'));
    expect(pane != null, true);
    expect(vertical.currentPane != null, true);
    expect(vertical.currentPane.name, 'v');
  });
  test('options object: vertical and key', async () => {
    expect(DG.TabControl.create({vertical: true}).root.classList.contains('d4-tab-vertical'), true);
    expect(DG.TabControl.create({}).root.classList.contains('d4-tab-vertical'), false);

    const key = `apitests-tabcontrol-${Date.now()}`;
    try {
      const tc = DG.TabControl.create({key: key});
      tc.addPane('alpha', () => ui.divText('a'));
      tc.currentPane = tc.addPane('beta', () => ui.divText('b'));
      expect(window.localStorage[`TabControl:${key}`], 'beta');
    }
    finally {
      window.localStorage.removeItem(`TabControl:${key}`);
    }
  });

  test('ui.tabControl: options object and deprecated positional args', async () => {
    const tabs = ui.tabControl({'first': ui.divText('1'), 'second': ui.divText('2')}, {vertical: true});
    expect(tabs.panes.length, 2);
    expect(tabs.root.classList.contains('d4-tab-vertical'), true);

    expect(ui.tabControl(null, true).root.classList.contains('d4-tab-vertical'), true);
    expect(ui.tabControl(null, false).root.classList.contains('d4-tab-vertical'), false);

    // untyped JS callers pass an explicit null to skip `vertical` and reach `key`
    expectNoThrow(() => {(ui.tabControl as any)(null, null);});
    expectNoThrow(() => {(DG.TabControl.create as any)(null);});

    const key = `apitests-tabcontrol-legacy-${Date.now()}`;
    try {
      const legacy = ui.tabControl({'alpha': ui.divText('a'), 'beta': ui.divText('b')}, true, key);
      legacy.currentPane = legacy.getPane('beta');
      expect(window.localStorage[`TabControl:${key}`], 'beta');
    }
    finally {
      window.localStorage.removeItem(`TabControl:${key}`);
    }
  });
}, {owner: 'agolovko@datagrok.ai'});

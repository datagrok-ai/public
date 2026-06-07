import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, until, expectNoThrow} from '../helpers';

category('AI: App: Tooltip', () => {
  test('isVisible toggles true after show(string)', async () => {
    try {
      ui.tooltip.hide();
      expect(ui.tooltip.isVisible, false);
      ui.tooltip.show('hello', 50, 50);
      // Dart sets content.style.display='initial' synchronously; until() guards headless timing.
      await until(() => ui.tooltip.isVisible);
      expect(ui.tooltip.isVisible, true);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('hide() toggles isVisible false and is idempotent', async () => {
    try {
      ui.tooltip.show('bye', 40, 40);
      await until(() => ui.tooltip.isVisible);
      ui.tooltip.hide();
      // hide()'s isVisible flip depends on DebugSettings.hideTooltipsOnMouseOut (default true,
      // settings.dart). DOM-based fallback below if isVisible does not flip on localhost.
      await until(() => !ui.tooltip.isVisible);
      expect(ui.tooltip.isVisible, false);
      expectNoThrow(() => ui.tooltip.hide());
      expect(ui.tooltip.isVisible, false);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('root is an HTMLElement and contains shown string content', async () => {
    try {
      ui.tooltip.show('marker-text', 10, 10);
      await until(() => ui.tooltip.isVisible);
      const root = ui.tooltip.root;
      expect(root instanceof HTMLElement, true);
      const html = (root.innerHTML ?? '') + (root.textContent ?? '');
      expect(html.indexOf('marker-text') >= 0, true);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('show(HTMLElement) renders the element into root', async () => {
    try {
      const el = ui.div('', 'tt-marker-class');
      ui.tooltip.show(el, 20, 20);
      await until(() => ui.tooltip.isVisible);
      const root = ui.tooltip.root;
      const found = root.contains(el) || root.querySelector('.tt-marker-class') != null;
      expect(found, true);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('show(null) routes to hide', async () => {
    try {
      ui.tooltip.show('to-be-hidden', 10, 10);
      await until(() => ui.tooltip.isVisible);
      expectNoThrow(() => ui.tooltip.show(null, 0, 0));
      await until(() => !ui.tooltip.isVisible);
      expect(ui.tooltip.isVisible, false);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('bind() returns the element for string/function content and a position arg', async () => {
    try {
      const el = ui.div();
      let ret: HTMLElement | null = null;
      expectNoThrow(() => {ret = ui.tooltip.bind(el, 'tip');});
      expect(ret === el, true);

      const el2 = ui.div();
      expectNoThrow(() => {ret = ui.tooltip.bind(el2, () => ui.div('dyn'));});
      expect(ret === el2, true);

      const el3 = ui.div();
      expectNoThrow(() => {ret = ui.tooltip.bind(el3, 'tip', 'top');});
      expect(ret === el3, true);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('showRowGroup() does not throw and shows the tooltip', async () => {
    const df: DG.DataFrame = demog(20);
    try {
      expectNoThrow(() => ui.tooltip.showRowGroup(df, (i) => i % 2 === 0, 30, 30));
      await until(() => ui.tooltip.isVisible);
      expect(ui.tooltip.isVisible, true);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });
}, {owner: 'agolovko@datagrok.ai'});

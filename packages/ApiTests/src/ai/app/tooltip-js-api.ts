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
      ui.tooltip.show(null, 0, 0);
      await until(() => !ui.tooltip.isVisible);
      expect(ui.tooltip.isVisible, false);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('bind() returns the same element for string/function content and a position arg', async () => {
    try {
      const el = ui.div();
      expect(ui.tooltip.bind(el, 'tip') === el, true);

      const el2 = ui.div();
      expect(ui.tooltip.bind(el2, () => ui.div('dyn')) === el2, true);

      const el3 = ui.div();
      expect(ui.tooltip.bind(el3, 'tip', 'top') === el3, true);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('showRowGroup() shows the tooltip', async () => {
    const df: DG.DataFrame = demog(20);
    try {
      ui.tooltip.showRowGroup(df, (i) => i % 2 === 0, 30, 30);
      await until(() => ui.tooltip.isVisible);
      expect(ui.tooltip.isVisible, true);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });
}, {owner: 'agolovko@datagrok.ai'});

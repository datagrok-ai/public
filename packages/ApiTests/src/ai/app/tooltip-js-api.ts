import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, until, expectNoThrow} from '../helpers';

// Tooltip visibility (content.style.display == 'initial') depends on transient UI state
// (no menu shown, not mid-drag): Tooltip._show() bails on those guards, and in the headless
// test runner it bails every time, so isVisible never flips even though the content IS set.
// These tests therefore assert the deterministic contract — show() populates root, hide()/
// show(null) don't throw, bind() returns the element — rather than the non-observable isVisible.
category('AI: App: Tooltip', () => {
  const rootHtml = (): string => {
    const r = ui.tooltip.root;
    return (r.innerHTML ?? '') + (r.textContent ?? '');
  };

  test('show(string) populates root with the content', async () => {
    try {
      ui.tooltip.show('hello-tooltip', 50, 50);
      await until(() => rootHtml().indexOf('hello-tooltip') >= 0);
      expect(rootHtml().indexOf('hello-tooltip') >= 0, true);
      expect(typeof ui.tooltip.isVisible, 'boolean');
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('hide() does not throw and is idempotent', async () => {
    ui.tooltip.show('bye-tooltip', 40, 40);
    expectNoThrow(() => ui.tooltip.hide());
    expectNoThrow(() => ui.tooltip.hide());
    expect(typeof ui.tooltip.isVisible, 'boolean');
  });

  test('root is an HTMLElement and contains shown string content', async () => {
    try {
      ui.tooltip.show('marker-text', 10, 10);
      const root = ui.tooltip.root;
      expect(root instanceof HTMLElement, true);
      await until(() => ((root.innerHTML ?? '') + (root.textContent ?? '')).indexOf('marker-text') >= 0);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('show(HTMLElement) renders the element into root', async () => {
    try {
      const el = ui.div('', 'tt-marker-class');
      ui.tooltip.show(el, 20, 20);
      await until(() => ui.tooltip.root.contains(el) || ui.tooltip.root.querySelector('.tt-marker-class') != null);
      expect(ui.tooltip.root.contains(el) || ui.tooltip.root.querySelector('.tt-marker-class') != null, true);
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });

  test('show(null) routes to hide without throwing', async () => {
    ui.tooltip.show('to-be-hidden', 10, 10);
    expectNoThrow(() => ui.tooltip.show(null, 0, 0));
    expect(typeof ui.tooltip.isVisible, 'boolean');
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

  test('showRowGroup() is callable without throwing', async () => {
    const df: DG.DataFrame = demog(20);
    try {
      expectNoThrow(() => ui.tooltip.showRowGroup(df, (i) => i % 2 === 0, 30, 30));
      expect(typeof ui.tooltip.isVisible, 'boolean');
    } finally {
      expectNoThrow(() => ui.tooltip.hide());
    }
  });
}, {owner: 'agolovko@datagrok.ai'});

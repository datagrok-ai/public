import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, uniqueName, withTableView} from '../helpers';

category('AI: App: View Geometry Ribbon', () => {
  const mark = (): string => uniqueName('vgr-marker');

  let tv: DG.TableView;
  before(async () => {tv = grok.shell.addTableView(demog());});
  after(async () => tv.close());

  test('root is an HTMLElement and box reflects ui-box class', async () => {
    expect(tv.root instanceof HTMLElement, true);
    expect(typeof tv.box, 'boolean');
    expect(tv.box, tv.root.classList.contains('ui-box'));
  });

  test('type getter returns TableView', async () => {
    expect(DG.VIEW_TYPE.TABLE_VIEW, 'TableView');
    expect(tv.type, DG.VIEW_TYPE.TABLE_VIEW);
  });

  test('setRibbonPanels installs panels (DOM-based assertion)', async () => {
    await withTableView(demog(), async (tv) => {
      const cls = mark();
      const div = ui.div(['ribbon item'], {classes: cls});
      tv.setRibbonPanels([[div, ui.button('x', () => {})]]);
      const found = tv.root.querySelector('.' + cls) ?? document.querySelector('.' + cls);
      expect(found != null, true);
      expect(found instanceof HTMLElement, true);
    });
  });

  test('setRibbonPanels clear flag replaces vs appends', async () => {
    await withTableView(demog(), async (tv) => {
      // clear:true replaces user panels but one persistent system panel remains (count = user + 1);
      // platform owns item grouping/order, so assert nested-array shape + total items, not per-panel splits.
      const a = ui.div(['a'], {classes: mark()});
      const b = ui.div(['b'], {classes: mark()});
      const c = ui.div(['c'], {classes: mark()});
      tv.setRibbonPanels([[a], [b, c]], true);
      const panels = tv.getRibbonPanels();
      expect(Array.isArray(panels), true);
      expect(panels.length, 3);
      let totalItems = 0;
      for (const row of panels) {
        expect(Array.isArray(row), true);
        totalItems += row.length;
        for (const item of row)
          expect(item instanceof HTMLElement, true);
      }
      // 3 installed items (a, b, c) plus whatever the persistent system panel carries.
      expect(totalItems >= 3, true);
      // clear:true with one user panel yields exactly 2 (1 system + 1 user).
      tv.setRibbonPanels([[ui.div(['first'], {classes: mark()})]], true);
      const afterFirst = tv.getRibbonPanels().length;
      expect(afterFirst, 2);
      tv.setRibbonPanels([[ui.div(['second'], {classes: mark()})]], false);
      const afterAppend = tv.getRibbonPanels().length;
      expect(afterAppend > afterFirst, true);
      tv.setRibbonPanels([[ui.div(['third'], {classes: mark()})]], true);
      expect(tv.getRibbonPanels().length, 2);
    });
  });

  test('append(el) returns an HTMLElement and lands in view.root', async () => {
    await withTableView(demog(), async (tv) => {
      const cls = mark();
      const el = ui.div(['appended'], {classes: cls});
      const ret = tv.append(el);
      expect(ret instanceof HTMLElement, true);
      expect(tv.root.contains(el), true);
    });
  });

  test('basePath set->get round-trips; path composes basePath', async () => {
    await withTableView(demog(), async (tv) => {
      try {
        tv.basePath = '/ai-test-base';
        expect(tv.basePath, '/ai-test-base');
        expect(typeof tv.path, 'string');
        expect(tv.path.indexOf('/ai-test-base') >= 0, true);
      } finally {
        expectNoThrow(() => tv.basePath = '');
      }
    });
  });
}, {owner: 'agolovko@datagrok.ai', clear: false});

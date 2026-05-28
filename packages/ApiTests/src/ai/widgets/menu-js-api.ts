// DG.Menu + DG.Balloon — core/client/d4/lib/src/widgets/menu/menu.dart (scenario: menu-js-api)
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolGetSet, expectNoThrow, subscribeAll, wait} from '../helpers';

category('AI: Widgets: Menu JS API', () => {
  test('popup builds fluent chain: item/separator/group/endGroup/colorPalette/fontEditor/col selectors', async () => {
    const df = demog();
    const m = DG.Menu.popup();
    expect(m instanceof DG.Menu, true);
    expect(m.item('A', () => {}) instanceof DG.Menu, true);
    expect(m.separator() instanceof DG.Menu, true);
    // group() returns the nested sub-menu; endGroup() returns its parent (root).
    const sub = m.group('Sub');
    expect(sub instanceof DG.Menu, true);
    expect(sub.item('inner', () => {}) instanceof DG.Menu, true);
    expect(sub.endGroup() instanceof DG.Menu, true);
    expect(m.colorPalette([[0xFF0000FF, 0xFF00FF00, 0xFFFF0000]], {onSelect: () => {}}) instanceof DG.Menu, true);
    expect(m.fontEditor('Arial', {fontSizeMin: 6, fontSizeMax: 72, fontSizeStep: 1,
      fontFamilies: ['Arial'], asGroup: '', onChange: () => {}}) instanceof DG.Menu, true);
    const single = m.singleColumnSelector(df, {columnFilter: (c) => c.isNumerical, onChange: () => {}});
    expect(single instanceof DG.Menu, true);
    expect(m.multiColumnSelector(df, {onChange: () => {}}) instanceof DG.Menu, true);
  });

  test('closeOnClick get/set round-trip', async () => {
    expectBoolGetSet(DG.Menu.popup(), 'closeOnClick');
  });

  test('show then hide lifecycle on a host element', async () => {
    const host = DG.Menu.popup();
    expect(host instanceof DG.Menu, true);
    const m = DG.Menu.popup().item('Action', () => {});
    const div = document.createElement('div');
    document.body.appendChild(div);
    try {
      expect(m.bind(div) instanceof DG.Menu, true);
      m.show({element: div});
      await wait(50);
      expect(m.root instanceof HTMLElement, true);
    } finally {
      expectNoThrow(() => m.hide());
      div.remove();
    }
  });

  test('onContextMenuItemClick is an rxjs Observable', async () => {
    const m = DG.Menu.popup();
    subscribeAll([m.onContextMenuItemClick])();
  });

  test('Balloon.warning does not throw', async () => {
    const b = new DG.Balloon();
    expectNoThrow(() => b.warning('AI test warning'));
    await wait(50);
    expectNoThrow(() => DG.Balloon.closeAll());
  });

  test('boundary: empty popup show/hide + double-hide no-op', async () => {
    const m = DG.Menu.popup();
    expectNoThrow(() => {
      m.show();
      m.hide();
      m.hide();
    });
  });
}, {owner: 'agolovko@datagrok.ai'});

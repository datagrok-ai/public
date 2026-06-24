import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

category('Flow: minimap', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('minimap draws a rect per node plus the viewport rectangle', async () => {
    const e = makeEditor();
    try {
      expect(e.container.querySelector('.ff-minimap') != null, true, 'minimap mounted');
      await addNode(e.flow, 'Constants/String', 0, 0);
      await addNode(e.flow, 'Utilities/ToString', 300, 0);
      // Redraw is rAF-coalesced — poll until both node rects are drawn.
      const drawn = await until(() => e.container.querySelectorAll('.ff-minimap-node').length === 2);
      expect(drawn, true, 'one minimap rect per node');
      expect(e.container.querySelector('.ff-minimap-viewport') != null, true, 'viewport rectangle drawn');
    } finally {
      destroyEditor(e);
    }
  });

  test('setMinimapCollapsed collapses and restores the minimap', async () => {
    const e = makeEditor();
    try {
      const mm = e.container.querySelector('.ff-minimap') as HTMLElement;
      expect(mm != null, true);
      expect(mm.dataset.collapsed, 'false', 'starts expanded');
      e.flow.setMinimapCollapsed(true);
      expect(mm.dataset.collapsed, 'true');
      e.flow.setMinimapCollapsed(false);
      expect(mm.dataset.collapsed, 'false');
      // The minimap is never hidden outright — only minimized.
      expect(mm.style.display !== 'none', true, 'minimap is always present');
    } finally {
      destroyEditor(e);
    }
  });

  test('clicking the minimap header toggles collapse', async () => {
    const e = makeEditor();
    try {
      await addNode(e.flow, 'Constants/String', 0, 0);
      const mm = e.container.querySelector('.ff-minimap') as HTMLElement;
      const header = mm.querySelector('.ff-minimap-header') as HTMLElement;
      header.click();
      expect(mm.dataset.collapsed, 'true', 'collapsed after clicking the header');
      header.click();
      expect(mm.dataset.collapsed, 'false', 'restored on second header click');
    } finally {
      destroyEditor(e);
    }
  });
});

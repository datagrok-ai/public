import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow, uniqueName} from '../helpers';

category('AI: App: Top Menu', () => {
  test('topMenu getter returns a Menu with a root element', async () => {
    const m = grok.shell.topMenu;
    expect(m instanceof DG.Menu, true);
    expect(m.root instanceof HTMLElement, true);
  });

  test('group/item added to topMenu is findable, then removable', async () => {
    const m = grok.shell.topMenu;
    expect(m.find(uniqueName('AI-TopMenu-missing')) == null, true);
    const name = uniqueName('AI-TopMenu');
    try {
      m.group(name).item('Ping', () => {});
      expect(m.find(name) instanceof DG.Menu, true);
    } finally {
      expectNoThrow(() => m.remove(name));
    }
    expect(m.find(name) == null, true);
  });
}, {owner: 'agolovko@datagrok.ai'});

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow} from '../helpers';

// DG.BrowsePanel — navigation tree accessors + item-tooltip binding, reached via grok.shell.browsePanel.
category('AI: App: BrowsePanel JS API', () => {
  test('grok.shell.browsePanel returns a BrowsePanel instance', async () => {
    const bp = grok.shell.browsePanel;
    expect(bp != null, true);
    expect(bp instanceof DG.BrowsePanel, true);
  });

  test('localTree returns a TreeViewGroup with an accessible children array', async () => {
    const tree = grok.shell.browsePanel.localTree;
    expect(tree != null, true);
    expect(tree instanceof DG.TreeViewGroup, true);
    let children: DG.TreeViewNode[] = [];
    expectNoThrow(() => {
      children = tree.children;
    });
    expect(Array.isArray(children), true);
  });

  test('mainTree returns a TreeViewGroup', async () => {
    const tree = grok.shell.browsePanel.mainTree;
    expect(tree != null, true);
    expect(tree instanceof DG.TreeViewGroup, true);
  });

  test('localTree and mainTree are distinct objects', async () => {
    const bp = grok.shell.browsePanel;
    // Wrappers are re-created per call, so === is unreliable; compare stable DOM roots.
    expect(bp.localTree.root !== bp.mainTree.root, true);
  });

  test('bindItemTooltip with string content does not throw', async () => {
    const el = ui.div('hover-target');
    try {
      expectNoThrow(() => grok.shell.browsePanel.bindItemTooltip('Tooltip text', el));
    } finally {
      el.remove();
    }
  });

  test('expandPath expands the groups along a path of node texts without selecting', async () => {
    const tree = grok.shell.browsePanel.mainTree;
    const before = tree.currentItem?.text;
    // Settings: always there, its sections loaded on the first expand (Domains is Beta-gated)
    await grok.shell.browsePanel.expandPath('Platform/Settings');
    const platform = tree.children.find((n) => n.text === 'Platform') as DG.TreeViewGroup;
    expect(platform?.expanded, true, 'Platform expanded');
    const settings = platform.children.find((n) => n.text === 'Settings') as DG.TreeViewGroup;
    expect(settings?.expanded, true, 'Settings expanded');
    for (let i = 0; i < 50 && settings.children.length === 0; i++)
      await delay(100);
    expect(settings.children.length > 0, true, 'Settings loaded its sections');
    expect(tree.currentItem?.text, before, 'the selection is untouched');
    let threw = false;
    try {
      await grok.shell.browsePanel.expandPath('Platform/No such node');
    } catch {
      threw = true;
    }
    expect(threw, false, 'an unknown path is ignored');
  });

  test('bindItemTooltip with a function content provider does not throw', async () => {
    const el = ui.div('hover-target');
    try {
      expectNoThrow(() => grok.shell.browsePanel.bindItemTooltip(() => ui.div('lazy tooltip'), el));
    } finally {
      el.remove();
    }
  });
}, {owner: 'agolovko@datagrok.ai'});

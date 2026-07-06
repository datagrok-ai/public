import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
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

  test('bindItemTooltip with a function content provider does not throw', async () => {
    const el = ui.div('hover-target');
    try {
      expectNoThrow(() => grok.shell.browsePanel.bindItemTooltip(() => ui.div('lazy tooltip'), el));
    } finally {
      el.remove();
    }
  });
}, {owner: 'agolovko@datagrok.ai'});

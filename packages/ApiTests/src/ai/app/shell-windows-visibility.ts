import * as grok from 'datagrok-api/grok';
import {category, test} from '@datagrok-libraries/test/src/test';
import {expectBoolGetSet} from '../helpers';

category('AI: App: Shell Windows Visibility', () => {
  test('window visibility flags round-trip', async () => {
    // If a single flag's getter doesn't mirror its setter headless, drop that flag — the
    // round-trip is the contract under test.
    const flags = ['showSidebar', 'showBrowse', 'showProjects', 'showFavorites',
      'showRibbon', 'autoShowToolbox', 'showHelp', 'showContextPanel',
      'presentationMode', 'simpleMode', 'hideTabsInPresentationMode'];
    for (const f of flags)
      expectBoolGetSet(grok.shell.windows, f);
  });
}, {owner: 'agolovko@datagrok.ai'});

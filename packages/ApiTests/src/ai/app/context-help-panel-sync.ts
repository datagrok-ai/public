import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolGetSet} from '../helpers';

// grok.shell.o get/set + setCurrentObject, plus help/context sub-panels: syncCurrentObject get/set
// and visible get/set (delegating to showHelp / showContextPanel).
category('AI: App: ContextPanel / HelpPanel Sync', () => {
  test('help.syncCurrentObject get/set round-trips', async () => {
    expect(typeof grok.shell.windows.help.syncCurrentObject, 'boolean');
    expectBoolGetSet(grok.shell.windows.help, 'syncCurrentObject');
  });

  test('context.syncCurrentObject get/set round-trips', async () => {
    expect(typeof grok.shell.windows.context.syncCurrentObject, 'boolean');
    expectBoolGetSet(grok.shell.windows.context, 'syncCurrentObject');
  });

  test('help.visible get/set round-trips (delegates to showHelp)', async () => {
    const original = grok.shell.windows.help.visible;
    try {
      expect(typeof grok.shell.windows.help.visible, 'boolean');
      // visible is a thin delegate to windows.showHelp, so they must stay in sync.
      grok.shell.windows.help.visible = true;
      expect(grok.shell.windows.help.visible, true);
      expect(grok.shell.windows.showHelp, true);
      grok.shell.windows.help.visible = false;
      expect(grok.shell.windows.help.visible, false);
      expect(grok.shell.windows.showHelp, false);
    } finally {
      grok.shell.windows.help.visible = original;
      expect(grok.shell.windows.help.visible, original);
    }
  });

  test('context.visible get/set round-trips (delegates to showContextPanel)', async () => {
    const original = grok.shell.windows.context.visible;
    try {
      expect(typeof grok.shell.windows.context.visible, 'boolean');
      // visible is a thin delegate to windows.showContextPanel, so they must stay in sync.
      grok.shell.windows.context.visible = true;
      expect(grok.shell.windows.context.visible, true);
      expect(grok.shell.windows.showContextPanel, true);
      grok.shell.windows.context.visible = false;
      expect(grok.shell.windows.context.visible, false);
      expect(grok.shell.windows.showContextPanel, false);
    } finally {
      grok.shell.windows.context.visible = original;
      expect(grok.shell.windows.context.visible, original);
    }
  });

  test('shell.o round-trip: set DataFrame, read back as DataFrame', async () => {
    const original = grok.shell.o;
    const df = demog(20);
    try {
      grok.shell.o = df;
      const read = grok.shell.o;
      expect(read instanceof DG.DataFrame, true);
      expect((read as DG.DataFrame).id, df.id);
    } finally {
      grok.shell.o = original;
    }
  });

  // Dropped: setting grok.shell.o to a plain string does not persist as the current object
  // (the platform keeps the prior entity/DataFrame); plain strings are not valid current objects.

  test('setCurrentObject(x, false) sets shell.o without locking', async () => {
    const original = grok.shell.o;
    const df = demog(20);
    try {
      // freeze=false: the object becomes current but is not pinned/frozen.
      grok.shell.setCurrentObject(df, false);
      const read = grok.shell.o;
      expect(read instanceof DG.DataFrame, true);
      expect((read as DG.DataFrame).id, df.id);
    } finally {
      grok.shell.o = original;
    }
  });
}, {owner: 'agolovko@datagrok.ai'});

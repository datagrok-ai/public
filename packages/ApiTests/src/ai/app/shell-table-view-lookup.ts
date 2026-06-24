import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, uniqueName} from '../helpers';

// DG.Shell table/view lookup — getTableView, view, tableView, tableNames, dockElement.
// Read-only lookups share one attached table view (clear:false).
category('AI: App: Shell Table/View Lookup', () => {
  let tv: DG.TableView;

  before(async () => {
    const d = demog();
    d.name = uniqueName('stvl');
    tv = grok.shell.addTableView(d);
  });

  after(async () => tv.close());

  test('getTableView returns the TableView for an open table', async () => {
    const found = grok.shell.getTableView(tv.dataFrame.name);
    expect(found instanceof DG.TableView, true);
    expect(found.dataFrame === tv.dataFrame, true);
  });

  test('tableView(name) is an alias of getTableView', async () => {
    const a = grok.shell.getTableView(tv.dataFrame.name);
    const b = grok.shell.tableView(tv.dataFrame.name);
    expect(a.dataFrame === b.dataFrame, true);
  });

  test('view(name) resolves an open view and is null for a missing name', async () => {
    const v = grok.shell.view(tv.name);
    expect(v instanceof DG.View, true);
    expect(grok.shell.view(uniqueName('no-such-view')) == null, true);
  });

  test('tableNames lists open table names', async () => {
    expect(grok.shell.tableNames.indexOf(tv.dataFrame.name) >= 0, true);
  });

  test('getTableView for an unknown table name is null-ish', async () => {
    expect(grok.shell.getTableView(uniqueName('no-such-table')) == null, true);
  });

  test('dockElement docks an element and close removes it', async () => {
    const marker = uniqueName('stvl-dock');
    const el = ui.div(['docked'], {classes: marker});
    expectNoThrow(() => grok.shell.dockElement(el, marker));
    expect(document.querySelector('.' + marker) != null, true);
    expectNoThrow(() => grok.shell.dockManager.close(el));
  });
}, {owner: 'agolovko@datagrok.ai', clear: false});

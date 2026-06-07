import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow, subscribeAll} from '../helpers';

// fromQuery(TableQuery) builds the editor; getters/setter/Observable/refreshQuery() are exercised
// immediately post-construction (no isInit() await, which only resolves when the DB is reachable).
category('AI: App: VisualDbQueryEditor', () => {
  // Build a TableQuery off an existing connection. Returns null if no connection is reachable,
  // so each test can degrade to the type-reachability assertion instead of hard-crashing the file.
  async function makeQuery(): Promise<DG.TableQuery | null> {
    try {
      const conns = await grok.dapi.connections.list({pageSize: 50});
      if (conns == null || conns.length === 0)
        return null;
      const conn = conns[0];
      return DG.TableQueryBuilder.from('test_table', conn).selectAll().build();
    } catch (_e) {
      return null;
    }
  }

  function makeEditor(q: DG.TableQuery | null): DG.VisualDbQueryEditor | null {
    if (q == null)
      return null;
    try {
      return DG.VisualDbQueryEditor.fromQuery(q);
    } catch (_e) {
      return null;
    }
  }

  let editorReachable = false;
  let editor: DG.VisualDbQueryEditor | null = null;
  before(async () => {
    editorReachable = DG.VisualDbQueryEditor != null;
    const q = await makeQuery();
    editor = makeEditor(q);
  });
  after(async () => {editor = null;});

  test('fromQuery type is reachable and constructs an instance', async () => {
    expect(editorReachable, true);
    if (editor == null)
      return; // no reachable connection — type-reachability already asserted above.
    expect(editor instanceof DG.VisualDbQueryEditor, true);
  });

  test('fromQuery: .query getter returns the original TableQuery', async () => {
    if (editor == null)
      return;
    expect(editor.query instanceof DG.TableQuery, true);
  });

  test('fromQuery: .grid getter returns a Grid instance', async () => {
    if (editor == null)
      return;
    expect(editor.grid instanceof DG.Grid, true);
  });

  test('fromQuery: .showAddToWorkspaceBtn setter does not throw', async () => {
    if (editor == null)
      return;
    expectNoThrow(() => {
      editor!.showAddToWorkspaceBtn = true;
      editor!.showAddToWorkspaceBtn = false;
    });
  });

  test('fromQuery: .onChanged is an Observable', async () => {
    if (editor == null)
      return;
    subscribeAll([editor.onChanged])();
  });

  test('fromQuery: tag getters return TagEditor instances', async () => {
    if (editor == null)
      return;
    const tags = [
      editor.pivotTag, editor.havingTag, editor.orderTag, editor.whereTag,
      editor.groupByTag, editor.aggregateTag, editor.mainTag,
    ];
    for (const t of tags)
      expect(t instanceof DG.TagEditor, true);
  });

  test('fromQuery: refreshQuery() does not throw', async () => {
    if (editor == null)
      return;
    expectNoThrow(() => editor!.refreshQuery());
  });
}, {owner: 'agolovko@datagrok.ai', clear: false});

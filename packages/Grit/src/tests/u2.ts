import * as grok from 'datagrok-api/grok';
import {before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import type {DomainSource} from '@datagrok-libraries/u2';
import {DomainApp} from '@datagrok-libraries/u2/src/dg/index.js';
import type {DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';
import {gritDb, IssueRow} from '../generated/db';
import type {GritDb} from '../generated/db-ui';
import {IssuesApp, openGrit} from '../package';

// The code tier over the generated typed handles: the Issues app opens over `getGritDb()`, the
// table handle is a `DomainTable<IssueRow>`, the actions show under their `when`, the validator
// refuses closing an unassigned issue. Nothing reaches the server: the rows probed are drafts.
category('Grit: issues app', () => {
  let db: GritDb;
  let issues: DomainTable<IssueRow>;

  /** The first load settled, or the reason it did not. */
  async function ready(src: DomainSource): Promise<void> {
    for (let i = 0; i < 200 && (src.state.value === 'idle' || src.state.value === 'loading'); i++)
      await delay(50);
    if (src.state.value !== 'ready')
      throw new Error(`source is ${src.state.value}: ${String(src.error.value)}`);
  }

  before(async () => {
    db = await openGrit();
    issues = db.tables.issues;
  });

  test('handles: one await, every table typed, the data clients beside them', async () => {
    expect(db.schema, 'grit');
    expect(db.data, gritDb);
    expect(issues.address, 'grit.issue');
    expect(db.tables.comments.address, 'grit.comment');
    expect(await openGrit(), db, 'cached per page');
    expect(issues.access.field('title'), 'editable');
    expect(issues.access.field('number'), 'readonly', 'autoNumber is engine-assigned');
  });

  test('app: opens over the handle as an IssuesApp — list page, presets, shortcuts, the entity page', async () => {
    const view = issues.app({name: 'Issues (test)', path: '/apps/Grit/IssuesTest', app: IssuesApp,
      children: {tables: ['comment']}});
    grok.shell.addView(view);
    try {
      const app = DomainApp.of(view)!;
      expect(app instanceof IssuesApp, true);
      expect(app.page.value, 'list');
      expect(app.shortcuts['m'], 'Assign to me');
      expect(app.shortcuts['c'], 'Close');
      const ribbon = app.ribbon();
      expect(ribbon.main.length, 5, 'New, Save, Discard, the ⋯ menu, Refresh');
      expect(ribbon.tools.length, 1, 'the search box');
      expect(ribbon.presets.length, 1, 'the Mine / Open switch');
      const presets = ribbon.presets[0];
      expect(('root' in presets ? presets.root : presets).dataset.u2, 'domain-presets');
      expect(app.ribbonGroups().length, 3, 'the positional shape appView takes');
      expect(app.ribbon(), ribbon, 'built once');
      expect(await app.open('?entity=new'), true);
      expect(app.page.value, 'entity');
      expect(app.entity.value, DomainApp.NEW);
      expect(await app.open(''), true, 'a pristine draft leaves without a prompt');
      expect(app.page.value, 'list');
    } finally {
      view.close();
    }
  });

  test('actions: Assign to me only for a row that is not mine; Close only for an open issue', async () => {
    const draft = issues.draft({title: 'u2 action probe'});
    try {
      await ready(draft);
      const row = draft.currentRow.value!;
      const names = () => issues.actions.for(row).map((a) => a.name);
      expect(names().includes('Assign to me'), true);
      issues.actions.for(row).find((a) => a.name === 'Assign to me')!.run();
      expect(row.assignee, grok.shell.user.id);
      expect(names().includes('Assign to me'), false, 'hidden once the row is mine');
      const closed = await gritDb.statuses.getByKey({name: 'closed'});
      if (closed === null)
        return;
      expect(names().includes('Close'), true);
      issues.actions.for(row).find((a) => a.name === 'Close')!.run();
      expect(row.status_id, closed.id);
      expect(names().includes('Close'), false, 'hidden once closed');
    } finally {
      draft.dispose();
    }
  });

  test('validator: closing an unassigned issue is refused', async () => {
    const closed = await gritDb.statuses.getByKey({name: 'closed'});
    if (closed === null)
      return;
    const draft = issues.draft({title: 'u2 validator probe'});
    try {
      await ready(draft);
      const row = draft.currentRow.value!;
      expect(issues.validators.check('status_id', closed.id, row), 'Assign before closing');
      row.assignee = grok.shell.user.id;
      expect(issues.validators.check('status_id', closed.id, row), null);
    } finally {
      draft.dispose();
    }
  });
});

import * as grok from 'datagrok-api/grok';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import type {DomainSource} from '@datagrok-libraries/u2';
import {domains, domainForm, domainList, DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';

const PROJECTS = 'grit.project';

// The u2 domain stack over Grit's own schema: a project created through `domainForm` over a
// pristine draft, listed by `domainList`, edited through `currentRow`, and the access the
// controls degrade by. One project per run, keyed uniquely and deleted afterwards.
category('Grit: u2 form and list', () => {
  const client = () => grok.dapi.domains.table(PROJECTS);
  const key = `U2${Date.now() % 1e10}${Math.floor(Math.random() * 1e3)}`;
  const query = `key = "${key}"`;
  let projects: DomainTable;
  let id: string;

  /** The first load settled, or the reason it did not. */
  async function ready(src: DomainSource): Promise<void> {
    for (let i = 0; i < 200 && (src.state.value === 'idle' || src.state.value === 'loading'); i++)
      await delay(50);
    if (src.state.value !== 'ready')
      throw new Error(`source is ${src.state.value}: ${String(src.error.value)}`);
  }

  before(async () => {
    projects = await domains.table(PROJECTS);
  });

  after(async () => {
    await client().deleteWhere(query);
  });

  test('create: a project through the form over a draft, inserted by save', async () => {
    const draft = projects.draft();
    const form = domainForm(draft);
    try {
      await ready(draft);
      expect(form.form !== null, true, 'the draft is the form\'s row');
      expect(draft.isDirty.value, false, 'pristine until touched');
      form.input('key')!.value.value = key;
      form.input('name')!.value.value = 'u2 project';
      await delay(10);
      expect(draft.isDirty.value, true);
      expect(await draft.save(), true, String(draft.error.value));
      const [row] = await client().query({filter: query});
      expect(row?.name, 'u2 project');
      id = row.id;
    } finally {
      form.dispose();
      draft.dispose();
    }
  });

  test('list: shows the row; edit through currentRow, save, discard', async () => {
    const src = projects.source({query});
    const list = domainList(src, {mode: 'cards'});
    try {
      await ready(src);
      expect(src.rows.items.value.length, 1);
      expect(src.rows.byKey(id)?.name, 'u2 project');
      list.list.selectedIndex.value = 0;
      await delay(10);
      expect(src.currentRow.value?.id, id, 'the selection is the current row');
      src.currentRow.value!.name = 'u2 project 2';
      expect(src.isDirty.value, true);
      expect(await src.save(), true, String(src.error.value));
      expect((await client().get(id)).name, 'u2 project 2');
      src.currentRow.value!.name = 'zzz';
      expect(src.isDirty.value, true);
      src.discard();
      expect(src.isDirty.value, false);
      expect(src.rows.byKey(id)!.name, 'u2 project 2', 'discard restores the cell');
    } finally {
      list.dispose();
      src.dispose();
    }
  });

  test('access: system and autoNumber columns are readonly, declared columns editable', async () => {
    expect(projects.access.can('insert'), true);
    expect(projects.access.field('name'), 'editable');
    expect(projects.access.field('id'), 'readonly');
    expect(projects.access.field('version'), 'readonly');
    expect(projects.access.field('nosuch'), 'hidden');
    const issues = await domains.table('grit.issue');
    expect(issues.access.field('title'), 'editable');
    expect(issues.access.field('created_on'), 'readonly');
    expect(issues.access.field('number'), 'readonly', 'autoNumber is engine-assigned');
  });
});

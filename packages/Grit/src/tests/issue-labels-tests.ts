import type * as _DG from 'datagrok-api/dg';
declare let DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {gritDb, IssueColumn, IssueUpdate} from '../generated/db';

// The many-to-many dogfood (GROK-20298): `issue.labels` is declared once in
// databases/grit/schema.json — {"via": "issue_label", "target": "label"}, with both
// junction FKs auto-resolved — and everything below comes from that one line, with no
// Grit code behind it. The engine itself is covered by the datlas suite; what these
// tests prove is the TYPED surface generated from Grit's real schema: expand, relation
// filter paths in both filter forms, inline set-semantics writes, the create-and-link
// transaction, and the batch refusal.
category('Grit: issue labels', () => {
  const marker = `LBL${Date.now() % 1e10}${Math.floor(Math.random() * 1e3)}`;
  const issueIds: string[] = [];
  const labelIds: string[] = [];
  // A project and three labels per TEST, torn down again: the filter assertions
  // are exact memberships within one project, so a leftover issue from another
  // test would make them order-dependent.
  let projectId: string | null = null;
  let bug: string;
  let urgent: string;
  let spare: string;

  // Display names of one issue's labels, as the expand returns them. The builder form
  // is what types `labels` onto the row — the point of generating the client at all.
  const linked = async (id: string): Promise<string[]> => {
    const rows = await gritDb.issues.query().where('id', '=', id).expand('labels');
    return (rows[0].labels ?? []).map((l) => l.name);
  };

  const insertIssue = async (number: number, title: string, labels?: string[]): Promise<string> => {
    const [row] = await gritDb.issues.insert(
      {project_id: projectId!, number: number, title: title, labels: labels});
    issueIds.push(row.id);
    return row.id;
  };

  // Both `project.key` and `label.name` are business keys, and a soft-deleted row
  // still holds its key — so every test's fixture gets its own suffix rather than
  // re-inserting (and silently deduping onto) the previous test's deleted rows.
  let run = 0;

  const setup = async (): Promise<void> => {
    const tag = `${marker}-${run++}`;
    const [project] = await gritDb.projects.insert({key: tag, name: `${tag} labels`});
    projectId = project.id;
    const ins = await gritDb.labels.insert([{name: `${tag} bug`, color: '#ff0000'},
      {name: `${tag} urgent`}, {name: `${tag} spare`}]);
    for (const r of ins)
      labelIds.push(r.id);
    [bug, urgent, spare] = labelIds;
  };

  // The display names of this run's fixture labels, for assertions.
  const labelName = (kind: string): string => `${marker}-${run - 1} ${kind}`;

  const cleanup = async (): Promise<void> => {
    for (const id of issueIds)
      await gritDb.issues.delete(id).catch(() => {}); // cascades its issue_label rows
    issueIds.length = 0;
    for (const id of labelIds)
      await gritDb.labels.delete(id).catch(() => {});
    labelIds.length = 0;
    if (projectId != null)
      await gritDb.projects.delete(projectId).catch(() => {});
    projectId = null;
  };

  test('insert with links, typed expand, and the d42 chips column', async () => {
    await setup();
    try {
      // The relation rides the ordinary insert payload as the full set of target ids.
      const id = await insertIssue(1, `${labelName('issue')} tagged`, [urgent, bug, urgent]);
      // Deduped (the junction's businessKey makes re-linking idempotent) and ordered
      // by display name, whatever order they were written in.
      expect((await linked(id)).join(','), `${labelName('bug')},${labelName('urgent')}`);
      // An issue with no labels expands to [], never null.
      const bare = await insertIssue(2, `${labelName('issue')} bare`);
      expect((await linked(bare)).length, 0);

      // d42 (what a grid binds to): display names joined by ', ' plus the hidden
      // '~labels.id' companion carrying the ids in the SAME order.
      const df = await gritDb.issues.queryDf({filter: [{property: 'id', operator: '=', value: id}],
        expand: ['labels']});
      const names = df.columns.byName('labels');
      const ids = df.columns.byName('~labels.id');
      expect(names != null && ids != null, true, `columns: ${df.columns.names().join(', ')}`);
      expect(names.get(0), `${labelName('bug')}, ${labelName('urgent')}`);
      expect(ids.get(0).split(', ')[0], bug, `id companion: ${ids.get(0)}`);
      // The header says 'Labels', not the raw relation name.
      expect(names.getTag(DG.TAGS.FRIENDLY_NAME), 'Labels');
    } finally {
      await cleanup();
    }
  });

  test('filter: relation paths in both filter forms', async () => {
    await setup();
    try {
      const tagged = await insertIssue(1, `${labelName('issue')} tagged`, [bug]);
      const other = await insertIssue(2, `${labelName('issue')} other`, [urgent]);
      const bare = await insertIssue(3, `${labelName('issue')} bare`);
      const mine: _DG.DomainConditionNode<IssueColumn> =
        {property: 'project_id', operator: '=', value: projectId!};

      // Tree form: a relation segment matches the owners that HAVE such a link.
      const byName = await gritDb.issues.query({filter: [mine, 'and',
        {property: 'labels.name', operator: '=', value: labelName('bug')}]});
      expect(byName.map((r) => r.id).join(','), tagged);

      // The smart-filter STRING form takes the same path (values quoted, as everywhere).
      const byString = await gritDb.issues.query(
        {filter: `project_id = "${projectId}" and labels.name = "${labelName('urgent')}"`});
      expect(byString.map((r) => r.id).join(','), other);

      // 'labels.id' with a list is the ANY-of shape the filter panel's checkbox list
      // emits — the reason the manifest declares {"column": "labels.id"} explicitly.
      const anyOf = await gritDb.issues.query({filter: [mine, 'and',
        {property: 'labels.id', operator: '=', value: [bug, urgent]}], sort: 'number'});
      expect(anyOf.map((r) => r.id).join(','), `${tagged},${other}`);

      // '!=' on that leaf is the uncheck gesture: EXCLUDE the owners having those
      // links, rather than "has some other label" (which would keep them).
      const excluded = await gritDb.issues.query({filter: [mine, 'and',
        {property: 'labels.id', operator: '!=', value: [bug]}], sort: 'number'});
      expect(excluded.map((r) => r.id).join(','), `${other},${bare}`);

      // '= null' is the facet's "(no value)" bucket: no visible link at all.
      const unlinked = await gritDb.issues.query({filter: [mine, 'and',
        {property: 'labels.id', operator: '=', value: null}]});
      expect(unlinked.map((r) => r.id).join(','), bare);

      // The facet behind the declared 'labels.id' filter: DISTINCT owners per label,
      // with the label's display name.
      const res = await gritDb.issues.facets({filter: [mine],
        facets: [{id: 'l', kind: 'categories', column: 'labels.id'}]});
      const byValue: {[v: string]: any} = {};
      for (const c of (res.facets['l'] as any).categories)
        byValue[c.value] = c;
      expect(byValue[bug]?.total, 1, JSON.stringify(res.facets['l']));
      expect(byValue[bug]?.display, labelName('bug'), JSON.stringify(res.facets['l']));
    } finally {
      await cleanup();
    }
  });

  test('update: the list replaces the link set, [] clears it, absent leaves it alone', async () => {
    await setup();
    try {
      const id = await insertIssue(1, `${labelName('issue')} edited`, [bug]);
      // The generated <Table>Update names the relation next to the columns.
      const values: IssueUpdate = {labels: [urgent, spare]};
      const upd = await gritDb.issues.update(id, values);
      expect((await linked(id)).join(','), `${labelName('spare')},${labelName('urgent')}`);
      // A relations-only change still bumps the row version, so expectedVersion covers
      // the links exactly like an ordinary column edit.
      expect(upd.version, 2, 'a relations-only update must bump the version');

      // An ordinary column rides along; an ABSENT relation key leaves the links alone.
      await gritDb.issues.update(id, {description: 'still labelled'});
      expect((await linked(id)).length, 2, 'an absent relation key must not clear the links');

      // [] clears the (visible) link set.
      await gritDb.issues.update(id, {labels: []});
      expect((await linked(id)).length, 0);

      // A missing or invisible target is one code — linking is no existence oracle —
      // and the whole write rolls back.
      let code = '';
      try {
        await gritDb.issues.update(id, {labels: ['3f2504e0-4f89-11d3-9a0c-0305e82c3301']});
      } catch (e: any) {
        code = e.rows?.[0]?.errors?.[0]?.code ?? `${e.message}`;
      }
      expect(code, 'not-visible-or-missing');
      expect((await linked(id)).length, 0, 'a failed link write must roll back');
    } finally {
      await cleanup();
    }
  });

  test('transaction: create the label and link it in one atomic operation', async () => {
    await setup();
    try {
      const name = `${labelName('fresh')}`;
      const [label, issue] = await gritDb.transaction([
        {op: 'insert', table: 'label', ref: 'l1', values: {name: name}},
        {op: 'insert', table: 'issue', values: {project_id: projectId!, number: 1,
          title: `${labelName('issue')} tx`, labels: ['$l1']}},
      ]);
      labelIds.push(label.id);
      issueIds.push(issue.id);
      expect((await linked(issue.id)).join(','), name);
    } finally {
      await cleanup();
    }
  });

  test('batch refuses relation keys and points at insert/update', async () => {
    await setup();
    try {
      let message = '';
      try {
        const report = await gritDb.issues.batch(
          [{project_id: projectId!, number: 1, title: `${labelName('issue')} B`, labels: [bug]}]);
        message = report.error ?? '';
      } catch (e: any) {
        message = e.message ?? `${e}`;
      }
      expect(message.includes('not supported in batch mode'), true, `unexpected outcome: '${message}'`);
      expect(message.includes('insert or patch'), true, `unexpected outcome: '${message}'`);
    } finally {
      await cleanup();
    }
  });
}, {owner: 'askalkin@datagrok.ai'});

import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {gritDb, CommentRow, IssueRow} from '../generated/db';

// End-to-end issue CRUD over the typed clients generated from databases/grit/schema.json.
category('Grit: issue CRUD', () => {
  const unique = (prefix: string) => `${prefix}${Date.now() % 1e10}${Math.floor(Math.random() * 1e3)}`;

  const throws = async (action: () => Promise<any>): Promise<boolean> => {
    try {
      await action();
    } catch (e) {
      return true;
    }
    return false;
  };

  // Upsert-by-name id of a lookup row (status/priority/issue_type): businessKey ["name"]
  // reports a duplicate with the existing id, so this is idempotent. The canonical rows
  // ('open', 'high', ...) are the schema's seed data shared with the ingestion script —
  // deliberately NOT cleaned up.
  const lookupId = async (client: {insert: (values: any) => Promise<any[]>},
    values: {name: string, [key: string]: any}): Promise<string> => {
    const [r] = await client.insert(values);
    return r.existingId ?? r.id;
  };

  test('project and issue lifecycle', async () => {
    const me = await grok.dapi.users.current();
    const open = await lookupId(gritDb.statuses, {name: 'open'});
    const inProgress = await lookupId(gritDb.statuses, {name: 'in progress'});
    const high = await lookupId(gritDb.priorities, {name: 'high', priority: 3});
    const bug = await lookupId(gritDb.issueTypes, {name: 'bug'});
    const [project] = await gritDb.projects.insert({key: unique('T'), name: 'Test project'});
    const [issue] = await gritDb.issues.insert({
      project_id: project.id, number: 1, title: 'First bug',
      status_id: open, priority_id: high, type_id: bug, reporter: me.id, assignee: me.id,
    });
    try {
      const row = await gritDb.issues.get(issue.id);
      expect(row.status_id, open, 'status ref not persisted');
      expect(row.reporter, me.id, 'user column not persisted');

      await gritDb.issues.update(issue.id, {status_id: inProgress}, {version: 1});

      const [comment] = await gritDb.comments.insert({issue_id: issue.id, text: 'On it'});
      expect(comment.created, true);
      // one expand query returns the issue with its comments embedded (what the app renders from)
      const expanded: (IssueRow & {comment?: CommentRow[]})[] = await gritDb.issues.query(
        {filter: `project_id = "${project.id}" and number = 1`, expand: ['details:comment']});
      expect(expanded.length, 1);
      expect(expanded[0].comment?.length, 1, 'expanded payload must embed the comment');
      expect(expanded[0].comment![0].text, 'On it');

      // ref expand resolves the lookup name in the same query (what badges/grids show)
      const [withNames] = await gritDb.issues.query(
        {filter: `id = "${issue.id}"`, expand: ['status_id']});
      expect((withNames as any)['status_id.name'], 'in progress');

      const audit = await gritDb.issues.audit(issue.id);
      expect(audit.map((a) => a.op).join(','), 'insert,update');
      expect(audit[1].after.status_id, inProgress);

      expect(await throws(() => gritDb.projects.delete(project.id)), true,
        'deleting a project with live issues must be restricted');
      expect(await throws(() => gritDb.statuses.delete(inProgress)), true,
        'deleting a status referenced by live issues must be restricted');
    } finally {
      await gritDb.issues.delete(issue.id);
      await gritDb.projects.delete(project.id);
    }
    expect((await gritDb.comments.query({filter: `issue_id = "${issue.id}"`})).length, 0,
      'comments must cascade with their issue');
  });

  test('per-project issue numbering is deduplicated', async () => {
    const [project] = await gritDb.projects.insert({key: unique('T'), name: 'Numbering'});
    const [first] = await gritDb.issues.insert({project_id: project.id, number: 1, title: 'A'});
    try {
      const [dup] = await gritDb.issues.insert({project_id: project.id, number: 1, title: 'B'});
      expect(dup.status, 'duplicate');
      expect(dup.existingId, first.id);
    } finally {
      await gritDb.issues.delete(first.id);
      await gritDb.projects.delete(project.id);
    }
  });

  test('labels: N:N join cascades with the issue', async () => {
    const [project] = await gritDb.projects.insert({key: unique('T'), name: 'Labels'});
    const [issue] = await gritDb.issues.insert({project_id: project.id, number: 1, title: 'Labeled'});
    const [label] = await gritDb.labels.insert({name: unique('bug'), color: '#ff0000'});
    try {
      await gritDb.issueLabels.insert({issue_id: issue.id, label_id: label.id});
      expect((await gritDb.issueLabels.query({filter: `issue_id = "${issue.id}"`})).length, 1);

      const [dup] = await gritDb.issueLabels.insert({issue_id: issue.id, label_id: label.id});
      expect(dup.status, 'duplicate', 'the same label attached twice');

      await gritDb.issues.delete(issue.id);
      expect((await gritDb.issueLabels.query({filter: `issue_id = "${issue.id}"`})).length, 0,
        'issue_label rows must cascade with their issue');
    } finally {
      await gritDb.labels.delete(label.id);
      await gritDb.projects.delete(project.id);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});

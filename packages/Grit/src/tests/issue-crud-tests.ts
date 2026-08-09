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

  test('project and issue lifecycle', async () => {
    const me = await grok.dapi.users.current();
    const [project] = await gritDb.projects.insert({key: unique('T'), name: 'Test project'});
    const [issue] = await gritDb.issues.insert({
      project_id: project.id, number: 1, title: 'First bug',
      priority: 'high', reporter: me.id, assignee: me.id,
    });
    try {
      const row = await gritDb.issues.get(issue.id);
      expect(row.status, 'open', 'choices default not applied');
      expect(row.reporter, me.id, 'user column not persisted');

      await gritDb.issues.update(issue.id, {status: 'in progress'}, {version: 1});

      const [comment] = await gritDb.comments.insert({issue_id: issue.id, text: 'On it'});
      expect(comment.created, true);
      // one expand query returns the issue with its comments embedded (what the app renders from)
      const expanded: (IssueRow & {comment?: CommentRow[]})[] = await gritDb.issues.query(
        {filter: `project_id = "${project.id}" and number = 1`, expand: ['details:comment']});
      expect(expanded.length, 1);
      expect(expanded[0].comment?.length, 1, 'expanded payload must embed the comment');
      expect(expanded[0].comment![0].text, 'On it');

      const audit = await gritDb.issues.audit(issue.id);
      expect(audit.map((a) => a.op).join(','), 'insert,update');
      expect(audit[1].after.status, 'in progress');

      expect(await throws(() => gritDb.projects.delete(project.id)), true,
        'deleting a project with live issues must be restricted');
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

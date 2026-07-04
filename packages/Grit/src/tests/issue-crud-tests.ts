import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// End-to-end issue CRUD over the 'grit' domain schema (databases/grit/schema.json).
category('Grit: issue CRUD', () => {
  const projects = () => grok.dapi.domains.table('grit.project');
  const issues = () => grok.dapi.domains.table('grit.issue');
  const comments = () => grok.dapi.domains.table('grit.comment');
  const labels = () => grok.dapi.domains.table('grit.label');
  const issueLabels = () => grok.dapi.domains.table('grit.issue_label');
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
    const [project] = await projects().insert({key: unique('T'), name: 'Test project'});
    const [issue] = await issues().insert({
      project_id: project.id, number: 1, title: 'First bug',
      priority: 'high', reporter: me.id, assignee: me.id,
    });
    try {
      const row = await issues().get(issue.id);
      expect(row.status, 'open', 'choices default not applied');
      expect(row.reporter, me.id, 'user column not persisted');

      await issues().update(issue.id, {status: 'in progress'}, {version: 1});

      const [comment] = await comments().insert({issue_id: issue.id, text: 'On it'});
      expect(comment.created, true);
      const commentRows = await comments().query({filter: `issue_id = "${issue.id}"`});
      expect(commentRows.length, 1);
      expect(commentRows[0].text, 'On it');

      const audit = await issues().audit(issue.id);
      expect(audit.map((a) => a.op).join(','), 'insert,update');
      expect(audit[1].after.status, 'in progress');

      expect(await throws(() => projects().delete(project.id)), true,
        'deleting a project with live issues must be restricted');
    } finally {
      await issues().delete(issue.id);
      await projects().delete(project.id);
    }
    expect((await comments().query({filter: `issue_id = "${issue.id}"`})).length, 0,
      'comments must cascade with their issue');
  });

  test('per-project issue numbering is deduplicated', async () => {
    const [project] = await projects().insert({key: unique('T'), name: 'Numbering'});
    const [first] = await issues().insert({project_id: project.id, number: 1, title: 'A'});
    try {
      const [dup] = await issues().insert({project_id: project.id, number: 1, title: 'B'});
      expect(dup.status, 'duplicate');
      expect(dup.existingId, first.id);
    } finally {
      await issues().delete(first.id);
      await projects().delete(project.id);
    }
  });

  test('labels: N:N join cascades with the issue', async () => {
    const [project] = await projects().insert({key: unique('T'), name: 'Labels'});
    const [issue] = await issues().insert({project_id: project.id, number: 1, title: 'Labeled'});
    const [label] = await labels().insert({name: unique('bug'), color: '#ff0000'});
    try {
      await issueLabels().insert({issue_id: issue.id, label_id: label.id});
      expect((await issueLabels().query({filter: `issue_id = "${issue.id}"`})).length, 1);

      const [dup] = await issueLabels().insert({issue_id: issue.id, label_id: label.id});
      expect(dup.status, 'duplicate', 'the same label attached twice');

      await issues().delete(issue.id);
      expect((await issueLabels().query({filter: `issue_id = "${issue.id}"`})).length, 0,
        'issue_label rows must cascade with their issue');
    } finally {
      await labels().delete(label.id);
      await projects().delete(project.id);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});

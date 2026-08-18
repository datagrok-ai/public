//name: IngestGitHubIssues
//description: Ingests files/github-issues.json (all tickets of github.com/datagrok-ai/public) into the grit domain schema through the server Domain API — a worked example of grok.dapi.domains: lookup-table seeding via business-key upserts, bulk issue upload via batch(mode: 'upsert'), and a query round-trip to wire up the N:N labels
//language: javascript

// The data file is fetched from the GitHub API ahead of time and shipped with the
// package; re-running this script is safe — every write below merges by business key.
const data = JSON.parse(await grok.dapi.files.readAsText('System:AppData/Grit/github-issues.json'));
grok.shell.info(`Ingesting ${data.issues.length} tickets of ${data.repo}...`);

const table = (name) => grok.dapi.domains.table('grit.' + name);
const idOf = (report) => report.existingId ?? report.id;

// 1. Lookup tables (status / priority / issue_type). insert() on a table with a
// business key reports {status: 'duplicate', existingId} instead of failing, so
// this doubles as a get-or-create and yields the name → id maps we join by.
const seedLookup = async (name, rows) => {
  const reports = await table(name).insert(rows);
  return new Map(rows.map((row, i) => [row.name, idOf(reports[i])]));
};

const statusIds = await seedLookup('status', [
  {name: 'open', description: 'Ready to be picked up'},
  {name: 'in progress', description: 'Actively being worked on'},
  {name: 'resolved', description: 'Fixed, awaiting verification'},
  {name: 'closed', description: 'Done, or will not be fixed'},
]);
const priorityIds = await seedLookup('priority', [
  {name: 'low', description: 'Nice to have', priority: 1},
  {name: 'medium', description: 'Normal course of business', priority: 2},
  {name: 'high', description: 'Should be addressed soon', priority: 3},
  {name: 'critical', description: 'Drop everything', priority: 4},
]);
const typeIds = await seedLookup('issue_type', [
  {name: 'bug', description: 'Something is broken'},
  {name: 'feature', description: 'New functionality or an enhancement'},
  {name: 'task', description: 'Chores, docs, and everything else'},
]);

// 2. Local users for the GitHub reporters/assignees, so the issues' user columns
// point at real platform users (rendered with the standard user pickers/avatars).
// A GitHub login is matched to an existing user by its normalized form ('alex-aprm'
// → the local 'alex.aprm'); the rest are created.
const normalize = (login) => login.toLowerCase().replace(/[^a-z0-9]/g, '');
const existingUsers = new Map(
  (await grok.dapi.users.list({pageSize: 1000})).map((u) => [normalize(u.login), u.id]));
const userIds = new Map();
let createdUsers = 0;
for (const gh of data.users) {
  let id = existingUsers.get(normalize(gh.login));
  if (id == null) {
    const user = DG.User.create();
    user.login = gh.login.toLowerCase();
    user.status = 'active';
    const parts = (gh.name ?? gh.login).trim().split(/\s+/);
    user.firstName = parts[0];
    user.lastName = parts.slice(1).join(' ');
    await grok.dapi.users.save(user);
    id = (await grok.dapi.users.filter(`login = "${user.login}"`).first()).id;
    createdUsers++;
  }
  userIds.set(gh.login, id);
}
grok.shell.info(`Users: ${userIds.size} mapped (${createdUsers} created)`);

// 3. The project the tickets land in (issue numbers are per-project, GitHub's fit as-is).
const projectId = idOf((await table('project').insert({
  key: 'PUB', name: data.repo,
  description: `Tickets imported from https://github.com/${data.repo}`,
}))[0]);

// 4. Labels, with their GitHub colors.
const labelReports = await table('label').insert(data.labels.map((l) => ({name: l.name, color: l.color})));
const labelIds = new Map(data.labels.map((l, i) => [l.name, idOf(labelReports[i])]));

// 5. The tickets — ONE bulk upsert, merged by the issue's (project_id, number)
// business key. GitHub has no priority/type fields, so both are derived from labels.
const priorityOf = (labels) =>
  labels.includes('critical') || labels.includes('blocker') ? 'critical' :
  labels.includes('high') ? 'high' : labels.includes('low') ? 'low' : 'medium';
const typeOf = (labels) =>
  labels.includes('bug') ? 'bug' :
  labels.includes('enhancement') || labels.includes('feature') ? 'feature' : 'task';

const report = await table('issue').batch(data.issues.map((i) => ({
  project_id: projectId,
  number: i.number,
  title: i.title,
  description: i.body,
  status_id: statusIds.get(i.state === 'closed' ? 'closed' : 'open'),
  priority_id: priorityIds.get(priorityOf(i.labels)),
  type_id: typeIds.get(typeOf(i.labels)),
  reporter: userIds.get(i.reporter) ?? null,
  assignee: userIds.get(i.assignee) ?? null,
})), {mode: 'upsert'});
grok.shell.info(`Issues: ${report.inserted} inserted, ${report.updated} updated, ${report.errorCount} errors`);

// 6. The N:N issue ↔ label links: query the issue numbers back (system columns like
// `id` always ride along), then bulk-upsert the junction rows — the
// (issue_id, label_id) key dedupes.
const issues = await table('issue').query(
  {filter: `project_id = "${projectId}"`, columns: ['number'], limit: 10000});
const issueIds = new Map(issues.map((r) => [r.number, r.id]));
const links = data.issues.flatMap((i) => i.labels
  .filter((l) => labelIds.has(l))
  .map((l) => ({issue_id: issueIds.get(i.number), label_id: labelIds.get(l)})));
const linkReport = await table('issue_label').batch(links, {mode: 'upsert'});
grok.shell.info(`Labels: ${linkReport.inserted} links created, ${linkReport.skipped} already there`);

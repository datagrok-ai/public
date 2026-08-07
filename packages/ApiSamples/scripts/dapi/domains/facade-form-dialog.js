// UI-FACADE §1.2 — an insert / edit form in a dialog, in TWO lines.
//
// `domains.table(...)` is the ONE await: it prefetches the typed client, the registry
// metadata and your capabilities, so every widget factory on the handle is synchronous.
// Free: inputs generated from the registry (types, labels, choices, nullability),
// reference pickers with async typeahead for FK / user columns, client validation with
// the server's own rules, server errors mapped onto the named fields, the standard 409
// dialog, writable-column gating. Cancel / Esc discards silently.
//
// In a package: `import {domains} from '@datagrok-libraries/domain-ui';`
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

const project = await domains.pick('grit.project');   // DomainRow | null (§1.4)
if (project == null)
  return;

// The example, verbatim:
const issues = await domains.table('grit.issue');
const saved = await issues.formDialog({values: {project_id: project.id}});

grok.shell.info(saved ? 'Issue created' : 'Cancelled — nothing was written');

// Editing an existing row is the same call with the row instead of the values.
const [values] = await issues.query({filter: DG.cond('project_id', '=', project.id), limit: 1});
if (values != null) {
  const row = new DG.DomainObjectHandler('grit.issue').rowFrom(values);
  await issues.formDialog({row});
}

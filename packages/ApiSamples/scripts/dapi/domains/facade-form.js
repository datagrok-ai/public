// The domain form, plain: a property form for ONE row of a domain table. The inputs
// (types, labels, choices, reference pickers), the validation and the one-transaction
// `save()` all come from the registry — nothing here is hand-wired.
//
// In a package: `import {domains} from '@datagrok-libraries/domain-ui';`
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

const issues = await domains.table('grit.issue');
const form = issues.form({values: {priority: 'critical'}});   // insert mode; {row} edits
grok.shell.newView('New issue', [form.root, ui.bigButton('Save', () => form.save())]);

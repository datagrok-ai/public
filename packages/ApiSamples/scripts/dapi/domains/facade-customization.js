// UI-FACADE §1.6c / §1.6d — override ONE thing, stay in the few-lines world:
// replace the generated input for one column, and add one validator. Everything else
// about the form (registry inputs, reference pickers, the editor behind it, the
// machine surface) stays exactly as it was.
//
// In a package: `import {domains} from '@datagrok-libraries/domain-ui';`
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

const issues = await domains.table('grit.issue');

// §1.6c — one form field:
const form = issues.form();
form.replaceInput('description', (p) => ui.input.textArea(p.caption));

// §1.6d — one validation, sync or async, rendered as the standard per-field marker:
form.addValidator('title', (s) => s.trim().length < 5 ? 'At least 5 characters' : null);

grok.shell.newView('New issue', [form.root]);

// Both are visible through the machine surface too — the customized form keeps it.
await form.ready;
form.props['title'] = 'x';
grok.shell.info(form.getWidgetStatus().inputs.find((i) => i.name === 'title').error);

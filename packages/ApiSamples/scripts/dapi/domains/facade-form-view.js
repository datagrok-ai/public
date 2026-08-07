// UI-FACADE §1.3 — a routed form page, in TWO lines.
//
// `domains.formView(w)` is `domains.view([w])`: a View — which IS a widget — hosting
// the form, with the app's URL identity, the form's functions (Save / Discard / Reset)
// on the ribbon, its pending-change count on the shell status bar, and the
// unsaved-changes gate + beforeunload while it is dirty.
//
// In a package: `import {domains} from '@datagrok-libraries/domain-ui';`
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

const me = await grok.dapi.users.current();

// The example, verbatim:
const issues = await domains.table('grit.issue');
grok.shell.preview = domains.formView(issues.form({values: {reporter: me.id}}));

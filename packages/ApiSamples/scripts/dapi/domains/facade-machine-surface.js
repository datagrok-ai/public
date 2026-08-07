// UI-FACADE §1.11 — what an AI agent (or a UI test) sees and does: the PLATFORM's own
// widget protocol, nothing domain-ui-specific.
//
// Discovery is `DG.Widget.getAll()`; description and live field state come from
// `getWidgetStatus()` (with the additive `inputs` array); named values are widget
// PROPERTIES, whose writes take the same single-writer path a keystroke takes; actions
// are REAL registered Funcs from `getFunctions()`.
//
// Free for the app developer: 0 lines.
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

const issues = await domains.table('grit.issue');
grok.shell.newView('New issue', [issues.form().root]);
await DG.delay(500);

// The example, verbatim:
const w = DG.Widget.getAll().find((w) => w.type === 'DomainForm');
w.getWidgetStatus();
// → {description: 'Creates a issue in grit.issue. Fields: ...',
//    error: 'Project: Value can\'t be empty', parts: {title: <el>, priority: <el>, ...},
//    inputs: [{name: 'title', type: 'string', value: null, required: true, valid: false},
//             {name: 'priority', type: 'string', choices: ['low', 'medium', 'high', 'critical']},
//             {name: 'assignee', type: 'ref', ref: 'User', value: null}, ...],
//    hitAreas: {}, shortcuts: {}, events: []}
w.props['priority'] = 'critical';       // the platform's property bag — same path as typing
w.props['assignee'] = 'Alex';           // a ref property resolves display text through the
                                        // SAME suggestion source its typeahead uses;
                                        // ambiguity → a validation error in the status
w.getWidgetStatus().error;              // re-read validation
const save = w.getFunctions().find((f) => f.name === 'Save');
await save.apply({widget: w});          // a real Func — typed, annotated, AI-discoverable

grok.shell.info(JSON.stringify(w.getWidgetStatus().inputs
  .filter((i) => ['priority', 'assignee'].includes(i.name))
  .map((i) => `${i.name} = ${i.value} (${i.valid ? 'valid' : i.error})`), null, 2));

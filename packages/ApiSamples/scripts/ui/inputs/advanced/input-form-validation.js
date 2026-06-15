// Manual form validation + reacting to validation events.

const sin = DG.Func.find({name: 'Sin'})[0];
const form = await DG.InputForm.forFuncCall(sin.prepare({x: 1}));
form.onValidationCompleted.subscribe(() => form.isValid ? grok.shell.info('Valid') : grok.shell.error('Invalid'));

const validateBtn = ui.button('Validate', () => form.validateInputs());
grok.shell.newView('Validation', [ui.divV([form.root, validateBtn])]);

// Pre-call and post-call actions
grok.functions.onBeforeRunAction.subscribe((_) => grok.shell.info('Runs before a call'));
grok.functions.onAfterRunAction.subscribe((_) => grok.shell.info('Runs after a call'));

// Filtering events to specific function calls
grok.functions
  .onAfterRunAction
  .pipe(rxjs.operators.filter((call) => call.func.name === 'Max'))
  .subscribe((call) => grok.shell.info(`Max of ${Array.from(call.inputs.values())}: ${call.getOutputParamValue()}`));

let x = await DG.Func.byName('Abs').apply({x: -5});
let y = await DG.Func.byName('Max').apply({nums: [-4, 9, 7]});

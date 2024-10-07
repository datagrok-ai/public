const intInput = ui.input.int('int', {value: 6, min: 2, max: 12, step: 1});
intInput.addValidator((v) => +v % 2 === 0 ? null : 'Value should be even');
const floatInput = ui.input.float('float', {value: 6.5758, min: 0, max: 10, step: 0.5});
floatInput.addValidator((v) => +v > 4.75 ? 'Value should be less than 4.75' : null);
const textInput = ui.input.string('string', {value: 'text'});
textInput.addValidator((v) => v.length > 5 ? null : 'Text should be longer than 5 characters');

ui.dialog({title: 'Inputs'})
  .add(intInput).add(floatInput).add(textInput)
  .onOK(() => grok.shell.info(`Int: ${intInput.value}, Float: ${floatInput.value}, Text: ${textInput.value}`))
  .show();
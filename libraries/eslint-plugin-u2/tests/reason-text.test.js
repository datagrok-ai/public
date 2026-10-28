'use strict';

const tester = require('./tester');
const rule = require('../rules/reason-text');

const empty = (code, ...texts) => ({code, errors: texts.map((text) => ({messageId: 'empty', data: {text}}))});

tester.run('reason-text', rule, {
  valid: [
    `button('Submit', run, {disabled: 'Name the campaign first'});`,
    `new TextInput({label: 'Name', disabled: computed(() => locked.value ? 'Only the owner can rename it' : null)});`,
    `input.disabled = file ? null : 'Pick a source file first';`,
    `button('Save', save, {disabled: true});`,
    `dialog.close({reason: 'cancel'});`,
    `const status = 'Disabled';`,
    `el.setAttribute('title', 'Disabled');`,
  ],
  invalid: [
    empty(`button('Merge', merge, {disabled: 'Disabled'});`, 'Disabled'),
    empty(`new ChoiceInput({label: 'Key', items, disabled: table ? null : \`Not available\`});`, 'Not available'),
    empty(`picker.disabled = !table && 'N/A';`, 'N/A'),
    empty(`const run = button('Run', go, {disabled: computed(() => {
      if (!ready.value)
        return 'Unavailable';
      return busy.value ? 'disabled.' : null;
    })});`, 'Unavailable', 'disabled.'),
    empty(`dialog.okDisabled('Not applicable');`, 'Not applicable'),
    empty(`pane.setAttribute('aria-description', 'Not available');`, 'Not available'),
  ],
});

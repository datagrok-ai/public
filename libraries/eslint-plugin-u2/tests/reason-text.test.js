'use strict';

const tester = require('./tester');
const rule = require('../rules/reason-text');

const empty = (code, ...texts) => ({code, errors: texts.map((text) => ({messageId: 'empty', data: {text}}))});
const U2 = `import {button, iconButton, TextInput} from '@datagrok-libraries/u2';
`;
const off = (code, what) => ({code: U2 + code, errors: [{messageId: 'boolean', data: {what}}]});

tester.run('reason-text', rule, {
  valid: [
    `button('Submit', run, {disabled: 'Name the campaign first'});`,
    `new TextInput({label: 'Name', disabled: computed(() => locked.value ? 'Only the owner can rename it' : null)});`,
    `input.disabled = file ? null : 'Pick a source file first';`,
    `button('Save', save, {disabled: true});`,
    `dialog.close({reason: 'cancel'});`,
    `const status = 'Disabled';`,
    `el.setAttribute('title', 'Disabled');`,
    `input.disabled = '';`,
    U2 + `const run = button('Run', go); run.disabled = ready ? null : 'Pick a source file first';`,
    U2 + `const native = document.createElement('button'); native.disabled = true;`,
    U2 + `const b = ui.button('Run', go); b.disabled = true; b.setAttribute('disabled', '');`,
    U2 + `function native() { const run = document.createElement('button'); run.disabled = true; }
    function toolbar() { const run = button('Run', go); return run; }`,
    U2 + `const run = button('Run', go);
    function busy(run) { run.disabled = true; }`,
    'button(\'Save\', save, {disabled: `${verb}…`});',
    `button('Save', save, {disabled: verb + '…'});`,
    'button(\'Save\', save, {disabled: ready ? null : `${missing} …`});',
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
    empty(`button('Save', save, {disabled: '...'});`, '...'),
    empty(`button('Save', save, {disabled: 'Not allowed'});`, 'Not allowed'),
    empty(`button('Save', save, {disabled: 'Cannot save'});`, 'Cannot save'),
    empty(`button('Export', run, {disabled: ' '});`, ' '),
    off(`const run = button('Run', go); run.disabled = true;`, 'run.disabled = true'),
    off(`const name = new TextInput({label: 'Name'}); name.disabled = !ready;`, 'name.disabled = !ready'),
    off(`const share = iconButton('share', go, {tooltip: 'Share'}); share.setAttribute('disabled', '');`,
      'share.setAttribute(\'disabled\', \'\')'),
    off(`function toolbar() { const run = button('Run', go); run.disabled = true; }
    function native() { const run = document.createElement('button'); run.disabled = true; }`, 'run.disabled = true'),
  ],
});

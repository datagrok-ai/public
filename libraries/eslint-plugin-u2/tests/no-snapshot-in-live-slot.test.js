'use strict';

const tester = require('./tester');
const rule = require('../rules/no-snapshot-in-live-slot');

const U2 = `import {span, div, button, computed, signal, Dialog, TextInput, textInput}
  from '@datagrok-libraries/u2';\n`;
const input = {ctor: {options: {0: ['label', 'tooltip']}}, methods: {setLabel: {at: [0]}}};
const slots = (valueAccessor) => ({
  functions: {span: {at: [0]}, div: {at: [0]}, button: {options: {2: ['disabled', 'tooltip']}},
    textInput: {options: {1: ['label']}, returns: 'TextInput'}},
  classes: {TextInput: {...input, valueAccessor}, Dialog: {methods: {onOK: {}, add: {at: [0]}}}},
});
const today = [slots(false)];
const accessor = [slots(true)];
const snapshot = (code, read, slot, fix, options = today) =>
  ({code: U2 + code, options, errors: [{messageId: 'snapshot', data: {read, slot, fix}}]});

tester.run('no-snapshot-in-live-slot', rule, {
  valid: [
    {code: U2 + `const name = new TextInput({label: 'Name'}); span(name.value);`, options: today},
    {code: U2 + `const name = new TextInput({label: 'Name'}); span(name.signal);`, options: accessor},
    {code: U2 + `const count = signal(0); span(computed(() => \`\${count.value} rows\`));`, options: today},
    {code: U2 + `const total = signal(0); button('Save', () => save(total.value), {});`, options: today},
    {code: U2 + `const caption = computed(() => 'Save'); const d = new Dialog(); d.onOK(run, {text: caption.value});`,
      options: today},
    {code: U2 + `const compound = signal('Aspirin'); new TextInput({label: 'Compound', value: compound.value});`,
      options: today},
    {code: U2 + `const n = signal(1); console.log(n.value); grok.shell.info(n.peek());`, options: today},
    {code: U2 + `const name = new TextInput({label: 'Name'});
      div([new TextInput({label: 'Copy', value: name.value}), span(name.signal)]);`, options: accessor},
  ],
  invalid: [
    snapshot(`const name = new TextInput({label: 'Name'}); span(name.value);`, 'name.value', '`span(…)`',
      'name.signal', accessor),
    snapshot(`const name = new TextInput({label: 'Name'}); span(name.value.value);`, 'name.value.value', '`span(…)`',
      'name.value'),
    snapshot(`const count = signal(0); span(\`\${count.value} rows\`);`, 'count.value', '`span(…)`', 'count'),
    snapshot(`const busy = computed(() => jobs.value > 0); button('Run', run, {disabled: busy.peek()});`,
      'busy.peek()', '`disabled` of `button(…)`', 'busy'),
    snapshot(`const ready = signal(false); div([span('Status'), ready.value ? 'ready' : 'waiting']);`,
      'ready.value', '`div(…)`', 'ready'),
    snapshot(`class Pane { private readonly count = signal(0); render() { return span(this.count.value); } }`,
      'this.count.value', '`span(…)`', 'this.count'),
    snapshot(`const name = textInput('Name', {}); const other = new TextInput({tooltip: name.value.peek()});`,
      'name.value.peek()', '`tooltip` of `new TextInput(…)`', 'name.value'),
    snapshot(`function show(name: TextInput) { name.setLabel(name.value); }`, 'name.value', '`name.setLabel(…)`',
      'name.signal', accessor),
    snapshot(`const table = new TextInput({label: 'Table'}); div([h3('Tables'), readout('table', table.value)]);`,
      'table.value', '`div(…)`', 'table.signal', accessor),
  ],
});

'use strict';

const tester = require('./tester');
const rule = require('../rules/no-raw-column-list');

const U2 = `import {ChoiceInput, MultiChoiceInput, Form, computed} from '@datagrok-libraries/u2';
import {columnInput, Platform} from '@datagrok-libraries/u2/src/dg/index.js';\n`;
const raw = (code, what, key = 'items') => ({code: U2 + code, errors: [{messageId: 'raw', data: {key, what}}]});

tester.run('no-raw-column-list', rule, {
  valid: [
    U2 + `const molecules = columnInput('Molecules', table, {filter: (c) => c.semType === 'Molecule'});`,
    U2 + `new ChoiceInput({label: 'Column', items: Platform.listedColumns(table).map((c) => c.name)});`,
    U2 + `new ChoiceInput({label: 'Mode', items: Object.values(MODE_LABELS)});`,
    U2 + `const names = df.columns.names(); expectSame(names, ['Molecule']);`,
    `ui.input.choice('Column', {items: df.columns.names()});`,
    U2 + `new ChoiceInput({label: 'Column', items: df.columns.length > 0 ? ['a'] : []});`,
  ],
  invalid: [
    raw(`new ChoiceInput({label: 'Molecule column', items: df.columns.names()});`, '.columns.names()'),
    raw(`let molecules;
      const columns = computed(() => table.value?.columns.toList().filter(isMolecule).map((c) => c.name) ?? []);
      molecules = new ChoiceInput({label: 'Molecules', items: columns});`, '.columns.toList()'),
    raw(`const columns = computed(() => [NONE_OPTION, ...incoming.value?.columns.names() ?? []]);
      const molCol = new ChoiceInput({label: 'Molecule column (incoming)', items: columns});`, '.columns.names()'),
    raw(`new MultiChoiceInput({label: 'Columns', items: Array.from(df.columns).map((c) => c.name)});`,
      'Array.from(….columns)'),
    raw(`new Form({fields: [{name: 'key', choices: [...df.columns].map((c) => c.name)}]});`, '[...columns]',
      'choices'),
  ],
});

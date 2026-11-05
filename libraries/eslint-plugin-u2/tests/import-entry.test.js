'use strict';

const tester = require('./tester');
const rule = require('../rules/import-entry');

const deep = (code, entry) => ({code, errors: [{messageId: 'deep', data: {source: /'([^']+)'/.exec(code)[1], entry}}]});

tester.run('import-entry', rule, {
  valid: [
    `import {Form, TextInput, ChoiceInput, Dialog} from '@datagrok-libraries/u2';`,
    `import {ViewHost, ObjectForm} from '@datagrok-libraries/u2/dg';`,
    `import {u2test} from '@datagrok-libraries/u2/testing';`,
    `import {parseFilter} from '@datagrok-libraries/u2/dg/filter';`,
    `import '@datagrok-libraries/u2/css/icons.css';`,
    `import {DomainApp} from '@datagrok-libraries/ems/src/dg';`,
    `import {div} from './elements.js';`,
    {code: `import {host} from '@datagrok-libraries/u2/src/dg';`, options: [{entries: ['', '/src/dg']}]},
  ],
  invalid: [
    deep(`import {host, appView, ObjectForm} from '@datagrok-libraries/u2/src/dg/index.js';`,
      '@datagrok-libraries/u2/dg'),
    deep(`import {ObjectForm} from '@datagrok-libraries/u2/src/dg/forms/object-form.js';`, '@datagrok-libraries/u2/dg'),
    deep(`import {host} from '@datagrok-libraries/u2/src/dg';`, '@datagrok-libraries/u2/dg'),
    deep(`export * from '@datagrok-libraries/u2/src/core/filter';`, '@datagrok-libraries/u2/filter'),
    deep(`export {FilterInput} from '@datagrok-libraries/u2/src/dg/filter';`, '@datagrok-libraries/u2/dg/filter'),
    deep(`import {installDom} from '@datagrok-libraries/u2/tests/dom-shim.js';`, '@datagrok-libraries/u2/testing'),
    deep(`const u2 = await import('@datagrok-libraries/u2/src/index.js');`, '@datagrok-libraries/u2'),
    deep(`import type {Owner} from '@datagrok-libraries/u2/src/core/owner.js';`, '@datagrok-libraries/u2'),
  ],
});

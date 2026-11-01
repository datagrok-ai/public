'use strict';

const tester = require('./tester');
const rule = require('../rules/comments');

tester.run('comments', rule, {
  valid: [
    `/** The campaign's name, as the catalog lists it. */\nexport function name(c) { return c.name; }`,
    `/* eslint-disable max-len */\nexport const x = 1;`,
    `/*! (c) Datagrok, MIT */\nexport const x = 1;`,
    `import {div} from '@datagrok-libraries/u2';\n/* a block after the code starts is no header */\nexport const x = div();`,
    `export class Tile {
      /** What the tile shows; read by the catalog. */
      label = '';
      // the platform fires twice on a reorder
      private _pending = false;
      /** Opens it. */
      open(): void {}
    }`,
    `// One bundle load shared by every test file.\nexport const x = 1;`,
  ],
  invalid: [
    {code: `/* Wave 2 adds the compute pane:\n   the tracker moved here from hit-design-app. */\nexport const x = 1;`,
      errors: [{messageId: 'header', line: 1}]},
    {code: `export class Session {
      /** Set once the view is closed. */
      private _closed = false;
      /** Re-reads the campaign. */
      protected reload(): void {}
      /** The live rows. */
      #rows = [];
      constructor(/** the owning view */ private readonly view: unknown) {}
    }`, errors: [
      {messageId: 'privateDoc', data: {access: 'private'}, line: 2},
      {messageId: 'privateDoc', data: {access: 'protected'}, line: 4},
      {messageId: 'privateDoc', data: {access: 'private'}, line: 6},
      {messageId: 'privateDoc', data: {access: 'private'}, line: 8},
    ]},
  ],
});

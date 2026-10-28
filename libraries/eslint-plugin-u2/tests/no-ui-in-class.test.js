'use strict';

const tester = require('./tester');
const rule = require('../rules/no-ui-in-class');

const options = [{controls: ['Control', 'Dialog', 'CatalogTable', 'TextInput'], factories: ['div', 'span']}];
const U2 = `import {AsyncValue, CatalogTable, Control, Dialog, TextInput, div, signal, span}
  from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';\n`;
const ui = (code, name, what, opts = options) => ({code: U2 + code, options: opts,
  errors: [{messageId: 'ui', data: {name, what}}]});

tester.run('no-ui-in-class', rule, {
  valid: [
    {code: U2 + `class Holder { static of(o) { return new CatalogTable(o); } }`, options},
    {code: U2 + `class Tracker { open() { return new Dialog('Progress Tracker'); } }`, options},
    {code: U2 + `class Tracker { readonly show = () => new Dialog('Progress Tracker').show(); }`, options},
    {code: U2 + `class Session {
      constructor(app, view) { this.listing = new AsyncValue(() => app.list()); this.moved = signal(0); }
    }`, options},
    {code: U2 + `class Session {
      constructor(o) { this.home = ViewHost.of(o.view).build(() => new CatalogTable(o.catalog)); }
    }`, options},
    {code: U2 + `class Panel {
      readonly picker = Control.build(() => [div([new TextInput({label: 'Name'})])]);
    }`, options},
    {code: U2 + `class Picker extends Control { readonly box = div(); }`, options},
    {code: `class Plain { readonly box = document.createElement('div'); }`, options},
  ],
  invalid: [
    ui(`export class ProgressTracker {
      private readonly dialog = new Dialog('Progress Tracker', {name: 'progress-tracker'});
      private readonly lanes = div([], 'hit-design-tracker-lanes');
    }`, 'ProgressTracker', 'new Dialog(…)'),
    ui(`export class CampaignCatalog {
      readonly table;
      constructor(o) {
        this.table = new CatalogTable({columns: o.columns, groupings: o.groupings});
      }
    }`, 'CampaignCatalog', 'new CatalogTable(…)'),
    ui(`class AbGrid {
      readonly root = div([span('Antibodies')], 'u2demo-ab-grid');
    }`, 'AbGrid', 'div(…)'),
    ui(`class Editor {
      constructor(scope) { this.form = Scope.runWith(scope, () => new TextInput({label: 'Name'})); }
    }`, 'Editor', 'new TextInput(…)'),
    ui(`import * as DG from 'datagrok-api/dg';
    class MonomerManagerView extends DG.ViewBase {
      constructor() { super(); this.editor = new TextInput({label: 'Symbol'}); }
    }`, 'MonomerManagerView', 'new TextInput(…)'),
    ui(`class Defaults { readonly dialog = new Dialog('Defaults'); }`, 'Defaults', 'new Dialog(…)', []),
  ],
});

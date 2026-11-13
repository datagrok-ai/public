'use strict';

const tester = require('./tester');
const rule = require('../rules/no-subclass');

const subclass = (code, name, lib, base) => ({code, errors: [{messageId: 'subclass', data: {name, lib, base}}]});

tester.run('no-subclass', rule, {
  valid: [
    `import * as DG from 'datagrok-api/dg';\nexport class MolstarViewer extends DG.JsViewer {}`,
    `import * as DG from 'datagrok-api/dg';\nclass Panel extends DG.Widget {}\nclass View extends DG.ViewBase {}`,
    `import {Control} from '@datagrok-libraries/u2';\nexport class ComputeEntries { static of() { return []; } }`,
    `import {Control} from '@datagrok-libraries/u2';\nclass Base {}\nclass Local extends Base {}`,
    `import {Control} from './control.js';\nexport class Card extends Control {}`,
  ],
  invalid: [
    subclass(`import {Control, signal} from '@datagrok-libraries/u2';
      export class TemplatePicker<T extends CampaignTemplate> extends Control {
        readonly names = signal<string[] | null>(null);
      }`, 'TemplatePicker', 'u2', 'Control'),
    subclass(`import {Control as U2Control} from '@datagrok-libraries/u2';
      export class MonomerEditor extends U2Control implements INewMonomerForm {}`, 'MonomerEditor', 'u2', 'Control'),
    subclass(`import {DomainApp} from '@datagrok-libraries/ems/src/dg';
      class IssuesApp extends DomainApp {}`, 'IssuesApp', 'ems', 'DomainApp'),
    subclass(`import * as DG from 'datagrok-api/dg';
      class Sample extends DG.U2.Control {}`, 'Sample', 'u2', 'Control'),
    subclass(`import * as u2 from '@datagrok-libraries/u2';
      const Home = class HitHome extends u2.Control {};`, 'HitHome', 'u2', 'Control'),
    subclass(`import {Component} from 'datagrok-api/u2core';
      export default class extends Component {}`, 'This class', 'u2', 'Component'),
  ],
});

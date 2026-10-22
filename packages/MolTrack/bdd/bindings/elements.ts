/* The MolTrack EMS app as a context: only its own pieces; the rest resolves from the u2/EMS contract. */
import {context} from '@datagrok-libraries/bdd';

export const ems = context('MolTrack EMS', {selector: '[data-u2="domain-app"]'});

ems.element('app ribbon', {selector: '.grok-view-container > .d4-ribbon',
  description: 'the view\'s ribbon (New, Save, Discard, ⋯), told apart from the New of the form\'s Batches tab'});
ems.element('compound list', {selector: '[data-u2="domain-list"] [data-u2="list"]',
  description: 'the virtual list inside the domain list: the control that reports every compound it holds'});
ems.element('compound count', {selector: '.u2-domain-status-count',
  description: 'the count slot of the status bar: "N", or "50 of N" while the list has fetched one page'});
ems.element('save problem', {selector: '.u2-domain-status-problem',
  description: 'the problem slot of the status bar: why the last save was refused'});
ems.element('molecule canvas', {selector: 'canvas.d4-input-molecule-canvas',
  description: 'the picture of a platform molecule input; a click opens the sketcher dialog'});
ems.element('filters pane', {selector: '.u2-domain-app-filters-pane',
  description: 'the condition builder on the list\'s left'});
ems.element('add condition', {selector: '[data-u2="filter-builder"] [data-u2-part="add"]',
  description: 'the + of the builder\'s conditions'});
ems.element('condition property', {selector: '[data-u2="filter-row"] [data-u2-part="prop"]',
  description: 'the column a condition is about'});
ems.element('condition operator', {selector: '[data-u2="filter-row"] [data-u2-part="op"]',
  description: 'what a condition asks of its column'});
ems.element('condition value', {selector: '[data-u2="filter-row"] [data-u2-part="value"]',
  description: 'the value a condition compares with: a molecule input on a structure column'});
ems.element('condition threshold', {selector: '[data-u2="filter-row"] [data-u2-part="parameter"]',
  description: 'the similarity a structure condition asks for at least'});
ems.element('condition remove', {selector: '[data-u2="filter-row"] [data-u2-part="remove"]',
  description: 'the − that removes a condition'});
ems.element('register mapping', {selector: '.moltrack-ems-register-mapping',
  description: 'the bulk registration\'s table of file columns and what each registers as'});
ems.element('register preview', {selector: '.moltrack-ems-register-preview',
  description: 'what each row of the file becomes, before anything is written'});
ems.element('register report', {selector: '.moltrack-ems-register-report',
  description: 'the outcome of the registration, row by row'});
ems.element('more actions', {selector: '[data-u2="actions-menu"]',
  description: 'the ⋯ of the ribbon: Import, Bulk edit and Trash on a list, the row\'s Delete and actions on a record'});

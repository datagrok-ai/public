/* The Access protocol (WO-5): the two constants answer every name, a built access answers what
   it lists and hides the rest, a row narrows edit/delete by its own `~can_*` columns — and the
   action surfaces drop what the access denies while keeping what is merely disabled. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Access} from '../src/core/access.js';
import {rowActions, actionsMenu} from '../src/components/actions/actions.js';

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const DATA = {can: {view: true, insert: true, edit: true, delete: false, approve: true},
  fields: {title: 'editable', number: 'readonly'}};

scoped('Access.full and Access.readOnly answer every capability and field, and are frozen', () => {
  assert.equal(Access.full.can('delete'), true);
  assert.equal(Access.full.can('approve'), true, 'a custom capability too');
  assert.equal(Access.full.field('anything'), 'editable');
  assert.equal(Access.readOnly.can('view'), true);
  assert.equal(Access.readOnly.can('edit'), false);
  assert.equal(Access.readOnly.can('approve'), false);
  assert.equal(Access.readOnly.field('anything'), 'readonly');
  assert.equal(Object.isFrozen(Access.full), true);
  assert.equal(Object.isFrozen(Access.readOnly), true);
});

scoped('Access.from: listed capabilities and fields as given, an unlisted field is hidden', () => {
  const access = Access.from(DATA);
  assert.equal(access.can('edit'), true);
  assert.equal(access.can('delete'), false);
  assert.equal(access.can('approve'), true);
  assert.equal(access.can('share'), false, 'unlisted is denied');
  assert.equal(access.field('title'), 'editable');
  assert.equal(access.field('number'), 'readonly');
  assert.equal(access.field('salary'), 'hidden', 'a restricted column is absent from fields');
  DATA.can.edit = false;
  assert.equal(access.can('edit'), true, 'a copy of the data, not a view over it');
  DATA.can.edit = true;
});

scoped('Access.row: a carried ~can_* column replaces the table\'s answer either way', () => {
  assert.deepEqual(Access.ROW_COLUMNS.map(([column]) => column), ['~can_edit', '~can_delete', '~can_share'],
    'the js-api DOMAIN_ACCESS_COLUMNS');
  const access = Access.from({...DATA, can: {...DATA.can, share: true}});
  assert.equal(access.row({title: 'x'}), access, 'a row without the columns answers the table');
  const denied = access.row({'~can_edit': false, '~can_delete': true});
  assert.equal(denied.can('edit'), false, 'the row narrows');
  assert.equal(denied.can('delete'), true,
    'the row widens: table-level delete is a false negative in row mode, the row is the server\'s truth');
  assert.equal(denied.can('share'), true, 'not carried: the table\'s answer');
  assert.equal(denied.field('title'), 'readonly', 'a row the caller may not edit has no editable field');
  assert.equal(denied.field('number'), 'readonly');
  assert.equal(denied.field('salary'), 'hidden', 'hidden stays hidden');
  assert.equal(access.row({'~can_edit': true}).field('title'), 'editable');
  const allowed = access.row({'~can_edit': true});
  assert.equal(allowed.can('edit'), true);
  assert.equal(allowed.can('delete'), false, 'not carried: the table said no');
  const granted = Access.from({can: {view: true, insert: false, edit: false, delete: false, share: false},
    fields: {title: 'editable'}}).row({'~can_edit': true, '~can_delete': false});
  assert.equal(granted.can('edit'), true, 'a row granted to the caller under a table with no grants');
  assert.equal(granted.field('title'), 'editable', 'and its editable field is editable');
  assert.equal(granted.can('delete'), false);
  assert.equal(granted.can('insert'), false, 'insert is table-level: never carried by a row');
  const unshared = access.row({'~can_share': false});
  assert.notEqual(unshared, access, 'a row carrying only ~can_share is narrowed too');
  assert.equal(unshared.can('share'), false);
  assert.equal(unshared.can('edit'), true);

  assert.equal(access.row({'~can_approve': false}).can('approve'), false,
    'generic: any ~can_<name> the table knows a capability for');
  assert.equal(access.row({'~can_share': null}), access, 'null (the server off row mode) is not carried');
  assert.equal(access.row({'~can_edit': 'yes'}), access, 'nor anything but a boolean');
  assert.equal(access.row(null), access);
  const fullRow = Access.full.row({'~can_delete': false});
  assert.equal(fullRow.can('delete'), false);
  assert.equal(fullRow.can('edit'), true);
  assert.equal(fullRow.can('approve'), true, 'full stays full for everything the row did not narrow');
  assert.equal(fullRow.field('anything'), 'editable');
});

scoped('field folds the write capability in: editable needs edit on a row, insert on a draft', () => {
  // the server lists every unrestricted column as editable even for a caller with no grants
  const none = Access.from({can: {view: true, insert: false, edit: false, delete: false, share: false},
    fields: {title: 'editable', number: 'readonly'}});
  assert.equal(none.field('title'), 'readonly', 'listed editable, but the caller may not edit');
  assert.equal(none.field('number'), 'readonly');
  assert.equal(none.field('salary'), 'hidden');
  assert.equal(none.forDraft().field('title'), 'readonly', 'nor insert');

  const inserter = Access.from({can: {view: true, insert: true, edit: false, delete: false, share: false},
    fields: {title: 'editable', number: 'readonly'}});
  assert.equal(inserter.field('title'), 'readonly', 'an existing row: edit is what counts');
  assert.equal(inserter.isDraft, false);
  const draft = inserter.forDraft();
  assert.equal(draft.isDraft, true);
  assert.equal(draft.field('title'), 'editable', 'a draft: insert is what counts');
  assert.equal(draft.field('number'), 'readonly');
  assert.equal(draft.forDraft(), draft, 'already the draft view');
  assert.equal(draft.row({'~can_edit': false}).field('title'), 'editable',
    'a row narrowing edit keeps the draft view, which gates on insert');

  // a draft reports itself through the editor's state column: its ~can_* cells are the frame's
  // defaults, not truth, so row() answers the draft view whatever they hold
  const draftRow = inserter.row({'~state': 'new', '~can_edit': false, '~can_delete': false});
  assert.equal(draftRow.isDraft, true);
  assert.equal(draftRow.field('title'), 'editable', 'insert is what counts on a draft');
  assert.equal(draftRow.can('delete'), false, 'the table\'s answer, not the cell\'s');
  assert.equal(Access.from(DATA).row({'~state': 'new', '~can_delete': true}).can('delete'), false);

  const editor = Access.from(DATA);
  assert.equal(editor.field('title'), 'editable');
  assert.equal(editor.row({'~can_edit': false}).forDraft().field('title'), 'editable');
  assert.equal(Access.full.forDraft().field('anything'), 'editable');
  assert.equal(Access.readOnly.forDraft().field('anything'), 'readonly');
});

function actions(log) {
  return [
    {name: 'Open', icon: 'external-link-alt', run: () => log.push('open')},
    {name: 'Escalate', icon: 'arrow-up', requires: 'edit', enabled: false, run: () => log.push('escalate')},
    {name: 'Delete', icon: 'trash', requires: 'delete', run: () => log.push('delete')},
    {name: 'Approve', requires: 'approve', run: () => log.push('approve')},
  ];
}

scoped('rowActions: a denied `requires` is not rendered, `enabled: false` is disabled; default full', () => {
  const log = [];
  const all = rowActions(actions(log));
  assert.deepEqual([...all.querySelectorAll('button')].map((b) => b.getAttribute('aria-label')),
    ['Open', 'Escalate', 'Delete'], 'Access.full by default: nothing outside EMS changes');

  const block = rowActions(actions(log), {access: Access.from(DATA)});
  const buttons = [...block.querySelectorAll('button')];
  assert.deepEqual(buttons.map((b) => b.getAttribute('aria-label')), ['Open', 'Escalate'],
    'permission ⇒ hidden');
  assert.equal(buttons[1].disabled, true, 'state ⇒ disabled');
  fire(buttons[0], 'click');
  assert.deepEqual(log, ['open']);

  const perRow = rowActions(actions(log), {access: Access.from(DATA), row: {'~can_edit': false}});
  assert.deepEqual([...perRow.querySelectorAll('button')].map((b) => b.getAttribute('aria-label')), ['Open'],
    'the row refines the table');
});

scoped('actionsMenu: the same rule over the full list', () => {
  const log = [];
  const menu = actionsMenu(actions(log), {access: Access.from(DATA)});
  menu.show({x: 10, y: 10});
  const items = [...document.querySelectorAll('[role="menuitem"]')];
  assert.deepEqual(items.map((el) => el.querySelector('.u2-menu-label').textContent),
    ['Open', 'Escalate', 'Approve'], 'Delete is gone; the custom permission the table grants stays');
  fire(items[2], 'click');
  assert.deepEqual(log, ['approve']);
  assert.equal(menu.isOpen.value, false);
});

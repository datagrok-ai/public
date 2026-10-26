/* EMS's kind catalog: one entry per `data-u2` value EMS stamps, in u2's shape (u2
   `src/core/kinds.ts`). The manifest generator turns it into this library's manifest `kinds`;
   tests/kinds.test.js fails while a stamped value has no entry. Contract:
   core/docs/features/ui2/CONTRACT.md §4, §5.8. */

import type {U2KindDef, U2StateKey} from '@datagrok-libraries/u2';

const INPUT_PARTS: Record<string, string> = {
  label: '[data-u2-part="label"]',
  editor: '[data-u2-part="editor"]',
  options: '[data-u2-part="options"]',
  error: '[data-u2-part="error"]',
};

function composite(kind: string, noun: string, state: U2StateKey[] = ['busy'],
  extra: Partial<U2KindDef> = {}): U2KindDef {
  return {kind, noun, family: 'composite', nameFrom: ['aria'], parts: {}, state, gestures: [], ...extra};
}

function button(kind: string): U2KindDef {
  return {kind, noun: 'button', family: 'action', role: 'button', nameFrom: ['text', 'aria', 'title'], parts: {},
    state: ['label', 'busy', 'pressed'], gestures: ['press']};
}

/** The type-ahead pickers and the search box are inputs: the shared parts, a label. */
function input(kind: string, noun: string, focusRole: string, state: U2StateKey[],
  extra: Partial<U2KindDef>): U2KindDef {
  return {kind, noun, family: 'input', focusRole, nameFrom: ['label', 'aria', 'placeholder'], parts: INPUT_PARTS,
    state: ['label', ...state], gestures: [], ...extra};
}

const PICK: Partial<U2KindDef> = {items: {role: 'option', noun: 'item', gestures: ['pick']},
  gestures: ['type', 'open', 'pick']};

export const EMS_KINDS: readonly U2KindDef[] = [
  composite('domain-app', 'domain app'),
  composite('domain-children', 'domain children'),
  composite('domain-data-table', 'domain data table', ['busy', 'items', 'selection'],
    {items: {role: 'row', noun: 'table row', gestures: ['select']}, gestures: ['select', 'scroll']}),
  composite('domain-filters', 'domain filters'),
  composite('domain-form', 'domain form', ['busy', 'valid', 'error']),
  composite('domain-grid', 'domain grid', ['busy', 'items', 'selection'],
    {items: {role: 'row', noun: 'table row', gestures: ['select']}, gestures: ['select', 'scroll']}),
  composite('domain-history', 'domain history', ['busy', 'items'],
    {items: {role: 'option', noun: 'item', gestures: ['select']}, gestures: ['select']}),
  composite('domain-list', 'domain list', ['busy', 'items', 'selection'],
    {items: {role: 'option', noun: 'item', gestures: ['select']}, gestures: ['select', 'scroll']}),
  input('domain-pick', 'domain pick', 'combobox', ['value', 'text', 'valid', 'error', 'open', 'placeholder'], PICK),
  composite('domain-presets', 'domain presets', ['value', 'current', 'items'],
    {role: 'group', items: {role: 'button', noun: 'button', gestures: ['press']}, gestures: ['pick']}),
  input('domain-search', 'domain search', 'textbox', ['value', 'text', 'placeholder'], {gestures: ['type', 'clear']}),
  composite('domain-tree', 'domain tree', ['busy', 'items', 'selection'],
    {items: {role: 'treeitem', noun: 'tree node', gestures: ['select', 'expand', 'collapse']},
      gestures: ['select', 'expand', 'collapse']}),
  button('save-button'),
  button('discard-button'),
  button('new-button'),
  button('refresh-button'),
  {kind: 'actions-menu', noun: 'dropdown button', family: 'action', role: 'button', nameFrom: ['aria', 'title'],
    parts: {}, state: ['label', 'open'], gestures: ['press', 'open', 'pick']},
  {kind: 'sort-menu', noun: 'dropdown button', aliases: ['sort'], family: 'action', role: 'button',
    nameFrom: ['aria', 'title'], parts: {}, state: ['label', 'open'], gestures: ['press', 'open', 'pick']},
  composite('binding-dialog', 'binding dialog'),
  composite('binding-edit-dialog', 'binding edit dialog'),
  composite('manifest-editor', 'manifest editor'),
  composite('manifest-panel', 'manifest panel', []),
  composite('manifest-tree', 'manifest tree', []),
];

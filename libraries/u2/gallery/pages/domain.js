import {Scope, VirtualList, rowActions, backends, DomainSource, MemoryDomainBackend} from '../../src/index.js';
import {divH, divV, span, button} from '../../src/core/elements.js';
import {propertyForm} from '../../src/dg/forms/object-form.js';
import {Editors} from '../../src/dg/forms/editors.js';
import {domainPick} from '../../src/dg/domain/pick.js';
import {saveButton, discardButton} from '../../src/dg/domain/buttons.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

for (const name of ['elements', 'inputs', 'choice', 'number', 'date', 'icons', 'buttons', 'tooltip', 'menu',
  'list', 'form', 'typeahead', 'async', 'notify', 'dialog', 'domain'])
  injectOnce(`u2-${name}-css`, `../../css/${name}.css`);

/** A `schema.json` as a package ships it — the memory backend answers the same seam the platform
 * fills over `grok.dapi.domains`, so everything below is what a Grit page does, without a server. */
const SCHEMA = {
  name: 'demo',
  tables: {
    person: {friendlyName: 'People', columns: {name: {type: 'string', required: true, isName: true}}},
    task: {
      friendlyName: 'Tasks',
      columns: {
        title: {type: 'string', required: true, isName: true},
        status: {type: 'string', choices: ['Open', 'Blocked', 'Done'], friendlyName: 'Status'},
        priority: {type: 'int', min: 1, max: 5, friendlyName: 'Priority'},
        assignee_id: {type: 'ref', ref: 'person', friendlyName: 'Assignee'},
        due: {type: 'datetime', friendlyName: 'Due'},
        notes: {type: 'string', editor: 'textarea'},
      },
    },
  },
};

const PEOPLE = [{id: 'u1', name: 'Ada'}, {id: 'u2', name: 'Grace'}, {id: 'u3', name: 'Edsger'}];
const day = 86400000;
const TASKS = Array.from({length: 40}, (_, i) => ({
  id: `t${i + 1}`, title: `Task #${i + 1}: something worth a row`, status: ['Open', 'Blocked', 'Done'][i % 3],
  priority: (i % 5) + 1, assignee_id: PEOPLE[i % 3].id, due: new Date(Date.now() + i * day).toISOString(),
  // row-level security as the server sends it: `withAccess` rows carry their own truth
  '~can_edit': i % 7 !== 3, '~can_delete': i % 2 === 0,
}));
const SYSTEM = ['id', 'version', 'created_on', 'updated_on'];

// the editor rule `src/dg/domain/form.ts` registers in the platform; here the page does, so a ref
// column is a picker over its table
Editors.register({
  match: (p) => /^\w+\.\w+$/.test(p.semType ?? ''),
  create: (p, options) => domainPick(p.semType, {...options, debounceMs: 0}),
});

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

/** The access a row is edited under: a draft under `insert`, an existing row as its own columns say. */
const viewOf = (src, row) => src.access.peek().row(row);

/** `u2.domain.list` in the platform: the same list, the handler's rendering, Open through the handler. */
function taskList(src) {
  const list = new VirtualList({
    itemHeight: 30,
    keyOf: (row) => row.id,
    contextActions: (row) => actionsOf(src, row),
    render: (row, _index, el) => {
      el.classList.toggle('u2-domain-list-deleted', row['~state'] === 'deleted');
      return divH([
        span(row.title || 'New task', row.title ? 'u2-domain-list-name' : 'u2-domain-list-name u2-domain-draft'),
        span(row.status, 'u2-p2'),
        rowActions(actionsOf(src, row)),
      ], 'u2-domain-list-item');
    },
  });
  list.name = 'tasks';
  list.root.style.height = '260px';
  list.setItems(src.rows.items);
  list.effect(() => {
    const row = src.rows.items.peek()[list.selectedIndex.value] ?? null;
    if (row?.id !== src.currentRow.peek()?.id)
      src.currentRow.value = row;
  });
  list.effect(() => {
    const row = src.currentRow.value;
    const at = row === null ? -1 : src.rows.items.value.findIndex((r) => r.id === row.id);
    if (at !== list.selectedIndex.peek())
      list.selectedIndex.value = at;
  });
  const onScroll = () => {
    const r = list.root;
    if (r.scrollTop + r.clientHeight > r.scrollHeight - 5 * 30)
      src.loadMore();
  };
  list.root.addEventListener('scroll', onScroll);
  list.own(() => list.root.removeEventListener('scroll', onScroll));
  return list;
}

/** Permission ⇒ hidden: Delete is offered only where the row's own `~can_delete` says so; a row
 * marked deleted stays, struck through, with Restore alone. */
function actionsOf(src, row) {
  if (row['~state'] === 'deleted')
    return [{name: 'Restore', icon: 'undo', run: () => src.edit.peek()?.unmarkDeleted(row.id)}];
  const actions = [{name: 'Delete', icon: 'trash-alt', requires: 'delete', run: () => src.edit.peek()?.markDeleted(row.id)}];
  const access = viewOf(src, row);
  return actions.filter((a) => a.requires === undefined || access.can(a.requires));
}

/** `u2.domain.form` in the platform: `propertyForm` over the current row under the row's access, a
 * fresh form per row, every write through the row proxy — which is the source's `EditState`. */
function taskForm(src, host) {
  let shown;
  let sub;
  host.own(() => {
    shown?.dispose();
    sub?.unsubscribe();
  });
  host.effect(() => {
    const row = src.currentRow.value;
    src.access.value;
    shown?.dispose();
    sub?.unsubscribe();
    shown = new Scope();
    host.root.replaceChildren(Scope.runWith(shown, () => {
      if (row === null)
        return span('Select a task to edit, or add one.', 'u2-domain-form-hint');
      const form = propertyForm(src.schema.properties, row, {access: viewOf(src, row), exclude: SYSTEM});
      sub = src.edit.peek()?.onChanged.subscribe(() => form.refresh());
      return form.root;
    }));
  });
}

export async function render(main) {
  main.append(el('h1', null, 'Domain source'));
  const intro = el('p');
  intro.innerHTML = 'A domain table as data, without a server: <code>MemoryDomainBackend</code> answers ' +
    'the seam the platform fills over <code>grok.dapi.domains</code>, <code>DomainSource</code> holds the rows, ' +
    'the current row and every pending change, and plain controls sit on top — a list over ' +
    '<code>rows</code>, <code>propertyForm</code> over <code>currentRow</code> (writes go through the row into the ' +
    'edit state), a picker over the referenced table, Save as one transaction. In the platform ' +
    '<code>u2.domain.list</code>, <code>u2.domain.form</code> and <code>u2.domain.pick</code> are this composition ' +
    'plus the table\'s handler; the access rules are the same: a field the caller may not see is absent, ' +
    'one they may not write is text, an action they may not run is gone — per row, from the ' +
    '<code>~can_edit</code> / <code>~can_delete</code> / <code>~can_share</code> columns the rows carry.';
  main.append(intro);

  backends.domain = new MemoryDomainBackend(SCHEMA, {rows: {person: PEOPLE, task: TASKS.map((t) => ({...t}))}});
  const src = new DomainSource({table: 'demo.task', pageSize: 15});
  src.name = 'source';
  src.start();

  const list = taskList(src);
  const formHost = divV([], 'u2-domain-form');
  formHost.style.cssText = 'flex:1;min-width:0;padding:0 var(--dg-space-m)';
  const form = new Scope();
  taskForm(src, {root: formHost, effect: (fn) => form.effect(fn), own: (fn) => form.own(fn)});

  const save = saveButton(src);
  const discard = discardButton(src);
  const add = button('New task', () => src.newRow({title: '', status: 'Open', priority: 3}, {pristine: true}));
  add.dataset.u2Name = 'newTask';
  const status = span('', 'u2-gallery-status');
  status.dataset.u2Name = 'summary';
  form.effect(() => status.textContent = `${src.summary.value} · state ${src.state.value} · ` +
    `${src.changeCount.value} change(s) · live scopes ${Scope.liveCount}`);

  const page = divV([
    divH([add, save.root, discard.root, status], 'u2-gallery-row'),
    divH([list.root, formHost], 'u2-gallery-row'),
  ]);
  page.querySelectorAll('.u2-gallery-row').forEach((r) => r.style.cssText = 'gap:var(--dg-space-m);align-items:flex-start');
  list.root.style.flex = '1';
  main.append(page);

  main.append(el('h2', null, 'What to try'));
  const hints = el('ul');
  for (const text of [
    'Pick a task: the form follows the selection; edit a field and Save enables — one transaction for every pending row.',
    'Task #4, #11, #18… carry ~can_edit = false: their form is text. Every other row carries ~can_delete only on even ' +
      'rows: Delete appears on hover for those, never for the others (permission ⇒ hidden).',
    'New task: a pristine draft — nothing to save until you type; its fields are editable under insert.',
    'Assignee is a ref column: the picker searches the people table by name and stores the id.',
    'Scroll the list to the bottom: the next page is appended into the same collection, pending edits kept.',
    'Delete marks a row — struck through, Restore on hover — and Save removes it; Discard drops every pending change.',
  ])
    hints.append(el('li', null, text));
  main.append(hints);

  main.append(el('h2', null, 'Disposal'));
  main.append(button('Dispose', () => {
    form.dispose();
    save.dispose();
    discard.dispose();
    list.dispose();
    src.dispose();
    status.textContent = `Disposed: live scopes = ${Scope.liveCount}`;
  }));
}

import {Scope, VirtualList, BasicTable, rowActions, backends, DomainSource, SharedSession, confirmDiscard,
  MemoryDomainBackend, Filters, FilterQueryInput, signal, computed} from '../../src/index.js';
import {divH, divV, span, button} from '../../src/core/elements.js';
import {propertyForm} from '../../src/dg/forms/object-form.js';
import {Editors} from '../../src/dg/forms/editors.js';
import {DomainPick} from '../../src/dg/domain/pick.js';
import {DomainSearch} from '../../src/dg/domain/search.js';
import {saveButton, discardButton} from '../../src/dg/domain/buttons.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

for (const name of ['elements', 'inputs', 'choice', 'number', 'icons', 'buttons', 'tooltip', 'menu', 'list', 'table',
  'form', 'typeahead', 'async', 'notify', 'dialog', 'filter', 'filter-query', 'domain'])
  injectOnce(`u2-${name}-css`, `../../css/${name}.css`);

/** A two-table slice of Grit's `schema.json` — a parent and the table that refers to it. */
const SCHEMA = {
  name: 'grit',
  tables: {
    project: {
      friendlyName: 'Projects', businessKey: ['key'],
      columns: {
        key: {type: 'string', required: true, friendlyName: 'Key'},
        name: {type: 'string', required: true, isName: true, searchable: true},
        description: {type: 'string', editor: 'textarea'},
      },
    },
    issue: {
      friendlyName: 'Issues',
      columns: {
        project_id: {type: 'ref', ref: 'project', required: true, friendlyName: 'Project'},
        title: {type: 'string', required: true, isName: true, searchable: true},
        status: {type: 'string', choices: ['Open', 'Blocked', 'Done'], friendlyName: 'Status'},
        priority: {type: 'int', min: 1, max: 5, friendlyName: 'Priority'},
      },
    },
  },
};

const PROJECTS = [
  {id: 'p1', key: 'GRIT', name: 'Grit', description: 'The issue tracker itself'},
  {id: 'p2', key: 'DG', name: 'Datagrok'},
  {id: 'p3', key: 'U2', name: 'u2'},
];
const ISSUES = Array.from({length: 12}, (_, i) => ({
  id: `i${i + 1}`, project_id: PROJECTS[i % 3].id, title: `Issue #${i + 1}: something to fix`,
  status: ['Open', 'Blocked', 'Done'][i % 3], priority: (i % 5) + 1,
}));
const SYSTEM = ['id', 'version', 'created_on', 'updated_on', 'author_id'];

// the editor rule `src/dg/domain/form.ts` registers in the platform: a ref column is a picker over
// its table — and a draft parent's `~new:` id shows as the draft's own name
Editors.register({
  match: (p) => /^\w+\.\w+$/.test(p.semType ?? ''),
  create: (p, options) => new DomainPick(p.semType, {...options, debounceMs: 0}),
});

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

const viewOf = (src, row) => src.access.peek().row(row);

/** The rows as a fresh array on every pending change, so a table over them repaints a typed title. */
function liveItems(host, src) {
  const tick = signal(0);
  let sub;
  host.own(() => sub?.unsubscribe());
  host.effect(() => {
    sub?.unsubscribe();
    sub = src.edit.value?.onChanged.subscribe(() => tick.value++);
  });
  return computed(() => (tick.value, [...src.rows.items.value]));
}

/** Delete where the row allows it; a row marked deleted stays, struck through, with Restore alone. */
function actionsOf(src, row) {
  if (row['~state'] === 'deleted')
    return [{name: 'Restore', icon: 'undo', run: () => src.edit.peek()?.unmarkDeleted(row.id)}];
  const actions = [{name: 'Delete', icon: 'trash-alt', requires: 'delete', run: () => src.edit.peek()?.markDeleted(row.id)}];
  return actions.filter((a) => a.requires === undefined || viewOf(src, row).can(a.requires));
}

/** `domains.form` in the platform: a fresh `propertyForm` per current row under the row's access. */
function formOf(src, host, empty) {
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
        return span(empty, 'u2-domain-form-hint');
      const form = propertyForm(src.schema.properties, row, {access: viewOf(src, row), exclude: SYSTEM});
      sub = src.edit.peek()?.onChanged.subscribe(() => form.refresh());
      return form.root;
    }));
  });
}

/** The gate every user-initiated move goes through while the session is dirty. */
function gated(session, action, apply, revert) {
  if (!session.isDirty.peek())
    return void apply();
  void confirmDiscard(session, {action}).then((ok) => ok ? apply() : revert());
}

function projectList(src, session) {
  const list = new VirtualList({
    itemHeight: 30,
    keyOf: (row) => row.id,
    contextActions: (row) => actionsOf(src, row),
    render: (row, _index, el) => {
      el.classList.toggle('u2-domain-list-deleted', row['~state'] === 'deleted');
      const draft = !row.name;
      return divH([
        span(draft ? 'New project' : row.name, draft ? 'u2-domain-list-name u2-domain-draft' : 'u2-domain-list-name'),
        span(row.key ?? '', 'u2-p2'),
        rowActions(actionsOf(src, row)),
      ], 'u2-domain-list-item');
    },
  });
  list.name = 'projects';
  list.root.style.height = '180px';
  list.setItems(src.rows.items);
  // a click while the session is dirty asks first; a cancel puts the selection back
  list.effect(() => {
    const at = list.selectedIndex.value;
    const row = src.rows.items.peek()[at] ?? null;
    const current = src.currentRow.peek();
    if (row?.id === current?.id)
      return;
    const back = src.rows.items.peek().findIndex((r) => r.id === current?.id);
    gated(session, 'switch projects', () => src.currentRow.value = row, () => list.selectedIndex.value = back);
  });
  list.effect(() => {
    const row = src.currentRow.value;
    const at = row === null ? -1 : src.rows.items.value.findIndex((r) => r.id === row.id);
    if (at !== list.selectedIndex.peek())
      list.selectedIndex.value = at;
  });
  return list;
}

function issueTable(src, host) {
  const table = new BasicTable({
    selectable: true,
    items: liveItems(host, src),
    onRowClick: (row) => src.currentRow.value = row,
    columns: [
      {header: 'Title', render: (row) => span(row.title || 'New issue',
        row.title ? 'u2-domain-list-name' : 'u2-domain-list-name u2-domain-draft')},
      {header: 'Status', render: (row) => row.status ?? '', width: '80px'},
      {header: 'Priority', render: (row) => String(row.priority ?? ''), width: '70px', align: 'right'},
      {header: '', render: (row) => rowActions(actionsOf(src, row)), width: '60px'},
    ],
  });
  table.name = 'issues';
  table.effect(() => {
    const row = src.currentRow.value;
    table.selectedIndex.value = row === null ? -1 : src.rows.items.value.findIndex((r) => r.id === row.id);
  });
  return table;
}

/** `domains.filters` over a non-platform backend: the query box over the source's own schema, the
 * tree AND-ed with the parent condition; a change while dirty goes through the gate. */
function issueFilters(src, session, host) {
  const tree = signal(Filters.group('and'));
  let input;
  host.effect(() => {
    src.state.value;
    if (input !== undefined || src.schema.properties.length === 0)
      return;
    input = host.runInScope(() => new FilterQueryInput({schema: src.schema, target: 'domain', inline: true,
      name: 'filters', placeholder: 'Filter issues…'}));
    let accepted = input.value.peek();
    input.effect(() => {
      const next = input.value.value;
      if (Filters.equals(next, accepted))
        return;
      gated(session, 'change the filter', () => tree.value = accepted = next, () => input.value.value = accepted);
    });
    host.root.append(input.root);
  });
  return tree;
}

export async function render(main) {
  main.append(el('h1', null, 'Domain session'));
  const intro = el('p');
  intro.innerHTML = 'Two tables, one unit of work: a <code>projects</code> source and an <code>issues</code> source ' +
    'share a <code>SharedSession</code>, so <code>saveButton(session)</code> writes every pending row of both as ' +
    'ONE transaction and <code>session.summary</code> counts them together. The issues follow the selected project; ' +
    'a project that is still a draft is referenced by its child through its <code>~new:</code> id and the two ' +
    'inserts land in order. <code>DomainSearch</code> writes the source\'s <code>search</code>, the filter box ' +
    'its <code>query</code>; the filter box, a project switch and New project ask before dropping unsaved ' +
    'changes. In the platform ' +
    'this page is <code>domains.form</code>, <code>domains.grid</code>, <code>domains.search</code>, ' +
    '<code>domains.filters</code> and <code>domains.children</code> over the same two sources — or ' +
    '<code>projects.app()</code>, which builds it.';
  main.append(intro);

  backends.domain = new MemoryDomainBackend(SCHEMA, {rows: {
    project: PROJECTS.map((p) => ({...p})), issue: ISSUES.map((i) => ({...i}))}});
  const session = new SharedSession();
  const projects = new DomainSource({table: 'grit.project', pageSize: 20, session});
  const issues = new DomainSource({table: 'grit.issue', pageSize: 50, session});
  projects.name = 'projects';
  issues.name = 'issues';
  projects.start();
  issues.start();

  const page = new Scope();
  // the list starts on its first project
  page.effect(() => {
    const rows = projects.rows.items.value;
    if (projects.state.value === 'ready' && rows.length > 0 && projects.currentRow.peek() === null)
      projects.currentRow.value = rows[0];
  });

  const list = projectList(projects, session);
  const projectFormHost = divV([], 'u2-domain-form');
  projectFormHost.style.cssText = 'flex:1;min-width:0;padding:0 var(--dg-space-m)';
  formOf(projects, {root: projectFormHost, effect: (fn) => page.effect(fn), own: (fn) => page.own(fn)},
    'Select a project, or add one.');

  const filterHost = divH([], 'u2-gallery-row');
  const filterHostControl = {root: filterHost, effect: (fn) => page.effect(fn), own: (fn) => page.own(fn),
    runInScope: (fn) => Scope.runWith(page, fn)};
  const tree = issueFilters(issues, session, filterHostControl);
  const search = new DomainSearch(issues, {placeholder: 'Search issues', debounceMs: 0});
  // the issues are the selected project's, narrowed by the filter box
  page.effect(() => {
    const project = projects.currentRow.value;
    const own = tree.value;
    issues.query.value = Filters.group('and', [Filters.cond('project_id', '=', project?.id ?? ''),
      ...(own.nodes.length === 0 ? [] : [own])]);
  });

  const tableHost = divV([], 'u2-gallery-table');
  const table = issueTable(issues, {effect: (fn) => page.effect(fn), own: (fn) => page.own(fn)});
  tableHost.append(table.root);
  tableHost.style.cssText = 'flex:1;min-width:0;max-height:260px;overflow:auto';
  const issueFormHost = divV([], 'u2-domain-form');
  issueFormHost.style.cssText = 'flex:1;min-width:0;padding:0 var(--dg-space-m)';
  formOf(issues, {root: issueFormHost, effect: (fn) => page.effect(fn), own: (fn) => page.own(fn)},
    'Select an issue, or add one.');

  const save = saveButton(session);
  const discard = discardButton(session);
  const newProject = button('New project', () => gated(session, 'add a project',
    () => projects.newRow({key: '', name: ''}, {pristine: true}), () => {}));
  newProject.dataset.u2Name = 'newProject';
  const newIssue = button('New issue', () => issues.newRow(
    {project_id: projects.currentRow.peek()?.id ?? null, title: '', status: 'Open', priority: 3}, {pristine: true}));
  newIssue.dataset.u2Name = 'newIssue';
  page.effect(() => newIssue.disabled = projects.currentRow.value === null);
  const status = span('', 'u2-gallery-status');
  status.dataset.u2Name = 'summary';
  page.effect(() => status.textContent = `${session.summary.value || 'nothing to save'} · projects ` +
    `${projects.summary.value} · issues ${issues.summary.value} · live scopes ${Scope.liveCount}`);

  const rows = [
    divH([newProject, newIssue, save.root, discard.root, status], 'u2-gallery-row'),
    divH([list.root, projectFormHost], 'u2-gallery-row'),
    divH([search.root, filterHost], 'u2-gallery-row'),
    divH([tableHost, issueFormHost], 'u2-gallery-row'),
  ];
  for (const r of rows)
    r.style.cssText = 'gap:var(--dg-space-m);align-items:flex-start';
  list.root.style.flex = '1';
  filterHost.style.flex = '1';
  main.append(divV(rows));

  main.append(el('h2', null, 'What to try'));
  const hints = el('ul');
  for (const text of [
    'Rename the project and an issue, then Save: "2 unsaved changes in 2 tables" lands as one transaction, ' +
      'one balloon.',
    'New project, then New issue under it: the issue\'s Project field shows the draft by name; fill the key and ' +
      'the name and Save — the project is inserted first, the issue\'s reference re-points to its real id.',
    'Edit an issue, then click another project, change the filter or press New project: the unsaved-changes ' +
      'dialog asks first; Cancel keeps the selection and the changes. A search while dirty is simply not run.',
    'Filter with the query box (status = "Open") or search by title: the rows narrow within the selected project.',
    'Delete marks a row — struck through, Restore on hover — and Save removes it; Discard drops every pending ' +
      'change in both tables.',
  ])
    hints.append(el('li', null, text));
  main.append(hints);

  main.append(el('h2', null, 'Disposal'));
  main.append(button('Dispose', () => {
    page.dispose();
    search.dispose();
    table.dispose();
    save.dispose();
    discard.dispose();
    list.dispose();
    issues.dispose();
    projects.dispose();
    status.textContent = `Disposed: live scopes = ${Scope.liveCount}`;
  }));
}

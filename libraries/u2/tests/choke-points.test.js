import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readdirSync, readFileSync} from 'node:fs';
import {join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {createRequire} from 'node:module';

const ts = createRequire(import.meta.url)('typescript');
const root = fileURLToPath(new URL('..', import.meta.url));

function walk(dir, keep, out = []) {
  for (const e of readdirSync(join(root, dir), {withFileTypes: true})) {
    const rel = `${dir}/${e.name}`;
    if (e.isDirectory())
      walk(rel, keep, out);
    else if (keep(rel))
      out.push(rel);
  }
  return out;
}

const SOURCES = walk('src', (f) => f.endsWith('.ts') && !/\.[dg]\.ts$/.test(f));
const SHEETS = walk('css', (f) => f.endsWith('.css'));
const read = (file) => readFileSync(join(root, file), 'utf8');

const blanked = new Map();

/** The text with its comments blanked and its line breaks kept, so a pattern sees code and strings only. */
function code(file, text) {
  if (blanked.get(file)?.text === text)
    return blanked.get(file).code;
  const out = text.split('');
  const blank = (from, to) => {
    for (let i = from; i < to; i++)
      if (out[i] !== '\n')
        out[i] = ' ';
  };
  if (file.endsWith('.css')) {
    for (const m of text.matchAll(/\/\*[\s\S]*?\*\//g))
      blank(m.index, m.index + m[0].length);
  }
  else {
    const sf = ts.createSourceFile(file, text, ts.ScriptTarget.Latest, false, ts.ScriptKind.TS);
    const visit = (node) => {
      for (const r of [...ts.getLeadingCommentRanges(text, node.pos) ?? [], ...ts.getTrailingCommentRanges(text, node.end) ?? []])
        blank(r.pos, r.end);
      for (const child of node.getChildren(sf))
        visit(child);
    };
    visit(sf);
  }
  const result = out.join('');
  blanked.set(file, {text, code: result});
  return result;
}

/** A declaration carrying `/* operable: <why> *\/` on its line or the line above. */
function operable(raw, at) {
  const line = raw.lastIndexOf('\n', at - 1);
  const end = raw.indexOf('\n', at);
  return /\/\*\s*operable:/.test(raw.slice(raw.lastIndexOf('\n', line - 1) + 1, end < 0 ? raw.length : end));
}

const DG_SOURCES = SOURCES.filter((f) => f.startsWith('src/dg/'));

const CHOKE_POINTS = [
  {
    id: 'owner policy',
    what: '`Scope.ambient` read(s)',
    files: SOURCES,
    pattern: /\bScope\.ambient\b/g,
    owners: ['src/core/owner.ts'],
    use: 'resolve the owner through Ownership.scopeOf(owner, what) or Ownership.adopt(owner, dispose) ' +
      '(src/core/owner.ts): one owner policy, not one more (CR4 C-1)',
    pinned: {
      'src/components/actions/buttons.ts': 4,
      'src/components/collections/table.ts': 1,
      'src/components/containers/dialog.ts': 2,
      'src/components/containers/tabs.ts': 1,
      'src/components/containers/wizard.ts': 1,
      'src/components/display/badge.ts': 1,
      'src/components/display/icon.ts': 1,
      'src/components/display/notify.ts': 1,
      'src/components/forms/property-grid.ts': 2,
      'src/components/inputs/range-slider.ts': 1,
      'src/components/navigation/menu.ts': 1,
      'src/core/elements.ts': 3,
      'src/core/input-base.ts': 1,
      'src/core/region.ts': 1,
      'src/dg/designer/prop-editors.ts': 1,
      'src/dg/funcs/param-sources.ts': 2,
      'src/dg/funcs/param-tables.ts': 1,
      'src/dg/viewers/viewer-control.ts': 1,
      'src/spec/spec.ts': 5,
    },
  },
  {
    id: 'thenable',
    what: 'duck-typed `.then` read(s)',
    files: SOURCES,
    pattern: /\.then\b(?!\s*\()|['"]then['"]\s+in\b/g,
    owners: [],
    use: 'ask Scope.thenable(x): reading `.then` on a platform proxy throws ("No such column: then" on a DG.Row, N-3)',
    pinned: {},
  },
  {
    id: 'aria-label',
    what: 'hand-written `aria-label`(s)',
    files: SOURCES,
    pattern: /setAttribute\(\s*['"]aria-label['"]|\.ariaLabel\s*=(?!=)/g,
    owners: ['src/components/actions/buttons.ts'],
    use: 'name a control through its tooltip (the `tooltip` of button/iconButton): src/components/actions/buttons.ts ' +
      'writes the name and follows a live tooltip, a hand-set aria-label goes stale (CR4 C-4)',
    pinned: {
      'src/components/actions/button-group.ts': 1,
      'src/components/collections/catalog-table.ts': 1,
      'src/components/collections/functions-browser.ts': 3,
      'src/components/collections/master-detail.ts': 1,
      'src/components/collections/tree.ts': 1,
      'src/components/containers/dialog.ts': 1,
      'src/components/display/badge.ts': 1,
      'src/components/display/icon.ts': 1,
      'src/components/filter/filter-builder.ts': 1,
      'src/components/filter/filter-row.ts': 1,
      'src/components/forms/access-grid.ts': 4,
      'src/components/inputs/choice-input.ts': 1,
      'src/components/inputs/color-input.ts': 1,
      'src/components/inputs/combobox.ts': 1,
      'src/components/inputs/date-input.ts': 4,
      'src/components/inputs/font-input.ts': 4,
      'src/components/inputs/list-input.ts': 2,
      'src/components/inputs/map-input.ts': 2,
      'src/components/inputs/multi-select.ts': 1,
      'src/components/inputs/number-input.ts': 1,
      'src/components/inputs/range-slider.ts': 3,
      'src/components/inputs/text-input.ts': 2,
      'src/components/inputs/typeahead.ts': 1,
      'src/components/navigation/breadcrumbs.ts': 2,
      'src/components/navigation/toolbar.ts': 2,
      'src/core/input-base.ts': 3,
      'src/dg/entities/func-call-history-browser.ts': 1,
      'src/dg/inputs/columns.ts': 3,
      'src/dg/shell/table-view-chrome.ts': 1,
    },
  },
  {
    id: 'aria-description',
    what: 'hand-written `aria-description`(s)',
    files: SOURCES,
    pattern: /setAttribute\(\s*['"]aria-description['"]|\.ariaDescription\s*=(?!=)/g,
    owners: ['src/core/reason.ts'],
    use: 'say why a control is off through `disabled: \'<reason>\'` or Reason.apply(el, reason) (src/core/reason.ts), ' +
      'which writes `aria-disabled` and the reason as the description together',
    pinned: {
      'src/components/inputs/choice-chip.ts': 1,
      'src/core/input-base.ts': 1,
    },
  },
  {
    id: 'overflow',
    what: 'overflow glyph or label literal(s)',
    files: SOURCES,
    pattern: /['"`](?:ellipsis-[hv](?:-alt)?|More actions|More|…|⋯|\\u2026|\\u22ef)['"`]/g,
    owners: ['src/components/navigation/toolbar.ts'],
    use: 'use Toolbar.OVERFLOW_ICON and Toolbar.OVERFLOW_LABEL (src/components/navigation/toolbar.ts): one overflow ' +
      'glyph and one label (CR4 C-3)',
    pinned: {
      'src/components/actions/buttons.ts': 1,
      'src/components/filter/filter-row.ts': 1,
      'src/components/navigation/breadcrumbs.ts': 1,
      'src/dg/designer/bind-picker.ts': 1,
    },
  },
  {
    id: 'sort comparator',
    what: 'text comparator(s)',
    files: SOURCES,
    pattern: /\.localeCompare\s*\(|\bIntl\.Collator\b/g,
    owners: ['src/components/collections/catalog-table.ts'],
    use: 'order text with the one comparator (COLLATOR, src/components/collections/catalog-table.ts): a second one ' +
      'puts "Item 10" before "Item 2" in one list and after it in the next (CR4 C-2)',
    pinned: {
      'src/dg/designer/func-picker.ts': 1,
      'src/dg/designer/palette.ts': 1,
      'src/dg/entities/functions-browser.ts': 1,
    },
  },
  {
    id: 'sort order',
    what: 'hand-parsed sort string(s)',
    files: SOURCES,
    pattern: /\.startsWith\(\s*['"]!['"]\s*\)|`!\$\{/g,
    owners: ['src/components/collections/sort-order.ts'],
    use: 'read and write a table\'s order through SortOrder.parse / of / next (src/components/collections/sort-order.ts): ' +
      'one sort contract for DataTable, BasicTable and CatalogTable (CR4 C-2)',
    pinned: {},
  },
  {
    id: 'handlers',
    what: 'raw listener(s)',
    files: SOURCES,
    pattern: new RegExp('\\.addEventListener\\s*\\(|\\.on(?:click|dblclick|input|change|keydown|keyup|keypress|' +
      'mousedown|mouseup|mousemove|mouseenter|mouseleave|mouseover|mouseout|pointerdown|pointerup|pointermove|' +
      'pointerenter|pointerleave|focus|blur|focusin|focusout|scroll|wheel|contextmenu|submit|dragstart|dragover|' +
      'dragenter|dragleave|dragend|drop|paste|copy|cut|resize|load|error)\\s*=(?!=)', 'g'),
    owners: ['src/core/elements.ts'],
    use: 'listen through scope.listen(target, type, handler) or an element factory\'s handler (src/core/elements.ts): ' +
      'the handler then runs through Scope.bind — a throw reported, a returned promise awaited, removed on dispose ' +
      '(CR4 C-7)',
    pinned: {
      'src/core/reason.ts': 1,
      'src/components/actions/buttons.ts': 1,
      'src/components/collections/data-table.ts': 2,
      'src/components/collections/functions-browser.ts': 4,
      'src/components/collections/tree.ts': 2,
      'src/components/containers/card.ts': 1,
      'src/components/display/badge.ts': 1,
      'src/components/display/notify.ts': 4,
      'src/components/filter/filter-dnd.ts': 5,
      'src/components/forms/form.ts': 1,
      'src/components/inputs/pick-input.ts': 2,
      'src/components/inputs/range-slider.ts': 2,
      'src/core/input-base.ts': 2,
      'src/core/tooltip.ts': 8,
      'src/dg/designer/bind-chip.ts': 1,
      'src/dg/designer/look-grid.ts': 1,
      'src/dg/entities/entity.ts': 2,
      'src/dg/entities/func-call-history-browser.ts': 3,
      'src/dg/funcs/func-form.ts': 2,
      'src/dg/funcs/param-tables.ts': 1,
      'src/dg/inputs/column-combo.ts': 6,
      'src/dg/inputs/columns.ts': 6,
      'src/dg/inputs/file-input.ts': 1,
      'src/dg/inputs/files-input.ts': 1,
      'src/dg/inputs/from-dart-input.ts': 1,
      'src/dg/inputs/func-call-input.ts': 2,
      'src/dg/inputs/pickers.ts': 2,
      'src/dg/shell/collection-selection.ts': 3,
      'src/dg/shell/view-sync.ts': 1,
    },
  },
  {
    id: 'table names',
    what: '`grok.shell.tableNames` read(s)',
    files: SOURCES,
    pattern: /\bshell\.tableNames\b/g,
    owners: ['src/dg/shell/platform.ts'],
    use: 'list tables through one normalizer that drops nameless ones (Platform.tableNames, src/dg/shell/platform.ts): ' +
      'a nameless table broke every funcForm (B-3, N-2)',
    pinned: {},
  },
  {
    id: 'column listing',
    what: 'raw column listing(s)',
    files: DG_SOURCES,
    pattern: /\.columns\s*\.\s*(?:names|toList)\s*\(/g,
    owners: ['src/dg/shell/platform.ts'],
    use: 'list a frame\'s columns by the platform\'s rule, no `~` and no empty names (Platform.listedColumns, ' +
      'src/dg/shell/platform.ts; core column_grid.dart:909-911): a picker offering `~canonical smiles` defaults to a ' +
      'hidden column (B-4)',
    pinned: {
      'src/dg/viewers/registrations.ts': 1,
    },
  },
  {
    id: 'formatter',
    what: 'ad-hoc value format(s)',
    files: SOURCES,
    pattern: /\.toLocale(?:Date|Time)?String\s*\(|['"](?:Yes|No)['"]/g,
    owners: ['src/core/format.ts'],
    use: 'format a value through one formatter (Format.value, src/core/format.ts): five formatters print one date five ' +
      'ways and a group header reads unlike its cells (N-5)',
    pinned: {},
  },
  {
    id: 'inert',
    what: '`inert` write(s)',
    files: SOURCES,
    pattern: /\.inert\s*=(?!=)|(?:set|toggle)Attribute\(\s*['"]inert['"]/g,
    owners: ['src/core/reason.ts'],
    use: 'make a region inert with its reason shown (Reason.region(el, reason), src/core/reason.ts): a greyed pane with ' +
      'no reason is dead and nothing says why',
    pinned: {},
  },
  {
    id: 'pointer-events',
    what: 'unexplained `pointer-events: none`',
    files: SHEETS,
    pattern: /pointer-events\s*:\s*none\b/g,
    allowed: operable,
    owners: [],
    use: 'a control under pointer-events: none is dead and nothing says so (B-1): keep it operable, or say why nobody ' +
      'needs to click it in a /* operable: … */ comment on the declaration',
    pinned: {
      'css/color.css': 1,
      'css/designer.css': 3,
      'css/elements.css': 1,
      'css/filter.css': 1,
      'css/inputs.css': 2,
      'css/list.css': 1,
      'css/message-input.css': 1,
      'css/notify.css': 1,
      'css/tour.css': 1,
      'css/typeahead.css': 1,
      'css/wizard.css': 1,
    },
  },
];

function findings(point, source) {
  const out = new Map();
  for (const file of point.files) {
    if (point.owners.includes(file))
      continue;
    const raw = source(file);
    const text = code(file, raw);
    const lines = [];
    for (const m of text.matchAll(point.pattern))
      if (!point.allowed?.(raw, m.index))
        lines.push(text.slice(0, m.index).split('\n').length);
    if (lines.length)
      out.set(file, lines);
  }
  return out;
}

function problems(point, source = read) {
  const found = findings(point, source);
  const out = [];
  for (const [file, lines] of found) {
    const pin = point.pinned[file] ?? 0;
    if (lines.length > pin)
      out.push(`${file}:${lines.join(',')}: ${lines.length} ${point.what}, ${pin} pinned — ${point.use}`);
  }
  for (const [file, pin] of Object.entries(point.pinned)) {
    const n = found.get(file)?.length ?? 0;
    if (n < pin) {
      out.push(`${file}: ${n} ${point.what}, ${pin} pinned — pins only shrink: ` +
        `${n === 0 ? 'delete its pin' : `pin it at ${n}`} in tests/choke-points.test.js ('${point.id}')`);
    }
  }
  return out;
}

const point = (id) => CHOKE_POINTS.find((p) => p.id === id);
const withText = (file, change) => (f) => f === file ? change(read(f)) : read(f);
const appended = (file, text) => withText(file, (t) => `${t}\n${text}\n`);

for (const p of CHOKE_POINTS) {
  test(`choke point "${p.id}": ${p.what} only in ${[...p.owners, 'the pinned sites'].join(', ')}`, () => {
    assert.deepEqual(problems(p), []);
  });
}

const MUTANTS = [
  ['owner policy', 'src/components/display/badge.ts', 'export const mutant = () => Scope.ambient;'],
  ['thenable', 'src/components/containers/dialog.ts', 'const mutant = (r: any) => typeof r?.then === \'function\';'],
  ['aria-label', 'src/components/display/badge.ts',
    'export const mutant = (el: HTMLElement, tip: string) => el.setAttribute(\'aria-label\', tip);'],
  ['aria-description', 'src/components/containers/section.ts',
    'export const mutant = (el: HTMLElement, why: string) => el.setAttribute(\'aria-description\', why);'],
  ['overflow', 'src/components/navigation/menu-bar.ts',
    'export const mutant = () => iconButton(\'ellipsis-v\', () => {}, {tooltip: \'More\'});'],
  ['sort comparator', 'src/components/collections/list.ts',
    'export const mutant = (a: string[]) => a.sort((x, y) => x.localeCompare(y));'],
  ['sort order', 'src/components/collections/data-table.ts',
    'export const mutant = (sort: string) => sort.startsWith(\'!\') ? sort.slice(1) : sort;'],
  ['handlers', 'src/components/display/badge.ts',
    'export const mutant = (el: HTMLElement, run: () => Promise<void>) => el.addEventListener(\'click\', run);'],
  ['table names', 'src/dg/inputs/columns.ts', 'export const mutant = () => grok.shell.tableNames;'],
  ['column listing', 'src/dg/inputs/file-input.ts', 'export const mutant = (df: DG.DataFrame) => df.columns.names();'],
  ['formatter', 'src/components/collections/table.ts', 'export const mutant = (d: Date) => d.toLocaleDateString();'],
  ['inert', 'src/components/containers/section.ts', 'export const mutant = (el: HTMLElement) => { el.inert = true; };'],
  ['pointer-events', 'css/inputs.css', '.u2-mutant > .u2-input-options { pointer-events: none; }'],
];

for (const [id, file, text] of MUTANTS) {
  test(`mutation: a new ${point(id).what.replace(/\(s\)$/, '')} in ${file} fails the "${id}" choke point`, () => {
    const found = problems(point(id), appended(file, text));
    assert.equal(found.length, 1, found.join('\n'));
    assert.ok(found[0].startsWith(`${file}:`), found[0]);
    assert.ok(found[0].endsWith(point(id).use), found[0]);
  });
}

test('mutation: a comment or a doc block naming a mechanism is not a use of it', () => {
  const file = 'src/components/display/badge.ts';
  const text = '// Scope.ambient, el.addEventListener(\'x\')\n/** `d.toLocaleDateString()` and \'More\' */';
  for (const id of ['owner policy', 'handlers', 'formatter', 'overflow'])
    assert.deepEqual(problems(point(id), appended(file, text)), [], id);
  const quoted = problems(point('owner policy'), appended(file, 'const q = \'/*\', s = () => Scope.ambient; // */'));
  assert.equal(quoted.length, 1, 'a comment opener inside a string hides nothing');
});

test('mutation: an owner may use its own mechanism', () => {
  assert.deepEqual(problems(point('owner policy'), appended('src/core/owner.ts', 'export const again = () => Scope.ambient;')), []);
  assert.deepEqual(problems(point('overflow'), appended('src/components/navigation/toolbar.ts', 'const glyph = \'…\';')), []);
});

test('mutation: a pointer-events rule that says why it is not operable passes', () => {
  assert.deepEqual(problems(point('pointer-events'),
    appended('css/inputs.css', '.u2-mutant::after { pointer-events: none; /* operable: a decoration over the field */ }')), []);
  assert.deepEqual(problems(point('pointer-events'),
    appended('css/inputs.css', '.u2-mutant::after {\n  /* operable: a decoration over the field */\n  pointer-events: none;\n}')), []);
});

test('mutation: a pinned site that goes away must take its pin with it', () => {
  const p = point('handlers');
  const [file, pin] = Object.entries(p.pinned).find(([, n]) => n > 1);
  const found = problems(p, withText(file, (t) => t.replace('.addEventListener(', '.removeEventListener(')));
  assert.deepEqual(found, [`${file}: ${pin - 1} ${p.what}, ${pin} pinned — pins only shrink: ` +
    `pin it at ${pin - 1} in tests/choke-points.test.js ('handlers')`]);
});

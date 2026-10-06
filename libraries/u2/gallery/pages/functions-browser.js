import {signal} from '../../src/index.js';
import {FunctionsBrowser} from '../../src/components/collections/functions-browser.js';

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function ensureCss() {
  const names = ['functions-browser', 'list', 'elements', 'inputs', 'buttons',
    'icons', 'icon-input', 'async', 'tooltip', 'menu'];
  for (const name of names) {
    if (document.getElementById(`u2-${name}-css`)) continue;
    const link = el('link');
    link.id = `u2-${name}-css`;
    link.rel = 'stylesheet';
    link.href = new URL(`../../css/${name}.css`, import.meta.url).href;
    document.head.append(link);
  }
}

const PACKAGES = ['Chem', 'Bio', 'ML', 'Math', 'Core', 'Charts'];
const VERBS = ['Find', 'Compute', 'Render', 'Convert', 'Predict', 'Align', 'Cluster', 'Export'];
const NOUNS = ['Similarity', 'Descriptors', 'Sequence', 'Model', 'Spectrum', 'Structure', 'Table', 'Report'];
const TAGS = ['chem', 'bio', 'ml', 'math', 'transform', 'viewers', 'widgets', 'search'];
const SIGNATURES = ['(table) : dataframe', '(x) : double', '(molecule) : string',
  '(sequence, cutoff) : dataframe', '() : void', '(a, b) : int'];
const ROLES = [
  {role: 'app', description: 'An application that gets shown in the app store.'},
  {role: 'panel', description: 'Context-specific widget on the context panel.'},
  {role: 'init', description: 'Runs when the containing package is initialized.'},
];

function mockItems(count) {
  const items = [];
  for (let i = 0; i < count; i++) {
    const pkg = PACKAGES[i % PACKAGES.length];
    const label = `${VERBS[i % VERBS.length]} ${NOUNS[(i * 3 + 1) % NOUNS.length]} ${i + 1}`;
    const name = `${pkg}:${label.replaceAll(' ', '')}`;
    const tags = [TAGS[i % TAGS.length]];
    if (i % 5 === 0) tags.push(TAGS[(i + 3) % TAGS.length]);
    if (i % 29 === 0) tags.push('internal');
    const roles = i % 7 === 0 ? [...tags, ROLES[i % ROLES.length].role] : tags;
    const description = `A mock ${tags[0]} function that ${label.toLowerCase()}s.`;
    items.push({
      name, label, tags, roles, description,
      signature: SIGNATURES[i % SIGNATURES.length],
      search: `${name} ${label} ${description}`.toLowerCase(),
    });
  }
  items.sort((a, b) => a.label.localeCompare(b.label));
  return items;
}

export function render(main) {
  ensureCss();
  main.append(el('h1', null, 'Functions browser'));
  const intro = el('p');
  intro.innerHTML = 'The u2 counterpart of the Dart function registry browser: plain terms AND ' +
    'over the precomputed search string, <code>#tag</code> / <code>@role</code> terms mirrored ' +
    'into the checkbox panes and back (a label click checks exclusively, clicking the sole ' +
    'checked one unchecks it), per-tag counts over the term-filtered set, and a virtualized ' +
    'list. 200 mock items here; the dg factory feeds it every registered <code>DG.Func</code>.';
  main.append(intro);

  const log = el('p', 'u2-gallery-status');
  const fb = new FunctionsBrowser({
    items: mockItems(200),
    search: '#chem',
    roles: ROLES,
    runAction: (item) => log.textContent = `Ran ${item.name}`,
    contextActions: (item) => [
      {name: 'Run', icon: 'play', run: () => log.textContent = `Ran ${item.name}`},
      {name: 'Copy name', icon: 'copy', run: () => log.textContent = `Copied "${item.name}"`},
    ],
    onActivate: (item) => log.textContent = `Activated ${item.name} (dblclick / Enter)`,
  });
  fb.root.style.height = '380px';
  fb.root.style.border = '1px solid var(--dg-border)';
  main.append(fb.root);

  const readout = el('p', 'u2-gallery-status');
  fb.effect(() => readout.textContent =
    `query = "${fb.query.value}" · checked tags = [${[...fb.checkedTags.value].join(', ')}]` +
    ` · selected = ${fb.selected.value?.name ?? 'none'}`);
  main.append(readout, log);

  main.append(el('h2', null, 'Live toggles'));
  const toggles = el('div');
  toggles.style.cssText = 'display:flex;gap:16px;flex-wrap:wrap';
  for (const name of ['showSearch', 'showTags', 'showSignature', 'showRunButton']) {
    const label = el('label');
    const box = el('input');
    box.type = 'checkbox';
    box.checked = fb[name].value;
    box.addEventListener('change', () => fb[name].value = box.checked);
    fb.effect(() => box.checked = fb[name].value);
    label.append(box, ` ${name}`);
    toggles.append(label);
  }
  main.append(toggles);
  const note = el('p', 'u2-gallery-status');
  note.textContent = 'The show* members are writable signals — the filter icon in the search ' +
    'row writes showTags too, so the checkbox above follows it.';
  main.append(note);

  main.append(el('h2', null, 'Widget mode'));
  const compact = new FunctionsBrowser({
    items: mockItems(200),
    showSignature: signal(false),
    showRunButton: false,
    visibleTags: ['chem', 'bio', 'ml'],
  });
  compact.root.style.height = '220px';
  compact.root.style.border = '1px solid var(--dg-border)';
  main.append(compact.root);
}

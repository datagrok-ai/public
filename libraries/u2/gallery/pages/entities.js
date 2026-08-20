import {Scope, bindText, computed} from '../../src/index.js';
import {TypeAhead} from '../../src/components/typeahead.js';
import {Tooltip} from '../../src/core/tooltip.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-typeahead-css', '../../css/typeahead.css');
injectOnce('u2-tooltip-css', '../../css/tooltip.css');
injectOnce('u2-entity-css', '../../css/entity.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

// A mock ObjectRenderer — the platform-free stand-in for u2/dg's handlerRenderer, which
// adapts the platform's ObjectHandler registry to the same interface.
const PROJECTS = [
  {name: 'Chem', owner: 'Ada Almeida', entities: 412},
  {name: 'Bio', owner: 'Bruno Bauer', entities: 267},
  {name: 'Curves', owner: 'Chen Costa', entities: 96},
  {name: 'HitTriage', owner: 'Dmitri Dubois', entities: 154},
  {name: 'PowerPack', owner: 'Elena Egan', entities: 320},
];

const projectRenderer = {
  caption: (p) => p.name,
  listItem: (p) => {
    const row = el('div', 'u2-typeahead-user');
    const text = el('div', 'u2-typeahead-user-text');
    text.append(el('div', 'u2-typeahead-user-name', p.name),
      el('div', 'u2-typeahead-user-secondary', `${p.owner} · ${p.entities} entities`));
    row.append(el('div', 'u2-avatar u2-avatar-2', p.name.substring(0, 2).toUpperCase()), text);
    return row;
  },
  markup: (p) => {
    const chip = el('span');
    chip.append(el('span', 'u2-avatar u2-avatar-2', p.name.substring(0, 2).toUpperCase()),
      document.createTextNode(` ${p.name}`));
    return chip;
  },
  tooltip: (p) => `${p.name} — owned by ${p.owner}, ${p.entities} entities`,
};

export async function render(main) {
  main.append(el('h1', null, 'Object renderers'));
  const intro = el('p');
  intro.innerHTML = 'The <code>ObjectRenderer&lt;T&gt;</code> protocol: caption / icon / ' +
    'list item / markup / tooltip / card as one object, accepted wherever a control takes ' +
    'per-item callbacks. In the platform, <code>u2/dg</code>\'s <code>handlerRenderer()</code> ' +
    'implements it over the ObjectHandler registry, so any handler-backed entity renders on ' +
    'u2 surfaces; this page uses a mock renderer — the core seam is platform-free.';
  main.append(intro);

  main.append(el('h2', null, 'TypeAhead with renderer — no itemText, no render callback'));
  const picker = new TypeAhead({
    source: PROJECTS,
    renderer: projectRenderer,
    placeholder: 'Project…',
  });
  picker.root.style.width = '280px';
  const selected = el('code');
  bindText(picker.scope, selected,
    computed(() => picker.selected.value ? picker.selected.value.name : 'null'));
  const status = el('p', 'u2-gallery-status');
  status.append('selected = ', selected);
  main.append(picker.root, status);

  main.append(el('h2', null, 'Inline markup + tooltip — what u2/dg\'s chip() builds on'));
  main.append(el('p', 'u2-gallery-status',
    'Hover for the renderer tooltip. In the platform the chip also sets grok.shell.o on ' +
    'click and stamps data-entity automation ids.'));
  const chips = el('p');
  const chipScope = new Scope();
  for (const p of PROJECTS) {
    const c = el('span', 'u2-chip');
    c.dataset.entity = 'Project';
    c.append(projectRenderer.markup(p));
    Tooltip.bind(c, () => projectRenderer.tooltip(p), chipScope);
    chips.append(c, ' ');
  }
  Scope.ambient?.own(() => chipScope.dispose());
  main.append(chips);
}

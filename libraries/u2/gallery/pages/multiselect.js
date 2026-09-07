import {computed, Scope, bindText} from '../../src/index.js';
import {MultiSelect} from '../../src/components/inputs/multi-select.js';
import {TextInput} from '../../src/components/inputs/text-input.js';
import {Form} from '../../src/components/forms/form.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-inputs-css', '../../css/inputs.css');
injectOnce('u2-form-css', '../../css/form.css');
injectOnce('u2-badge-css', '../../css/badge.css');
injectOnce('u2-multi-select-css', '../../css/multi-select.css');

const DESCRIPTORS = [
  ['mw', 'Molecular weight'], ['exactMass', 'Exact mass'], ['logP', 'LogP'], ['logD', 'LogD'],
  ['tpsa', 'TPSA'], ['hba', 'H-bond acceptors'], ['hbd', 'H-bond donors'],
  ['rotB', 'Rotatable bonds'], ['rings', 'Rings'], ['aromaticRings', 'Aromatic rings'],
  ['heavyAtoms', 'Heavy atoms'], ['charge', 'Formal charge'], ['csp3', 'Fraction Csp3'],
  ['mr', 'Molar refractivity'], ['qed', 'QED'], ['sa', 'Synthetic accessibility'],
  ['lipinski', 'Lipinski violations'], ['stereo', 'Stereocenters'], ['halogens', 'Halogen count'],
  ['basicN', 'Basic nitrogens'], ['acidic', 'Acidic groups'], ['assemblies', 'Ring assemblies'],
  ['spiro', 'Spiro atoms'], ['bridgeheads', 'Bridgehead atoms'], ['hetRatio', 'Heteroatom ratio'],
].map(([value, label]) => ({value, label}));

const SERIES = ['Series A', 'Series B', 'Series C', 'Series D', 'Series E'];

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function button(text, onClick) {
  const b = el('button', null, text);
  b.addEventListener('click', onClick);
  return b;
}

function readout(ms) {
  const value = el('code');
  const state = el('code');
  bindText(ms.scope, value, computed(() => JSON.stringify(ms.value.value)));
  bindText(ms.scope, state, computed(() =>
    `${ms.machineState.value}${ms.isOpen.value ? ', open' : ''}, ${ms.value.value.length} selected`));
  const row = el('p', 'u2-gallery-status');
  row.append('value = ', value, '  ·  machine = ', state);
  return row;
}

export async function render(main) {
  main.append(el('h1', null, 'Multi-select'));
  const intro = el('p');
  intro.innerHTML = 'Popup multi-select: the closed field shows the picks as tags with per-tag ✕, ' +
    'the popup is a checkbox list. Toggling never closes the popup — outside click, Esc or Enter do. ' +
    'Keys: ↓/Enter/Space open, ↑↓ move the active row, Space toggles it, Enter closes, Tab closes and ' +
    'keeps the state. <code>MultiChoiceInput</code> stays the inline always-open variant.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  main.append(el('h2', null, '25 descriptors — filter and select all/none'));
  main.append(el('p', 'u2-gallery-status',
    'The filter box appears by default above eight items; "Select all" acts on the rows the filter ' +
    'leaves visible. The value is always a fresh array in item order, never in click order.'));
  const descriptors = new MultiSelect({
    label: 'Descriptors',
    items: DESCRIPTORS,
    value: ['mw', 'logP'],
    placeholder: 'Pick descriptors…',
  });
  main.append(descriptors.root, readout(descriptors));

  main.append(el('h2', null, 'Five items — no filter, maxTags: 2'));
  main.append(el('p', 'u2-gallery-status',
    'Below the threshold the popup has no filter box. Tags past maxTags collapse into "+N more".'));
  const series = new MultiSelect({
    label: 'Series',
    items: SERIES,
    maxTags: 2,
    value: ['Series A', 'Series B', 'Series C', 'Series D'],
    placeholder: 'Pick series…',
  });
  main.append(series.root, readout(series));

  main.append(el('h2', null, 'Per-tag removal'));
  main.append(el('p', 'u2-gallery-status',
    'The ✕ on a tag drops that one value without opening the popup — the fastest way to undo one pick.'));
  const assays = new MultiSelect({
    label: 'Assays',
    items: ['IC50', 'EC50', 'Ki', 'Solubility', 'Permeability', 'hERG'],
    value: ['IC50', 'Ki', 'hERG'],
    maxTags: 6,
    placeholder: 'Pick assays…',
  });
  main.append(assays.root, readout(assays));

  main.append(el('h2', null, 'Free tags (allowCustom)'));
  main.append(el('p', 'u2-gallery-status',
    'With allowCustom the popup always carries the filter box, and Enter turns whatever is typed ' +
    'into a tag — the parity path for Dart\'s TagsInput. A typed word that already exists selects ' +
    'that item instead of adding a second one, ignoring case.'));
  const tags = new MultiSelect({
    label: 'Tags',
    items: ['acid', 'base', 'zwitterion'],
    allowCustom: true,
    maxTags: 6,
    placeholder: 'Pick or type a tag…',
  });
  main.append(tags.root, readout(tags));

  main.append(el('h2', null, 'Validation'));
  main.append(el('p', 'u2-gallery-status',
    'Validators run through the input base — the same machinery every u2 input uses.'));
  const axes = new MultiSelect({
    label: 'Plot axes',
    items: DESCRIPTORS.slice(0, 8),
    placeholder: 'Pick at least two…',
  });
  axes.addValidator((v) => v.length < 2 ? 'Pick at least 2 descriptors' : null);
  main.append(axes.root, readout(axes));

  main.append(el('h2', null, 'Programmatic value and items'));
  main.append(el('p', 'u2-gallery-status',
    'Writing the value re-renders the tags and the checkboxes (open the popup and watch both). ' +
    'setItems prunes values that no longer exist, without counting as a user edit.'));
  const buttons = el('p');
  buttons.append(
    button('Set Lipinski set', () => descriptors.value.value = ['mw', 'logP', 'hba', 'hbd']),
    ' ',
    button('Clear', () => descriptors.value.value = []),
    ' ',
    button('setItems: first 6 only', () => descriptors.setItems(DESCRIPTORS.slice(0, 6))),
    ' ',
    button('setItems: all 25', () => descriptors.setItems(DESCRIPTORS)));
  main.append(buttons);

  main.append(el('h2', null, 'Inside a form'));
  main.append(el('p', 'u2-gallery-status',
    'The field is an ordinary editor slot, so the label column aligns with its neighbours and the ' +
    'form aggregates its validity.'));
  const form = new Form();
  const name = new TextInput({label: 'Report name', value: 'Q3 profiling'});
  const columns = new MultiSelect({
    label: 'Columns',
    items: DESCRIPTORS,
    value: ['mw', 'tpsa', 'qed'],
    placeholder: 'Pick columns…',
  });
  columns.addValidator((v) => v.length === 0 ? 'At least one column is required' : null);
  form.addAll([name, columns]);
  const formStatus = el('code');
  bindText(form.scope, formStatus, computed(() =>
    `{"Report name": ${JSON.stringify(name.value.value)}, "Columns": ` +
    `${JSON.stringify(columns.value.value)}} · validity = ${form.validity.value ?? 'ok'}`));
  const formLine = el('p', 'u2-gallery-status');
  formLine.append('form = ', formStatus);
  main.append(form.root, formLine);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Dispose closes any open popup and releases every effect and listener.'));
  main.append(button('Dispose multi-selects', () => {
    for (const ms of [descriptors, series, assays, tags, axes, columns, name, form])
      ms.dispose();
    refresh();
  }));
  refresh();
}

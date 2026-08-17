import {computed, Scope, bindText} from '../../src/index.js';

// The library imports @floating-ui/dom by bare specifier (correct for plugin bundlers); the
// gallery runs as raw ESM, so it needs an import map. Belongs in index.html once its owner
// adds it — kept here so this page is self-sufficient.
const vendor = (p) => new URL(`../../node_modules/${p}`, import.meta.url).href;

function injectOnce(id, make) {
  if (document.getElementById(id)) return;
  const e = make();
  e.id = id;
  document.head.append(e);
}

injectOnce('u2-floating-ui-importmap', () => {
  const s = document.createElement('script');
  s.type = 'importmap';
  s.textContent = JSON.stringify({imports: {
    '@floating-ui/dom': vendor('@floating-ui/dom/dist/floating-ui.dom.mjs'),
    '@floating-ui/core': vendor('@floating-ui/core/dist/floating-ui.core.mjs'),
    '@floating-ui/utils': vendor('@floating-ui/utils/dist/floating-ui.utils.mjs'),
    '@floating-ui/utils/dom': vendor('@floating-ui/utils/dist/floating-ui.utils.dom.mjs'),
  }});
  return s;
});

injectOnce('u2-combobox-css', () => {
  const l = document.createElement('link');
  l.rel = 'stylesheet';
  l.href = new URL('../../css/combobox.css', import.meta.url).href;
  return l;
});

const COMPOUNDS = [
  'Acetic acid', 'Acetone', 'Acetonitrile', 'Adenine', 'Aniline', 'Anthracene', 'Ascorbic acid',
  'Benzaldehyde', 'Benzene', 'Benzoic acid', 'Butanol', 'Caffeine', 'Chlorobenzene', 'Chloroform',
  'Cholesterol', 'Citric acid', 'Cyclohexane', 'Cyclohexanone', 'Dichloromethane', 'Diethyl ether',
  'Dimethylformamide', 'Dimethyl sulfoxide', 'Dopamine', 'Ethanol', 'Ethyl acetate',
  'Ethylene glycol', 'Formaldehyde', 'Formic acid', 'Glucose', 'Glycerol', 'Glycine', 'Heptane',
  'Hexane', 'Histamine', 'Ibuprofen', 'Indole', 'Isopropanol', 'Methanol', 'Morphine',
  'Naphthalene', 'Nicotine', 'Nitrobenzene', 'Octanol', 'Paracetamol', 'Phenol', 'Pyridine',
  'Pyrrole', 'Quinoline', 'Salicylic acid', 'Serotonin', 'Styrene', 'Testosterone',
  'Tetrahydrofuran', 'Toluene', 'Urea', 'Xylene',
];

const TARGETS = [
  'ABL1 kinase', 'ACE2', 'Adenosine A2A receptor', 'ALK', 'Androgen receptor', 'AURKA', 'BCL-2',
  'BRAF V600E', 'BRD4', 'BTK', 'Carbonic anhydrase II', 'CDK4', 'CDK9', 'COX-2', 'DHFR', 'DPP-4',
  'EGFR', 'ERK2', 'Estrogen receptor alpha', 'FGFR1', 'GSK-3 beta', 'HDAC6', 'HER2', 'HIV protease',
  'IDH1', 'JAK2', 'KRAS G12C', 'MEK1', 'mTOR', 'PARP1', 'PD-L1', 'PI3K alpha', 'PPAR gamma',
  'Proteasome 20S', 'ROCK1', 'SGLT2', 'SRC kinase', 'Thrombin', 'Trypsin', 'VEGFR2',
];

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

function fakeProvider(query, abort) {
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => {
      if (query.trim().toLowerCase() === 'err') {
        reject(new Error('Provider unavailable (HTTP 503)'));
        return;
      }
      const q = query.toLowerCase();
      resolve(TARGETS.filter((t) => t.toLowerCase().includes(q)));
    }, 300);
    abort.addEventListener('abort', () => {
      clearTimeout(timer);
      reject(new DOMException('Aborted', 'AbortError'));
    });
  });
}

function readout(combo) {
  const value = el('b');
  const state = el('code');
  bindText(combo.scope, value, computed(() => combo.value.value || '(empty)'));
  bindText(combo.scope, state, computed(() => `${combo.machineState.value}${combo.isOpen.value ? ', open' : ''}`));
  const row = el('p', 'u2-gallery-status');
  row.append('value = ', value, '  ·  machine = ', state);
  return row;
}

export async function render(main) {
  const {Combobox} = await import('../../src/components/combobox.js');

  main.append(el('h1', null, 'Combobox'));
  const intro = el('p');
  intro.innerHTML = 'Hand-written WAI-ARIA 1.2 combobox machine (see the header comment in ' +
    '<code>src/components/combobox.ts</code> for the Zag.js evaluation), popup via the shared ' +
    '<code>Overlay</code> layer, async data via <code>AsyncSource</code>. ' +
    'Keys: ↓ opens/moves, ↑ moves, Enter selects, Esc closes then clears, Tab closes and keeps the value.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  main.append(el('h2', null, 'Sync — 56 compounds, filtered as you type'));
  const sync = new Combobox({items: COMPOUNDS, placeholder: 'Compound…'});
  sync.root.style.width = '260px';
  main.append(sync.root, readout(sync));

  main.append(el('h2', null, 'Async — 300 ms provider, loading / empty / error+Retry'));
  const hint = el('p', 'u2-gallery-status',
    'Debounced 150 ms, each keystroke aborts the in-flight request. Type "err" to make the provider fail.');
  const remote = new Combobox({items: fakeProvider, placeholder: 'Target…'});
  remote.root.style.width = '260px';
  main.append(hint, remote.root, readout(remote));

  main.append(el('h2', null, 'Disposal'));
  const note = el('p', 'u2-gallery-status',
    'Dispose closes any open popup, drops the overlay listeners and both component scopes.');
  main.append(note, button('Dispose comboboxes', () => {
    sync.dispose();
    remote.dispose();
    refresh();
  }));
  refresh();
}

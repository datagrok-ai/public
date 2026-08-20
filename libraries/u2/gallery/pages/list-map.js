import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {ListInput} from '../../src/components/list-input.js';
import {MapInput} from '../../src/components/map-input.js';

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
injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-buttons-css', '../../css/buttons.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

export async function render(main) {
  main.append(el('h1', null, 'List & map inputs'));
  const intro = el('p');
  intro.innerHTML = 'The two collection editors ported from Dart. <code>ListInput</code> is a ' +
    'comma-separated field with the same quote-aware split — <code>a,"b,c",d</code> is three ' +
    'items — and the same expand toggle, which swaps the one-line field for a textarea without ' +
    'losing what is typed. <code>MapInput</code> edits key/value rows; two divergences from Dart: ' +
    'duplicate keys report through the standard validity signal instead of a silent ' +
    '<code>validate()</code> override, and a row whose key is empty never reaches the value.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const parts = [];
  const section = (title, builder) => {
    main.append(el('h2', null, title));
    const component = Control.build(builder);
    parts.push(component);
    main.append(component.root);
    return component;
  };

  const synonyms = signal(['ASA', '2-acetoxybenzoic acid']);
  section('List', () => {
    const input = new ListInput({label: 'Synonyms', bind: synonyms,
      tooltipText: 'Comma-separated; items with a comma are quoted'});
    return [input, readout('synonyms', computed(() => JSON.stringify(synonyms.value)))];
  });
  main.append(el('p', 'u2-gallery-status',
    'Type a quoted item — ASA,"2-acetoxybenzoic acid, USP" — and watch it stay one item. The ' +
    'expand icon on the right swaps in a textarea for longer lists.'));

  section('List: validation', () => {
    const input = new ListInput({label: 'Assays', value: ['IC50']});
    input.addValidator((v) => v.length < 2 ? 'Name at least two assays' : null);
    return [input, readout('validity', computed(() => input.validity.value ?? 'valid'))];
  });

  const params = signal({solvent: 'DMSO', temperature: '25'});
  section('Map', () => {
    const input = new MapInput({label: 'Params', bind: params});
    return [input, readout('params', computed(() => JSON.stringify(params.value)))];
  });
  main.append(el('p', 'u2-gallery-status',
    'The + and − on a row add below it and remove it; the last row is emptied, never removed. ' +
    'Give two rows the same key and both key cells turn invalid, and so does the input.'));

  main.append(el('h2', null, 'Programmatic value'));
  const buttons = el('p');
  buttons.append(
    button('Set synonyms', () => synonyms.value = ['ASA', 'aspirin', '2-acetoxybenzoic acid']),
    ' ',
    button('Clear synonyms', () => synonyms.value = []),
    ' ',
    button('Set params', () => params.value = {solvent: 'water', pH: '7.4'}),
    ' ',
    button('Clear params', () => params.value = {}));
  main.append(buttons);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Control.build(...) owner: disposing it releases the field listeners, the ' +
    'map\'s delegated row listeners and the value effects, and live scopes drop back to the page ' +
    'baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {RadioInput} from '../../src/components/inputs/radio-input.js';

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
injectOnce('u2-radio-css', '../../css/radio.css');

const STAGES = ['Discovery', 'Preclinical', 'Clinical'];

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

function shown(value) {
  return value === null ? 'null' : String(value);
}

export async function render(main) {
  main.append(el('h1', null, 'Radio input'));
  const intro = el('p');
  intro.innerHTML = 'Native <code>input type=radio</code> sharing one group name, so arrow keys, ' +
    'label clicks and screen readers work without machinery of our own. Use it for three to five ' +
    'options that should all stay visible; above that <code>ChoiceInput</code> reads better. ' +
    'The <code>buttons</code> variant is Dart\'s <code>buttonStyle</code> — the same group, ' +
    'rendered as a segmented row.';
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

  const stage = signal('Preclinical');
  section('Radios', () => {
    const input = new RadioInput({label: 'Stage', items: STAGES, bind: stage});
    return [input, readout('stage', computed(() => shown(stage.value)))];
  });

  section('Buttons variant, with per-item tooltips', () => {
    const input = new RadioInput({label: 'Stage', items: STAGES, buttons: true, bind: stage,
      itemTooltips: {
        Discovery: 'Hit finding and lead identification',
        Preclinical: 'Tox and PK before first-in-human',
        Clinical: 'Phase I onwards',
      }});
    return [input, readout('same signal', computed(() => shown(stage.value)))];
  });

  section('Validation', () => {
    const input = new RadioInput({label: 'Route', items: ['Oral', 'IV', 'Topical']});
    input.addValidator((v) => v === null ? 'Pick a route' : null);
    return [input, readout('validity', computed(() => input.validity.value ?? 'valid'))];
  });

  let series;
  section('Changing the items', () => {
    series = new RadioInput({label: 'Series', items: ['Series A', 'Series B', 'Series C'],
      value: 'Series B'});
    return [series, readout('series', computed(() => shown(series.value.value)))];
  });
  main.append(el('p', 'u2-gallery-status',
    'setItems keeps the value while its item survives, and clears it when the item is gone.'));
  const buttons = el('p');
  buttons.append(
    button('Keep Series B', () => series.setItems(['Series B', 'Series C'])),
    ' ',
    button('Drop Series B', () => series.setItems(['Series A', 'Series C'])),
    ' ',
    button('All three', () => series.setItems(['Series A', 'Series B', 'Series C'])));
  main.append(buttons);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Control.build(...) owner: disposing it releases the group listener and ' +
    'the value effect, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

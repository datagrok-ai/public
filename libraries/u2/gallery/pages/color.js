import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {ColorInput} from '../../src/components/color-input.js';

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
injectOnce('u2-color-css', '../../css/color.css');

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
  main.append(el('h1', null, 'Color input'));
  const intro = el('p');
  intro.innerHTML = 'A hex field plus a swatch that fronts the browser\'s own ' +
    '<code>input type=color</code> picker — the platform palette lives in Dart, and u2\'s L0-L2 ' +
    'stay platform-free. Text that is not a hex color (<code>#rgb</code> or <code>#rrggbb</code>) ' +
    'stays on screen and only marks the input invalid: the value signal and the swatch keep the ' +
    'last good color, exactly as the date editor does with unparseable dates. The layout is the ' +
    'platform\'s: a 100px hex field with the swatch on the options rail right of it, and ' +
    '<code>swatchOnly</code> drops both the field and the rail underline.';
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

  const fill = signal('#1f77b4');
  section('Hex field and picker', () => {
    const input = new ColorInput({label: 'Fill', bind: fill,
      tooltipText: 'Type a hex color, or click the swatch for the native picker'});
    return [input, readout('fill', computed(() => fill.value))];
  });

  section('Swatch only', () => {
    const input = new ColorInput({label: 'Marker', value: '#ff7f0e', swatchOnly: true});
    return [input, readout('marker', computed(() => input.value.value))];
  });

  section('Validation', () => {
    const input = new ColorInput({label: 'Accent', value: '#2ca02c'});
    return [input, readout('validity', computed(() => input.validity.value ?? 'valid'))];
  });

  section('Two inputs, one signal', () => [
    new ColorInput({label: 'Series color', bind: fill}),
    new ColorInput({label: 'Same color', bind: fill}),
  ]);

  main.append(el('h2', null, 'Programmatic value'));
  const buttons = el('p');
  buttons.append(
    button('Blue', () => fill.value = '#1f77b4'),
    ' ',
    button('Orange', () => fill.value = '#ff7f0e'),
    ' ',
    button('Short form (#0f0)', () => fill.value = '#0f0'));
  main.append(buttons);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Control.build(...) owner: disposing it releases the text, picker and ' +
    'swatch listeners, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

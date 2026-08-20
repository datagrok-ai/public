import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {SliderInput} from '../../src/components/inputs/slider-input.js';

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
injectOnce('u2-slider-css', '../../css/slider.css');

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
  main.append(el('h1', null, 'Slider input'));
  const intro = el('p');
  intro.innerHTML = 'A bare <code>input type=range</code>, the counterpart of Dart\'s ' +
    '<code>SliderInput</code>. The step defaults to <code>(max - min) / 100</code>; without ' +
    'bounds the track runs 0…100, as a bare range element does. A naked track shows no number — ' +
    'the platform\'s quirk, ported: the value surfaces only in the tooltip, on hover and while ' +
    'dragging. When the exact number matters as much as the range, reach for ' +
    '<code>NumberInput</code> with its slider option instead.';
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

  const opacity = signal(50);
  section('Bounded', () => {
    const input = new SliderInput({label: 'Opacity', min: 0, max: 100, bind: opacity,
      tooltipText: 'Step defaults to (100 - 0) / 100 = 1'});
    return [input, readout('opacity', computed(() => shown(opacity.value)))];
  });

  const threshold = signal(0.5);
  section('Formatted tooltip', () => {
    const input = new SliderInput({label: 'Threshold', min: 0, max: 1, step: 0.01, bind: threshold,
      format: (v) => `${Math.round(v * 100)}%`});
    return [input, readout('threshold', computed(() => shown(threshold.value)))];
  });

  section('Coarse step', () => {
    const input = new SliderInput({label: 'Bins', min: 0, max: 100, step: 10, value: 30});
    return [input, readout('bins', computed(() => shown(input.value.value)))];
  });

  const shared = signal(3);
  section('Two sliders, one signal', () => [
    new SliderInput({label: 'From', min: 0, max: 10, bind: shared}),
    new SliderInput({label: 'Also from', min: 0, max: 10, bind: shared}),
    readout('shared', computed(() => shown(shared.value))),
  ]);

  section('Inline / compact', () => {
    const zoom = new SliderInput({inline: true, min: 10, max: 400, step: 10, value: 100});
    const toolbar = divH([span('Zoom'), zoom]);
    toolbar.style.gap = 'var(--dg-space-m)';
    toolbar.style.alignItems = 'center';
    toolbar.style.padding = 'var(--dg-space-s)';
    toolbar.style.border = 'var(--dg-border-width) solid var(--dg-border)';
    toolbar.style.borderRadius = 'var(--dg-radius)';
    return [toolbar];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Control.build(...) owner: disposing it releases the tracks\' effects and ' +
    'listeners, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

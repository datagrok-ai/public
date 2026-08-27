import {signal, computed, Scope, Control, RangeSlider} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-tooltip-css', '../../css/tooltip.css');
injectOnce('u2-range-slider-css', '../../css/range-slider.css');

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
  main.append(el('h1', null, 'Range slider'));
  const intro = el('p');
  intro.innerHTML = 'Two-handle range selector, the counterpart of Dart\'s ' +
    '<code>RangeSlider</code>. Drag a handle to move one end, the band between them to move ' +
    'the whole window, click the track to pull the nearest handle. Each handle is a keyboard ' +
    'slider (arrows, PageUp/Down, Home/End); the value shows in the tooltip while dragging. ' +
    'Vertical runs bottom-to-top.';
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

  section('Basic', () => {
    const slider = new RangeSlider({min: 0, max: 100, lo: 20, hi: 60});
    return [slider, readout('range', computed(() =>
      `${slider.lo.value.toFixed(0)} … ${slider.hi.value.toFixed(0)}`))];
  });

  section('Step + minRange', () => {
    const slider = new RangeSlider({min: 0, max: 100, step: 5, minRange: 10, lo: 30, hi: 70});
    return [slider, readout('range', computed(() =>
      `${slider.lo.value} … ${slider.hi.value} (snaps to 5, at least 10 apart)`))];
  });

  section('Formatted', () => {
    const slider = new RangeSlider({min: 0, max: 1, lo: 0.2, hi: 0.8,
      format: (v) => `${Math.round(v * 100)}%`});
    return [slider, readout('range', computed(() =>
      `${Math.round(slider.lo.value * 100)}% … ${Math.round(slider.hi.value * 100)}%`))];
  });

  section('Vertical', () => {
    const slider = new RangeSlider({vertical: true, min: 0, max: 10, lo: 3, hi: 7});
    slider.root.style.height = '160px';
    const row = divH([slider, readout('range', computed(() =>
      `${slider.lo.value.toFixed(1)} … ${slider.hi.value.toFixed(1)}`))]);
    row.style.gap = 'var(--dg-space-xl)';
    row.style.alignItems = 'center';
    return [row];
  });

  section('Two sliders, shared signals', () => {
    const lo = signal(25);
    const hi = signal(75);
    const first = new RangeSlider({min: 0, max: 100, lo, hi});
    const second = new RangeSlider({min: 0, max: 100, lo, hi});
    return [first, second, readout('shared', computed(() => `${lo.value} … ${hi.value}`))];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Control.build(...) owner: disposing it releases the sliders\' effects ' +
    'and listeners, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

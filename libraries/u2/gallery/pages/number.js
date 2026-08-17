import {signal, computed, Scope, Component} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {NumberInput} from '../../src/components/number-input.js';

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
injectOnce('u2-number-css', '../../css/number.css');

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
  main.append(el('h1', null, 'Number input'));
  const intro = el('p');
  intro.innerHTML = 'A text editor, not <code>input type=number</code>: the native spinner and ' +
    'validation fight the token chrome. Text that does not parse stays on screen and marks the ' +
    'input invalid — the value signal keeps the last good number until it parses again. ' +
    'Min/max are applied on commit (blur, Enter, spinner), never mid-keystroke. ' +
    'Hover or focus an editor for the spinner; ArrowUp/ArrowDown step, Shift steps by ten, and ' +
    'the wheel steps only while the editor is focused.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const parts = [];
  const section = (title, builder) => {
    main.append(el('h2', null, title));
    const component = Component.build(builder);
    parts.push(component);
    main.append(component.root);
    return component;
  };

  const count = signal(10);
  section('Int', () => {
    const input = new NumberInput({label: 'Count', mode: 'int', min: 0, max: 100, bind: count,
      tooltipText: 'Whole numbers, clamped to 0…100 on commit'});
    return [input, readout('count', computed(() => shown(count.value)))];
  });

  const threshold = signal(0.5);
  section('Float', () => {
    const input = new NumberInput({label: 'Threshold', step: 0.1, bind: threshold,
      tooltipText: 'Decimals and scientific notation (1e-3)'});
    return [input, readout('threshold', computed(() => shown(threshold.value)))];
  });

  const shared = signal(42);
  section('Two editors, one signal', () => [
    new NumberInput({label: 'Dose', bind: shared, step: 5}),
    new NumberInput({label: 'Same dose', bind: shared, step: 5}),
    readout('dose', computed(() => shown(shared.value))),
  ]);

  section('Validation', () => {
    const even = new NumberInput({label: 'Even', mode: 'int', value: 2, step: 2});
    even.addValidator((v) => v !== null && v % 2 !== 0 ? 'Must be even' : null);
    return [even, readout('validity', computed(() => even.validity.value ?? 'valid'))];
  });

  section('Inline / compact', () => {
    const top = new NumberInput({inline: true, mode: 'int', value: 100, min: 1, max: 1000});
    const toolbar = divH([span('Top rows'), top]);
    toolbar.style.gap = 'var(--dg-space-m)';
    toolbar.style.alignItems = 'center';
    toolbar.style.padding = 'var(--dg-space-s)';
    toolbar.style.border = 'var(--dg-border-width) solid var(--dg-border)';
    toolbar.style.borderRadius = 'var(--dg-radius)';
    return [toolbar];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Component.build(...) owner: disposing it releases the editors\' effects, ' +
    'their key/wheel/spinner listeners and the readout bindings, and live scopes drop back to ' +
    'the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

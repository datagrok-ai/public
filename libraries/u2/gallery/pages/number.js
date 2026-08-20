import {signal, computed, Scope, Control} from '../../src/index.js';
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
    '<b>Min/max validate, they do not clamp</b> (platform parity): an out-of-range number you type ' +
    'stays visible, reaches the value signal, and shows the range message, so a form can block on ' +
    'it. Only the bounded controls stay inside the range — the slider inherently, the spinner and ' +
    'the clicker deliberately. <b>The chrome is the platform\'s, ported as it stands</b>: at rest ' +
    'a bare field, and hovering anywhere on the input (the label included) reveals the − + pair on ' +
    'the options rail, left of the units, and the slider under the field at its full width — ' +
    'overlapping the row below by 9px, exactly as <code>d4.css</code> positions it. u2\'s own ' +
    'hover spinner has no platform counterpart and is opt-in (<code>spinner: true</code>). ' +
    'ArrowUp/ArrowDown step, Shift steps by ten, and the wheel steps only while the editor is ' +
    'focused.';
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

  const count = signal(10);
  section('Int', () => {
    const input = new NumberInput({label: 'Count', mode: 'int', min: 0, max: 100, bind: count,
      tooltipText: 'Whole numbers; 0…100 is validated, not enforced'});
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

  const dose = signal(250);
  section('Slider + postfix (what a float property gets)', () => {
    const input = new NumberInput({label: 'Dose', min: 0, max: 1000, slider: true, clicker: true,
      postfix: 'mg', format: (v) => v.toFixed(1), bind: dose,
      tooltipText: 'Hover the row: the slider appears under the field, the − + pair left of the units'});
    return [input, readout('dose', computed(() => shown(dose.value)))];
  });

  const replicates = signal(3);
  section('Clicker (what a bounded int property gets)', () => {
    const input = new NumberInput({label: 'Replicates', mode: 'int', min: 1, max: 10, clicker: true,
      bind: replicates, tooltipText: 'showPlusMinus: − / + step by `step ?? 1` and stop at the bounds'});
    return [input, readout('replicates', computed(() => shown(replicates.value)))];
  });

  const spun = signal(7);
  section('Hover spinner (u2 only, opt-in)', () => {
    const input = new NumberInput({label: 'Rows', mode: 'int', min: 0, max: 100, spinner: true,
      bind: spun, tooltipText: 'No platform counterpart: off unless spinner: true'});
    return [input, readout('rows', computed(() => shown(spun.value)))];
  });

  const outOfRange = signal(5);
  section('Out of range', () => {
    const input = new NumberInput({label: 'Percent', min: 0, max: 100, bind: outOfRange,
      tooltipText: 'Type 1000: the text stays, the value is 1000, the message shows'});
    return [input, readout('percent', computed(() => shown(outOfRange.value))),
      readout('validity', computed(() => input.validity.value ?? 'valid'))];
  });

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
    'Each section is a Control.build(...) owner: disposing it releases the editors\' effects, ' +
    'their key/wheel/spinner listeners and the readout bindings, and live scopes drop back to ' +
    'the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

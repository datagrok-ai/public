import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {QNum} from '../../src/core/qnum.js';
import {QNumInput} from '../../src/components/inputs/qnum-input.js';
import {BigIntInput} from '../../src/components/inputs/bigint-input.js';

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
  main.append(el('h1', null, 'Big & qualified numbers'));
  const intro = el('p');
  intro.innerHTML = 'Two ports of platform inputs that a plain number editor cannot carry. ' +
    '<code>QNumInput</code> edits a <em>qualified</em> number — <code>&lt;5.2</code>, ' +
    '<code>10</code>, <code>&gt;1e3</code> — whose value is a single double with the qualifier ' +
    'packed into its two lowest mantissa bits, exactly as the platform stores a ' +
    '<code>qnum</code> column. <code>BigIntInput</code> edits an integer of any size, kept as ' +
    'text end to end so no digit is ever routed through a double. Text that does not parse stays ' +
    'on screen and marks the input invalid, as everywhere else in u2.';
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

  const ic50 = signal(QNum.less(5.2));
  section('Qualified number', () => {
    const input = new QNumInput({label: 'IC50', bind: ic50, postfix: 'µM',
      placeholder: '<, > or a plain number'});
    return [input, readout('packed', computed(() => {
      const v = ic50.value;
      return v === null ? 'null' : `${v} (q=${QNum.getQ(v)}, value=${QNum.getValue(v)})`;
    }))];
  });

  section('Qualified number, formatted', () => {
    const input = new QNumInput({label: 'IC50', value: QNum.less(5.2), format: (v) => v.toFixed(2)});
    return [input, readout('shown', computed(() => QNum.format(input.value.value, (v) => v.toFixed(2))))];
  });

  const molecules = signal(123456789012345678901234567890n);
  section('Integer of any size', () => {
    const input = new BigIntInput({label: 'Molecules', bind: molecules});
    return [input, readout('digits', computed(() => {
      const v = molecules.value;
      return v === null ? 'null' : `${String(v).length} digits, exact`;
    }))];
  });

  section('The same number through a double', () => {
    const text = String(molecules.peek());
    return [span(`Number('${text}') = ${Number(text)} — which is why the editor never routes ` +
      'digits through one.', 'u2-gallery-status')];
  });

  section('Invalid text', () => {
    const qnum = new QNumInput({label: 'IC50', value: QNum.create(1)});
    const big = new BigIntInput({label: 'Molecules', value: 7n});
    return [qnum, big, span('Type "<abc" or "1.5": the text stays, the input turns red, and the ' +
      'value keeps the last good number.', 'u2-gallery-status')];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Control.build(...) owner: disposing it releases the editors\' effects ' +
    'and listeners, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

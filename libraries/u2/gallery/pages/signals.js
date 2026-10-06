import {signal, computed, Scope, Control, bindText, bindValue} from '../../src/index.js';

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

export async function render(main) {
  main.append(el('h1', null, 'Signals & lifecycle'));
  const intro = el('p');
  intro.innerHTML = 'L1 reactive core: vendored signals engine behind <code>u2.signal</code>, ' +
    'effects owned by <code>Scope</code>/<code>Control</code>, disposal releases everything. ' +
    'In Datagrok the same disposal rides <code>DG.Widget.kill</code> via <code>u2/dg host()</code>.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes (leak detector ground truth): ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  main.append(el('h2', null, 'Derived state'));
  const count = signal(0);
  const doubled = computed(() => count.value * 2);
  const counter = new Control();
  const countText = el('b');
  const doubledText = el('b');
  bindText(counter.scope, countText, count);
  bindText(counter.scope, doubledText, doubled);
  const row1 = el('div');
  row1.append(button('+1', () => { count.value++; refresh(); }), ' count = ', countText,
    ', computed ×2 = ', doubledText);
  counter.root.append(row1);
  main.append(counter.root);

  main.append(el('h2', null, 'Two-way binding, echo-suppressed'));
  const name = signal('world');
  const form = new Control();
  const a = el('input');
  const b = el('input');
  bindValue(form.scope, a, name);
  bindValue(form.scope, b, name);
  const mirror = el('span');
  bindText(form.scope, mirror, computed(() => `Hello, ${name.value}!`));
  const row2 = el('div');
  row2.append(a, ' ⇄ ', b, '  →  ', mirror);
  form.root.append(row2);
  main.append(form.root);

  main.append(el('h2', null, 'Disposal'));
  const note = el('p', 'u2-gallery-status',
    'Dispose both components above, then press +1 or type: nothing updates, live scopes drop.');
  const dispose = button('Dispose components', () => {
    counter.dispose();
    form.dispose();
    refresh();
  });
  main.append(note, dispose);
}

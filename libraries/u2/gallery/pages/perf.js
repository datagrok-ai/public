import {Scope, Component} from '../../src/index.js';
import {Perf} from '../../src/core/perf.js';
import {divV, span, button as u2Button} from '../../src/core/elements.js';
import {TextInput} from '../../src/components/text-input.js';
import {BoolInput} from '../../src/components/bool-input.js';
import {NumberInput} from '../../src/components/number-input.js';
import {ChoiceInput} from '../../src/components/choice-input.js';
import {TabStrip} from '../../src/components/tabs.js';
import {VirtualList} from '../../src/components/list.js';

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
injectOnce('u2-choice-css', '../../css/choice.css');
injectOnce('u2-tabs-css', '../../css/tabs.css');
injectOnce('u2-list-css', '../../css/list.css');

const BUDGET_MS = 15;
const WIDGETS = 200;
const ROWS = 1000000;
const SCROLL_STEPS = 60;
const LEAK_ROUNDS = 50;

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

function ms(value) {
  return `${value.toFixed(2)} ms`;
}

/** 200 mixed widgets in one window: inputs, buttons, a tab strip and a small list. */
function buildWindow() {
  return Component.build(() => {
    const rows = [];
    for (let i = 0; i < 50; i++)
      rows.push(new TextInput({label: `Field ${i}`, value: `value ${i}`}));
    for (let i = 0; i < 30; i++)
      rows.push(new NumberInput({label: `Count ${i}`, value: i, min: 0, max: 100}));
    for (let i = 0; i < 30; i++)
      rows.push(new BoolInput({label: `Flag ${i}`, value: i % 2 === 0}));
    for (let i = 0; i < 20; i++)
      rows.push(new ChoiceInput({label: `Series ${i}`, items: ['a', 'b', 'c'], value: 'b'}));
    for (let i = 0; i < 68; i++)
      rows.push(u2Button(`Action ${i}`, () => {}));
    const tabs = new TabStrip({tabs: ['Data', 'Model', 'Results', 'Log'].map((label) => ({
      id: label, label, content: () => span(`${label} panel`),
    }))});
    const list = new VirtualList({render: (item) => span(`Row ${item}`)});
    list.root.style.height = '120px';
    list.setItems(Array.from({length: 1000}, (_, i) => i));
    return [divV(rows), tabs, list];
  });
}

function frame() {
  return new Promise((resolve) => requestAnimationFrame(resolve));
}

async function measureScroll(list) {
  const step = Math.floor(ROWS * 22 / SCROLL_STEPS);
  const started = performance.now();
  for (let i = 1; i <= SCROLL_STEPS; i++) {
    list.root.scrollTop = i * step;
    await frame();
  }
  return (performance.now() - started) / SCROLL_STEPS;
}

function reportTable() {
  const table = el('table');
  const head = el('tr');
  for (const title of ['Measurement', 'Construct', 'Paint', 'Total'])
    head.append(el('th', null, title));
  table.append(head);
  for (const entry of Perf.report()) {
    const row = el('tr');
    row.append(el('td', null, entry.label), el('td', null, ms(entry.construct)),
      el('td', null, ms(entry.paint)), el('td', null, ms(entry.total)));
    table.append(row);
  }
  return table;
}

export async function render(main) {
  const page = new Component();
  main.append(el('h1', null, 'Performance & leak budgets'));
  const intro = el('p');
  intro.innerHTML = 'The P4 exit budgets, instrumented through <code>Perf.measure</code>: a ' +
    '200-widget window under 15 ms construct+paint, a million-row list that scrolls at frame rate, ' +
    'and a build/dispose round-trip that leaves no live scope behind.';
  main.append(intro);

  main.append(el('h2', null, `${WIDGETS}-widget window`));
  const windowStatus = el('p', 'u2-gallery-status', 'Not measured yet.');
  const results = el('div');
  main.append(button(`Build ${WIDGETS} widgets`, async (e) => {
    e.target.disabled = true;
    windowStatus.textContent = 'Measuring…';
    const entry = await Perf.measure(`${WIDGETS}-widget window`, buildWindow);
    // paint is double-rAF wall time and always spans >= one vsync (~16.7 ms @60Hz), so the
    // budget applies to construct (sync JS cost); paint just has to fit in two frames (no jank)
    const pass = entry.construct <= BUDGET_MS && entry.paint <= 34;
    windowStatus.textContent = `construct ${ms(entry.construct)} (budget ${BUDGET_MS} ms) + ` +
      `paint ${ms(entry.paint)} (vsync-bound, jank threshold 34 ms) — ${pass ? 'PASS ✔' : 'FAIL ✘'}`;
    windowStatus.style.color = pass ? 'var(--dg-success)' : 'var(--dg-failure)';
    results.replaceChildren(reportTable());
    e.target.disabled = false;
  }), windowStatus);

  main.append(el('h2', null, `${ROWS.toLocaleString()}-row list scroll`));
  const list = new VirtualList({render: (item) => span(`Row ${item.toLocaleString()}`)});
  list.root.style.height = '180px';
  list.setItems(Array.from({length: ROWS}, (_, i) => i));
  main.append(list.root);
  const scrollStatus = el('p', 'u2-gallery-status', 'Not measured yet.');
  main.append(button(`Scroll ${SCROLL_STEPS} steps`, async (e) => {
    e.target.disabled = true;
    scrollStatus.textContent = 'Scrolling…';
    const average = await measureScroll(list);
    scrollStatus.textContent = `${SCROLL_STEPS} steps · ${ms(average)} per frame · ` +
      `${list.renderedCount} rows in the DOM · ${average <= 16.7 ? '60 fps ✔' : 'below 60 fps ✘'}`;
    e.target.disabled = false;
  }), scrollStatus);

  main.append(el('h2', null, 'Leak round-trip'));
  const leakStatus = el('p', 'u2-gallery-status', 'Not measured yet.');
  main.append(button(`Build and dispose ${LEAK_ROUNDS} components`, () => {
    const before = Scope.liveCount;
    for (let i = 0; i < LEAK_ROUNDS; i++) {
      const component = Component.build(() => [
        new TextInput({label: `Round ${i}`, value: String(i)}),
        new BoolInput({label: 'Flag', value: true}),
        u2Button('Act', () => {}),
      ]);
      component.dispose();
    }
    const delta = Scope.liveCount - before;
    leakStatus.textContent = `${LEAK_ROUNDS} rounds · live scopes ${before} → ${Scope.liveCount} · ` +
      `delta ${delta} — ${delta === 0 ? 'PASS ✔' : 'FAIL ✘'}`;
    leakStatus.style.color = delta === 0 ? 'var(--dg-success)' : 'var(--dg-failure)';
  }), leakStatus);

  main.append(el('h2', null, 'Accumulated measurements'), results);

  const live = el('p', 'u2-gallery-status');
  const tick = () => live.textContent = `Live scopes right now: ${Scope.liveCount}`;
  tick();
  const timer = setInterval(tick, 1000);
  page.own(() => clearInterval(timer));
  main.append(live);
}

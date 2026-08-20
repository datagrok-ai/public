import {signal, computed, Scope, Control, bindText} from '../../src/index.js';
import {Accordion} from '../../src/components/accordion.js';

function injectOnce(id, make) {
  if (document.getElementById(id)) return;
  const e = make();
  e.id = id;
  document.head.append(e);
}

injectOnce('u2-accordion-css', () => {
  const l = document.createElement('link');
  l.rel = 'stylesheet';
  l.href = new URL('../../css/accordion.css', import.meta.url).href;
  return l;
});

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

function paragraph(text) {
  const p = el('div');
  p.style.padding = '2px 0';
  p.textContent = text;
  return p;
}

const TITLES = ['Overview', 'Lazy details', 'Live values', 'Removable'];

export async function render(main) {
  main.append(el('h1', null, 'Accordion'));
  const intro = el('p');
  intro.innerHTML = 'Panes expand independently (no exclusive mode, like the platform accordion). ' +
    'Headers are a roving-tabindex group: ↑/↓/Home/End move between them, Enter or Space toggles. ' +
    'Function content is built on first expansion and kept afterwards — collapsing sets ' +
    '<code>display:none</code>, it never detaches. Each pane builds under its own scope, so ' +
    '<code>removePane</code> disposes whatever the builder created.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const log = el('div', 'u2-gallery-status');
  const logLine = (text) => log.append(el('div', null, text));

  const ticks = signal(0);

  const eager = el('div');
  eager.append(
    paragraph('Built eagerly and passed in as an element — present in the DOM from the start.'),
    paragraph('Type below, collapse the pane, expand it again — the value survives.')
  );
  const memo = el('input');
  memo.placeholder = 'Scratch…';
  eager.append(memo);

  const lazy = () => {
    const started = performance.now();
    const box = el('div');
    box.append(paragraph('Built lazily on first expansion — never rebuilt afterwards.'));
    for (let i = 0; i < 200; i++)
      box.append(paragraph(`row ${i}: ${'█'.repeat(1 + (i % 10))}`));
    logLine(`"Lazy details" content built in ${(performance.now() - started).toFixed(1)} ms ` +
      `at ${new Date().toLocaleTimeString()}`);
    return box;
  };

  const live = () => {
    const content = Control.build(() => {
      const value = el('b');
      bindText(Scope.ambient, value, computed(() => String(ticks.value)));
      const line = el('div');
      line.append('ticks = ', value);
      return [
        paragraph('A Control built inside the builder — adopted by the pane scope via Scope.ambient.'),
        line,
        button('Tick', () => ticks.value++),
      ];
    });
    logLine(`"Live values" content built (Scope.liveCount is now ${Scope.liveCount})`);
    return content.root;
  };

  const acc = new Accordion();
  acc.root.style.maxWidth = '520px';
  acc.root.style.border = 'var(--dg-border-width) solid var(--dg-border)';
  acc.addPane('Overview', eager, true);
  acc.addPane('Lazy details', lazy);
  acc.addPane('Live values', live);
  acc.addPane('Removable', paragraph('Remove this pane below and watch the live scope count drop.'));
  main.append(acc.root);

  const state = el('code');
  bindText(acc.scope, state, computed(() => TITLES
    .map((t) => `${t}: ${acc.getPane(t) ? acc.getPane(t).expanded.value : '(gone)'}`)
    .join('   ')));
  const readout = el('p', 'u2-gallery-status');
  readout.append(state);
  main.append(readout);

  main.append(el('h2', null, 'Programmatic expand / collapse'));
  const hint = el('p', 'u2-gallery-status',
    'The pane signal is two-way: writing it moves the chevron and the content, exactly as a click does.');
  const row = el('div');
  for (const title of TITLES) {
    row.append(button(`Toggle "${title}"`, () => {
      const pane = acc.getPane(title);
      if (pane)
        pane.expanded.value = !pane.expanded.value;
      else
        logLine(`getPane("${title}") → undefined`);
    }), ' ');
  }
  const all = el('div');
  all.style.paddingTop = 'var(--dg-space-m)';
  const setAll = (on) => {
    for (const title of TITLES) {
      const pane = acc.getPane(title);
      if (pane)
        pane.expanded.value = on;
    }
  };
  all.append(button('Expand all', () => setAll(true)), ' ', button('Collapse all', () => setAll(false)));
  main.append(hint, row, all);

  main.append(el('h2', null, 'removePane'));
  const removeNote = el('p', 'u2-gallery-status',
    'Removing a pane disposes its scope — the effect, and any component the lazy builder created, go with it.');
  main.append(removeNote, button('removePane("Live values")', () => {
    acc.removePane('Live values');
    logLine(`removePane("Live values") → getPane is ${acc.getPane('Live values')}`);
    refresh();
  }), ' ', button('removePane("Removable")', () => {
    acc.removePane('Removable');
    refresh();
  }));

  main.append(el('h2', null, 'Build log'));
  main.append(log);

  main.append(el('h2', null, 'Disposal'));
  const note = el('p', 'u2-gallery-status',
    'Dispose drops the header listeners, every pane effect and every pane scope — live scopes drop back.');
  main.append(note, button('Dispose accordion', () => {
    acc.dispose();
    refresh();
  }));
  refresh();
}

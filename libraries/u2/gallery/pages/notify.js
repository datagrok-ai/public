import {notify} from '../../src/index.js';
import {divH, button} from '../../src/core/elements.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-buttons-css', '../../css/buttons.css');
injectOnce('u2-notify-css', '../../css/notify.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

export async function render(main) {
  main.append(el('h1', null, 'Notifications'));
  const intro = el('p');
  intro.innerHTML = 'The u2 counterpart of the platform\'s <code>Balloon</code>: ' +
    '<code>notify.info / warning / error</code> stack in the top-right corner. Info balloons ' +
    'auto-hide after 5 s — hovering pauses the countdown; warnings and errors stay until ' +
    'closed. Hover a balloon for the copy and close icons; a click anywhere dismisses it. ' +
    'String content always renders as text, never as markup.';
  main.append(intro);

  main.append(el('h2', null, 'Types'));
  const row = divH([
    button('Info', () => notify.info('Project saved.')),
    button('Warning', () => notify.warning('The query returned 10,000 rows; only the first 1,000 are shown.')),
    button('Error', () => notify.error('Connection "chembl" refused: timeout after 30 s.',
      {copyText: 'stack: QueryRunner.run (query_runner.dart:118)'})),
    button('Close all', () => notify.closeAll()),
  ]);
  row.style.gap = 'var(--dg-space-m)';
  main.append(row);

  main.append(el('h2', null, 'Keys'));
  const keys = divH([
    button('singleKey', () => {
      if (notify.info('Only one of me at a time.', {singleKey: 'demo', autoHide: false}) === null)
        console.log('suppressed: a "demo" balloon is already open');
    }),
    button('oneTimeKey', () => {
      if (notify.info('You will only ever see me once (per browser).', {oneTimeKey: 'gallery-demo'}) === null)
        console.log('suppressed: already shown once — clear localStorage["u2-notify-gallery-demo"] to reset');
    }),
  ]);
  keys.style.gap = 'var(--dg-space-m)';
  main.append(keys);
  main.append(el('p', 'u2-gallery-status',
    'A suppressed show returns null and logs to the console here.'));

  main.append(el('h2', null, 'Element content'));
  main.append(button('Rich balloon', () => {
    const content = document.createElement('div');
    const strong = document.createElement('b');
    strong.textContent = '3 tests failed';
    content.append(strong, document.createElement('br'), 'chem: substructure-search, to-mol, render');
    notify.warning(content);
  }));
}

import {computed, Scope, bindText} from '../../src/index.js';
import {TypeAhead} from '../../src/components/typeahead.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-typeahead-css', '../../css/typeahead.css');

// platform-free stand-ins for DG.User — the same shape src/dg/user-input.ts renders
const FIRST = ['Ada', 'Bruno', 'Chen', 'Dmitri', 'Elena', 'Farid', 'Grace', 'Hugo', 'Iris', 'Jonas',
  'Klara', 'Liam', 'Maya', 'Nikolai', 'Olga', 'Pavel', 'Quinn', 'Rosa', 'Sofia', 'Tomas',
  'Ulrich', 'Vera', 'Wanda', 'Xenia', 'Yusuf', 'Zara', 'Anton', 'Bianca', 'Cyril', 'Dana'];
const LAST = ['Almeida', 'Bauer', 'Costa', 'Dubois', 'Egan', 'Farrell', 'Gruber', 'Haas', 'Ivanov',
  'Jansen', 'Kowalski', 'Lindqvist', 'Moreau', 'Novak', 'Olsen', 'Petrov', 'Quintana', 'Rossi',
  'Sorensen', 'Tanaka', 'Ueda', 'Vargas', 'Weber', 'Xu', 'Yilmaz', 'Zeller', 'Adler', 'Blum',
  'Conti', 'Duarte'];

const USERS = FIRST.map((first, i) => {
  const last = LAST[i];
  const login = `${first}${last}`.toLowerCase();
  return {name: `${first} ${last}`, login, email: `${login}@datagrok.ai`};
});

function hue(login) {
  let hash = 0;
  for (let i = 0; i < login.length; i++)
    hash = (hash * 31 + login.charCodeAt(i)) % 9973;
  return hash % 6 + 1;
}

function initials(name) {
  const capitals = name.match(/[A-Z]/g);
  return capitals ? capitals.slice(0, 2).join('') : name.substring(0, 2).toUpperCase();
}

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

function renderUser(user) {
  const row = el('div', 'u2-typeahead-user');
  const text = el('div', 'u2-typeahead-user-text');
  text.append(el('div', 'u2-typeahead-user-name', user.name),
    el('div', 'u2-typeahead-user-secondary', user.email));
  row.append(el('div', `u2-avatar u2-avatar-${hue(user.login)}`, initials(user.name)), text);
  return row;
}

function fakeProvider(query, abort) {
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => {
      if (query.trim().toLowerCase() === 'err') {
        reject(new Error('Directory unavailable (HTTP 503)'));
        return;
      }
      const q = query.toLowerCase();
      resolve(USERS.filter((u) => u.name.toLowerCase().includes(q) || u.login.includes(q)));
    }, 300);
    abort.addEventListener('abort', () => {
      clearTimeout(timer);
      reject(new DOMException('Aborted', 'AbortError'));
    });
  });
}

function readout(typeahead) {
  const selected = el('code');
  const state = el('code');
  bindText(typeahead.scope, selected, computed(() => JSON.stringify(typeahead.selected.value)));
  bindText(typeahead.scope, state, computed(() =>
    `${typeahead.machineState.value}${typeahead.isOpen.value ? ', open' : ''}, text "${typeahead.text.value}"`));
  const row = el('p', 'u2-gallery-status');
  row.append('selected = ', selected, '  ·  machine = ', state);
  return row;
}

export async function render(main) {
  main.append(el('h1', null, 'Type-ahead'));
  const intro = el('p');
  intro.innerHTML = 'Generic <code>TypeAhead&lt;T&gt;</code>: the combobox machine over items of ' +
    'any type, with a caller-supplied row renderer. The machine owns the row (hover, active, ' +
    'ARIA); <code>render(item)</code> only supplies its content. ' +
    'Keys: ↓ opens/moves, ↑ moves, Enter picks, Esc closes then clears, Tab closes and keeps the pick. ' +
    'Typing after a pick clears <code>selected</code> — a stale object never survives an edit.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  main.append(el('h2', null, 'Sync — 30 users, complex rows'));
  const sync = new TypeAhead({
    source: USERS,
    itemText: (u) => u.name,
    render: renderUser,
    placeholder: 'Assignee…',
  });
  sync.root.style.width = '280px';
  main.append(sync.root, readout(sync));

  main.append(el('h2', null, 'Async — 300 ms directory, loading / empty / error+Retry'));
  main.append(el('p', 'u2-gallery-status',
    'Debounced 150 ms, each keystroke aborts the in-flight request. Type "err" to make it fail.'));
  const remote = new TypeAhead({
    source: fakeProvider,
    itemText: (u) => u.name,
    render: renderUser,
    placeholder: 'Search the directory…',
  });
  remote.root.style.width = '280px';
  main.append(remote.root, readout(remote));

  main.append(el('h2', null, 'minChars: 2 — the popup shows a hint until it is worth searching'));
  const gated = new TypeAhead({
    source: fakeProvider,
    itemText: (u) => u.name,
    render: renderUser,
    minChars: 2,
    debounceMs: 300,
    placeholder: 'At least two characters…',
  });
  gated.root.style.width = '280px';
  main.append(gated.root, readout(gated));

  main.append(el('h2', null, 'Plain strings — the default renderer'));
  const logins = new TypeAhead({
    source: USERS.map((u) => u.login),
    itemText: (login) => login,
    placeholder: 'Login…',
  });
  logins.root.style.width = '280px';
  main.append(logins.root, readout(logins));

  main.append(el('h2', null, 'Programmatic selection'));
  main.append(el('p', 'u2-gallery-status',
    'Setting selected updates the input text; clearing it empties the field.'));
  const buttons = el('p');
  buttons.append(
    button('Select Grace Gruber', () => sync.selected.value = USERS[6]),
    ' ',
    button('Clear', () => sync.selected.value = null));
  main.append(buttons);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Dispose closes any open popup and drops every option scope, listener and effect.'));
  main.append(button('Dispose type-aheads', () => {
    for (const t of [sync, remote, gated, logins])
      t.dispose();
    refresh();
  }));
  refresh();
}

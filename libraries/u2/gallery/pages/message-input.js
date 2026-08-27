import {computed, signal, Scope, bindText} from '../../src/index.js';
import {MessageInput} from '../../src/components/inputs/message-input.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-buttons-css', '../../css/buttons.css');
injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-tooltip-css', '../../css/tooltip.css');
injectOnce('u2-typeahead-css', '../../css/typeahead.css');
injectOnce('u2-message-input-css', '../../css/message-input.css');

// platform-free stand-ins for DG.User — the ids look like the real thing so the composed token is
// byte-identical in shape to what src/dg/inputs/message-input.ts emits
const PEOPLE = [
  {id: '3f2a1c40-1111-4a01-9c11-0a1b2c3d4e01', name: 'Ada Almeida', login: 'adaalmeida'},
  {id: '3f2a1c40-1111-4a01-9c11-0a1b2c3d4e02', name: 'Bruno Bauer', login: 'brunobauer'},
  {id: '3f2a1c40-1111-4a01-9c11-0a1b2c3d4e03', name: 'Chen Costa', login: 'chencosta'},
  {id: '3f2a1c40-1111-4a01-9c11-0a1b2c3d4e04', name: 'Dmitri Dubois', login: 'dmitridubois'},
  {id: '3f2a1c40-1111-4a01-9c11-0a1b2c3d4e05', name: 'Elena Egan', login: 'elenaegan'},
  {id: '3f2a1c40-1111-4a01-9c11-0a1b2c3d4e06', name: 'Farid Farrell', login: 'faridfarrell'},
];

const TOKEN = '<span>#\\{x\\.([^."]+)\\."(.*?)"\\}<\\/span>';

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

function personRow(person) {
  const row = el('div', 'u2-typeahead-user');
  const text = el('div', 'u2-typeahead-user-text');
  text.append(el('div', 'u2-typeahead-user-name', person.name),
    el('div', 'u2-typeahead-user-secondary', `${person.login}@datagrok.ai`));
  row.append(text);
  return row;
}

function peopleProvider() {
  return {
    trigger: '@',
    minChars: 0,
    debounceMs: 0,
    fetch: (query) => {
      const q = query.trim().toLowerCase();
      return Promise.resolve(q ? PEOPLE.filter((p) => p.name.toLowerCase().includes(q) ||
        p.login.includes(q)) : PEOPLE);
    },
    renderItem: personRow,
    caption: (p) => p.name,
    toMarkup: (p) => `<span>#{x.${p.id}."${p.name}"}</span>`,
    tokenPattern: TOKEN,
    renderToken: (m) => el('span', null, m[2]),
  };
}

function valueReadout(compose) {
  const code = el('code');
  const value = signal(compose.value);
  const update = () => value.value = compose.value;
  compose.root.addEventListener('input', update);
  compose.root.addEventListener('click', update);
  compose.root.addEventListener('keyup', update);
  bindText(compose.scope, code, computed(() => value.value === '' ? '(empty)' : value.value));
  const row = el('p', 'u2-gallery-status');
  row.append('value = ', code);
  return row;
}

export async function render(main) {
  main.append(el('h1', null, 'Message input'));
  const intro = el('p');
  intro.innerHTML = 'A contenteditable compose box whose value is <b>Datagrok markup text</b>. ' +
    'Type <code>@</code> (or press the @ button) to open a caret-anchored people popup; picking a ' +
    'row inserts an <b>atomic chip</b> that serializes back to the exact token ' +
    '<code>&lt;span&gt;#{x.&lt;id&gt;."&lt;name&gt;"}&lt;/span&gt;</code>. ' +
    'Backspace next to a chip deletes it whole. ' +
    'Keys in the popup: ↓/↑ move, Enter picks, Esc closes and gives the typed text back. ' +
    'The first box carries a <code>draftKey</code>, so whatever is left unsent in it — chips ' +
    'included — comes back after a reload. ' +
    'The provider here is a static fixture — the platform one lives in ' +
    '<code>src/dg/inputs/message-input.ts</code>.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  main.append(el('h2', null, 'Compose — Ctrl+Enter sends, and the unsent text outlives a reload'));
  const sentLog = el('code', null, '(nothing sent yet)');
  const compose = new MessageInput({
    placeholder: 'Message… (@ to mention)',
    mentions: [peopleProvider()],
    sendOn: 'ctrlEnter',
    draftKey: 'e2e',
    onSend: (markup) => sentLog.textContent = markup,
  });
  compose.root.style.width = '420px';
  main.append(compose.root, valueReadout(compose));
  const sentRow = el('p', 'u2-gallery-status');
  sentRow.append('sent = ', sentLog);
  main.append(sentRow);

  main.append(el('h2', null, 'Restore a draft — tokens come back as chips'));
  main.append(el('p', 'u2-gallery-status',
    'Assigning markup with tokens in it re-creates the chips; reading the value back is byte-exact.'));
  const restore = el('p');
  restore.append(
    button('Restore a two-mention draft', () => {
      compose.value = 'Hi ' +
        `<span>#{x.${PEOPLE[0].id}."${PEOPLE[0].name}"}</span> and ` +
        `<span>#{x.${PEOPLE[2].id}."${PEOPLE[2].name}"}</span>, please review this.`;
      compose.root.dispatchEvent(new Event('input', {bubbles: true}));
    }),
    ' ',
    button('Clear', () => {
      compose.clear();
      compose.root.dispatchEvent(new Event('input', {bubbles: true}));
    }));
  main.append(restore);

  main.append(el('h2', null, 'Enter sends, Shift+Enter breaks the line'));
  const quickLog = el('code', null, '(nothing sent yet)');
  const quick = new MessageInput({
    placeholder: 'Enter sends…',
    mentions: [peopleProvider()],
    sendOn: 'enter',
    onSend: (markup) => quickLog.textContent = markup,
  });
  quick.root.style.width = '420px';
  const quickRow = el('p', 'u2-gallery-status');
  quickRow.append('sent = ', quickLog);
  main.append(quick.root, quickRow);

  main.append(el('h2', null, 'Auto-grow — maxHeight 120px'));
  const tall = new MessageInput({
    placeholder: 'Keep typing — the box grows to 120px, then scrolls',
    mentions: [peopleProvider()],
    maxHeight: 120,
  });
  tall.root.style.width = '420px';
  main.append(tall.root);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Dispose closes any open mention popup and drops every session scope, listener and effect.'));
  main.append(button('Dispose message inputs', () => {
    for (const c of [compose, quick, tall])
      c.dispose();
    refresh();
  }));
  refresh();
}

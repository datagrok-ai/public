import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {TagsInput} from '../../src/components/inputs/tags-input.js';

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
injectOnce('u2-badge-css', '../../css/badge.css');
injectOnce('u2-tags-css', '../../css/tags.css');

const KEYWORDS = ['acid', 'base', 'ester', 'amine', 'ketone', 'alcohol', 'aromatic', 'chiral'];
const USERS = [
  {login: 'adaalmeida', name: 'Ada Almeida', team: 'Chemistry'},
  {login: 'brunobauer', name: 'Bruno Bauer', team: 'Biology'},
  {login: 'chencosta', name: 'Chen Costa', team: 'Chemistry'},
  {login: 'dianadias', name: 'Diana Dias', team: 'Data science'},
  {login: 'erikevans', name: 'Erik Evans', team: 'Biology'},
];

/** The ObjectRenderer protocol: caption for text, markup for chips, listItem for popup rows. */
const USER_RENDERER = {
  caption: (u) => u.name,
  markup: (u) => {
    const el = document.createElement('span');
    el.textContent = u.name;
    return el;
  },
  listItem: (u) => {
    const row = document.createElement('div');
    const name = document.createElement('div');
    name.textContent = u.name;
    const team = document.createElement('div');
    team.className = 'u2-gallery-status';
    team.textContent = `${u.login} · ${u.team}`;
    row.append(name, team);
    return row;
  },
};

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

/** A source that answers after a beat, so the loading and error rows are visible. */
function userSource(failFirst) {
  let failed = false;
  return async (query) => {
    await new Promise((resolve) => setTimeout(resolve, 400));
    if (failFirst && !failed) {
      failed = true;
      throw new Error('The directory did not answer');
    }
    const q = query.toLowerCase();
    return USERS.filter((u) => u.name.toLowerCase().includes(q) || u.login.includes(q));
  };
}

export async function render(main) {
  main.append(el('h1', null, 'Tags input'));
  const intro = el('p');
  intro.innerHTML = 'Chips over items of any type — the port of the platform\'s ' +
    '<code>TagsInput</code>. Suggestions come from a static list or an async source (the same ' +
    '<code>AsyncSource</code> every u2 async control uses, so stale answers are dropped), chips ' +
    'and rows are drawn by an <code>ObjectRenderer</code>, and <code>allowNew</code> turns typed ' +
    'text into an item. For a plain string list with a checkbox popup, use ' +
    '<code>MultiSelect</code> instead.';
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

  const keywords = signal(['acid']);
  section('Strings, with suggestions', () => {
    const input = new TagsInput({label: 'Keywords', items: KEYWORDS, bind: keywords,
      placeholder: 'Type to search…'});
    return [input, readout('keywords', computed(() => JSON.stringify(keywords.value)))];
  });
  main.append(el('p', 'u2-gallery-status',
    'Arrow keys and Home/End browse the popup, Enter accepts, Esc closes it, and Backspace on an ' +
    'empty box drops the last chip. Items already picked leave the suggestion list.'));

  section('Free entry (allowNew)', () => {
    const input = new TagsInput({label: 'Keywords', items: KEYWORDS, allowNew: true,
      placeholder: 'Type and press Enter'});
    return [input, readout('value', computed(() => JSON.stringify(input.value.value)))];
  });

  section('Objects, rendered', () => {
    const input = new TagsInput({label: 'Reviewers', items: USERS, renderer: USER_RENDERER,
      value: [USERS[0]], placeholder: 'Type a name…'});
    return [input, readout('logins',
      computed(() => input.value.value.map((u) => u.login).join(', ') || '—'))];
  });

  section('Async source', () => {
    const input = new TagsInput({label: 'Reviewers', source: userSource(false),
      renderer: USER_RENDERER, minChars: 2, placeholder: 'At least two characters…'});
    return [input, readout('picked', computed(() => String(input.value.value.length)))];
  });

  section('Async source that fails once', () => {
    const input = new TagsInput({label: 'Reviewers', source: userSource(true),
      renderer: USER_RENDERER, placeholder: 'Type anything…'});
    return [input, el('p', 'u2-gallery-status',
      'The first query fails into the error row; Retry runs the same query again.')];
  });

  section('Protected chips (deleteEnabled)', () => {
    const input = new TagsInput({label: 'Roles', items: ['Author', 'Reviewer', 'Owner'],
      value: ['Owner', 'Reviewer'], deleteEnabled: (item) => item !== 'Owner'});
    return [input, el('p', 'u2-gallery-status',
      'Owner keeps its chip but loses the ✕, and Backspace will not take it either.')];
  });

  section('Validation', () => {
    const input = new TagsInput({label: 'Keywords', items: KEYWORDS, allowNew: true});
    input.addValidator((v) => v.length < 2 ? 'Pick at least two' : null);
    return [input, readout('validity', computed(() => input.validity.value ?? 'valid'))];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Control.build(...) owner: disposing it closes the popup and releases the ' +
    'per-render chip and row scopes, so live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}

/* ObjectRenderer as a TypeAhead option: caption supplies the item text, listItem the rows,
   and explicit callbacks win over the renderer's members. The dg handlerRenderer implements
   the same interface over ObjectHandler and is exercised in the platform (U2Demo). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {TypeAhead} from '../src/components/typeahead.js';

const PROJECTS = [
  {name: 'Chem', owner: 'ada'},
  {name: 'Bio', owner: 'bruno'},
  {name: 'Curves', owner: 'chen'},
];

const renderer = {
  caption: (p) => p.name,
  listItem: (p) => {
    const row = document.createElement('div');
    row.className = 'test-project';
    row.textContent = `${p.name} (${p.owner})`;
    return row;
  },
  markup: (p) => {
    const chipEl = document.createElement('span');
    chipEl.textContent = `#${p.name}`;
    return chipEl;
  },
};

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function mount(component) {
  document.body.append(component.root);
  return component;
}

function open(typeahead) {
  const input = typeahead.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  return input;
}

smoke('renderer alone: rows from listItem, filtering and text from caption', async () => {
  const picker = mount(new TypeAhead({source: PROJECTS, renderer}));
  const input = open(picker);
  await flush();

  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.deepEqual(popup.querySelectorAll('.test-project').map((el) => el.textContent),
    ['Chem (ada)', 'Bio (bruno)', 'Curves (chen)']);

  input.value = 'cur';
  fire(input, 'input');
  await flush();
  assert.deepEqual(popup.querySelectorAll('.test-project').map((el) => el.textContent),
    ['Curves (chen)']);

  fire(input, 'keydown', {key: 'ArrowDown'});
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(picker.selected.value, PROJECTS[2]);
  assert.equal(picker.text.value, 'Curves', 'input text comes from renderer.caption');
  picker.dispose();
});

smoke('renderer without listItem falls back to markup', async () => {
  const picker = mount(new TypeAhead({source: PROJECTS, renderer: {caption: renderer.caption, markup: renderer.markup}}));
  open(picker);
  await flush();
  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.deepEqual(popup.querySelectorAll('.u2-typeahead-option span').map((el) => el.textContent),
    ['#Chem', '#Bio', '#Curves']);
  picker.dispose();
});

smoke('caption-only renderer renders the default text rows', async () => {
  const picker = mount(new TypeAhead({source: PROJECTS, renderer: {caption: (p) => p.name}}));
  open(picker);
  await flush();
  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.deepEqual(popup.querySelectorAll('.u2-typeahead-text').map((el) => el.textContent),
    ['Chem', 'Bio', 'Curves']);
  picker.dispose();
});

smoke('explicit itemText and render win over the renderer members', async () => {
  const picker = mount(new TypeAhead({
    source: PROJECTS,
    renderer,
    itemText: (p) => p.owner,
    render: (p) => {
      const row = document.createElement('div');
      row.className = 'test-override';
      row.textContent = p.owner;
      return row;
    },
  }));
  const input = open(picker);
  await flush();
  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.deepEqual(popup.querySelectorAll('.test-override').map((el) => el.textContent),
    ['ada', 'bruno', 'chen']);
  assert.equal(popup.querySelectorAll('.test-project').length, 0);

  input.value = 'bru';
  fire(input, 'input');
  await flush();
  assert.equal(popup.querySelectorAll('.test-override').length, 1, 'filter uses explicit itemText');
  picker.dispose();
});

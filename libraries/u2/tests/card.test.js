/* Card: sections that appear only when given, signal-bound title slots, and the clickable /
   selectable variants. Same contract as the other UI suites — every test leaves the live-scope
   count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, Scope} from '../src/index.js';
import {Card} from '../src/components/containers/card.js';

function ui(name, body) {
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

function el(text) {
  const node = document.createElement('div');
  node.textContent = text;
  return node;
}

const text = (card, cls) => card.root.querySelector(cls)?.textContent;

ui('sections render only when they are given, in header/media/body/footer order', () => {
  const bare = new Card();
  assert.equal(bare.root.dataset.u2, 'card');
  assert.equal(bare.root.classList.contains('u2-card'), true);
  assert.equal(bare.root.querySelector('.u2-card-header'), null, 'no header without header content');
  assert.equal(bare.root.querySelector('.u2-card-media'), null);
  assert.equal(bare.root.querySelector('.u2-card-footer'), null);
  assert.notEqual(bare.root.querySelector('.u2-card-body'), null, 'the body always exists');
  assert.equal(bare.root.getAttribute('role'), null, 'a plain card is not a button');
  bare.dispose();

  const full = new Card({
    title: 'Aspirin', subtitle: 'NSAID', icon: 'capsules',
    actions: [el('★')],
    media: 'https://example.org/aspirin.png',
    children: [el('180.16 g/mol')],
    footer: [el('ChEMBL')],
  });
  assert.deepEqual([...full.root.children].map((c) => c.className.split(' ').pop()),
    ['u2-card-header', 'u2-card-media', 'u2-card-body', 'u2-card-footer']);
  assert.equal(text(full, '.u2-card-title'), 'Aspirin');
  assert.equal(text(full, '.u2-card-subtitle'), 'NSAID');
  assert.equal(full.root.querySelector('.u2-card-icon').classList.contains('fa-capsules'), true);
  assert.equal(text(full, '.u2-card-actions'), '★');
  assert.equal(full.root.querySelector('.u2-card-media-img').src, 'https://example.org/aspirin.png');
  assert.equal(text(full, '.u2-card-body'), '180.16 g/mol');
  assert.equal(text(full, '.u2-card-footer'), 'ChEMBL');
  full.dispose();

  const titleOnly = new Card({title: 'Only'});
  assert.notEqual(titleOnly.root.querySelector('.u2-card-header'), null);
  assert.equal(titleOnly.root.querySelector('.u2-card-subtitle'), null);
  assert.equal(titleOnly.root.querySelector('.u2-card-actions'), null);
  titleOnly.dispose();

  const element = document.createElement('canvas');
  const custom = new Card({media: element});
  assert.equal(custom.root.querySelector('.u2-card-media').firstChild, element);
  custom.dispose();
});

ui('title and subtitle follow their signals, and stop when the card is disposed', () => {
  const title = signal('Draft');
  const subtitle = signal('v1');
  const card = new Card({title, subtitle});
  assert.equal(text(card, '.u2-card-title'), 'Draft');
  assert.equal(text(card, '.u2-card-subtitle'), 'v1');

  title.value = 'Final';
  subtitle.value = 'v2';
  assert.equal(text(card, '.u2-card-title'), 'Final');
  assert.equal(text(card, '.u2-card-subtitle'), 'v2');

  card.dispose();
  title.value = 'Gone';
  assert.equal(text(card, '.u2-card-title'), 'Final', 'the binding died with the card');
});

ui('clickable: role, tab stop, and click / Enter / Space run onClick', () => {
  let clicks = 0;
  const card = new Card({onClick: () => clicks++});
  assert.equal(card.root.classList.contains('u2-card-clickable'), true);
  assert.equal(card.root.getAttribute('role'), 'button');
  assert.equal(card.root.tabIndex, 0);

  fire(card.root, 'click');
  fire(card.root, 'keydown', {key: 'Enter'});
  fire(card.root, 'keydown', {key: ' '});
  fire(card.root, 'keydown', {key: 'a'});
  assert.equal(clicks, 3, 'other keys do nothing');

  card.dispose();
  fire(card.root, 'click');
  assert.equal(clicks, 3, 'the listeners died with the card');

  const affordance = new Card({clickable: true});
  assert.equal(affordance.root.classList.contains('u2-card-clickable'), true);
  fire(affordance.root, 'click');
  affordance.dispose();
});

ui('selectable: a click toggles the selection, and class and aria follow it', () => {
  const card = new Card({selectable: true});
  assert.equal(card.root.classList.contains('u2-card-selectable'), true);
  assert.equal(card.root.classList.contains('u2-card-clickable'), false, 'selection is its own variant');
  assert.equal(card.root.getAttribute('aria-pressed'), 'false');
  assert.equal(card.selected.value, false);

  fire(card.root, 'click');
  assert.equal(card.selected.value, true);
  assert.equal(card.root.classList.contains('u2-card-selected'), true);
  assert.equal(card.root.getAttribute('aria-pressed'), 'true');

  fire(card.root, 'keydown', {key: ' '});
  assert.equal(card.selected.value, false);
  assert.equal(card.root.classList.contains('u2-card-selected'), false);
  card.dispose();

  const preselected = new Card({selectable: true, selected: true});
  assert.equal(preselected.selected.value, true);
  assert.equal(preselected.root.classList.contains('u2-card-selected'), true);
  preselected.dispose();
});

ui('a Signal handed as selected is adopted, and writes flow both ways', () => {
  const picked = signal(false);
  const card = new Card({selectable: true, selected: picked});
  assert.equal(card.selected, picked, 'the signal is adopted, not copied');

  picked.value = true;
  assert.equal(card.root.classList.contains('u2-card-selected'), true, 'an outside write shows');

  fire(card.root, 'click');
  assert.equal(picked.value, false, 'a click writes the outside signal');
  card.dispose();
});

ui('selectable and onClick both run on one click', () => {
  let clicks = 0;
  const card = new Card({selectable: true, onClick: () => clicks++});
  fire(card.root, 'click');
  assert.equal(clicks, 1);
  assert.equal(card.selected.value, true);
  card.dispose();
});

ui('activation belongs to the innermost control the event landed on', () => {
  let clicks = 0;
  const input = document.createElement('input');
  const action = document.createElement('button');
  const heading = document.createElement('h3');
  heading.textContent = 'Aspirin';
  const card = new Card({selectable: true, onClick: () => clicks++,
    actions: [action], children: [input, heading]});

  assert.equal(fire(input, 'keydown', {key: ' '}), true, 'Space in a child input is not swallowed');
  assert.equal(card.selected.value, false, 'and it never activates the card');

  fire(action, 'click');
  assert.equal(clicks, 0, 'a click on a child button is the button\'s own');
  assert.equal(card.selected.value, false);

  fire(heading, 'click');
  assert.equal(clicks, 1, 'a click on plain content still activates');
  assert.equal(card.selected.value, true);

  const inner = new Card({selectable: true});
  card.body.append(inner.root);
  fire(inner.root, 'click');
  assert.equal(inner.selected.value, true, 'the inner card takes its own click');
  assert.equal(clicks, 1, 'and it never reaches the card around it');
  assert.equal(card.selected.value, true);
  inner.dispose();
  card.dispose();
});

ui('the body is exposed, so content can be appended after construction', () => {
  const card = new Card({children: [el('first')]});
  assert.equal(card.body.classList.contains('u2-card-body'), true);
  card.body.append(el('second'));
  assert.equal(card.body.children.length, 2);
  assert.equal(card.body.textContent, 'firstsecond');
  card.dispose();
});

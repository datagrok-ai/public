/* IconInput: the control, the anchored popup (search + VirtualGrid), pick/cancel semantics and
   the brand face. The shim lays nothing out, so the popup's grid is given a viewport. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, signal} from '../src/index.js';
import {IconInput} from '../src/components/inputs/icon-input.js';
import {ICON_NAMES, BRAND_ICON_NAMES} from '../src/components/display/icon-names.js';

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

/** The popup is anchored to the control, and an overlay whose anchor left the document closes
 * itself — so every input that opens one is mounted. */
function make(options = {}) {
  const input = new IconInput({label: 'Icon', ...options});
  document.body.append(input.root);
  return input;
}

function control(input) {
  return input.root.querySelector('.u2-icon-input');
}

function popup() {
  return document.body.querySelector('.u2-icon-input-popup');
}

function open(input) {
  fire(control(input), 'click');
  const box = popup();
  const grid = box.querySelector('.u2-grid');
  grid.clientWidth = 280;
  grid.clientHeight = 224;
  fire(grid, 'scroll');
  return box;
}

function shown(box) {
  return box.querySelectorAll('.u2-grid-cell').map((c) => c.title);
}

function type(box, text) {
  const field = box.querySelector('input');
  field.value = text;
  fire(field, 'input');
}

ui('the control shows the current glyph and name, or a dash when empty', () => {
  const input = make({value: 'star'});
  assert.equal(input.root.dataset.u2, 'icon-input');
  assert.equal(control(input).getAttribute('role'), 'button');
  assert.equal(control(input).getAttribute('aria-expanded'), 'false');
  assert.notEqual(input.root.querySelector('.u2-icon-input-preview .fa-star'), null);
  assert.equal(input.root.querySelector('.u2-icon-input-name').textContent, 'star');
  input.value.value = '';
  assert.equal(input.root.querySelector('.u2-icon-input-preview').children.length, 0);
  assert.equal(input.root.querySelector('.u2-icon-input-name').textContent, '—');
  input.dispose();
});

ui('click opens the popup over every icon; typing filters by name with spaces as hyphens', () => {
  const input = make();
  const box = open(input);
  assert.equal(control(input).getAttribute('aria-expanded'), 'true');
  assert.equal(box.querySelector('.u2-grid-content').style.height,
    `${Math.ceil((ICON_NAMES.length + BRAND_ICON_NAMES.length) / 10) * 28}px`, 'all names, 10 per row');
  type(box, 'Arrow to');
  const names = shown(box);
  assert.equal(names.includes('arrow-to-bottom'), true);
  assert.equal(names.every((n) => n.includes('arrow-to')), true);
  type(box, 'zzzz');
  assert.equal(shown(box).length, 0);
  assert.equal(input.value.value, '', 'filtering is not a pick');
  input.dispose();
  assert.equal(popup(), null, 'disposing the input takes the popup with it');
});

ui('a click on a cell commits the name and closes; Esc and an outside click close without a write', () => {
  const changes = [];
  const input = make({value: 'star', onChanged: (v) => changes.push(v)});
  let box = open(input);
  type(box, 'flask');
  fire(box.querySelectorAll('.u2-grid-cell').find((c) => c.title === 'flask').firstChild, 'click');
  assert.equal(input.value.value, 'flask');
  assert.deepEqual(changes, ['flask']);
  assert.equal(popup(), null);
  assert.equal(control(input).getAttribute('aria-expanded'), 'false');

  box = open(input);
  fire(document.body, 'pointerdown');
  assert.equal(popup(), null, 'outside pointerdown closes');
  box = open(input);
  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(popup(), null, 'Escape closes');
  assert.equal(input.value.value, 'flask', 'neither dismissal wrote');
  assert.deepEqual(changes, ['flask']);
  input.dispose();
});

ui('arrows typed in the search box drive the grid, and Enter picks the lead', () => {
  const input = make();
  const box = open(input);
  type(box, 'flask');
  const field = box.querySelector('input');
  fire(field, 'keydown', {key: 'ArrowRight'});
  fire(field, 'keydown', {key: 'ArrowRight'});
  const grid = box.querySelector('.u2-grid');
  assert.equal(grid.getAttribute('aria-activedescendant') !== null, true);
  const lead = shown(box)[1];
  fire(field, 'keydown', {key: 'Enter'});
  assert.equal(input.value.value, lead);
  assert.equal(popup(), null);
  input.dispose();
});

ui('Enter and Space on the control toggle the popup; Backspace clears a nullable value', () => {
  const input = make({value: 'star'});
  fire(control(input), 'keydown', {key: 'Enter'});
  assert.notEqual(popup(), null);
  fire(control(input), 'keydown', {key: ' '});
  assert.equal(popup(), null);
  fire(control(input), 'keydown', {key: 'Backspace'});
  assert.equal(input.value.value, '');
  const required = make({value: 'star', nullable: false});
  fire(control(required), 'keydown', {key: 'Delete'});
  assert.equal(required.value.value, 'star');
  input.dispose();
  required.dispose();
});

ui('a disabled input does not open; a bound signal is adopted; a custom name list is honoured', () => {
  const bound = signal('cog');
  const input = make({bind: bound, names: ['cog', 'star', 'github']});
  input.enabled = false;
  fire(control(input), 'click');
  assert.equal(popup(), null);
  input.enabled = true;
  const box = open(input);
  assert.deepEqual(shown(box), ['cog', 'star', 'github']);
  assert.notEqual(box.querySelectorAll('.u2-grid-cell')[2].querySelector('.fab.fa-github'), null,
    'a brand renders in its one face');
  fire(box.querySelectorAll('.u2-grid-cell')[1].firstChild, 'click');
  assert.equal(bound.value, 'star');
  input.dispose();
});

ui('opening and closing repeatedly leaves no scope behind', () => {
  const input = make();
  const live = Scope.liveCount;
  for (let i = 0; i < 5; i++) {
    open(input);
    fire(document.body, 'pointerdown');
  }
  assert.equal(Scope.liveCount, live);
  input.dispose();
});

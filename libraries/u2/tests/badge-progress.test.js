/* Badges and progress on the node DOM shim: variant skins, signal-bound counts, tag removal and
   the bar's clamp / width / aria contract. Same contract as smoke.test.js — every test leaves the
   live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, computed, Scope, Control} from '../src/index.js';
import {badge, countBadge, dot, tag} from '../src/components/badge.js';
import {ProgressBar} from '../src/components/progress.js';

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

function track(bar) {
  return bar.root.querySelector('[role="progressbar"]');
}

function fill(bar) {
  return bar.root.querySelector('.u2-progress-fill');
}

ui('badge: every variant gets its own class, and text follows a signal', () => {
  const plain = badge('New');
  assert.equal(plain.textContent, 'New');
  assert.equal(plain.dataset.u2, 'badge');
  assert.equal(plain.classList.contains('u2-badge'), true);
  assert.equal(plain.classList.contains('u2-badge-default'), true);

  for (const variant of ['accent', 'success', 'warning', 'error'])
    assert.equal(badge('x', {variant}).classList.contains(`u2-badge-${variant}`), true);

  const status = signal('Draft');
  const owner = Control.build(() => badge(status, {variant: 'warning'}));
  assert.equal(owner.root.textContent, 'Draft');
  status.value = 'Published';
  assert.equal(owner.root.textContent, 'Published');

  owner.dispose();
  status.value = 'Archived';
  assert.equal(owner.root.textContent, 'Published', 'the binding died with its owner');
});

ui('countBadge: caps at max, hides at zero and follows its signal', () => {
  assert.equal(countBadge(5).textContent, '5');
  assert.equal(countBadge(5).dataset.u2, 'count-badge');
  assert.equal(countBadge(120).textContent, '99+');
  assert.equal(countBadge(120, {max: 9}).textContent, '9+');
  assert.equal(countBadge(0).style.display, 'none');
  assert.equal(countBadge(1).style.display, '');

  const unread = signal(0);
  const owner = Control.build(() => countBadge(unread));
  const el = owner.root;
  assert.equal(el.style.display, 'none', 'zero hides the badge');

  unread.value = 3;
  assert.equal(el.textContent, '3');
  assert.equal(el.style.display, '');

  unread.value = 250;
  assert.equal(el.textContent, '99+');

  owner.dispose();
  unread.value = 7;
  assert.equal(el.textContent, '99+', 'the effect died with its owner');
});

ui('badge and countBadge: a signal argument needs an ambient owner', () => {
  assert.throws(() => badge(signal('x')), /needs an owner/);
  assert.throws(() => countBadge(signal(1)), /needs an owner/);
});

ui('dot: variant class and data-u2', () => {
  const plain = dot();
  assert.equal(plain.dataset.u2, 'dot');
  assert.equal(plain.classList.contains('u2-dot-default'), true);
  assert.equal(dot('success').classList.contains('u2-dot-success'), true);
  assert.equal(plain.getAttribute('aria-hidden'), 'true');
});

ui('tag: reports removal, and has no ✕ without a callback', () => {
  const plain = tag('Alpha');
  assert.equal(plain.dataset.u2, 'tag');
  assert.equal(plain.textContent, 'Alpha');
  assert.equal(plain.querySelector('.u2-tag-remove'), null, 'no callback, no remove button');

  const removed = [];
  const closable = tag('Beta', {onRemove: () => removed.push('Beta')});
  const remove = closable.querySelector('.u2-tag-remove');
  assert.equal(remove.getAttribute('aria-label'), 'Remove Beta');
  fire(remove, 'click');
  assert.deepEqual(removed, ['Beta']);
});

ui('ProgressBar: clamps the value, drives the fill width and aria-valuenow', () => {
  const bar = new ProgressBar({showPercent: true});
  const percent = bar.root.querySelector('.u2-progress-percent');
  assert.equal(bar.root.dataset.u2, 'progress-bar');
  assert.equal(track(bar).getAttribute('aria-valuemin'), '0');
  assert.equal(track(bar).getAttribute('aria-valuemax'), '1');
  assert.equal(fill(bar).style.width, '0%');

  bar.value.value = 0.4;
  assert.equal(fill(bar).style.width, '40%');
  assert.equal(track(bar).getAttribute('aria-valuenow'), '0.4');
  assert.equal(percent.textContent, '40%');

  bar.value.value = 1.5;
  assert.equal(bar.value.value, 1, 'above the range clamps to 1');
  assert.equal(fill(bar).style.width, '100%');
  assert.equal(percent.textContent, '100%');

  bar.value.value = -2;
  assert.equal(bar.value.value, 0);
  assert.equal(fill(bar).style.width, '0%');
  assert.equal(track(bar).getAttribute('aria-valuenow'), '0');
  bar.dispose();
});

ui('ProgressBar: indeterminate drops aria-valuenow, and switching back restores it', () => {
  const bar = new ProgressBar({indeterminate: true, showPercent: true});
  const percent = bar.root.querySelector('.u2-progress-percent');
  assert.equal(bar.root.classList.contains('u2-progress-indeterminate'), true);
  assert.equal(track(bar).hasAttribute('aria-valuenow'), false);

  bar.value.value = 0.5;
  assert.equal(track(bar).hasAttribute('aria-valuenow'), false, 'a value means nothing while indeterminate');
  assert.equal(fill(bar).style.width, '', 'the sliding fill is sized by CSS');
  assert.equal(percent.textContent, '');

  bar.indeterminate = false;
  assert.equal(bar.root.classList.contains('u2-progress-indeterminate'), false);
  assert.equal(track(bar).getAttribute('aria-valuenow'), '0.5');
  assert.equal(fill(bar).style.width, '50%');
  assert.equal(percent.textContent, '50%');

  bar.indeterminate = true;
  assert.equal(track(bar).hasAttribute('aria-valuenow'), false);
  bar.dispose();
});

ui('ProgressBar: the description line takes a signal and dies with the bar', () => {
  const value = signal(0);
  const bar = new ProgressBar({description: computed(() => `${Math.round(value.value * 100)} of 100 rows`)});
  bar.effect(() => bar.value.value = value.value);
  const description = bar.root.querySelector('.u2-progress-description');
  assert.equal(description.textContent, '0 of 100 rows');

  value.value = 0.25;
  assert.equal(description.textContent, '25 of 100 rows');
  assert.equal(fill(bar).style.width, '25%');

  bar.dispose();
  value.value = 1;
  assert.equal(description.textContent, '25 of 100 rows');
});

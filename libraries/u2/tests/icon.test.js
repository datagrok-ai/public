/* Icons on the node DOM shim: the platform class convention, variant/size/cls mapping and the
   tooltip fallback. Same contract as smoke.test.js — every test leaves the live-scope count
   where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope, Component, Tooltip} from '../src/index.js';
import {icon} from '../src/components/icon.js';

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

ui('icon: platform class convention, in order, with no ambient scope', () => {
  assert.equal(Scope.ambient, undefined);
  const el = icon('search');
  assert.equal(el.tagName, 'I');
  assert.equal(el.className, 'grok-icon u2-icon fal fa-search');
  assert.equal(el.dataset.u2, 'icon');
});

ui('icon: variants map to the Font Awesome weight classes', () => {
  assert.equal(icon('star', {variant: 'light'}).classList.contains('fal'), true);
  assert.equal(icon('star', {variant: 'regular'}).classList.contains('far'), true);
  assert.equal(icon('star', {variant: 'solid'}).classList.contains('fas'), true);
  assert.equal(icon('star', {variant: 'solid'}).classList.contains('fal'), false);
});

ui('icon: size adds a sizing class, and no size adds none', () => {
  for (const size of ['small', 'medium', 'large'])
    assert.equal(icon('cog', {size}).classList.contains(`u2-icon-${size}`), true);
  assert.equal(icon('cog').className, 'grok-icon u2-icon fal fa-cog');
});

ui('icon: cls appends, including several classes at once', () => {
  const el = icon('sync', {variant: 'solid', size: 'small', cls: 'u2-icon-action my-icon'});
  assert.equal(el.className, 'grok-icon u2-icon fas fa-sync u2-icon-small u2-icon-action my-icon');
});

ui('icon: without an owner the tooltip falls back to the title attribute', () => {
  const el = icon('info-circle', {tooltip: 'What this does'});
  assert.equal(el.title, 'What this does');
  assert.equal(el.getAttribute('aria-label'), 'What this does');
  assert.equal(icon('info-circle').hasAttribute('title'), false);
});

ui('icon: with an owner the tooltip binds to the scope and dies with it', () => {
  const owner = Component.build(() => icon('info-circle', {tooltip: 'Bound'}));
  const el = owner.root;
  assert.equal(el.hasAttribute('title'), false, 'the bound tooltip replaces the title attribute');
  assert.equal(el.getAttribute('aria-label'), 'Bound');

  Tooltip.show('Bound', el);
  assert.equal(Tooltip.isVisible, true);
  owner.dispose();
  assert.equal(Tooltip.isVisible, false, 'disposing the owner hides its tooltip');
});

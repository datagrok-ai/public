/* EntityCard (src/dg/entities/entity.ts) over tests/dg-stub.mjs: the platform card's gallery title
   is the card's title part and its label, read live; a card without one is named by the object's
   caption. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';

register('./dg-stub.mjs', import.meta.url);
const {entityCard} = await import('../src/dg/entities/entity.js');

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

function galleryCard(name) {
  const card = document.createElement('div');
  card.className = 'grok-gallery-grid-item';
  const title = document.createElement('div');
  title.className = 'grok-gallery-grid-item-title';
  title.textContent = name;
  const date = document.createElement('div');
  date.textContent = 'updated today';
  card.append(title, date);
  return card;
}

smoke('the gallery title is the title part and the label, and follows a rename', async () => {
  const card = entityCard({name: 'Lab space'},
    {renderer: {caption: (x) => x.name, card: (x) => galleryCard(x.name)}, currentOnClick: false});
  document.body.append(card.root);
  const title = card.root.querySelector('[data-u2-part="title"]');
  assert.equal(title?.textContent, 'Lab space');
  assert.equal(card.getState().label, 'Lab space');

  title.textContent = 'Renamed space';
  assert.equal(card.getState().label, 'Renamed space');
  card.dispose();
});

smoke('a card without a gallery title is named by the caption', async () => {
  const markup = (x) => {
    const el = document.createElement('span');
    el.textContent = `#${x.name}`;
    return el;
  };
  const card = entityCard({name: 'Chem'}, {renderer: {caption: (x) => x.name, markup}, currentOnClick: false});
  document.body.append(card.root);
  assert.equal(card.root.querySelector('[data-u2-part="title"]'), null);
  assert.equal(card.getState().label, 'Chem');
  card.dispose();
});

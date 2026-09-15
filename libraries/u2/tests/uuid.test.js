/* `uuid4` — the v4 the draft ids and the memory backend's stamps are made of, without
   `crypto.randomUUID` (secure contexts only). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {uuid4} from '../src/core/uuid.js';
import {Rows} from '../src/sources/rows-like.js';

function scoped(name, body) {
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

scoped('the shape is a v4 uuid, and every one is its own', () => {
  const made = new Set();
  for (let i = 0; i < 1000; i++) {
    const id = uuid4();
    assert.match(id, /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/);
    made.add(id);
  }
  assert.equal(made.size, 1000);
});

scoped('a draft id is the prefix and one of them', () => {
  const id = Rows.draftId();
  assert.equal(Rows.isDraft(id), true);
  assert.match(id.slice(Rows.DRAFT_PREFIX.length), /^[0-9a-f-]{36}$/);
});

/* The sample gallery's contract (P3.5 F12): every sample parses and renders with zero broken
   nodes against a bare `SpecContext` — no bind, no on, no dependence on any hosting context —
   and the localStorage wrapper behind 'Save to gallery…' never throws or answers garbage. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {SpecContext, parseSpec, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {SAMPLES} from '../src/dg/designer/samples.js';
import {brokenCount} from '../src/dg/designer/node-ref.js';
import {GALLERY_KEY, listGallery, loadGallery, saveToGallery} from '../src/dg/designer/gallery.js';

registerAll();

/** Every test runs against a clean document and must leave the live-scope count where it was. */
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

scoped('every sample parses and renders zero-broken against a bare SpecContext', () => {
  assert.deepEqual(SAMPLES.map((s) => s.name),
    ['Blank', 'Data entry', 'Master–detail', 'Settings', 'Dashboard']);
  for (const sample of SAMPLES) {
    const spec = parseSpec(JSON.parse(JSON.stringify(sample.spec)));
    const instance = renderSpec(spec, new SpecContext());
    assert.equal(brokenCount(instance), 0, `${sample.name} renders with no broken nodes`);
    assert.ok(instance.nodes().size >= 1, `${sample.name} rendered something`);
    instance.dispose();
  }
});

scoped('no sample carries a bind or an event — the bare-context contract holds structurally', () => {
  const walk = (node, visit) => {
    visit(node);
    for (const child of node.children ?? [])
      walk(child, visit);
  };
  for (const sample of SAMPLES) {
    walk(sample.spec.root, (node) => {
      assert.equal(node.bind, undefined, `${sample.name}: ${node.tag} binds nothing`);
      assert.equal(node.on, undefined, `${sample.name}: ${node.tag} wires no command`);
    });
  }
});

function stubStorage() {
  const map = new Map();
  globalThis.localStorage = {
    getItem: (k) => map.has(k) ? map.get(k) : null,
    setItem: (k, v) => map.set(k, String(v)),
    removeItem: (k) => map.delete(k),
  };
  return map;
}

test('gallery wrapper: save, load, list — an existing name overwrites', () => {
  const map = stubStorage();
  try {
    assert.deepEqual(loadGallery(), {});
    assert.deepEqual(listGallery(), []);
    assert.equal(saveToGallery('My form', {a: 1}), true, 'a save that landed says so');
    saveToGallery('Other', {b: 2});
    assert.deepEqual(listGallery(), ['My form', 'Other']);
    assert.deepEqual(loadGallery()['My form'], {a: 1});
    saveToGallery('My form', {a: 2});
    assert.deepEqual(loadGallery()['My form'], {a: 2}, 'the same name overwrites');
    assert.equal(listGallery().length, 2);
    assert.ok(map.has(GALLERY_KEY), 'everything lives under the one key');
  } finally {
    delete globalThis.localStorage;
  }
});

test('gallery wrapper: corrupt JSON, a non-record payload and missing storage all answer empty', () => {
  stubStorage();
  try {
    globalThis.localStorage.setItem(GALLERY_KEY, '{nope');
    assert.deepEqual(loadGallery(), {});
    assert.deepEqual(listGallery(), []);
    globalThis.localStorage.setItem(GALLERY_KEY, '[1, 2]');
    assert.deepEqual(loadGallery(), {}, 'a non-record payload is refused');
    globalThis.localStorage.setItem(GALLERY_KEY, 'null');
    assert.deepEqual(loadGallery(), {});
  } finally {
    delete globalThis.localStorage;
  }
  assert.deepEqual(loadGallery(), {}, 'no storage at all is just empty');
  assert.equal(saveToGallery('x', {a: 1}), false, 'and a save without storage reports the failure');
});

test('gallery wrapper: a storage whose getItem throws answers empty instead of throwing', () => {
  globalThis.localStorage = {
    getItem: () => {
      throw new Error('SecurityError');
    },
    setItem: () => {
      throw new Error('SecurityError');
    },
  };
  try {
    assert.deepEqual(loadGallery(), {}, 'a SecurityError never escapes to the caller');
    assert.deepEqual(listGallery(), []);
    assert.equal(saveToGallery('x', {a: 1}), false, 'a refused write reports the failure');
  } finally {
    delete globalThis.localStorage;
  }
});

// emitFiltered — the fence-aware chunk filter that turns SDK text_deltas into client `chunk`
// messages. It must never lose or duplicate a character: what the browser accumulates has to
// equal what the model emitted. Fence state is keyed per session id, so each test uses its own.

const test = require('node:test');
const assert = require('node:assert');
const {emitFiltered, flushFenceState} = require('../dist/stream-filter.js');

let seq = 0;
// Collects emitted chunks for one throwaway session, and reassembles what the browser would show.
function sink() {
  const sid = `t${seq++}`;
  const chunks = [];
  const ws = {send: (s) => {
    const m = JSON.parse(s);
    if (m.type === 'chunk')
      chunks.push(m.content);
  }};
  return {
    sid, ws, chunks,
    feed: (...deltas) => deltas.forEach((d) => emitFiltered(ws, sid, d)),
    flush: () => flushFenceState(ws, sid),
    get text() { return chunks.join(''); },
  };
}

test('plain prose passes through unchanged', () => {
  const s = sink();
  s.feed('Hello, ', 'world.\n');
  s.flush();
  assert.equal(s.text, 'Hello, world.\n');
});

test('text split mid-word across deltas reassembles exactly', () => {
  const s = sink();
  s.feed('The quick ', 'brown fo', 'x jumps.\n');
  s.flush();
  assert.equal(s.text, 'The quick brown fox jumps.\n');
});

test('a fenced code block survives intact', () => {
  const s = sink();
  const src = 'Here:\n```js\nconst a = 1;\n```\nDone.\n';
  s.feed(src);
  s.flush();
  assert.equal(s.text, src);
});

test('a fence arriving one character at a time is not mangled', () => {
  const s = sink();
  const src = 'a\n```js\nx\n```\nb\n';
  for (const ch of src)
    s.feed(ch);
  s.flush();
  assert.equal(s.text, src, 'char-by-char streaming must reassemble byte-identically');
});

test('a partial trailing line starting with a backtick is held back until resolved', () => {
  const s = sink();
  s.feed('text\n', '``');
  assert.ok(!s.text.includes('``'), 'a possible fence marker must not be emitted early');
  s.feed('`js\n');
  s.flush();
  assert.equal(s.text, 'text\n```js\n');
});

test('inline backticks mid-line are not mistaken for a fence', () => {
  const s = sink();
  const src = 'Use the `grok.shell.v` global here.\n';
  s.feed(src);
  s.flush();
  assert.equal(s.text, src);
});

test('a line that merely starts with a backtick but is not a fence still emits', () => {
  const s = sink();
  s.feed('`not a fence` really\n');
  s.flush();
  assert.equal(s.text, '`not a fence` really\n');
});

test('multi-line deltas batch into few emits rather than one per line', () => {
  const s = sink();
  s.feed('a\nb\nc\nd\ne\nf\n');
  s.flush();
  assert.equal(s.text, 'a\nb\nc\nd\ne\nf\n');
  assert.ok(s.chunks.length <= 2, `expected batching, got ${s.chunks.length} chunks`);
});

test('flush emits a held-back partial so nothing is lost at end of turn', () => {
  const s = sink();
  s.feed('done\n', '`lonely');
  assert.ok(!s.text.includes('lonely'), 'held back while it could still become a fence');
  s.flush();
  assert.equal(s.text, 'done\n`lonely', 'flush must release it');
});

test('flushing twice does not duplicate the carry', () => {
  const s = sink();
  s.feed('x\n', '`held');
  s.flush();
  s.flush();
  assert.equal(s.text, 'x\n`held');
});

test('sessions do not leak fence state into each other', () => {
  const a = sink();
  const b = sink();
  a.feed('```js\n');       // a is now inside a fence
  b.feed('plain text\n');  // b must be unaffected
  a.feed('code\n```\n');
  a.flush();
  b.flush();
  assert.equal(a.text, '```js\ncode\n```\n');
  assert.equal(b.text, 'plain text\n');
});

test('consecutive fenced blocks toggle mode correctly', () => {
  const s = sink();
  const src = '```js\none\n```\nmid\n```py\ntwo\n```\nend\n';
  s.feed(src);
  s.flush();
  assert.equal(s.text, src);
});

test('a fence marker with no language is handled', () => {
  const s = sink();
  const src = 'a\n```\nraw\n```\nb\n';
  s.feed(src);
  s.flush();
  assert.equal(s.text, src);
});

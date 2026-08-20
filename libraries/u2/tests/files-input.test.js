/* filesInput — per-file rows, progress, cancellation and the emit-when-settled contract
   (plan WO-15.3). `DG.FileInfo` and `grok.shell` come from tests/dg-stub.mjs; the files are plain
   objects carrying the three members the input reads (name, size, stream), which is what lets the
   progress and cancel paths be driven exactly. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';

register('./dg-stub.mjs', import.meta.url);
const {filesInput} = await import('../src/dg/files-input.js');
const grok = await import('datagrok-api/grok');

function filesTest(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const messages = [];
    grok.shell.error = (m) => messages.push(m);
    grok.shell.warning = (m) => messages.push(m);
    try {
      await body(messages);
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function file(name, bytes = [1, 2, 3]) {
  const data = new Uint8Array(bytes);
  return {name, size: data.length, stream: () => new ReadableStream({
    start(c) {
      c.enqueue(data);
      c.close();
    },
  })};
}

/** A file whose second half arrives only when the test says so. */
function pausedFile(name) {
  let resume = () => {};
  const gate = new Promise((r) => resume = r);
  const f = {name, size: 4, stream: () => new ReadableStream({
    async start(c) {
      c.enqueue(new Uint8Array([1, 2]));
      await gate;
      c.enqueue(new Uint8Array([3, 4]));
      c.close();
    },
  })};
  return {file: f, resume: () => resume()};
}

function failingFile(name) {
  return {name, size: 4, stream: () => new ReadableStream({
    start(c) {
      c.error(new Error('read failed'));
    },
  })};
}

async function settle() {
  await flush();
  await flush();
}

const rows = (input) => [...input.root.querySelectorAll('.u2-files-row')];
const names = (input) => rows(input).map((r) => r.querySelector('.u2-files-name').textContent);
const summary = (input) => input.root.querySelector('.u2-files-summary').textContent;

filesTest('the picker sits on the shared options rail', async () => {
  const input = filesInput('Files');
  const rail = input.root.querySelector('.u2-input-options');
  assert.equal(rail.parentNode, input.box, 'so it crosses the platform bridge with the editor');
  assert.equal(rail.querySelectorAll('input').length, 1, 'the hidden file picker');
  assert.equal(rail.querySelectorAll('.u2-btn').length, 1, 'and the folder-open icon');
  input.dispose();
});

filesTest('several files load into rows and one value', async () => {
  const input = filesInput('Files');
  input.add([file('a.csv'), file('b.csv'), file('c.csv')]);
  assert.deepEqual(names(input), ['a.csv', 'b.csv', 'c.csv']);

  await settle();
  assert.equal(input.value.peek().length, 3);
  assert.deepEqual(input.value.peek().map((f) => f.fileName), ['a.csv', 'b.csv', 'c.csv']);
  assert.equal(summary(input), '3 loaded');
  assert.equal(input.validity.peek(), null);
  input.dispose();
});

filesTest('the value waits for the last read, and progress shows meanwhile', async () => {
  const input = filesInput('Files');
  const paused = pausedFile('big.csv');
  input.add([file('small.csv'), paused.file]);
  await settle();

  assert.equal(input.validity.peek(), 'Files are still loading.');
  assert.deepEqual(input.value.peek(), [], 'nothing is published while a read is running');
  assert.equal(rows(input)[1].querySelector('.u2-files-status').textContent, '50%');

  paused.resume();
  await settle();
  assert.equal(input.value.peek().length, 2);
  assert.equal(rows(input)[1].querySelector('.u2-files-status').textContent, '');
  assert.equal(input.validity.peek(), null);
  input.dispose();
});

filesTest('adding a file again replaces it instead of duplicating it', async () => {
  const input = filesInput('Files');
  input.add([file('a.csv', [1, 2, 3])]);
  await settle();
  input.add([file('a.csv', [1, 2, 3, 4, 5])]);
  await settle();

  assert.deepEqual(names(input), ['a.csv']);
  assert.equal(input.value.peek().length, 1);
  assert.equal(input.value.peek()[0].data.length, 5, 'the newer file won');
  input.dispose();
});

filesTest('✕ drops the file from the rows and from the value', async () => {
  const input = filesInput('Files');
  input.add([file('a.csv'), file('b.csv')]);
  await settle();

  fire(rows(input)[0].querySelector('button'), 'click');
  await settle();
  assert.deepEqual(names(input), ['b.csv']);
  assert.deepEqual(input.value.peek().map((f) => f.fileName), ['b.csv']);
  input.dispose();
});

filesTest('a file that cannot be read is reported, and blocks the input', async (messages) => {
  const input = filesInput('Files');
  input.add([file('a.csv'), failingFile('broken.csv')]);
  await settle();

  assert.deepEqual(messages, ['Failed to load broken.csv']);
  assert.equal(input.validity.peek(), 'Some files failed to load.');
  assert.equal(summary(input), '1 loaded, 1 failed');
  assert.deepEqual(input.value.peek().map((f) => f.fileName), ['a.csv'],
    'what did load is still published');
  input.dispose();
});

filesTest('accept refuses what it does not name', async (messages) => {
  const input = filesInput('Files', {accept: ['.csv', '.sdf']});
  input.add([file('a.CSV'), file('notes.txt')]);
  await settle();

  assert.deepEqual(names(input), ['a.CSV'], 'the extension test is case-insensitive');
  assert.deepEqual(messages, ['notes.txt is not .csv or .sdf.']);
  input.dispose();
});

filesTest('an empty input is invalid only when it is not nullable', async () => {
  const optional = filesInput('Files');
  assert.equal(optional.validity.peek(), null);
  optional.dispose();

  const required = filesInput('Files', {nullable: false});
  assert.equal(required.validity.peek(), 'At least one file is required.');
  required.add([file('a.csv')]);
  await settle();
  assert.equal(required.validity.peek(), null);
  required.dispose();
});

filesTest('a value written from outside becomes the rows', async () => {
  const input = filesInput('Files');
  input.value.value = [{fileName: 'remote.csv', path: 'geo/remote.csv'}];
  assert.deepEqual(names(input), ['remote.csv']);
  assert.equal(summary(input), '1 loaded');
  assert.equal(input.validity.peek(), null);
  input.dispose();
});

filesTest('a read still running when the input goes away is cancelled, not settled', async () => {
  const input = filesInput('Files');
  const paused = pausedFile('big.csv');
  input.add([paused.file]);
  await settle();
  assert.deepEqual(input.value.peek(), []);

  input.dispose();
  paused.resume();
  await settle();
  assert.deepEqual(input.value.peek(), [], 'nothing is written into a disposed component');
  assert.equal(rows(input).length, 1, 'and the rows are not rebuilt');
});

/* filesInput — per-file rows, progress, cancellation and the emit-when-settled contract
   (plan WO-15.3). `DG.FileInfo` and `grok.shell` come from tests/dg-stub.mjs; the files are plain
   objects carrying the three members the input reads (name, size, stream), which is what lets the
   progress and cancel paths be driven exactly. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {FileInfo} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {filesInput} = await import('../src/dg/inputs/files-input.js');
const grok = await import('datagrok-api/grok');
const ui = await import('datagrok-api/ui');

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
  assert.equal(input.peek().length, 3);
  assert.deepEqual(input.peek().map((f) => f.fileName), ['a.csv', 'b.csv', 'c.csv']);
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
  assert.deepEqual(input.peek(), [], 'nothing is published while a read is running');
  assert.equal(rows(input)[1].querySelector('.u2-files-status').textContent, '50%');

  paused.resume();
  await settle();
  assert.equal(input.peek().length, 2);
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
  assert.equal(input.peek().length, 1);
  assert.equal(input.peek()[0].data.length, 5, 'the newer file won');
  input.dispose();
});

filesTest('✕ drops the file from the rows and from the value', async () => {
  const input = filesInput('Files');
  input.add([file('a.csv'), file('b.csv')]);
  await settle();

  fire(rows(input)[0].querySelector('button'), 'click');
  await settle();
  assert.deepEqual(names(input), ['b.csv']);
  assert.deepEqual(input.peek().map((f) => f.fileName), ['b.csv']);
  input.dispose();
});

filesTest('a file that cannot be read is reported, and blocks the input', async (messages) => {
  const input = filesInput('Files');
  input.add([file('a.csv'), failingFile('broken.csv')]);
  await settle();

  assert.deepEqual(messages, ['Failed to load broken.csv']);
  assert.equal(input.validity.peek(), 'Some files failed to load.');
  assert.equal(summary(input), '1 loaded, 1 failed');
  assert.deepEqual(input.peek().map((f) => f.fileName), ['a.csv'],
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
  input.value = [{fileName: 'remote.csv', path: 'geo/remote.csv'}];
  assert.deepEqual(names(input), ['remote.csv']);
  assert.equal(summary(input), '1 loaded');
  assert.equal(input.validity.peek(), null);
  input.dispose();
});

/** The platform's "Select a file" door (`ui.pickFileFromFiles`) for the body, answering `pick(options)`;
 * removed after. */
async function withFilesDialog(pick, body) {
  ui.setPickFileFromFiles(pick);
  try {
    await body();
  } finally {
    ui.setPickFileFromFiles(undefined);
  }
}

function pickFromFiles(input) {
  fire(input.root.querySelector('.u2-input-options [aria-label="Add file from Files"]'), 'click');
}

filesTest('where the platform has its Files dialog, its icon follows Open files on the rail', () =>
  withFilesDialog(async () => null, async () => {
    const input = filesInput('Files');
    const icons = [...input.root.querySelectorAll('.u2-input-options [data-u2="icon-button"]')];
    assert.deepEqual(icons.map((el) => el.getAttribute('aria-label')), ['Open files', 'Add file from Files'],
      'the platform order, its tooltips verbatim');
    assert.ok(icons[1].querySelector('.fa-folder-tree'));
    input.dispose();
  }));

filesTest('files picked from Files join the rows as loaded, keyed by their share path', () => {
  const asked = [];
  return withFilesDialog(async (options) => {
    asked.push(options);
    return [new FileInfo('Admin:Home/a.csv'), new FileInfo('Demo:Files/a.csv')];
  }, async () => {
    const changes = [];
    const input = filesInput('Files', {onChanged: (v) => changes.push(v.map((f) => f.fullPath))});
    input.add([file('local.csv')]);
    await settle();
    pickFromFiles(input);
    await settle();
    assert.deepEqual(asked, [{multiple: true}]);
    assert.deepEqual(names(input), ['local.csv', 'Admin:Home/a.csv', 'Demo:Files/a.csv'],
      'two shares\' files of one name are two rows');
    assert.equal(summary(input), '3 loaded');
    assert.deepEqual(changes.at(-1), ['local.csv', 'Admin:Home/a.csv', 'Demo:Files/a.csv']);

    pickFromFiles(input);
    await settle();
    assert.equal(rows(input).length, 3, 'picking a file again replaces its row');
    input.dispose();
  });
});

filesTest('a folder, or a file accept does not name, picked from Files is refused', (messages) =>
  withFilesDialog(async () => [new FileInfo('Admin:Home/geo', {directory: true}),
    new FileInfo('Admin:Home/notes.txt'), new FileInfo('Admin:Home/b.csv')], async () => {
    const input = filesInput('Files', {accept: ['.csv']});
    pickFromFiles(input);
    await settle();
    assert.deepEqual(names(input), ['Admin:Home/b.csv']);
    assert.deepEqual(messages, ['Admin:Home/geo is a folder.', 'Admin:Home/notes.txt is not .csv.']);
    input.dispose();
  }));

filesTest('a pick that lands after the input went away writes nothing', () => {
  let answer;
  return withFilesDialog(() => new Promise((resolve) => answer = resolve), async () => {
    const changes = [];
    const input = filesInput('Files', {onChanged: (v) => changes.push(v)});
    pickFromFiles(input);
    input.dispose();
    answer([new FileInfo('Admin:Home/late.csv')]);
    await settle();
    assert.deepEqual(changes, []);
    assert.equal(rows(input).length, 0);
  });
});

filesTest('a read still running when the input goes away is cancelled, not settled', async () => {
  const input = filesInput('Files');
  const paused = pausedFile('big.csv');
  input.add([paused.file]);
  await settle();
  assert.deepEqual(input.peek(), []);

  input.dispose();
  paused.resume();
  await settle();
  assert.deepEqual(input.peek(), [], 'nothing is written into a disposed component');
  assert.equal(rows(input).length, 1, 'and the rows are not rebuilt');
});

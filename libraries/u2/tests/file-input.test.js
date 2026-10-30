/* fileInput — the path/validity state machine and the local-file paths (plan WO-15.1/15.2).
   `grok.dapi.files`, `grok.shell` and `DG.FileInfo` come from tests/dg-stub.mjs; everything under
   test is the shipped module. The existence check debounces 250 ms, as the platform's does, so the
   waits below are real. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Entity, FileInfo} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {FileInput, fileInput} = await import('../src/dg/inputs/file-input.js');
const grok = await import('datagrok-api/grok');
const ui = await import('datagrok-api/ui');

const DEBOUNCE = 300;

function wait(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

/** Every test runs against a clean document, a clean stub, and must leave the scope count where
 * it was. */
function fileTest(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const errors = [];
    grok.dapi.files.exists = async () => false;
    grok.dapi.files.list = async () => [];
    grok.shell.error = (m) => errors.push(m);
    grok.shell.warning = (m) => errors.push(m);
    try {
      await body(errors);
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function type(input, text) {
  const field = input.root.querySelector('input.u2-file-path');
  field.value = text;
  fire(field, 'input');
}

function drop(input, files) {
  fire(input.root.querySelector('.u2-file'), 'drop', {dataTransfer: {files}});
}

test('qualify: the path the platform is asked about', () => {
  assert.equal(FileInput.qualify(''), null);
  assert.equal(FileInput.qualify('geo/cities.csv'), 'geo/cities.csv');
  assert.equal(FileInput.qualify('geo/'), 'geo', 'a trailing slash is dropped');
  assert.equal(FileInput.qualify('geo/cities.csv', 'Demo:Files'), 'Demo:Files/geo/cities.csv');
  assert.equal(FileInput.qualify('/geo', 'Demo:Files'), 'Demo:Files/geo', 'one separator, not two');
  assert.equal(FileInput.qualify('Demo:Files/geo', 'Demo:Files'), 'Demo:Files/geo',
    'a path that already names the connection is left alone');
  assert.equal(FileInput.qualify('/', 'Demo:Files'), 'Demo:Files/', 'the root of the connection');
});

fileTest('an existing path resolves to the file it names', async () => {
  const input = fileInput('Path');
  grok.dapi.files.exists = async () => true;
  grok.dapi.files.list = async () => [new FileInfo('geo/cities.csv'), new FileInfo('geo/other.csv')];

  type(input, 'geo/cities.csv');
  assert.equal(input.validity.peek(), 'File geo/cities.csv is loading.', 'typing invalidates at once');
  await wait(DEBOUNCE);
  assert.equal(input.validity.peek(), null);
  assert.equal(input.peek().fullPath, 'geo/cities.csv');
  input.dispose();
});

fileTest('a path that is not there says so, and keeps the value it had', async () => {
  const input = fileInput('Path');
  grok.dapi.files.exists = async () => true;
  grok.dapi.files.list = async () => [new FileInfo('geo/cities.csv')];
  type(input, 'geo/cities.csv');
  await wait(DEBOUNCE);

  grok.dapi.files.exists = async () => false;
  type(input, 'geo/missing.csv');
  await wait(DEBOUNCE);
  assert.equal(input.validity.peek(), 'File does not exist.');
  assert.equal(input.peek().fullPath, 'geo/cities.csv', 'the last resolved file stands');
  input.dispose();
});

fileTest('an answer to an abandoned path never overwrites a newer verdict', async () => {
  const input = fileInput('Path');
  const seen = [];
  grok.dapi.files.exists = async (path) => {
    seen.push(path);
    if (path === 'slow.csv') {
      await wait(200);
      return true;
    }
    return false;
  };
  grok.dapi.files.list = async () => [new FileInfo('slow.csv')];

  type(input, 'slow.csv');
  await wait(DEBOUNCE);              // the slow check is now in flight
  type(input, 'fast.csv');
  await wait(DEBOUNCE + 200);
  assert.deepEqual(seen, ['slow.csv', 'fast.csv']);
  assert.equal(input.validity.peek(), 'File does not exist.');
  assert.equal(input.peek(), null, 'the abandoned answer is dropped');
  input.dispose();
});

fileTest('the local-file picker sits on the shared options rail', async () => {
  const input = fileInput('Path');
  const rail = input.root.querySelector('.u2-input-options');
  assert.equal(rail.parentNode, input.box, 'so it crosses the platform bridge with the editor');
  assert.equal(rail.querySelectorAll('.u2-btn').length, 1, 'the folder-open icon, and only it');
  assert.equal(rail.querySelectorAll('input').length, 0, 'the hidden picker is not rail furniture');
  assert.equal(input.root.querySelector('.u2-input-editor').querySelectorAll('input[type=file]')
    .length, 1, 'it sits in the editor instead');
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

fileTest('where the platform has its Files dialog, its icon follows Open file on the rail', () =>
  withFilesDialog(async () => null, async () => {
    const input = fileInput('Path');
    const icons = [...input.root.querySelectorAll('.u2-input-options [data-u2="icon-button"]')];
    assert.deepEqual(icons.map((el) => el.getAttribute('aria-label')), ['Open file', 'Add file from Files'],
      'the platform order, its tooltips verbatim');
    assert.ok(icons[1].querySelector('.fa-folder-tree'));
    input.dispose();
  }));

fileTest('a file picked from Files is the value at once, its share path in the box', (errors) =>
  withFilesDialog(async () => new FileInfo('Admin:Home/merge.csv'), async () => {
    let asked = 0;
    grok.dapi.files.exists = async () => {
      asked++;
      return true;
    };
    const changes = [];
    const input = fileInput('Path', {onChanged: (v) => changes.push(v?.fullPath)});
    pickFromFiles(input);
    await flush();
    assert.equal(input.peek().fullPath, 'Admin:Home/merge.csv');
    assert.equal(input.root.querySelector('input.u2-file-path').value, 'Admin:Home/merge.csv');
    assert.equal(input.validity.peek(), null);
    assert.deepEqual(changes, ['Admin:Home/merge.csv']);
    await wait(DEBOUNCE);
    assert.equal(asked, 0, 'the dialog\'s file is not checked again');
    assert.deepEqual(errors, []);
    input.dispose();
  }));

fileTest('closing the Files dialog without a pick leaves the value', () =>
  withFilesDialog(async () => null, async () => {
    const input = fileInput('Path', {value: new FileInfo('Admin:Home/kept.csv')});
    pickFromFiles(input);
    await flush();
    assert.equal(input.peek().fullPath, 'Admin:Home/kept.csv');
    input.dispose();
  }));

fileTest('a folder picked from Files in file mode is judged as its typed path is', () =>
  withFilesDialog(async () => new FileInfo('Admin:Home/geo', {directory: true}), async () => {
    grok.dapi.files.exists = async () => true;
    grok.dapi.files.list = async () => [new FileInfo('Admin:Home/geo', {directory: true})];
    const input = fileInput('Path');
    pickFromFiles(input);
    await flush();
    assert.equal(input.root.querySelector('input.u2-file-path').value, 'Admin:Home/geo');
    await wait(DEBOUNCE);
    assert.equal(input.validity.peek(), 'File does not exist.');
    assert.equal(input.peek(), null);
    input.dispose();
  }));

fileTest('a pick that lands after the input went away writes nothing', () => {
  let answer;
  return withFilesDialog(() => new Promise((resolve) => answer = resolve), async () => {
    const input = fileInput('Path');
    pickFromFiles(input);
    input.dispose();
    answer(new FileInfo('Admin:Home/late.csv'));
    await flush();
    assert.equal(input.peek(), null);
  });
});

fileTest('the Files dialog is asked for what the mode wants, over the shares', () => {
  const asked = [];
  return withFilesDialog(async (options) => {
    asked.push(options);
    return null;
  }, async () => {
    const input = fileInput('Folder', {mode: 'folder'});
    pickFromFiles(input);
    await flush();
    assert.deepEqual(asked, [{connection: undefined, mode: 'folder'}]);
    input.dispose();
  });
});

const demoFiles = () => new Entity('Files', {nqName: 'Demo:Files'});

function chooseInConnection(input) {
  fire(input.root.querySelector('.u2-input-options [aria-label="Choose directory or file"]'), 'click');
}

fileTest('an input on one connection offers the platform\'s chooser over that share alone', () =>
  withFilesDialog(async () => null, async () => {
    const input = fileInput('Path', {connection: demoFiles()});
    const icons = [...input.root.querySelectorAll('.u2-input-options [data-u2="icon-button"]')];
    assert.deepEqual(icons.map((el) => el.getAttribute('aria-label')), ['Choose directory or file'],
      'the platform\'s tooltip verbatim, and no local picker');
    assert.ok(icons[0].querySelector('.fa-folder-tree'));
    input.dispose();
  }));

fileTest('a pick in the connection\'s chooser is the value, its path relative to the share', (errors) => {
  const connection = demoFiles();
  const asked = [];
  return withFilesDialog(async (options) => {
    asked.push(options);
    return new FileInfo('Demo:Files/geo/cities.csv');
  }, async () => {
    const changes = [];
    const input = fileInput('Path', {connection, mode: 'any', onChanged: (v) => changes.push(v?.fullPath)});
    chooseInConnection(input);
    await flush();
    assert.equal(asked.length, 1);
    assert.equal(asked[0].connection, connection, 'the dialog browses the input\'s own share');
    assert.equal(asked[0].mode, 'all', 'the platform\'s word for files and folders');
    assert.deepEqual(changes, ['Demo:Files/geo/cities.csv']);
    assert.equal(input.root.querySelector('input.u2-file-path').value, 'geo/cities.csv');
    assert.deepEqual(errors, []);
    input.dispose();
  });
});

fileTest('the connection\'s root picked in folder mode shows as /', () =>
  withFilesDialog(async () => new FileInfo('Demo:Files', {directory: true}), async () => {
    const changes = [];
    const input = fileInput('Folder', {connection: demoFiles(), mode: 'folder', onChanged: (v) => changes.push(v?.fullPath)});
    chooseInConnection(input);
    await flush();
    assert.deepEqual(changes, ['Demo:Files']);
    assert.equal(input.root.querySelector('input.u2-file-path').value, '/');
    input.dispose();
  }));

fileTest('an empty box is invalid only when the input is not nullable', async () => {
  const optional = fileInput('Path');
  assert.equal(optional.validity.peek(), null);
  optional.dispose();

  const required = fileInput('Path', {nullable: false});
  assert.equal(required.validity.peek(), 'Can\'t be empty.');
  required.dispose();
});

fileTest('folder mode asks for a folder, and says folder', async () => {
  const input = fileInput('Folder', {mode: 'folder'});
  assert.equal(input.root.querySelector('.u2-input-options'), null,
    'no local-file affordance in folder mode');
  grok.dapi.files.exists = async () => true;
  grok.dapi.files.list = async () => [new FileInfo('geo/cities.csv')];

  type(input, 'geo/cities.csv');
  await wait(DEBOUNCE);
  assert.equal(input.validity.peek(), 'Folder does not exist.');

  grok.dapi.files.list = async () => [new FileInfo('geo/maps', {directory: true})];
  type(input, 'geo/maps');
  await wait(DEBOUNCE);
  assert.equal(input.validity.peek(), null);
  assert.equal(input.peek().fullPath, 'geo/maps');
  input.dispose();
});

fileTest('a connection namespaces the path and hides the local picker', async () => {
  const asked = [];
  grok.dapi.files.exists = async (path) => {
    asked.push(path);
    return true;
  };
  grok.dapi.files.list = async () =>
    [new FileInfo('Demo:Files/geo/cities.csv')];
  const input = fileInput('Path', {connection: new Entity('Files', {nqName: 'Demo:Files'})});
  assert.equal(input.root.querySelector('.u2-input-options'), null);

  type(input, 'geo/cities.csv');
  await wait(DEBOUNCE);
  assert.deepEqual(asked, ['Demo:Files/geo/cities.csv']);
  assert.equal(input.validity.peek(), null);
  assert.equal(input.root.querySelector('input.u2-file-path').value, 'geo/cities.csv',
    'the box shows the path relative to the connection');
  input.dispose();
});

fileTest('a path the listing does not carry does not resolve', async () => {
  const input = fileInput('Path');
  grok.dapi.files.exists = async () => true;
  grok.dapi.files.list = async () => [];
  type(input, 'geo/cities.csv');
  await wait(DEBOUNCE);
  assert.equal(input.validity.peek(), 'File does not exist.',
    'a share that answers every call with nothing leaves no silent verdict');
  assert.equal(input.peek(), null);
  input.dispose();
});

fileTest('the root of a connection resolves with no FileInfo to show for it', async () => {
  grok.dapi.files.exists = async () => true;
  grok.dapi.files.list = async () => [];
  const input = fileInput('Path', {connection: new Entity('Files', {nqName: 'Demo:Files'})});
  type(input, '/');
  await wait(DEBOUNCE);
  assert.equal(input.validity.peek(), null, 'js-api cannot build a directory FileInfo for it');
  assert.equal(input.peek(), null);
  input.dispose();
});

fileTest('the extension filter is for file mode alone — a dotted folder name is not one', async () => {
  const asked = [];
  grok.dapi.files.exists = async () => true;
  grok.dapi.files.list = async (path, recursive, ext) => {
    asked.push(ext);
    return [new FileInfo('geo/v1.2', {directory: true})];
  };
  const input = fileInput('Path', {mode: 'folder'});
  type(input, 'geo/v1.2');
  await wait(DEBOUNCE);
  assert.deepEqual(asked, [null]);
  assert.equal(input.validity.peek(), null);
  assert.equal(input.peek().fullPath, 'geo/v1.2');
  input.dispose();
});

fileTest('a value written from outside drops the verdict on the path it replaces', async () => {
  const input = fileInput('Path');
  type(input, 'geo/missing.csv');
  await wait(DEBOUNCE);
  assert.equal(input.validity.peek(), 'File does not exist.');

  input.value = new FileInfo('geo/cities.csv');
  assert.equal(input.root.querySelector('input.u2-file-path').value, 'geo/cities.csv');
  assert.equal(input.validity.peek(), null, 'the platform setter clears and revalidates');

  type(input, 'geo/other.csv');
  input.value = new FileInfo('geo/cities.csv');
  await wait(DEBOUNCE);
  assert.equal(input.validity.peek(), null, 'a check in flight cannot re-assert a verdict either');
  input.dispose();
});

fileTest('a dropped file becomes the value; several are refused', async (errors) => {
  const input = fileInput('Path');
  drop(input, [new File(['a,b\n1,2'], 'local.csv'), new File(['x'], 'other.csv')]);
  await flush();
  assert.deepEqual(errors, ['Multiple files are not supported.']);
  assert.equal(input.peek(), null);

  drop(input, [new File(['a,b\n1,2'], 'local.csv')]);
  await flush();
  assert.equal(input.peek().fileName, 'local.csv');
  assert.equal(input.peek().data.length, 7);
  assert.equal(input.validity.peek(), null);
  input.dispose();
});

fileTest('a file dropped while a path check runs is what the input keeps', async () => {
  const input = fileInput('Path');
  grok.dapi.files.exists = async () => true;
  grok.dapi.files.list = async () => [new FileInfo('geo/cities.csv')];

  type(input, 'geo/cities.csv');
  drop(input, [new File(['x'], 'local.csv')]);
  await wait(DEBOUNCE);
  assert.equal(input.peek().fileName, 'local.csv', 'the check that was overtaken is ignored');
  assert.equal(input.validity.peek(), null);
  input.dispose();
});

fileTest('typing over a local file in flight abandons its read', async () => {
  const input = fileInput('Path');
  drop(input, [new File(['a,b\n1,2'], 'local.csv')]);
  type(input, 'geo/cities.csv');
  await wait(DEBOUNCE);
  assert.equal(input.peek(), null, 'the file that is no longer in the box never lands');
  assert.equal(input.validity.peek(), 'File does not exist.');
  input.dispose();
});

fileTest('a local read still running when the input goes away lands nowhere', async () => {
  const input = fileInput('Path');
  drop(input, [new File(['a,b\n1,2'], 'local.csv')]);
  input.dispose();
  await flush();
  assert.equal(input.peek(), null, 'nothing is written into a disposed input');
});

fileTest('emptying the box clears the file', async () => {
  const input = fileInput('Path');
  grok.dapi.files.exists = async () => true;
  grok.dapi.files.list = async () => [new FileInfo('geo/cities.csv')];
  type(input, 'geo/cities.csv');
  await wait(DEBOUNCE);
  assert.equal(input.peek().fullPath, 'geo/cities.csv');

  type(input, '');
  assert.equal(input.peek(), null, 'no path, no file');
  await wait(DEBOUNCE);
  assert.equal(input.peek(), null, 'and no late answer brings it back');
  input.dispose();
});

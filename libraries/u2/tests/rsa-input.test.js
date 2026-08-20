/* rsaInput — the PEM-vs-BASE64 sniff and the upload/masked swap (plan WO-16). `grok`, `DG` and
   `ui.makeDroppable` come from tests/dg-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';

register('./dg-stub.mjs', import.meta.url);
const {RsaInput, rsaInput} = await import('../src/dg/rsa-input.js');
const grok = await import('datagrok-api/grok');
const DG = await import('datagrok-api/dg');
const {drops} = await import('datagrok-api/ui');

const PEM = '-----BEGIN RSA PRIVATE KEY-----\nMIIBOgIBAAJBAK\n-----END RSA PRIVATE KEY-----';

function rsaTest(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const errors = [];
    grok.shell.error = (m) => errors.push(m);
    try {
      await body(errors);
    } finally {
      grok.dapi.files.readAsText = async () => '';
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** A key sitting in a file share, as the Files browser hands it to the drag channel. */
function shared(name) {
  return new DG.FileInfo(`keys/${name}`, null, {nqName: 'Home'});
}

/** Drops `payload` on the input's own registration, in the shape the running platform uses. */
function drop(payload) {
  return drops[drops.length - 1].doDrop(payload);
}

const upload = (input) => input.root.querySelector('.u2-rsa-upload');
const masked = (input) => input.root.querySelector('.u2-rsa-masked');

test('content: PEM text as it stands, anything else BASE64', () => {
  assert.equal(RsaInput.content(new TextEncoder().encode(`\n${PEM}\n  `)), PEM,
    'the header is looked for after trimming');
  const der = new Uint8Array([0x30, 0x82, 0x01, 0x3a, 0x02, 0x01, 0x00]);
  assert.equal(RsaInput.content(der), 'MIIBOgIBAA==');
  assert.equal(RsaInput.content(new Uint8Array(0)), '');
});

rsaTest('a key file fills the input and swaps it to the masked view', async () => {
  const input = rsaInput('Key');
  assert.equal(upload(input).hidden, false);
  assert.equal(masked(input).hidden, true);

  await input.accept(new File([PEM], 'server.pem'));
  assert.equal(input.value.peek(), PEM);
  assert.equal(upload(input).hidden, true);
  assert.equal(masked(input).hidden, false);
  assert.equal(masked(input).querySelector('input').value, RsaInput.PLACEHOLDER,
    'the key itself never reaches the DOM');
  input.dispose();
});

rsaTest('a binary key is kept as BASE64', async () => {
  const input = rsaInput('Key');
  await input.accept(new File([new Uint8Array([0, 1, 2, 3])], 'server.der'));
  assert.equal(input.value.peek(), 'AAECAw==');
  input.dispose();
});

rsaTest('a file the accept list does not name is refused', async (errors) => {
  const input = rsaInput('Key');
  await input.accept(new File([PEM], 'notes.txt'));
  assert.equal(input.value.peek(), null);
  assert.deepEqual(errors, ['Expected a .pem,.der file.']);

  const relaxed = rsaInput('Key', {accept: '.txt'});
  await relaxed.accept(new File([PEM], 'notes.txt'));
  assert.equal(relaxed.value.peek(), PEM);
  input.dispose();
  relaxed.dispose();
});

rsaTest('Change clears the key and brings the upload back', async () => {
  const input = rsaInput('Key');
  await input.accept(new File([PEM], 'server.pem'));

  fire(masked(input).querySelector('button.u2-btn'), 'click');
  assert.equal(input.value.peek(), null);
  assert.equal(upload(input).hidden, false);
  assert.equal(masked(input).hidden, true);
  input.dispose();
});

rsaTest('a required input is invalid until a key is in', async () => {
  const optional = rsaInput('Key');
  assert.equal(optional.validity.peek(), null);
  optional.dispose();

  const input = rsaInput('Key', {nullable: false});
  assert.equal(input.validity.peek(), 'Can\'t be empty.');
  await input.accept(new File([PEM], 'server.pem'));
  assert.equal(input.validity.peek(), null);
  input.dispose();
});

rsaTest('dropping a key file on the input stores it', async () => {
  const input = rsaInput('Key');
  fire(input.root.querySelector('.u2-rsa'), 'drop',
    {dataTransfer: {files: [new File([PEM], 'server.pem')]}});
  await flush();
  assert.equal(input.value.peek(), PEM);
  input.dispose();
});

rsaTest('several files at once are refused, as the single-file input does', async (errors) => {
  const input = rsaInput('Key');
  fire(input.root.querySelector('.u2-rsa'), 'drop',
    {dataTransfer: {files: [new File([PEM], 'a.pem'), new File([PEM], 'b.pem')]}});
  await flush();
  assert.deepEqual(errors, ['Multiple files are not supported.']);
  assert.equal(input.value.peek(), null);
  input.dispose();
});

rsaTest('a key dragged out of the Files browser is read and stored', async () => {
  grok.dapi.files.readAsText = async (f) => f === key ? PEM : '';
  const key = shared('server.pem');

  const input = rsaInput('Key');
  assert.equal(drops[drops.length - 1].acceptDrop(key), true);
  await drop({dragObject: key});
  assert.equal(input.value.peek(), PEM);
  assert.equal(masked(input).hidden, false);
  input.dispose();
});

rsaTest('the dragged object is taken however the platform hands it over', async () => {
  grok.dapi.files.readAsText = async () => PEM;

  const input = rsaInput('Key');
  // platform builds before the DragDropArgs wrapper pass the dragged object itself
  await drop(shared('server.pem'));
  assert.equal(input.value.peek(), PEM);
  input.dispose();
});

rsaTest('a drop that is not a key file of a share is left alone', async (errors) => {
  grok.dapi.files.readAsText = async () => PEM;
  const input = rsaInput('Key');

  for (const payload of [null, {}, 'Demog', DG.FileInfo.fromBytes('server.pem', null),
    shared('notes.txt')]) {
    assert.equal(drops[drops.length - 1].acceptDrop(payload), false);
    await drop({dragObject: payload});
  }
  assert.equal(input.value.peek(), null);
  assert.deepEqual(errors, [], 'somebody else\'s drag is not the input\'s to complain about');
  input.dispose();
});

rsaTest('a shared file that is not a key says so', async (errors) => {
  grok.dapi.files.readAsText = async () => 'ssh-rsa AAAAB3Nza';

  const input = rsaInput('Key');
  await drop({dragObject: shared('server.pem')});
  assert.equal(input.value.peek(), null);
  assert.deepEqual(errors, ['Invalid key file']);
  input.dispose();
});

/** A deployed js-api bundle passes the wrapped handle to the accessor (`table-info.ts:73`), so
 * `connection` throws for every FileInfo. Rejecting the drop over that would be a silent no-op. */
rsaTest('a key whose connection cannot be read is still taken', async () => {
  grok.dapi.files.readAsText = async () => PEM;
  const key = shared('server.pem');
  Object.defineProperty(key, 'connection', {get() {
    throw new RangeError('Invalid argument');
  }});

  const input = rsaInput('Key');
  assert.equal(drops[drops.length - 1].acceptDrop(key), true);
  await drop({dragObject: key});
  assert.equal(input.value.peek(), PEM);
  input.dispose();
});

rsaTest('a key file that cannot be read says so instead of nothing', async (errors) => {
  grok.dapi.files.readAsText = async () => {
    throw new Error('no such connection');
  };

  const input = rsaInput('Key');
  await drop({dragObject: shared('server.pem')});
  assert.equal(input.value.peek(), null);
  assert.deepEqual(errors, ['Could not read the key file.']);
  input.dispose();
});

rsaTest('a drop that lands after dispose is read and thrown away', async () => {
  grok.dapi.files.readAsText = async () => PEM;

  const input = rsaInput('Key');
  input.dispose();
  // the platform's drop registration has no counterpart to undo it (rsa-input.ts:_makeDroppable)
  await drop({dragObject: shared('server.pem')});
  assert.equal(input.value.peek(), null, 'nothing is written into a disposed input');
});

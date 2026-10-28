import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Dialog} from '../src/index.js';
import {Wizard} from '../src/components/containers/wizard.js';

const rowLike = () => new Proxy({}, {get: (_, name) => {
  throw new Error(`No such column: ${String(name)}`);
}});

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = Scope.onUnhandledError;
    const reported = [];
    Scope.onUnhandledError = (e) => reported.push(String(e));
    try {
      await body();
    } finally {
      Scope.onUnhandledError = saved;
      resetDom();
      await flush();
    }
    assert.deepEqual(reported, [], 'nothing reported');
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

smoke('a dialog whose OK answers a row-like proxy (a DG.Row) closes as OK', () => {
  const closed = [];
  const dialog = Dialog.create('Add row').onOK(() => rowLike()).onClose((reason) => closed.push(reason)).show();
  fire(dialog.getButton('OK'), 'click');
  assert.equal(dialog.isOpen.value, false);
  assert.deepEqual(closed, ['ok']);
});

smoke('a dialog whose OK answers a thenable that is no Promise waits for it', async () => {
  let settle;
  const dialog = Dialog.create('Save').onOK(() => ({then: (ok) => settle = ok})).show();
  fire(dialog.getButton('OK'), 'click');
  await flush();
  assert.equal(dialog.isOpen.value, true, 'waiting on the answer');
  settle(true);
  await flush();
  assert.equal(dialog.isOpen.value, false);
});

smoke('a wizard step whose commit answers a row-like proxy moves on; a finish answering one completes', () => {
  const w = new Wizard({
    steps: [
      {id: 'one', title: 'One', content: document.createElement('div'), commit: () => rowLike()},
      {id: 'two', title: 'Two', content: document.createElement('div')},
    ],
    onFinish: () => rowLike(),
  });
  document.body.append(w.root);
  w.next();
  assert.equal(w.currentStep.value, 'two');
  w.next();
  assert.equal(w.completed.value, true);
  assert.equal(w.busy.value, false);
  w.dispose();
});

smoke('a wizard commit answering a thenable that is no Promise waits; its rejection keeps the step', async () => {
  let fail;
  const w = new Wizard({steps: [
    {id: 'one', title: 'One', content: document.createElement('div'), commit: () => ({then: (_ok, no) => fail = no})},
    {id: 'two', title: 'Two', content: document.createElement('div')},
  ]});
  document.body.append(w.root);
  w.next();
  assert.equal(w.busy.value, true, 'waiting on the commit');
  fail(new Error('Server refused'));
  await flush();
  assert.equal(w.currentStep.value, 'one');
  assert.equal(w.error.value, 'Server refused');
  w.dispose();
});

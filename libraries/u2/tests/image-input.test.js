/* ImageInput: the value is the source, REMOVE clears it, and an empty value hides both. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {ImageInput} from '../src/components/image-input.js';

const LOGO = 'https://datagrok.ai/img/logo.svg';

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

function mount(component) {
  document.body.append(component.root);
  return component;
}

function image(input) {
  return input.root.querySelector('.u2-image-preview');
}

function remove(input) {
  return input.root.querySelector('.u2-image-remove');
}

smoke('the value drives the preview; REMOVE clears it', async () => {
  const changes = [];
  const input = mount(new ImageInput({label: 'Logo', value: LOGO, alt: 'Datagrok',
    onChanged: (v) => changes.push(v)}));
  assert.equal(input.root.dataset.u2, 'image-input');
  assert.equal(image(input).src, LOGO);
  assert.equal(image(input).alt, 'Datagrok');
  assert.equal(image(input).hidden, false);
  assert.equal(remove(input).hidden, false);

  fire(remove(input), 'click');
  assert.equal(input.value.value, '');
  assert.deepEqual(changes, ['']);
  assert.equal(image(input).src, '');
  assert.equal(image(input).hidden, true, 'nothing to preview, nothing shown');
  assert.equal(remove(input).hidden, true, 'and nothing to remove');
  input.dispose();
});

smoke('an empty input starts hidden and shows up on a programmatic write', async () => {
  const input = mount(new ImageInput({label: 'Logo'}));
  assert.equal(input.value.value, '');
  assert.equal(remove(input).hidden, true);

  input.value.value = LOGO;
  assert.equal(image(input).src, LOGO);
  assert.equal(image(input).hidden, false);
  assert.equal(remove(input).hidden, false);
  input.dispose();
});

smoke('validators come from the input base; dispose kills the listeners', async () => {
  const before = Scope.liveCount;
  const input = mount(new ImageInput({label: 'Logo'}));
  input.addValidator((v) => v === '' ? 'Pick an image' : null);
  assert.equal(input.validity.value, 'Pick an image');
  assert.equal(input.root.querySelector('.u2-input-error').textContent, 'Pick an image');

  input.value.value = LOGO;
  assert.equal(input.validity.value, null);

  input.dispose();
  assert.equal(Scope.liveCount, before);
  fire(remove(input), 'click');
  assert.equal(input.value.value, LOGO, 'listeners died with the scope');
});

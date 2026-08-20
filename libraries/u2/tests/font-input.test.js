/* FontInput: the strict format round-trip, the controls that compose it, the size preset popup,
   and the hidden option that keeps an unknown family alive. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {FontInput} from '../src/components/font-input.js';

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

function size(input) {
  return input.root.querySelector('.u2-font-size');
}

function family(input) {
  return input.root.querySelector('.u2-font-family');
}

function bold(input) {
  return input.root.querySelector('.u2-font-bold');
}

function italic(input) {
  return input.root.querySelector('.u2-font-italic');
}

function popup() {
  return document.body.querySelector('.u2-font-sizes');
}

smoke('parse/compose: the strict format, quoted families and the generic ones', async () => {
  assert.deepEqual(FontInput.parse('bold italic 14px "Roboto"'),
    {bold: true, italic: true, size: 14, family: 'Roboto'});
  assert.deepEqual(FontInput.parse('normal normal 12px monospace'),
    {bold: false, italic: false, size: 12, family: 'monospace'});
  assert.equal(FontInput.compose({bold: true, italic: true, size: 14, family: 'Roboto'}),
    'bold italic 14px "Roboto"');
  assert.equal(FontInput.compose({bold: false, italic: false, size: 12, family: 'Monospace'}),
    'normal normal 12px monospace', 'a generic family is written lowercase and unquoted');

  for (const bad of ['', '14px Roboto', 'heavy normal 14px "Roboto"', 'bold oblique 14px "Roboto"',
    'bold normal 14 "Roboto"', 'bold normal px "Roboto"', 'normal normal 12px fantasy'])
    assert.equal(FontInput.parse(bad), null, `"${bad}" is not the strict format`);
});

smoke('the value reaches the controls, and every control writes it back', async () => {
  const input = mount(new FontInput({label: 'Font', value: 'bold italic 14px "Roboto"'}));
  assert.equal(input.root.dataset.u2, 'font-input');
  assert.equal(size(input).value, '14');
  assert.equal(family(input).value, 'Roboto');
  assert.equal(bold(input).classList.contains('u2-font-toggle-active'), true);
  assert.equal(italic(input).getAttribute('aria-pressed'), 'true');

  fire(bold(input), 'click');
  assert.equal(input.value.value, 'normal italic 14px "Roboto"');
  assert.equal(bold(input).classList.contains('u2-font-toggle-active'), false);

  size(input).value = '24';
  fire(size(input), 'input');
  assert.equal(input.value.value, 'normal italic 24px "Roboto"');

  family(input).value = 'Arial';
  fire(family(input), 'change');
  assert.equal(input.value.value, 'normal italic 24px "Arial"');

  input.value.value = 'bold normal 9px "Open Sans"';
  assert.equal(size(input).value, '9', 'a programmatic write reaches the controls');
  assert.equal(family(input).value, 'Open Sans');
  assert.equal(bold(input).classList.contains('u2-font-toggle-active'), true);
  assert.equal(italic(input).classList.contains('u2-font-toggle-active'), false);
  input.dispose();
});

smoke('an unknown family joins the select as a hidden option and survives the round-trip',
  async () => {
    const input = mount(new FontInput({label: 'Font', value: 'normal normal 11px "Comic Sans MS"'}));
    const options = family(input).querySelectorAll('option');
    const unknown = options.find((o) => o.value === 'Comic Sans MS');
    assert.ok(unknown, 'the family the select does not offer is added');
    assert.equal(unknown.hidden, true);
    assert.equal(family(input).value, 'Comic Sans MS');

    fire(italic(input), 'click');
    assert.equal(input.value.value, 'normal italic 11px "Comic Sans MS"', 'and stays the family');
    input.dispose();
  });

smoke('unparseable text leaves the defaults; the next edit normalizes the value', async () => {
  const input = mount(new FontInput({label: 'Font', value: 'nonsense'}));
  assert.equal(size(input).value, '');
  assert.equal(family(input).value, 'Roboto');
  assert.equal(bold(input).classList.contains('u2-font-toggle-active'), false);

  fire(bold(input), 'click');
  assert.equal(input.value.value, 'bold normal 10px "Roboto"');
  input.dispose();
});

smoke('the size popup opens under the box, picks a size, and closes with it', async () => {
  const input = mount(new FontInput({label: 'Font', value: 'normal normal 12px "Roboto"',
    sizes: [10, 12, 14]}));
  assert.equal(popup(), null);

  fire(size(input), 'mousedown');
  await flush();
  assert.deepEqual(popup().querySelectorAll('.u2-font-size-option').map((o) => o.textContent),
    ['10', '12', '14']);

  fire(popup().querySelectorAll('.u2-font-size-option')[2], 'pointerdown');
  assert.equal(input.value.value, 'normal normal 14px "Roboto"');
  assert.equal(size(input).value, '14');
  assert.equal(popup(), null, 'picking closes the popup');

  fire(size(input), 'mousedown');
  await flush();
  assert.ok(popup());
  fire(size(input), 'blur');
  assert.equal(popup(), null, 'so does leaving the box');
  input.dispose();
});

smoke('clicking the open box closes the popup, and the next click reopens it', async () => {
  const input = mount(new FontInput({label: 'Font', sizes: [10, 12]}));
  fire(size(input), 'mousedown');
  await flush();
  assert.ok(popup());

  fire(size(input), 'mousedown');
  assert.equal(popup(), null, 'a click on the box it is already under toggles it shut');

  fire(size(input), 'mousedown');
  await flush();
  assert.ok(popup(), 'and the click after that brings it back');
  input.dispose();
});

smoke('dispose closes the popup and kills the listeners', async () => {
  const before = Scope.liveCount;
  const input = mount(new FontInput({label: 'Font'}));
  assert.equal(input.value.value, 'normal normal 10px "Roboto"', 'the default value is the format');
  fire(size(input), 'mousedown');
  await flush();
  const open = popup();

  input.dispose();
  assert.equal(open.isConnected, false);
  assert.equal(Scope.liveCount, before);

  fire(bold(input), 'click');
  assert.equal(input.value.value, 'normal normal 10px "Roboto"', 'listeners died with the scope');
});

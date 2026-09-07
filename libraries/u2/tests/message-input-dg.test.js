/* The platform mention layer: the smart-filter the users source builds, the exact markup token a
   chip serializes to, and the placeholder-then-resolve chip of a restored token. `grok.dapi` and
   `DG.ObjectHandler` come from tests/dg-stub.mjs — per its contract everything on `dapi` is a field
   a test replaces, and the handler registry gets its lookup installed here. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {User} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const DG = await import('datagrok-api/dg');
const {messageInput, userMentionProvider, USER_TOKEN} = await import('../src/dg/inputs/message-input.js');

/** The handler surface `HandlerRenderer` asks for; the caption doubles as the markup text so an
 * assertion can tell handler output from the plain-caption fallback. */
DG.ObjectHandler.forEntity = () => ({
  getCaption: (user) => user.friendlyName,
  renderMarkup: (user) => {
    const el = document.createElement('span');
    el.className = 'test-user-markup';
    el.textContent = user.friendlyName;
    return el;
  },
  renderListItem: (user) => {
    const el = document.createElement('div');
    el.className = 'test-user-row';
    el.textContent = user.friendlyName;
    return el;
  },
});

const ADA = new User('ada', {id: 'u-1', friendlyName: 'Ada Almeida'});
const listed = [];
let found = async () => ADA;

grok.dapi.users = {
  list: async (options) => {
    listed.push(options);
    return [ADA];
  },
  find: (id) => found(id),
};

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      listed.length = 0;
      found = async () => ADA;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

smoke('the users source builds the four-field like filter, and none for an empty query', async () => {
  const provider = userMentionProvider();
  await provider.fetch('ada', new AbortController().signal);
  assert.deepEqual(listed[0], {pageSize: 20,
    filter: 'firstName like "ada" or lastName like "ada" or login like "ada" or email like "ada"'});

  await provider.fetch('', new AbortController().signal);
  assert.deepEqual(listed[1], {pageSize: 20});

  await provider.fetch('a"b\\c', new AbortController().signal);
  assert.ok(!listed[2].filter.includes('"b'), `quotes are sanitized away: ${listed[2].filter}`);
});

smoke('toMarkup composes the platform token, collapsing whitespace in the name', async () => {
  const provider = userMentionProvider();
  assert.equal(provider.toMarkup(ADA), '<span>#{x.u-1."Ada Almeida"}</span>');
  assert.equal(provider.caption(ADA), 'Ada Almeida');

  const messy = new User('messy', {id: 'u-2', friendlyName: 'Ada\tde\nAlmeida\r'});
  assert.equal(provider.toMarkup(messy), '<span>#{x.u-2."Ada de Almeida "}</span>');

  const own = new User('own', {id: 'u-3', friendlyName: 'Ignored'});
  own.toMarkup = () => '<span>#{x.real."From the platform"}</span>';
  assert.equal(provider.toMarkup(own), '<span>#{x.real."From the platform"}</span>',
    'the platform\'s own toMarkup always wins');
});

smoke('USER_TOKEN matches what toMarkup emits and captures the id and the name', async () => {
  const match = new RegExp(USER_TOKEN).exec(`hi ${userMentionProvider().toMarkup(ADA)}!`);
  assert.ok(match);
  assert.equal(match[0], '<span>#{x.u-1."Ada Almeida"}</span>');
  assert.equal(match[1], 'u-1');
  assert.equal(match[2], 'Ada Almeida');
});

smoke('renderToken shows the name at once and swaps in the handler markup on resolve', async () => {
  const provider = userMentionProvider();
  const match = new RegExp(USER_TOKEN).exec(provider.toMarkup(ADA));
  const el = provider.renderToken(match);
  assert.equal(el.textContent, 'Ada Almeida');
  assert.equal(el.querySelector('.test-user-markup'), null, 'the name stands in until the entity lands');

  await flush();
  assert.ok(el.querySelector('.test-user-markup'), 'resolved content replaces the placeholder');
  assert.equal(el.textContent, 'Ada Almeida');
});

smoke('a failing lookup keeps the name text', async () => {
  found = async () => {
    throw new Error('offline');
  };
  const provider = userMentionProvider();
  const el = provider.renderToken(new RegExp(USER_TOKEN).exec(provider.toMarkup(ADA)));
  await flush();
  assert.equal(el.textContent, 'Ada Almeida');
  assert.equal(el.querySelector('.test-user-markup'), null);
});

smoke('messageInput: an inserted user round-trips through the value as a chip', async () => {
  const control = messageInput({placeholder: 'Message…'});
  document.body.append(control.root);
  control.insertMention(ADA);

  const chip = control.root.querySelector('.u2-msg-chip');
  assert.ok(chip);
  assert.equal(chip.dataset.token, '<span>#{x.u-1."Ada Almeida"}</span>');
  assert.ok(chip.querySelector('.test-user-markup'), 'the chip renders through the object handler');
  assert.equal(control.value, '<span>#{x.u-1."Ada Almeida"}</span> ');

  const markup = control.value;
  control.value = markup;
  await flush();
  const restored = control.root.querySelectorAll('.u2-msg-chip');
  assert.equal(restored.length, 1);
  assert.equal(restored[0].dataset.token, '<span>#{x.u-1."Ada Almeida"}</span>');
  assert.equal(control.value, markup, 'restore then serialize is byte-exact');
  control.dispose();
});

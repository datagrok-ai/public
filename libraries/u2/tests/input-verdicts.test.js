/* An input's verdict follows everything its validators read, from the moment they are added — a rule over another
   field ("To is before From") re-judges when that field changes, before its own value ever did. The required
   asterisk says `required: true` and nothing else: a choice that can never be empty (`nullable: false`) refuses an
   empty value without wearing the marker. Every test leaves the live-scope count where it found it. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {signal, Scope, Form, TextInput, NumberInput, ChoiceInput, MultiChoiceInput, ChipsInput} from '../src/index.js';
import {TabStrip} from '../src/components/containers/tabs.js';
import {EMPTY} from '../src/core/text.js';

function verdict(name, body) {
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

verdict('a validator over another input re-judges when that input changes, from the start', () => {
  const from = new NumberInput({label: 'From', value: 1});
  const to = new NumberInput({label: 'To', value: 5});
  from.addValidator((v) => v !== null && to.value !== null && v > to.value ? 'From is after To' : null);
  assert.equal(from.validity.peek(), null);
  to.value = 0;
  assert.equal(from.validity.peek(), 'From is after To', 'To moved under From: From is judged again');
  to.value = 9;
  assert.equal(from.validity.peek(), null);
  from.dispose();
  to.dispose();
});

verdict('a validator over a plain signal follows it; one added later follows its own reads too', () => {
  const taken = signal(['ada']);
  const login = new TextInput({label: 'Login', value: 'bob'});
  login.addValidator((v) => taken.value.includes(v) ? 'Taken' : null);
  taken.value = ['ada', 'bob'];
  assert.equal(login.validity.peek(), 'Taken');
  const banned = signal(false);
  login.addValidator(() => banned.value ? 'Banned' : null);
  taken.value = [];
  banned.value = true;
  assert.equal(login.validity.peek(), 'Banned');
  login.dispose();
});

verdict('adding a validator inside an effect does not make that effect follow the validator', () => {
  const other = signal(0);
  let runs = 0;
  const input = new TextInput({label: 'Name', value: 'x'});
  input.effect(() => {
    runs++;
    if (runs === 1)
      input.addValidator(() => other.value > 0 ? 'Too many' : null);
  });
  other.value = 1;
  assert.equal(runs, 1, 'the effect that added it did not re-run');
  assert.equal(input.validity.peek(), 'Too many');
  input.dispose();
});

verdict('a tab marked by its Form names the field it shows in error', () => {
  const form = new Form();
  const name = new TextInput({label: 'Name', required: true});
  form.add(name);
  const strip = new TabStrip({tabs: [{id: 'details', title: 'Details', content: form}]});
  document.body.append(strip.root);
  form.validate();
  assert.equal(strip.root.querySelector('[data-u2-item="details"]').title, `Name: ${EMPTY}`);
  strip.dispose();
  name.dispose();
  form.dispose();
});

verdict('the asterisk says required: true; nullable: false refuses an empty value without it', () => {
  const pick = new ChoiceInput({label: 'Merge by', items: ['Name', 'ID'], value: 'Name', nullable: false});
  assert.equal(pick.root.classList.contains('u2-input-required'), false, 'a choice that is never empty: no asterisk');
  assert.equal(pick.getState().required, false);
  const key = new TextInput({label: 'Key', nullable: false});
  assert.equal(key.root.classList.contains('u2-input-required'), false);
  assert.equal(key.validity.peek(), EMPTY, 'still refuses an empty value');
  const name = new TextInput({label: 'Name', required: true});
  assert.equal(name.root.classList.contains('u2-input-required'), true);
  assert.equal(name.getState().required, true);
  for (const input of [pick, key, name])
    input.dispose();
});

verdict('an empty list input says the empty verdict, not its list placeholder (emptyText)', () => {
  const inputs = [new MultiChoiceInput({label: 'Steps', items: ['a', 'b'], required: true}),
    new MultiChoiceInput({label: 'Steps', items: ['a'], nullable: false, emptyText: 'No steps'}),
    new ChipsInput({label: 'Species', items: ['Human', 'Mouse'], required: true})];
  for (const input of inputs) {
    assert.equal(input.validity.peek(), EMPTY);
    input.dispose();
  }
});

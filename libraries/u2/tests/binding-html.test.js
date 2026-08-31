/* The HTML live tier (UB-4): every prop a plain element declares — text, cls, href/src and the
   appearance group — binds through an engine-owned effect writing the DOM in place. The element
   the node built is the element it keeps: identity is asserted on every follow. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, htmlProps, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {APPEARANCE_CATEGORY} from '../src/spec/appearance.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function spec(name, body) {
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

function registry() {
  const reg = new Registry();
  registerAll(reg);
  return reg;
}

spec('bound text follows the source with element identity unchanged', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {s: 'Hello'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'h3', name: 'title', bind: {text: '$.s'}}]},
  }, ctx, reg);
  const el = instance.node('title');
  assert.equal(el.textContent, 'Hello');
  ctx.data.s.value = 'World';
  assert.equal(el.textContent, 'World', 'in place, synchronously — no re-render involved');
  assert.equal(instance.node('title'), el, 'the element the node built is the element it keeps');
  ctx.data.s.value = null;
  assert.equal(el.textContent, '', 'nullish clears the text');
  instance.dispose();
  ctx.data.s.value = 'after';
  assert.equal(el.textContent, '', 'the effect died with the instance');
});

spec('bound cls tracks exactly the current set, never accumulating', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {c: 'a b'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'span', name: 'chip', bind: {cls: '$.c'}},
  }, ctx, reg);
  const el = instance.node('chip');
  assert.equal(el.className, 'a b');
  ctx.data.c.value = 'c';
  assert.equal(el.className, 'c', 'the previous set was removed');
  ctx.data.c.value = 'a b';
  ctx.data.c.value = 'c';
  assert.equal(el.className, 'c', 'alternation never accumulates');
  ctx.data.c.value = '';
  assert.equal(el.className, '', 'empty clears everything the binding added');
  instance.dispose();
});

spec('a bind supersedes a literal on the same key, uniformly — text and cls alike (F-1)', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {t: 'Bound', c: 'a b'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'h3', name: 'title', props: {text: 'Literal'}, bind: {text: '$.t'}},
      {tag: 'span', name: 'chip', props: {cls: 'keepme'}, bind: {cls: '$.c'}},
    ]},
  }, ctx, reg);
  assert.equal(instance.node('title').textContent, 'Bound', 'the bound text wins over the literal');
  const chip = instance.node('chip');
  assert.equal(chip.className, 'a b', 'the literal cls is replaced, not merged');
  ctx.data.c.value = 'c';
  assert.equal(chip.className, 'c', 'and alternation still tracks exactly the bound set');
  instance.dispose();
});

spec('bound href sets the attribute and a nullish value removes it', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {link: 'https://datagrok.ai'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'a', name: 'anchor', props: {text: 'Home'}, bind: {href: '$.link'}},
  }, ctx, reg);
  const el = instance.node('anchor');
  assert.equal(el.getAttribute('href'), 'https://datagrok.ai');
  ctx.data.link.value = null;
  assert.equal(el.getAttribute('href'), null, 'nullish removes the attribute');
  ctx.data.link.value = '/help';
  assert.equal(el.getAttribute('href'), '/help');
  instance.dispose();
});

spec('a bound appearance prop restyles a plain element live', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {tone: 'red'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', name: 'card', props: {padding: '8px'}, bind: {color: '$.tone'}},
  }, ctx, reg);
  const el = instance.node('card');
  assert.equal(el.style.color, 'red');
  assert.equal(el.style.padding, '8px', 'literal appearance still applies alongside');
  ctx.data.tone.value = 'green';
  assert.equal(el.style.color, 'green', 'the mount-scope effect follows the source');
  assert.equal(instance.node('card'), el, 'identity stable');
  ctx.data.tone.value = null;
  assert.equal(el.style.color, '', 'nullish clears the style');
  instance.dispose();
});

spec('a bind on an undeclared HTML prop is the standard placeholder', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'span', bind: {ghost: '$.s'}}]},
  }, new SpecContext({data: {s: 'x'}}), reg);
  const errors = instance.root.querySelectorAll('.u2-spec-error');
  assert.equal(errors.length, 1);
  assert.match(errors[0].textContent, /span: has no prop "ghost" to bind/);
  instance.dispose();
});

spec('src on a non-img, href on a non-a: undeclared, refused', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'h3', bind: {href: '$.s'}}]},
  }, new SpecContext({data: {s: 'x'}}), reg);
  assert.match(instance.root.querySelectorAll('.u2-spec-error')[0].textContent,
    /has no prop "href" to bind/);
  instance.dispose();
});

spec('a forward-ref HTML bind resolves after the render pass', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'h3', name: 'title', bind: {text: '$.who'}},
      {tag: 'u2-text-input', name: 'who', props: {value: 'world'}},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  assert.equal(instance.node('title').textContent, 'world');
  const editor = instance.node('who').root.querySelector('input');
  editor.value = 'typed';
  fire(editor, 'input');
  assert.equal(instance.node('title').textContent, 'typed', 'and follows live after the flush');
  instance.dispose();
});

spec('a forward-ref bind under an adopting root wires in place — no re-render loop', async () => {
  const reg = registry();
  const counts = [];
  reg.register({
    tag: 'u2-fake-target',
    description: 'Counted bind-source target, capped so a runaway loop fails instead of hanging',
    props: [{name: 'value', type: 'string', bindable: true, twoWay: true}],
    create: (props) => {
      counts.push({...props});
      if (counts.length > 20)
        throw new Error('runaway rebuild');
      return new TextInput({value: typeof props.value === 'string' ? props.value : undefined});
    },
    example: {tag: 'u2-fake-target'},
  });
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'h3', name: 'title', bind: {text: '$.who'}},
      {tag: 'u2-fake-target', name: 'who', props: {value: 'world'}},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0, 'built, never broken');
  assert.equal(counts.length, 1, 'the flush wired the element in place — nothing re-rendered');
  const el = instance.node('title');
  assert.equal(el.textContent, 'world', 'the value resolved');

  instance.node('who').value.value = 'typed';
  await flush();
  assert.equal(el.textContent, 'typed', 'and follows the live source afterwards');
  assert.equal(instance.node('title'), el, 'the element the node built is the element it keeps');
  assert.equal(counts.length, 1, 'a source change never re-renders anything');
  instance.dispose();
});

spec('htmlProps advertises the live tier: own props and the appearance group bindable', () => {
  const props = htmlProps('a');
  for (const name of ['text', 'cls', 'href'])
    assert.equal(props.find((p) => p.name === name)?.bindable, true, `${name} is live-bindable`);
  assert.equal(props.find((p) => p.name === 'src'), undefined, 'src stays img-only');
  assert.equal(htmlProps('img').find((p) => p.name === 'src')?.bindable, true);
  for (const p of props.filter((p) => p.category === APPEARANCE_CATEGORY))
    assert.equal(p.bindable, true, `${p.name} keeps its shared-group bindable flag`);
});

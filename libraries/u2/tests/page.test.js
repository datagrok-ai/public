import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush} from './dom-shim.js';
import {layout} from './layout.mjs';
import {Form, LabelColumn, Page, TextInput, button, iconButton, signal, span, statusLine} from '../src/index.js';
import {u2test, mount, byKind, byName, state, expectOperable, user} from '../src/testing/index.js';

/** A captioned text input whose caption measures `width`. */
function field(label, width) {
  const input = new TextInput({label, name: label});
  layout(input.root.querySelector('[data-u2-part="label"]'), {offsetWidth: width});
  return input;
}

const part = (page, name) => [...page.root.children].find((c) => c.dataset.u2Part === name) ?? null;

test('page: a header, a body that scrolls alone, a pinned action row; sections with a heading and actions',
  u2test(() => {
    let saved = 0;
    let refreshed = 0;
    const page = mount(() => new Page({name: 'template-page', header: span('Hit Design', 'my-header'),
      sections: [
        {title: 'Continue campaign', actions: [iconButton('sync', () => refreshed++, {tooltip: 'Refresh'})],
          content: span('the campaigns'), name: 'campaigns'},
        {title: 'New campaign', content: span('the picker'), maxWidth: 720},
      ],
      footer: [statusLine(signal(null), 'Creates KHL-3'), button('Cancel', () => {}), button('Create', () => saved++,
        {primary: true})]}));
    assert.deepEqual(byKind(document, 'page'), [page.root]);
    assert.equal(state(page.root).name, 'template-page');
    assert.deepEqual([...page.root.children].map((c) => c.dataset.u2Part), ['header', 'body', 'footer'],
      'the body between a header and a footer, neither of which scrolls with it');
    assert.equal(page.header.textContent, 'Hit Design');
    const body = part(page, 'body');
    assert.equal(body, page.body);
    const sections = [...body.children].filter((c) => c.dataset.u2Part === 'section');
    assert.deepEqual(sections.map((s) => s.querySelector('[data-u2-part="section-title"]').textContent),
      ['Continue campaign', 'New campaign']);
    assert.equal(sections[0].dataset.u2Name, 'campaigns');
    assert.equal(sections[1].querySelector('[data-u2-part="section-content"]').style.maxWidth, '720px');
    user.click(byName(page, 'Refresh'));
    assert.equal(refreshed, 1, 'the heading\'s actions act');
    user.click(byName(page.footer, 'Create'));
    assert.equal(saved, 1);
    assert.deepEqual([...page.footer.children].map((c) => c.textContent.trim()), ['Creates KHL-3', 'Cancel', 'Create'],
      'the line first, then the buttons: the row lays the buttons out on the right');
    expectOperable(page.root);
  }));

test('page: maxWidth caps the body and lets its forms go as wide; a form keeps its own cap', u2test(() => {
  const page = mount(() => new Page({maxWidth: 820, content: new Form()}));
  assert.equal(page.root.style.getPropertyValue('--u2-page-max-width'), '820px');
  const capped = mount(() => new Form({maxWidth: 480}));
  assert.equal(capped.root.style.getPropertyValue('--u2-form-max-width'), '480px');
  const compact = mount(() => new Form({density: 'compact'}));
  assert.equal(compact.root.dataset.dgDensity, 'compact', 'the compact token set, opted into');
}));

test('a wide form keeps its maxWidth and its page\'s; with neither it fills the row',
  {skip: typeof CSSStyleSheet !== 'function' && 'needs a DOM with the CSS cascade'}, u2test(() => {
    const capped = mount(() => new Form({layout: 'wide', maxWidth: 640}));
    assert.equal(getComputedStyle(capped.root).maxWidth, '640px');
    const page = mount(() => new Page({maxWidth: 820, content: new Form({layout: 'wide'})}));
    assert.equal(getComputedStyle(page.root.querySelector('.u2-form')).maxWidth, '820px');
    const free = mount(() => new Form({layout: 'wide'}));
    assert.equal(getComputedStyle(free.root).maxWidth, 'none');
  }));

test('page: its forms share one label column, unless a form keeps its own', u2test(async () => {
  const top = new Form({layout: 'normal'}).add(field('Name', 60));
  const bottom = new Form({layout: 'normal'}).add(field('Natural analog', 97));
  const own = new Form({layout: 'normal', labelColumn: false}).add(field('Id', 20));
  const page = mount(() => new Page({sections: [{title: 'Identity', content: top},
    {title: 'Details', content: [bottom, own]}]}));
  await flush();
  assert.ok(page.labelColumn instanceof LabelColumn);
  assert.equal(page.labelColumn.width.value, 97);
  assert.equal(top.root.style.getPropertyValue('--u2-form-label-width'), '97px', 'the narrower form follows');
  assert.equal(bottom.root.style.getPropertyValue('--u2-form-label-width'), '97px');
  assert.equal(own.root.style.getPropertyValue('--u2-form-label-width'), '20px', 'labelColumn: false keeps its own');
}));

test('page: labelColumn false leaves each form its own; a footer asked for is made on first ask', u2test(async () => {
  const a = new Form({layout: 'normal'}).add(field('Name', 60));
  const b = new Form({layout: 'normal'}).add(field('Natural analog', 97));
  const page = mount(() => new Page({content: [a, b], labelColumn: false}));
  await flush();
  assert.equal(page.labelColumn, undefined);
  assert.equal(a.root.style.getPropertyValue('--u2-form-label-width'), '60px');
  assert.equal(part(page, 'footer'), null, 'no footer until one is asked for');
  a.addButtons((row) => row.append(button('Save', () => {})), {host: page.footer});
  assert.equal(part(page, 'footer'), page.footer);
  assert.notEqual(byName(page.footer, 'Save'), null, 'a form\'s button row placed in the pinned row');
}));

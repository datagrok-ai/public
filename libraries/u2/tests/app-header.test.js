/* AppHeader: the icon in its three forms, items as text / element / link / markdown bullets, the
   titled variant, Learn more as a URL or an action, and the state block bdd reads. Same contract as
   the other UI suites — every test leaves the live-scope count where it found it. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {layout} from './layout.mjs';
import {signal, Scope, Control, span} from '../src/index.js';
import {AppHeader} from '../src/components/display/app-header.js';

function ui(name, body) {
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

const part = (header, name) => header.root.querySelector(`[data-u2-part="${name}"]`);
const items = (header) => [...header.root.querySelectorAll('[data-u2-part="item"]')];
const texts = (header) => items(header).map((li) => li.textContent);

ui('an image icon, the bullets in order and an external Learn more', () => {
  const header = new AppHeader({icon: '/images/hit-design.png',
    items: ['Sketch molecules', 'Calculate properties'], learnMore: 'https://example.com/README.md'});
  document.body.append(header.root);
  assert.equal(header.root.dataset.u2, 'app-header');
  assert.equal(header.root.classList.contains('u2-app-header'), true);
  assert.equal(header.root.style.getPropertyValue('--u2-app-header-icon-size'), '150px', 'the default icon size');
  const icon = part(header, 'icon');
  assert.equal(icon.tagName, 'DIV', 'a background, as the platform header drew it: a shrunk icon resamples the same, ' +
    'an image that fails to load leaves its square empty');
  assert.equal(icon.style.backgroundImage, 'url("/images/hit-design.png")');
  assert.equal(header.root.querySelector('ul').classList.contains('u2-app-header-items'), true);
  assert.deepEqual(texts(header), ['Sketch molecules', 'Calculate properties']);
  const more = part(header, 'learn-more');
  assert.equal(more.tagName, 'A');
  assert.equal(more.href, 'https://example.com/README.md');
  assert.equal(more.target, '_blank');
  assert.equal(more.textContent, 'Learn more');
  assert.ok(more.querySelector('.u2-app-header-external.fa-external-link-square'), 'a URL carries the external mark');
  assert.equal(part(header, 'title'), null);

  const state = header.getState();
  assert.equal(state.kind, 'app-header');
  assert.equal(state.label, undefined, 'no title, no label');
  assert.deepEqual(state.items.list.map((i) => [i.key, i.label]), [['0', 'Sketch molecules'], ['1', 'Calculate properties']]);
  assert.equal(state.items.count, 2);
  assert.deepEqual(state.extra, {learnMore: 'Learn more', learnMoreHref: 'https://example.com/README.md'});
  header.dispose();
});

ui('description: one item per markdown bullet, after the items, links made links', () => {
  const header = new AppHeader({icon: 'x.png', items: ['First'],
    description: '-  Configure your workflow\n- Read [the guide](https://example.com/guide) first\n\n* Starred\nPlain line\n '});
  assert.deepEqual(texts(header), ['First', 'Configure your workflow', 'Read the guide first', 'Starred', 'Plain line']);
  const a = items(header)[2].querySelector('a');
  assert.equal(a.textContent, 'the guide');
  assert.equal(a.href, 'https://example.com/guide');
  assert.equal(part(header, 'learn-more'), null, 'no learnMore, no link');
  assert.equal(header.getState().extra, undefined);
  header.dispose();
});

ui('an item is text, an element, a control, a link spec or a live text', async () => {
  const clicks = [];
  const live = signal('Rows: 1');
  const control = Control.build(() => span('a control'));
  const element = document.createElement('b');
  element.textContent = 'an element';
  const header = new AppHeader({icon: 'x.png', items: [
    'plain',
    element,
    control,
    {text: 'Docs', href: 'https://example.com/docs'},
    {text: 'Open help', onClick: () => clicks.push('help')},
    live,
  ]});
  assert.deepEqual(texts(header), ['plain', 'an element', 'a control', 'Docs', 'Open help', 'Rows: 1']);
  assert.equal(items(header)[1].firstChild, element, 'an element is placed as it is');
  const [docs, help] = [items(header)[3].querySelector('a'), items(header)[4].querySelector('a')];
  assert.equal(docs.target, '_blank');
  assert.equal(docs.querySelector('.u2-app-header-external'), null, 'only Learn more carries the mark');
  fire(help, 'click');
  assert.deepEqual(clicks, ['help']);
  live.value = 'Rows: 2';
  await flush();
  assert.equal(texts(header)[5], 'Rows: 2', 'a signal item follows its signal');
  header.dispose();
  control.dispose();
  fire(help, 'click');
  assert.deepEqual(clicks, ['help'], 'disposed with the header');
});

ui('a title: a 90px icon, the name and subtitle beside it, the title as the label', () => {
  const header = new AppHeader({icon: 'moltrack.png', title: 'MolTrack', subtitle: 'Track chemical data',
    description: '- Register compounds', bottomLine: true, name: 'moltrack-home'});
  assert.equal(header.root.classList.contains('u2-app-header-named'), true);
  assert.equal(header.root.classList.contains('u2-app-header-rule'), true);
  assert.equal(header.root.style.getPropertyValue('--u2-app-header-icon-size'), '90px');
  assert.equal(part(header, 'title').textContent, 'MolTrack');
  assert.equal(part(header, 'subtitle').textContent, 'Track chemical data');
  assert.equal(part(header, 'title').parentElement, part(header, 'subtitle').parentElement);
  assert.equal(header.root.dataset.u2Name, 'moltrack-home');
  const state = header.getState();
  assert.equal(state.label, 'MolTrack');
  assert.equal(state.name, 'moltrack-home');
  header.dispose();

  const sized = new AppHeader({icon: 'case.png', title: 'Case', iconSize: 120});
  assert.equal(sized.root.style.getPropertyValue('--u2-app-header-icon-size'), '120px', 'iconSize wins over the title default');
  assert.equal(sized.root.classList.contains('u2-app-header-rule'), false);
  sized.dispose();
});

ui('Learn more as an action: its own text, no external mark, reported by its text and no page', () => {
  let opened = 0;
  const header = new AppHeader({icon: 'flask', learnMore: {text: 'How it works', onClick: () => opened++}});
  const icon = part(header, 'icon');
  assert.equal(icon.tagName, 'I', 'a bare name is a Font Awesome icon');
  assert.equal(icon.classList.contains('fa-flask'), true);
  const more = part(header, 'learn-more');
  assert.equal(more.textContent, 'How it works');
  assert.equal(more.querySelector('.u2-app-header-external'), null);
  fire(more, 'click');
  assert.equal(opened, 1);
  assert.deepEqual(header.getState().extra, {learnMore: 'How it works'});
  assert.equal(header.getState().items.count, 0);
  assert.equal(header.getItems().truncated, false);
  header.dispose();
});

ui('an element icon is used as it is', () => {
  const svg = document.createElement('span');
  const header = new AppHeader({icon: svg, iconSize: 64});
  assert.equal(part(header, 'icon'), svg);
  assert.equal(svg.classList.contains('u2-app-header-icon'), true);
  assert.equal(header.root.style.getPropertyValue('--u2-app-header-icon-size'), '64px');
  header.dispose();
});

ui('its own width lays it out, whatever the host: under 480px a smaller icon, under 400px the bullets go below', async () => {
  const header = new AppHeader({icon: 'atom', title: 'MolTrack', items: ['Register compounds']});
  document.body.append(header.root);
  const observer = ResizeObserver.instances.find((o) => o.targets.includes(header.root) && !o.disconnected);
  assert.ok(observer, 'the header watches its own width');
  const look = () => ['u2-app-header-compact', 'u2-app-header-narrow'].filter((c) => header.root.classList.contains(c));
  const resize = async (width) => {
    layout(header.root, {clientWidth: width});
    observer.callback([]);
    assert.deepEqual(look(), before, 'nothing changes inside the observer: other observers would loop');
    await flush();
  };
  let before = [];
  await resize(1000);
  assert.deepEqual(look(), []);
  // the platform landings at 1366 px: Chem MPO 569 px, the Hit Triage family 578 px — the old header's look
  await resize(569);
  assert.deepEqual(look(), []);
  await resize(440);
  assert.deepEqual(look(), ['u2-app-header-compact']);
  before = look();
  await resize(320);
  assert.deepEqual(look(), ['u2-app-header-compact', 'u2-app-header-narrow']);
  before = look();
  await resize(0);
  assert.equal(header.root.classList.contains('u2-app-header-narrow'), true, 'hidden: the last layout stays');
  header.dispose();
  assert.equal(observer.disconnected, true);
});

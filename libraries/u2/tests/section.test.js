/* Section: one collapsible header + body — toggle by click/keyboard/signal, content kept when
   collapsed, the collapsible: false inert variant. Every test leaves the live-scope count where
   it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, Scope} from '../src/index.js';
import {Section} from '../src/components/containers/section.js';

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

function el(text) {
  const node = document.createElement('div');
  node.textContent = text;
  return node;
}

ui('structure: header with chevron and title, body, role/tabindex/aria on the header', () => {
  const s = new Section({title: 'Study'});
  assert.equal(s.root.dataset.u2, 'section');
  assert.equal(s.root.classList.contains('u2-section'), true);
  assert.equal(s.header.classList.contains('u2-section-header'), true);
  assert.equal(s.body.classList.contains('u2-section-body'), true);
  assert.deepEqual([...s.root.children], [s.header, s.body]);
  assert.notEqual(s.header.querySelector('.u2-section-chevron'), null);
  assert.equal(s.header.querySelector('.u2-section-chevron').getAttribute('aria-hidden'), 'true',
    'the chevron never reaches the accessible name');
  assert.equal(s.header.querySelector('.u2-section-title').textContent, 'Study');
  assert.equal(s.header.textContent, 'Study', 'the chevron adds no text');
  assert.equal(s.header.getAttribute('role'), 'button');
  assert.equal(s.header.tabIndex, 0);
  assert.equal(s.header.getAttribute('aria-expanded'), 'true');
  assert.notEqual(s.body.id, '', 'the body carries an id');
  assert.equal(s.header.getAttribute('aria-controls'), s.body.id);
  assert.equal(s.body.getAttribute('role'), 'region');
  assert.equal(s.body.getAttribute('aria-labelledby'), s.header.id);
  assert.equal(s.expanded.value, true);
  assert.notEqual(s.body.style.display, 'none', 'expanded by default');
  s.dispose();
});

ui('click toggles: body hides but content is kept, aria follows, and a second click restores', () => {
  const s = new Section({title: 'Study'});
  const content = el('kept');
  s.add(content);
  fire(s.header, 'click');
  assert.equal(s.expanded.value, false);
  assert.equal(s.header.getAttribute('aria-expanded'), 'false');
  assert.equal(s.body.style.display, 'none');
  assert.equal(content.parentNode, s.body, 'collapsed content is hidden, never detached');

  fire(s.header, 'click');
  assert.equal(s.expanded.value, true);
  assert.notEqual(s.body.style.display, 'none');
  s.dispose();
});

ui('a click on the chevron or the title bubbles to the header and toggles', () => {
  const s = new Section({title: 'Study'});
  fire(s.header.querySelector('.u2-section-chevron'), 'click');
  assert.equal(s.expanded.value, false);
  fire(s.header.querySelector('.u2-section-title'), 'click');
  assert.equal(s.expanded.value, true);
  s.dispose();
});

ui('Enter and Space toggle; other keys do nothing', () => {
  const s = new Section({title: 'Study'});
  fire(s.header, 'keydown', {key: 'Enter'});
  assert.equal(s.expanded.value, false);
  fire(s.header, 'keydown', {key: ' '});
  assert.equal(s.expanded.value, true);
  fire(s.header, 'keydown', {key: 'a'});
  assert.equal(s.expanded.value, true);
  s.dispose();
});

ui('programmatic writes and construction-time expanded: false drive the body', () => {
  const s = new Section({title: 'Study'});
  s.expanded.value = false;
  assert.equal(s.body.style.display, 'none');
  s.expanded.value = true;
  assert.notEqual(s.body.style.display, 'none');
  s.dispose();

  const closed = new Section({title: 'Closed', expanded: false});
  assert.equal(closed.body.style.display, 'none');
  assert.equal(closed.header.getAttribute('aria-expanded'), 'false');
  closed.dispose();
});

ui('a Signal handed as expanded is adopted, and writes flow both ways', () => {
  const outside = signal(true);
  const s = new Section({title: 'Study', expanded: outside});
  assert.equal(s.expanded, outside, 'the signal is adopted, not copied');

  outside.value = false;
  assert.equal(s.body.style.display, 'none', 'an outside write collapses');

  fire(s.header, 'click');
  assert.equal(outside.value, true, 'a click writes the outside signal');
  s.dispose();
});

ui('collapsible: false is a plain heading — no chevron, no role, clicks and writes inert', () => {
  const s = new Section({title: 'Plain', collapsible: false});
  assert.equal(s.header.querySelector('.u2-section-chevron'), null);
  assert.equal(s.header.getAttribute('role'), null);
  assert.equal(s.header.tabIndex, -1);
  assert.equal(s.header.getAttribute('aria-expanded'), null);

  fire(s.header, 'click');
  assert.equal(s.expanded.value, true, 'a click changes nothing');
  s.expanded.value = false;
  assert.notEqual(s.body.style.display, 'none', 'the display effect is not wired');
  s.dispose();
});

ui('add appends controls and elements into the body and chains', () => {
  const s = new Section({title: 'Study'});
  const inner = new Section({title: 'Inner'});
  const plain = el('row');
  assert.equal(s.add(inner, plain), s);
  assert.deepEqual([...s.body.children], [inner.root, plain]);
  inner.dispose();
  s.dispose();
});

ui('a signal title follows live and stops on dispose', () => {
  const title = signal('Draft');
  const s = new Section({title});
  const label = s.header.querySelector('.u2-section-title');
  assert.equal(label.textContent, 'Draft');
  title.value = 'Final';
  assert.equal(label.textContent, 'Final');
  s.dispose();
  title.value = 'Gone';
  assert.equal(label.textContent, 'Final', 'the binding died with the section');
});

ui('dispose drops the header listeners', () => {
  const s = new Section({title: 'Study'});
  s.dispose();
  fire(s.header, 'click');
  assert.equal(s.expanded.value, true, 'the listeners died with the section');
});

// the shim loads no stylesheets, so the skin rules the audit hangs on are pinned at the source
test('section.css keeps the collapsed chevron persistent (hover-only is expanded-only)', () => {
  const css = readFileSync(new URL('../css/section.css', import.meta.url), 'utf8');
  assert.match(css, /\.u2-section-header\[aria-expanded='false'\] \.u2-section-chevron/,
    'a collapsed section must keep its chevron visible');
  assert.match(css, /\.u2-section-header\[hidden\]/,
    'the display: flex header needs its own [hidden] rule');
});

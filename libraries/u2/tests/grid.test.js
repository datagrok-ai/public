/* VirtualGrid on the DOM shim: the shim lays nothing out, so each test gives the root a viewport
   (clientWidth → columns, clientHeight → rows) and fires `scroll` to re-render. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, signal} from '../src/index.js';
import {VirtualGrid} from '../src/components/collections/grid.js';

function grid(name, body) {
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

const NUMBERS = Array.from({length: 1000}, (_, i) => i);

/** 10 columns × 5 rows of 20px cells; the grid is told its size the way a layout would. */
function make(options = {}, width = 200, height = 100) {
  const g = new VirtualGrid({cellWidth: 20, cellHeight: 20, render: (n) => {
    const s = document.createElement('span');
    s.textContent = String(n);
    return s;
  }, ...options});
  document.body.append(g.root);
  g.root.clientWidth = width;
  g.root.clientHeight = height;
  fire(g.root, 'scroll');
  return g;
}

function cells(g) {
  return g.root.querySelectorAll('.u2-grid-cell');
}

function rendered(g) {
  return cells(g).map((c) => Number(c.dataset.index)).sort((a, b) => a - b);
}

grid('columns follow the width; the visible window plus overscan is what gets rendered', () => {
  const g = make();
  g.setItems(NUMBERS);
  fire(g.root, 'scroll');
  assert.equal(g.columns.value, 10);
  assert.equal(g.root.getAttribute('role'), 'listbox');
  // 5 visible rows + 2 overscan rows below (none above at the top), 10 per row
  assert.equal(g.renderedCount, 70);
  assert.deepEqual(rendered(g).slice(0, 3), [0, 1, 2]);
  const last = cells(g).find((c) => c.dataset.index === '69');
  assert.equal(last.style.left, '180px');
  assert.equal(last.style.top, '120px');
  assert.equal(g.root.querySelector('.u2-grid-content').style.height, '2000px');
  g.dispose();
});

grid('scrolling moves the window and recycles cells through the pool', () => {
  const g = make();
  g.setItems(NUMBERS);
  fire(g.root, 'scroll');
  g.root.scrollTop = 400;                       // row 20
  fire(g.root, 'scroll');
  const shown = rendered(g);
  assert.equal(shown[0], 180, 'two overscan rows above row 20');
  assert.equal(shown[shown.length - 1], 269, 'two overscan rows below the viewport');
  assert.equal(g.renderedCount, 90);
  assert.equal(cells(g).find((c) => c.dataset.index === '200').textContent, '200');
  g.dispose();
});

grid('a narrower width re-flows the same items into more rows', () => {
  const g = make();
  g.setItems(NUMBERS);
  g.root.clientWidth = 100;
  fire(g.root, 'scroll');
  assert.equal(g.columns.value, 5);
  assert.equal(g.root.querySelector('.u2-grid-content').style.height, '4000px');
  assert.equal(cells(g).find((c) => c.dataset.index === '7').style.left, '40px');
  g.dispose();
});

grid('keys move by cell, row and page, clamp at both ends, and Enter activates', () => {
  const picked = [];
  const g = make({onActivate: (n, i) => picked.push([n, i])});
  g.setItems(NUMBERS);
  fire(g.root, 'scroll');
  const key = (k) => fire(g.root, 'keydown', {key: k});
  key('ArrowRight');
  assert.equal(g.selectedIndex.value, 0);
  key('ArrowRight');
  assert.equal(g.selectedIndex.value, 1);
  key('ArrowDown');
  assert.equal(g.selectedIndex.value, 11);
  key('ArrowUp');
  assert.equal(g.selectedIndex.value, 1);
  key('PageDown');
  assert.equal(g.selectedIndex.value, 41, 'a page is visible rows - 1, in cells');
  key('End');
  assert.equal(g.selectedIndex.value, 999);
  key('ArrowRight');
  assert.equal(g.selectedIndex.value, 999, 'clamped');
  assert.equal(g.root.scrollTop > 0, true, 'the selection scrolled into view');
  assert.equal(g.root.getAttribute('aria-activedescendant'), cells(g).find((c) => c.dataset.index === '999').id);
  key('Home');
  assert.equal(g.selectedIndex.value, 0);
  key('Enter');
  key(' ');
  assert.deepEqual(picked, [[0, 0], [0, 0]]);
  g.dispose();
});

grid('a click selects and activates the cell; the selection class follows the lead', () => {
  const picked = [];
  const g = make({onActivate: (n) => picked.push(n)});
  g.setItems(NUMBERS);
  fire(g.root, 'scroll');
  const cell = cells(g).find((c) => c.dataset.index === '23');
  fire(cell.firstChild, 'click');
  assert.deepEqual(picked, [23]);
  assert.equal(g.selectedIndex.value, 23);
  assert.equal(cell.classList.contains('u2-grid-cell-selected'), true);
  assert.equal(cell.getAttribute('aria-selected'), 'true');
  g.selectedIndex.value = 24;
  assert.equal(cell.classList.contains('u2-grid-cell-selected'), false);
  g.dispose();
});

grid('a signal source keeps the grid live, and keyOf keeps the selection on the item', () => {
  const items = signal(['a', 'b', 'c', 'd']);
  const g = make({keyOf: (s) => s});
  g.setItems(items);
  fire(g.root, 'scroll');
  assert.equal(g.renderedCount, 4);
  g.selectedIndex.value = 2;
  items.value = ['x', 'c', 'y'];
  assert.equal(g.renderedCount, 3);
  assert.equal(g.selectedIndex.value, 1, 'c moved to index 1');
  items.value = ['x', 'y'];
  assert.equal(g.selectedIndex.value, -1, 'c is gone');
  g.dispose();
});

grid('an empty grid renders nothing and ignores keys', () => {
  const g = make();
  fire(g.root, 'keydown', {key: 'ArrowRight'});
  assert.equal(g.renderedCount, 0);
  assert.equal(g.selectedIndex.value, -1);
  g.dispose();
});

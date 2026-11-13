const test = require('node:test');
const assert = require('node:assert/strict');
const path = require('path');
const {pathToFileURL} = require('url');

/** What the probe reads of a page and no more: elements with attributes, an inline style read back as the computed
 * one (empty where unset, as happy-dom answers) unless a test sets what the cascade gives (`computed`), the style
 * rules a test lists in `document.styleSheets`, the rectangles a test sets and a hit test the test answers. */
class Element {
  constructor(tag, attrs = {}, children = []) {
    this.tagName = tag.toUpperCase();
    this.attrs = new Map(Object.entries(attrs));
    this.style = {};
    this.rect = rect(0, 0, 0, 0);
    this.parentElement = null;
    this.children = [];
    for (const child of children)
      this.append(child);
  }

  append(child) {
    child.parentElement?.children.splice(child.parentElement.children.indexOf(child), 1);
    child.parentElement = this;
    this.children.push(child);
    return this;
  }

  get className() { return this.attrs.get('class') ?? ''; }
  get textContent() { return (this.attrs.get('text') ?? '') + this.children.map((c) => c.textContent).join(''); }
  getAttribute(name) { return this.attrs.get(name) ?? null; }
  hasAttribute(name) { return this.attrs.has(name); }
  getBoundingClientRect() { return this.rect; }
  getClientRects() { return [this.rect]; }

  contains(other) {
    for (let at = other; at !== null; at = at.parentElement) {
      if (at === this)
        return true;
    }
    return false;
  }

  matches(selector) {
    return selector.split(',').some((one) => Element.complex(this, one.trim().split(/\s+/)));
  }

  closest(selector) {
    for (let at = this; at !== null; at = at.parentElement) {
      if (at.matches(selector))
        return at;
    }
    return null;
  }

  querySelectorAll(selector) {
    return this.children.flatMap((c) => [...c.matches(selector) ? [c] : [], ...c.querySelectorAll(selector)]);
  }

  /** Compounds and `>` read right to left: the last is this element, the ones before it an ancestor. */
  static complex(el, parts) {
    if (!Element.compound(el, parts[parts.length - 1]))
      return false;
    const rest = parts.slice(0, -1);
    if (rest.length === 0)
      return true;
    if (rest[rest.length - 1] === '>')
      return el.parentElement !== null && Element.complex(el.parentElement, rest.slice(0, -1));
    for (let at = el.parentElement; at !== null; at = at.parentElement) {
      if (Element.complex(at, rest))
        return true;
    }
    return false;
  }

  /** A tag, `.class` and `[attr]`/`[attr="v"]`; a pseudo-class never matches: nothing here is hovered or focused. */
  static compound(el, selector) {
    const [, tag, rest] = /^([a-z*]*)(.*)$/.exec(selector);
    if (tag !== '' && tag !== '*' && tag.toUpperCase() !== el.tagName)
      return false;
    const bare = rest.replace(/\[[^\]]*\]/g, '');
    const classes = el.className.split(/\s+/);
    return !bare.includes(':') && [...bare.matchAll(/\.([\w-]+)/g)].every(([, name]) => classes.includes(name)) &&
      [...rest.matchAll(/\[([\w-]+)(?:="([^"]*)")?\]/g)].every(([, name, value]) =>
        el.hasAttribute(name) && (value === undefined || el.getAttribute(name) === value));
  }
}

function rect(left, top, width, height) {
  return {left, top, width, height, right: left + width, bottom: top + height};
}

const el = (tag, attrs, children) => new Element(tag, attrs, children);
const CONTROL = Symbol.for('u2.control');

/** `root` as a u2 control's root (the engine reads the control off it, whichever copy built it). */
function asControl(root) {
  root[CONTROL] = {root};
  return root;
}

let hit = () => null;
Object.assign(globalThis, {
  Element, innerWidth: 1024, innerHeight: 768,
  getComputedStyle: (e) => new Proxy(e.style, {get: (style, name) => e.computed?.[name] ?? style[name] ?? ''}),
  document: {elementFromPoint: (x, y) => hit(x, y), getElementById: () => null, styleSheets: []},
});
test.afterEach(() => document.styleSheets = []);

let U2;
test.before(async () => {
  U2 = await import(pathToFileURL(path.resolve(__dirname, '../../dist/src/u2core/index.js')).href);
});

const findings = (root, options) => U2.operability(root, options).map((f) => [f.rule, f.detail]);
const centre = (e) => [e.rect.left + e.rect.width / 2, e.rect.top + e.rect.height / 2];
const hits = (target, by) => (x, y) => x === centre(target)[0] && y === centre(target)[1] ? by : null;

test('covered: an ancestor whose computed overflow is empty (no initial values, as on happy-dom) clips nothing', () => {
  const edit = el('button', {'aria-label': 'Edit'});
  edit.rect = rect(100, 10, 16, 28);
  const cover = el('div', {class: 'u2-input-box'});
  const page = el('div', {}, [el('div', {}, [edit]), cover]);
  hit = hits(edit, cover);
  assert.deepEqual(findings(page, {hitTest: true}), [['covered', 'button "Edit" is covered by div.u2-input-box']]);
  edit.parentElement.style.overflowX = 'hidden';
  assert.deepEqual(findings(page, {hitTest: true}), [], 'a centre outside a clipping ancestor is out of view');
});

test('an aria-hidden element is not offered: the native colour input under its swatch', () => {
  const picker = el('input', {class: 'u2-color-picker', tabindex: '-1', 'aria-hidden': 'true'});
  picker.style.pointerEvents = 'none';
  const box = el('div', {}, [el('button', {'aria-label': 'Pick Color', text: '■'}), picker]);
  assert.deepEqual(findings(box), []);
  picker.attrs.delete('aria-hidden');
  assert.deepEqual(findings(box), [['no-pointer', 'input.u2-color-picker takes no pointer: pointer-events: none']]);
});

test('covered: a header scrolled under a sticky column of its own scroller is out of view, not covered', () => {
  const header = el('th', {tabindex: '0', text: 'Last Modified'});
  header.rect = rect(240, 0, 100, 28);
  const pinned = el('th', {class: 'u2-table-pin-end'});
  pinned.style.position = 'sticky';
  const scroller = el('div', {class: 'u2-table-scroll'}, [el('tr', {}, [header, pinned])]);
  scroller.style.overflowX = 'auto';
  scroller.rect = rect(0, 0, 300, 100);
  const page = el('div', {}, [scroller]);
  hit = hits(header, pinned);
  assert.deepEqual(findings(page, {hitTest: true}), []);

  pinned.style.position = 'absolute';
  assert.deepEqual(findings(page, {hitTest: true}),
    [['covered', 'th "Last Modified" is covered by th.u2-table-pin-end']], 'a layer over the header is a cover');
  const outside = el('div', {class: 'page-header'});
  outside.style.position = 'sticky';
  page.append(outside);
  hit = hits(header, outside);
  assert.deepEqual(findings(page, {hitTest: true}), [['covered', 'th "Last Modified" is covered by div.page-header']],
    'a sticky element of another scroller is a cover');
  const bar = el('div', {class: 'toolbar'}, [header, el('span', {class: 'badge'})]);
  bar.style.position = 'sticky';
  scroller.append(bar);
  hit = hits(header, bar.children[1]);
  assert.deepEqual(findings(page, {hitTest: true}), [['covered', 'th "Last Modified" is covered by span.badge']],
    'a sticky bar the header rides in moves with it');
});

test('what a platform widget hosted below the probed root holds is the platform\'s to judge', () => {
  const canvas = el('canvas', {tabindex: '1000'});
  canvas.rect = rect(10, 10, 200, 100);
  const overlay = el('canvas');
  const editor = el('input', {class: 'ui-input-editor', disabled: ''});
  const off = el('button', {disabled: '', text: 'Save'});
  const inner = asControl(el('div', {'data-u2': 'button-bar'}, [off]));
  const grid = el('div', {class: 'd4-grid', 'data-widget': 'true'}, [canvas, overlay, editor, inner]);
  const host = asControl(el('div', {'data-u2': 'columns-input'}, [grid]));
  hit = hits(canvas, overlay);
  assert.deepEqual(findings(host, {hitTest: true}),
    [['disabled-without-reason', 'button "Save" is disabled and says no reason']],
    'a u2 control inside the platform widget is u2\'s again');

  assert.deepEqual(findings(grid, {hitTest: true}).map(([rule]) => rule),
    ['covered', 'disabled-without-reason', 'disabled-without-reason'], 'probed itself, the widget is judged');
  grid.attrs.delete('data-widget');
  grid[CONTROL] = {root: grid, toDart() {}};
  assert.deepEqual(findings(host).map(([rule]) => rule), ['disabled-without-reason'],
    'a platform wrapper on the root marks it too');
});

test('inside an ancestor that is off and says why, a nested control is off with it: no finding for what it holds', () => {
  const and = el('button', {disabled: '', text: 'and'});
  const advanced = el('a', {tabindex: '0', text: 'Advanced'});
  const header = el('div', {class: 'u2-fb-header'}, [and, advanced]);
  const rows = asControl(el('div', {class: 'u2-fb-root'}, [header]));
  advanced.computed = header.computed = rows.computed = {pointerEvents: 'none'};
  const panel = el('div', {class: 'u2-fb'}, [rows]);
  panel.style.pointerEvents = 'none';
  const builder = asControl(el('div', {'data-u2': 'filter-builder', 'aria-disabled': 'true'}, [panel]));
  const cell = el('div', {}, [builder]);
  const BUILDER = [['disabled-without-reason', 'button "and" is disabled and says no reason'],
    ['no-pointer', 'a "Advanced" takes no pointer: pointer-events: none on div.u2-fb-root']];
  assert.deepEqual(findings(cell), BUILDER, 'an off builder that says nothing');
  builder.attrs.set('aria-description', ' ');
  assert.deepEqual(findings(cell), BUILDER, 'a blank reason is none');
  builder.attrs.set('aria-description', 'Fill in the fields above first');
  assert.deepEqual(findings(cell), [], 'the builder\'s reason is theirs');

  const input = el('input', {class: 'u2-input-editor'});
  input.rect = rect(40, 40, 120, 24);
  const surface = el('div', {class: 'u2-designer-surface'}, [asControl(el('div', {class: 'u2-input-root'}, [input]))]);
  const pane = el('div', {class: 'u2-designer-pane'});
  const designer = asControl(el('div', {'data-u2': 'designer'}, [surface, pane]));
  hit = hits(input, pane);
  assert.deepEqual(findings(designer, {hitTest: true}),
    [['covered', 'input.u2-input-editor is covered by div.u2-designer-pane']], 'a canvas that says nothing');
  surface.attrs.set('aria-disabled', 'true');
  surface.attrs.set('aria-description', 'Switch to Run to use the controls');
  assert.deepEqual(findings(designer, {hitTest: true}), [], 'design mode: the canvas is off, and says why');
});

/** u2's rules for `rowActions` (css/list.css): no pointer until the row is hovered, or a button in it has focus. */
const ROW_ACTIONS = {cssRules: [
  {selectorText: '.u2-row-actions', style: {pointerEvents: 'none'}},
  {selectorText: '.u2-list-row:hover .u2-row-actions, .u2-row:hover .u2-row-actions, .u2-row-actions:focus-within',
    style: {pointerEvents: 'auto'}},
]};
const CROSS_ORIGIN = {get cssRules() { throw new Error('SecurityError: cannot access rules'); }};
const DEAD = ['no-pointer', 'button "Open" takes no pointer: pointer-events: none on div.u2-row-actions'];

/** A row with its action block, which at rest takes no pointer. */
function row(cls = 'u2-list-row') {
  const open = el('button', {'aria-label': 'Open'});
  open.rect = rect(260, 10, 24, 24);
  const actions = el('div', {class: 'u2-row-actions', role: 'toolbar'}, [open]);
  actions.computed = {pointerEvents: 'none'};
  const item = el('div', {class: cls}, [el('div', {class: 'u2-domain-list-item'}, [el('span', {text: 'Aspirin'}),
    actions])]);
  item.rect = rect(0, 0, 300, 44);
  return {open, actions, item};
}

test('no-pointer: a row action the row\'s hover reveals takes the pointer there, and is not hit-tested at rest', () => {
  const {open, actions, item} = row();
  const list = asControl(el('div', {'data-u2': 'domain-list'}, [item]));
  hit = hits(open, item);
  assert.deepEqual(findings(list, {hitTest: true}), [DEAD], 'no rule reveals it');
  document.styleSheets = [CROSS_ORIGIN, ROW_ACTIONS];
  assert.deepEqual(findings(list, {hitTest: true}), [], 'revealed on the list row\'s hover');
  item.attrs.set('class', 'u2-table-row u2-row');
  assert.deepEqual(findings(list, {hitTest: true}), [], 'and on any .u2-row\'s (a catalog table row)');
  list.append(actions);
  assert.deepEqual(findings(list, {hitTest: true}), [DEAD], 'outside a row nothing reveals it');
});

test('no-pointer: a row action still dead with its row hovered is a finding', () => {
  document.styleSheets = [ROW_ACTIONS,
    {cssRules: [{selectorText: '.u2-row-actions .u2-busy', style: {pointerEvents: 'none'}}]}];
  const busy = row();
  busy.open.attrs.set('class', 'u2-busy');
  const inline = row();
  inline.open.style.pointerEvents = 'none';
  const fenced = row();
  const steps = el('div', {class: 'u2-wizard-steps'}, [fenced.item]);
  steps.computed = fenced.item.computed = {pointerEvents: 'none'};
  assert.deepEqual(findings(asControl(el('div', {}, [busy.item, inline.item, steps]))), [DEAD, DEAD,
    ['no-pointer', 'button "Open" takes no pointer: pointer-events: none on div.u2-wizard-steps']],
  'its own none by a rule or inline, and a row inside a region that takes no pointer cannot be hovered');

  document.styleSheets = [{cssRules: [ROW_ACTIONS.cssRules[0],
    {selectorText: '.u2-row-actions:hover, .u2-row-actions:focus-within', style: {pointerEvents: 'auto'}},
    {selectorText: '.u2-list-row:hover .u2-row-actions', style: {pointerEvents: 'none'}}]}];
  assert.deepEqual(findings(asControl(el('div', {}, [row().item]))), [DEAD],
    'a block revealed on its own hover or focus only, or kept off on the row\'s: the pointer never gets in');
});

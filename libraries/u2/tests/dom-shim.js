/* Minimal DOM for node: enough to construct u2 components, dispatch events and read back
   attributes, classes and text. Layout metrics (rects, offsetWidth, clientWidth) are settable
   fakes — nothing here lays anything out. Installed on import; no dependencies. */

const DEFAULT_STYLE = {
  position: 'static', display: 'block', direction: 'ltr',
  overflow: 'visible', overflowX: 'visible', overflowY: 'visible',
  transform: 'none', translate: 'none', scale: 'none', rotate: 'none', perspective: 'none',
  filter: 'none', backdropFilter: 'none', willChange: 'auto', contain: 'none',
  width: '0px', height: '0px', columnGap: '0px', scrollbarGutter: 'auto',
  paddingLeft: '0px', paddingRight: '0px', paddingTop: '0px', paddingBottom: '0px',
  marginLeft: '0px', marginRight: '0px',
};

const METRICS = ['offsetWidth', 'offsetHeight', 'clientWidth', 'clientHeight',
  'scrollWidth', 'scrollHeight', 'scrollTop', 'scrollLeft', 'clientTop', 'clientLeft'];

function camel(name) {
  return name.replace(/-([a-z])/g, (_, c) => c.toUpperCase());
}

function dashed(name) {
  return name.replace(/[A-Z]/g, (c) => '-' + c.toLowerCase());
}

function createStyle() {
  const style = {};
  const hide = (name, value) => Object.defineProperty(style, name, {value});
  hide('setProperty', (name, value) => style[name] = String(value));
  hide('removeProperty', (name) => delete style[name]);
  hide('getPropertyValue', (name) => style[name] ?? '');
  Object.defineProperty(style, 'cssText', {
    get: () => Object.keys(style).map((k) => `${dashed(k)}: ${style[k]}`).join('; '),
    set: (text) => {
      for (const part of String(text).split(';')) {
        const i = part.indexOf(':');
        if (i > 0)
          style[camel(part.slice(0, i).trim())] = part.slice(i + 1).trim();
      }
    },
  });
  return style;
}

class DOMRect {
  constructor(x = 0, y = 0, width = 0, height = 0) {
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }

  get left() {
    return this.x;
  }

  get top() {
    return this.y;
  }

  get right() {
    return this.x + this.width;
  }

  get bottom() {
    return this.y + this.height;
  }
}

class DomEvent {
  constructor(type, init = {}) {
    this.type = type;
    this.bubbles = false;
    this.cancelable = true;
    this.defaultPrevented = false;
    this.target = null;
    this.currentTarget = null;
    Object.assign(this, init);
    this._stopped = false;
    this._done = false;
  }

  preventDefault() {
    this.defaultPrevented = true;
  }

  stopPropagation() {
    this._stopped = true;
  }

  stopImmediatePropagation() {
    this._stopped = true;
    this._done = true;
  }
}

class CustomEvent extends DomEvent {
  constructor(type, init = {}) {
    super(type, init);
    this.detail = init.detail ?? null;
  }
}

class KeyboardEvent extends DomEvent {
  constructor(type, init = {}) {
    super(type, {key: '', shiftKey: false, ctrlKey: false, altKey: false, metaKey: false, ...init});
  }
}

class MouseEvent extends DomEvent {
  constructor(type, init = {}) {
    super(type, {button: 0, clientX: 0, clientY: 0, ...init});
  }
}

class PointerEvent extends MouseEvent {
  constructor(type, init = {}) {
    super(type, {pointerId: 1, ...init});
  }
}

class WheelEvent extends MouseEvent {
  constructor(type, init = {}) {
    super(type, {deltaY: 0, ...init});
  }
}

class EventTargetBase {
  constructor() {
    this._listeners = new Map();
  }

  addEventListener(type, handler, options) {
    const capture = options === true || (!!options && options.capture === true);
    const list = this._listeners.get(type) ?? [];
    if (!list.some((l) => l.handler === handler && l.capture === capture))
      list.push({handler, capture});
    this._listeners.set(type, list);
  }

  removeEventListener(type, handler, options) {
    const capture = options === true || (!!options && options.capture === true);
    const list = this._listeners.get(type);
    if (!list)
      return;
    const index = list.findIndex((l) => l.handler === handler && l.capture === capture);
    if (index >= 0)
      list.splice(index, 1);
  }

  _invoke(event, capture) {
    const list = this._listeners.get(event.type);
    if (!list)
      return;
    event.currentTarget = this;
    for (const entry of list.slice()) {
      if (entry.capture !== capture)
        continue;
      entry.handler.call(this, event);
      if (event._done)
        return;
    }
  }

  dispatchEvent(event) {
    event.target = this;
    const path = [];
    for (let node = this; node; node = node.parentNode)
      path.push(node);
    for (let i = path.length - 1; i > 0 && !event._stopped; i--)
      path[i]._invoke(event, true);
    if (!event._stopped)
      this._invoke(event, true);
    if (!event._done)
      this._invoke(event, false);
    if (event.bubbles) {
      for (let i = 1; i < path.length && !event._stopped; i++)
        path[i]._invoke(event, false);
    }
    event.currentTarget = null;
    return !event.defaultPrevented;
  }
}

const observers = new Set();

function notifyMutation(parent) {
  if (observers.size === 0)
    return;
  for (const observer of observers) {
    if (observer._queued || !observer._watches(parent))
      continue;
    observer._queued = true;
    queueMicrotask(() => {
      observer._queued = false;
      observer._callback([], observer);
    });
  }
}

class MutationObserver {
  constructor(callback) {
    this._callback = callback;
    this._targets = [];
    this._queued = false;
  }

  observe(target, options = {}) {
    this._targets.push({target, subtree: !!options.subtree});
    observers.add(this);
  }

  disconnect() {
    this._targets.length = 0;
    observers.delete(this);
  }

  takeRecords() {
    return [];
  }

  _watches(node) {
    return this._targets.some((t) => t.target === node || (t.subtree && t.target.contains(node)));
  }
}

class ResizeObserver {
  constructor(callback) {
    this.callback = callback;
    this.targets = [];
    this.disconnected = false;
    ResizeObserver.instances.push(this);
  }

  observe(target) {
    if (!this.targets.includes(target))
      this.targets.push(target);
  }

  unobserve(target) {
    const i = this.targets.indexOf(target);
    if (i >= 0)
      this.targets.splice(i, 1);
  }

  disconnect() {
    this.disconnected = true;
    this.targets.length = 0;
  }
}
ResizeObserver.instances = [];

class ShadowRoot {}

class Node extends EventTargetBase {
  constructor(nodeName, nodeType) {
    super();
    this.nodeName = nodeName;
    this.nodeType = nodeType;
    this.parentNode = null;
    this.childNodes = [];
    this.ownerDocument = document;
  }

  get isConnected() {
    let node = this;
    while (node.parentNode)
      node = node.parentNode;
    return node === document;
  }

  get parentElement() {
    return this.parentNode instanceof Element ? this.parentNode : null;
  }

  get nextSibling() {
    const siblings = this.parentNode ? this.parentNode.childNodes : [];
    return siblings[siblings.indexOf(this) + 1] ?? null;
  }

  contains(node) {
    for (let n = node; n; n = n.parentNode) {
      if (n === this)
        return true;
    }
    return false;
  }

  remove() {
    if (this.parentNode)
      this.parentNode.removeChild(this);
  }

  replaceWith(...nodes) {
    const parent = this.parentNode;
    if (!parent)
      return;
    const index = parent.childNodes.indexOf(this);
    parent.removeChild(this);
    for (let i = 0; i < nodes.length; i++)
      parent._insertAt(toNode(nodes[i]), index + i);
    notifyMutation(parent);
  }
}

class Text extends Node {
  constructor(data) {
    super('#text', 3);
    this.data = String(data);
  }

  get textContent() {
    return this.data;
  }

  set textContent(value) {
    this.data = String(value);
  }
}

function toNode(x) {
  return x instanceof Node ? x : new Text(x);
}

class ClassList {
  constructor(el) {
    this._el = el;
  }

  add(...names) {
    for (const name of names)
      this._el._classes.add(name);
  }

  remove(...names) {
    for (const name of names)
      this._el._classes.delete(name);
  }

  contains(name) {
    return this._el._classes.has(name);
  }

  toggle(name, force) {
    const on = force === undefined ? !this.contains(name) : force;
    if (on)
      this.add(name);
    else
      this.remove(name);
    return on;
  }

  get length() {
    return this._el._classes.size;
  }
}

class Element extends Node {
  constructor(tag) {
    super(tag.toUpperCase(), 1);
    this.tagName = tag.toUpperCase();
    this.style = createStyle();
    this.rect = new DOMRect();
    this._attributes = new Map();
    this._classes = new Set();
    for (const metric of METRICS)
      this[metric] = 0;
  }

  get classList() {
    if (!this._classList)
      this._classList = new ClassList(this);
    return this._classList;
  }

  get className() {
    return [...this._classes].join(' ');
  }

  set className(value) {
    this._classes = new Set(String(value).split(/\s+/).filter((c) => c));
  }

  get dataset() {
    if (!this._dataset) {
      const el = this;
      this._dataset = new Proxy({}, {
        get: (_, key) => typeof key === 'string' ? el.getAttribute('data-' + dashed(key)) ?? undefined : undefined,
        set: (_, key, value) => {
          el.setAttribute('data-' + dashed(key), value);
          return true;
        },
        has: (_, key) => el.hasAttribute('data-' + dashed(key)),
        deleteProperty: (_, key) => {
          el.removeAttribute('data-' + dashed(key));
          return true;
        },
      });
    }
    return this._dataset;
  }

  get children() {
    return this.childNodes.filter((c) => c instanceof Element);
  }

  get firstChild() {
    return this.childNodes[0] ?? null;
  }

  get textContent() {
    return this.childNodes.map((c) => c.textContent).join('');
  }

  set textContent(value) {
    for (const child of this.childNodes)
      child.parentNode = null;
    this.childNodes.length = 0;
    if (value !== '' && value !== null && value !== undefined)
      this.appendChild(new Text(value));
    notifyMutation(this);
  }

  get id() {
    return this.getAttribute('id') ?? '';
  }

  set id(value) {
    this.setAttribute('id', value);
  }

  get title() {
    return this.getAttribute('title') ?? '';
  }

  set title(value) {
    this.setAttribute('title', value);
  }

  get type() {
    return this.getAttribute('type') ?? '';
  }

  set type(value) {
    this.setAttribute('type', value);
  }

  get tabIndex() {
    return this.hasAttribute('tabindex') ? Number(this.getAttribute('tabindex')) : -1;
  }

  set tabIndex(value) {
    this.setAttribute('tabindex', String(value));
  }

  get hidden() {
    return this.hasAttribute('hidden');
  }

  set hidden(value) {
    if (value)
      this.setAttribute('hidden', '');
    else
      this.removeAttribute('hidden');
  }

  get disabled() {
    return this.hasAttribute('disabled');
  }

  set disabled(value) {
    if (value)
      this.setAttribute('disabled', '');
    else
      this.removeAttribute('disabled');
  }

  get offsetParent() {
    return this.isConnected ? this.parentNode : null;
  }

  getAttribute(name) {
    return name === 'class' ? this.className : this._attributes.get(name) ?? null;
  }

  setAttribute(name, value) {
    if (name === 'class')
      this.className = value;
    else
      this._attributes.set(name, String(value));
  }

  removeAttribute(name) {
    if (name === 'class')
      this._classes.clear();
    else
      this._attributes.delete(name);
  }

  getAttributeNames() {
    const names = [...this._attributes.keys()];
    return this._classes.size > 0 ? ['class', ...names] : names;
  }

  hasAttribute(name) {
    return name === 'class' ? this._classes.size > 0 : this._attributes.has(name);
  }

  appendChild(node) {
    return this._insertAt(node, this.childNodes.length, true);
  }

  append(...nodes) {
    for (const node of nodes)
      this._insertAt(toNode(node), this.childNodes.length);
    notifyMutation(this);
  }

  prepend(...nodes) {
    for (let i = 0; i < nodes.length; i++)
      this._insertAt(toNode(nodes[i]), i);
    notifyMutation(this);
  }

  insertBefore(node, reference) {
    const index = reference ? this.childNodes.indexOf(reference) : this.childNodes.length;
    return this._insertAt(node, index < 0 ? this.childNodes.length : index, true);
  }

  replaceChildren(...nodes) {
    for (const child of this.childNodes)
      child.parentNode = null;
    this.childNodes.length = 0;
    for (const node of nodes)
      this._insertAt(toNode(node), this.childNodes.length);
    notifyMutation(this);
  }

  removeChild(node) {
    const index = this.childNodes.indexOf(node);
    if (index >= 0) {
      this.childNodes.splice(index, 1);
      node.parentNode = null;
      notifyMutation(this);
    }
    return node;
  }

  matches(selector) {
    return matchesSelector(this, selector);
  }

  closest(selector) {
    for (let el = this; el instanceof Element; el = el.parentNode) {
      if (el.matches(selector))
        return el;
    }
    return null;
  }

  querySelector(selector) {
    return this.querySelectorAll(selector)[0] ?? null;
  }

  querySelectorAll(selector) {
    const found = [];
    collect(this, (el) => {
      if (matchesSelector(el, selector))
        found.push(el);
    });
    return found;
  }

  getBoundingClientRect() {
    return this.rect;
  }

  getClientRects() {
    return [this.rect];
  }

  focus() {
    if (document.activeElement === this)
      return;
    const previous = document.activeElement;
    document.activeElement = this;
    if (previous) {
      previous.dispatchEvent(new DomEvent('blur'));
      previous.dispatchEvent(new DomEvent('focusout', {bubbles: true}));
    }
    this.dispatchEvent(new DomEvent('focus'));
    this.dispatchEvent(new DomEvent('focusin', {bubbles: true}));
  }

  blur() {
    if (document.activeElement !== this)
      return;
    document.activeElement = document.body;
    this.dispatchEvent(new DomEvent('blur'));
  }

  scrollIntoView() {}

  setPointerCapture() {}

  releasePointerCapture() {}

  select() {}

  click() {
    this.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  }

  _insertAt(node, index, notify = false) {
    if (node.parentNode)
      node.parentNode.removeChild(node);
    node.parentNode = this;
    this.childNodes.splice(index, 0, node);
    if (notify)
      notifyMutation(this);
    return node;
  }
}

class HTMLElement extends Element {}
class HTMLDivElement extends HTMLElement {}
class HTMLButtonElement extends HTMLElement {}
class HTMLAnchorElement extends HTMLElement {}
class HTMLOptionElement extends HTMLElement {}

class HTMLInputElement extends HTMLElement {
  constructor(tag = 'input') {
    super(tag);
    this.value = '';
    this.checked = false;
    this.type = 'text';
  }
}

class HTMLTextAreaElement extends HTMLElement {
  constructor(tag = 'textarea') {
    super(tag);
    this.value = '';
  }
}

/** `value` is a plain field, not a projection of the selected option — enough for change tests. */
class HTMLSelectElement extends HTMLElement {
  constructor(tag = 'select') {
    super(tag);
    this.value = '';
  }
}

function Option(text = '', value = '') {
  const option = new HTMLOptionElement('option');
  option.textContent = text;
  option.value = value;
  return option;
}

const TAGS = {
  input: HTMLInputElement, textarea: HTMLTextAreaElement, select: HTMLSelectElement,
  option: HTMLOptionElement, button: HTMLButtonElement, a: HTMLAnchorElement, div: HTMLDivElement,
};

function collect(root, visit) {
  for (const child of root.childNodes) {
    if (!(child instanceof Element))
      continue;
    visit(child);
    collect(child, visit);
  }
}

/* Selector subset: tag, .class, #id, [attr], [attr="value"], :not(...), comma lists and
   descendant combinators. Any other pseudo-class never matches. */
const selectorCache = new Map();

function splitTop(text, separator) {
  const parts = [];
  let depth = 0;
  let current = '';
  for (const ch of text) {
    if (ch === '(' || ch === '[')
      depth++;
    else if (ch === ')' || ch === ']')
      depth--;
    if (depth === 0 && (separator === ',' ? ch === ',' : /\s/.test(ch))) {
      if (current.trim())
        parts.push(current.trim());
      current = '';
      continue;
    }
    current += ch;
  }
  if (current.trim())
    parts.push(current.trim());
  return parts;
}

function parseCompound(text) {
  const compound = {tag: null, id: null, classes: [], attrs: [], nots: [], never: false};
  const token = /:not\(([^()]*)\)|:[\w-]+|\[([^\]]*)\]|\.([\w-]+)|#([\w-]+)|([\w-]+)|\*/g;
  let match;
  while ((match = token.exec(text)) !== null) {
    if (match[1] !== undefined)
      compound.nots.push(parseCompound(match[1]));
    else if (match[0].startsWith(':'))
      compound.never = true;
    else if (match[2] !== undefined) {
      const eq = match[2].indexOf('=');
      const name = (eq < 0 ? match[2] : match[2].slice(0, eq)).trim();
      const value = eq < 0 ? null : match[2].slice(eq + 1).trim().replace(/^["']|["']$/g, '');
      compound.attrs.push({name, value});
    } else if (match[3] !== undefined)
      compound.classes.push(match[3]);
    else if (match[4] !== undefined)
      compound.id = match[4];
    else if (match[5] !== undefined)
      compound.tag = match[5].toUpperCase();
  }
  return compound;
}

function parseSelector(selector) {
  let parsed = selectorCache.get(selector);
  if (!parsed) {
    parsed = splitTop(selector, ',').map((part) => splitTop(part, ' ').map(parseCompound));
    selectorCache.set(selector, parsed);
  }
  return parsed;
}

function matchesCompound(el, compound) {
  if (compound.never || (compound.tag && el.tagName !== compound.tag))
    return false;
  if (compound.id !== null && el.id !== compound.id)
    return false;
  for (const cls of compound.classes) {
    if (!el.classList.contains(cls))
      return false;
  }
  for (const attr of compound.attrs) {
    if (!el.hasAttribute(attr.name) || (attr.value !== null && el.getAttribute(attr.name) !== attr.value))
      return false;
  }
  for (const not of compound.nots) {
    if (matchesCompound(el, not))
      return false;
  }
  return true;
}

function matchesSelector(el, selector) {
  for (const chain of parseSelector(selector)) {
    if (!matchesCompound(el, chain[chain.length - 1]))
      continue;
    let index = chain.length - 2;
    let ancestor = el.parentNode;
    while (index >= 0 && ancestor instanceof Element) {
      if (matchesCompound(ancestor, chain[index]))
        index--;
      ancestor = ancestor.parentNode;
    }
    if (index < 0)
      return true;
  }
  return false;
}

class Document extends EventTargetBase {
  constructor() {
    super();
    this.nodeName = '#document';
    this.nodeType = 9;
    this.compatMode = 'CSS1Compat';
    this.defaultView = globalThis;
    this.parentNode = null;
  }

  createElement(tag) {
    return new (TAGS[tag] ?? HTMLElement)(tag);
  }

  createTextNode(data) {
    return new Text(data);
  }

  querySelector(selector) {
    return this.documentElement.querySelector(selector);
  }

  querySelectorAll(selector) {
    return this.documentElement.querySelectorAll(selector);
  }

  getElementById(id) {
    return this.documentElement.querySelector(`#${id}`);
  }
}

let document = new Document();
document.documentElement = new HTMLElement('html');
document.documentElement.parentNode = document;
document.head = new HTMLElement('head');
document.body = new HTMLElement('body');
document.documentElement.append(document.head, document.body);
document.documentElement.clientWidth = 1024;
document.documentElement.clientHeight = 768;
document.activeElement = document.body;

const frames = new Map();
let frameSeq = 0;

function requestAnimationFrame(fn) {
  const id = ++frameSeq;
  frames.set(id, setTimeout(() => {
    frames.delete(id);
    fn(performance.now());
  }, 0));
  return id;
}

function cancelAnimationFrame(id) {
  clearTimeout(frames.get(id));
  frames.delete(id);
}

function getComputedStyle(el) {
  const computed = Object.assign({}, DEFAULT_STYLE, el.style);
  // custom properties are stored (and looked up) under their raw -- name, never camelized
  computed.getPropertyValue = (name) =>
    computed[name.startsWith('--') ? name : camel(name)] ?? '';
  return computed;
}

const windowTarget = new EventTargetBase();

Object.assign(globalThis, {
  window: globalThis, document, Node, Text, Element, HTMLElement, HTMLDivElement, HTMLInputElement,
  HTMLTextAreaElement, HTMLSelectElement, HTMLButtonElement, HTMLAnchorElement, HTMLOptionElement,
  ShadowRoot, DOMRect, Option, Event: DomEvent, CustomEvent, KeyboardEvent, MouseEvent, PointerEvent,
  WheelEvent, MutationObserver, ResizeObserver, getComputedStyle, requestAnimationFrame,
  cancelAnimationFrame, innerWidth: 1024, innerHeight: 768, scrollX: 0, scrollY: 0,
  addEventListener: (...args) => windowTarget.addEventListener(...args),
  removeEventListener: (...args) => windowTarget.removeEventListener(...args),
  dispatchEvent: (event) => windowTarget.dispatchEvent(event),
});

/** Dispatches a bubbling event; `init` fields land on the event object. */
export function fire(el, type, init = {}) {
  const types = {keydown: KeyboardEvent, keyup: KeyboardEvent, click: MouseEvent, dblclick: MouseEvent};
  const pointer = type.startsWith('pointer') || type === 'mousedown' || type === 'contextmenu';
  const Type = pointer ? PointerEvent : types[type] ?? DomEvent;
  return el.dispatchEvent(new Type(type, {bubbles: true, ...init}));
}

/** Lets queued microtasks and animation frames run. */
export function flush() {
  return new Promise((resolve) => setTimeout(resolve, 1));
}

/** Empties the document between tests so leftover overlays never cross-contaminate. */
export function resetDom() {
  document.body.replaceChildren();
  document.activeElement = document.body;
}

/**
 * Minimal inert DOM for Node.js — installed by `startDatagrok` when no real
 * `document` exists. It is just enough for the pure-DOM parts of the js-api
 * (`ui.div()`-style builders, cash-dom, widgets' element bookkeeping) to run
 * without a browser: elements are plain containers, nothing ever renders.
 * Dart-backed UI (dialogs, viewers, grid) still degrades via the window Proxy.
 */

class FakeClassList {
  private _set = new Set<string>();
  add(...names: string[]) { for (const n of names) if (n) this._set.add(n); }
  remove(...names: string[]) { for (const n of names) this._set.delete(n); }
  contains(name: string): boolean { return this._set.has(name); }
  toggle(name: string, force?: boolean): boolean {
    const on = force ?? !this._set.has(name);
    on ? this._set.add(name) : this._set.delete(name);
    return on;
  }
  get value(): string { return [...this._set].join(' '); }
  set value(v: string) { this._set = new Set(v.split(/\s+/).filter((s) => s)); }
  forEach(f: (name: string) => void) { this._set.forEach(f); }
  get length(): number { return this._set.size; }
  toString(): string { return this.value; }
}

class FakeNode {
  childNodes: FakeNode[] = [];
  parentNode: FakeNode | null = null;
  parentElement: FakeNode | null = null;
  nodeType = 1;
  nodeName = '';
  textContent: any = '';

  appendChild<T extends FakeNode>(child: T): T {
    child.remove();
    this.childNodes.push(child);
    child.parentNode = child.parentElement = this;
    return child;
  }
  append(...children: any[]) { for (const c of children) this.appendChild(typeof c === 'object' && c !== null ? c : new FakeText(String(c))); }
  prepend(...children: any[]) { for (const c of children.reverse()) this.insertBefore(typeof c === 'object' && c !== null ? c : new FakeText(String(c)), this.firstChild); }
  insertBefore<T extends FakeNode>(child: T, ref: FakeNode | null): T {
    child.remove();
    const i = ref ? this.childNodes.indexOf(ref) : -1;
    i >= 0 ? this.childNodes.splice(i, 0, child) : this.childNodes.push(child);
    child.parentNode = child.parentElement = this;
    return child;
  }
  removeChild<T extends FakeNode>(child: T): T {
    const i = this.childNodes.indexOf(child);
    if (i >= 0) this.childNodes.splice(i, 1);
    child.parentNode = child.parentElement = null;
    return child;
  }
  replaceChild(newChild: FakeNode, oldChild: FakeNode): FakeNode {
    this.insertBefore(newChild, oldChild);
    return this.removeChild(oldChild);
  }
  remove() { this.parentNode?.removeChild(this); }
  cloneNode(_deep?: boolean): FakeNode { return new (this.constructor as any)(); }
  contains(other: FakeNode | null): boolean {
    for (let n = other; n != null; n = n.parentNode)
      if (n === this) return true;
    return false;
  }
  get firstChild(): FakeNode | null { return this.childNodes[0] ?? null; }
  get lastChild(): FakeNode | null { return this.childNodes[this.childNodes.length - 1] ?? null; }
  get nextSibling(): FakeNode | null {
    const sib = this.parentNode?.childNodes;
    return sib ? sib[sib.indexOf(this) + 1] ?? null : null;
  }
  get previousSibling(): FakeNode | null {
    const sib = this.parentNode?.childNodes;
    return sib ? sib[sib.indexOf(this) - 1] ?? null : null;
  }
  hasChildNodes(): boolean { return this.childNodes.length > 0; }
  addEventListener() {}
  removeEventListener() {}
  dispatchEvent(): boolean { return true; }
}

class FakeText extends FakeNode {
  nodeType = 3;
  nodeName = '#text';
  constructor(text: string = '') { super(); this.textContent = text; }
  get data() { return this.textContent; }
  set data(v: any) { this.textContent = v; }
}

class FakeDocumentFragment extends FakeNode {
  nodeType = 11;
  nodeName = '#document-fragment';
}

function makeStyle(): any {
  const store: Record<string, string> = {};
  return new Proxy(store, {
    get(t, prop: string) {
      if (prop === 'setProperty') return (k: string, v: string) => { t[k] = v; };
      if (prop === 'getPropertyValue') return (k: string) => t[k] ?? '';
      if (prop === 'removeProperty') return (k: string) => { const v = t[k]; delete t[k]; return v ?? ''; };
      if (prop === 'cssText') return Object.entries(t).map(([k, v]) => `${k}: ${v}`).join('; ');
      return t[prop] ?? '';
    },
    set(t, prop: string, value) { t[prop] = value; return true; },
  });
}

class FakeElement extends FakeNode {
  tagName: string;
  namespaceURI: string | null = null;
  id = '';
  classList = new FakeClassList();
  style = makeStyle();
  dataset: Record<string, string> = {};
  attributes: Record<string, string> = {};
  innerHTML = '';
  innerText: any = '';
  value: any = '';
  checked = false;
  disabled = false;
  tabIndex = -1;
  scrollTop = 0;
  scrollLeft = 0;
  offsetWidth = 0;
  offsetHeight = 0;
  clientWidth = 0;
  clientHeight = 0;
  onclick: any = null;
  ownerDocument: any;

  constructor(tagName: string = 'div') {
    super();
    this.tagName = tagName.toUpperCase();
    this.nodeName = this.tagName;
  }
  get className(): string { return this.classList.value; }
  set className(v: string) { this.classList.value = v; }
  get children(): FakeElement[] { return this.childNodes.filter((n) => n instanceof FakeElement) as FakeElement[]; }
  get firstElementChild(): FakeElement | null { return this.children[0] ?? null; }
  get lastElementChild(): FakeElement | null { return this.children[this.children.length - 1] ?? null; }
  setAttribute(name: string, value: any) { this.attributes[name] = String(value); }
  getAttribute(name: string): string | null { return this.attributes[name] ?? null; }
  removeAttribute(name: string) { delete this.attributes[name]; }
  hasAttribute(name: string): boolean { return name in this.attributes; }
  querySelector(): any { return null; }
  querySelectorAll(): any[] { return []; }
  getElementsByClassName(): any[] { return []; }
  getElementsByTagName(): any[] { return []; }
  closest(): any { return null; }
  matches(): boolean { return false; }
  getBoundingClientRect() { return {x: 0, y: 0, left: 0, top: 0, right: 0, bottom: 0, width: 0, height: 0}; }
  getContext(): any { return null; }  // canvas — nothing renders under Node
  focus() {}
  blur() {}
  click() {}
  scrollIntoView() {}
  insertAdjacentElement(_where: string, el: any) { this.appendChild(el); return el; }
  insertAdjacentHTML() {}
  cloneNode(_deep?: boolean): FakeElement { return new FakeElement(this.tagName); }
  toString(): string { return `[FakeElement ${this.tagName}]`; }
}

class FakeHTMLCollection extends Array {}

class FakeEvent {
  type: string;
  target: any = null;
  defaultPrevented = false;
  constructor(type: string = '', _init?: any) { this.type = type; }
  preventDefault() { this.defaultPrevented = true; }
  stopPropagation() {}
  stopImmediatePropagation() {}
}

class FakeObserver {
  observe() {}
  unobserve() {}
  disconnect() {}
  takeRecords(): any[] { return []; }
}

/** Installs `document` + DOM classes on globalThis when running outside a browser.
 *  Call AFTER the dart2js bundle import: its environment detection must keep
 *  seeing `typeof document === 'undefined'` (the proven non-browser path). */
export function installDomStub(): void {
  const g = globalThis as any;
  if (typeof g.document !== 'undefined' && g.document?.createElement && !g.document.__dgStub)
    return;   // real DOM present — never override

  const document: any = {
    __dgStub: true,
    readyState: 'complete',
    currentScript: null,
    cookie: '',
    createElement: (tag: string) => {
      const el = new FakeElement(tag);
      el.ownerDocument = document;
      return el;
    },
    createElementNS: (ns: string, tag: string) => {
      const el = new FakeElement(tag);
      el.namespaceURI = ns;
      el.ownerDocument = document;
      return el;
    },
    createTextNode: (text: string) => new FakeText(text),
    createDocumentFragment: () => new FakeDocumentFragment(),
    createEvent: () => new FakeEvent(),
    querySelector: () => null,
    querySelectorAll: () => [],
    getElementById: () => null,
    getElementsByClassName: () => [],
    getElementsByTagName: () => [],
    addEventListener: () => {},
    removeEventListener: () => {},
    dispatchEvent: () => true,
    hasFocus: () => false,
    execCommand: () => false,
  };
  document.documentElement = document.createElement('html');
  document.head = document.createElement('head');
  document.body = document.createElement('body');
  document.defaultView = g.window;

  g.document = document;
  g.Node = FakeNode;
  g.Element = FakeElement;
  g.HTMLElement = FakeElement;
  g.HTMLCollection = FakeHTMLCollection;
  g.Text = FakeText;
  g.DocumentFragment = FakeDocumentFragment;
  g.Event = FakeEvent;
  g.CustomEvent = FakeEvent;
  g.KeyboardEvent = FakeEvent;
  g.MouseEvent = FakeEvent;
  g.MutationObserver = FakeObserver;
  g.ResizeObserver = FakeObserver;
  g.IntersectionObserver = FakeObserver;
  // canvas API classes: js-api extends their prototypes at import time
  g.CanvasRenderingContext2D ??= class CanvasRenderingContext2D {};
  g.CanvasGradient ??= class CanvasGradient {};
  g.CanvasPattern ??= class CanvasPattern {};
  g.Image ??= class Image extends FakeElement { constructor() { super('img'); } };
  g.Path2D ??= class Path2D {};
  for (const tag of ['Div', 'Span', 'Input', 'TextArea', 'Select', 'Button', 'Anchor', 'Image',
    'Canvas', 'Table', 'TableRow', 'TableCell', 'Label', 'Paragraph', 'Heading', 'UList', 'LI', 'IFrame'])
    g[`HTML${tag}Element`] ??= FakeElement;
  g.requestAnimationFrame ??= (f: Function) => setTimeout(f, 0);
  g.cancelAnimationFrame ??= (id: any) => clearTimeout(id);
  g.getComputedStyle ??= () => makeStyle();
}

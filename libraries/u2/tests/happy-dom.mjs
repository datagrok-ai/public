import {register} from 'node:module';
import assert from 'node:assert';
import strict from 'node:assert/strict';
import {realpathSync} from 'node:fs';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {isDeepStrictEqual} from 'node:util';
import {GlobalWindow, PropertySymbol} from 'happy-dom';

const window = new GlobalWindow({url: 'http://localhost/', width: 1024, height: 768, console: globalThis.console});
for (const key of [...Object.getOwnPropertyNames(window), ...Object.getOwnPropertySymbols(window)]) {
  if (['constructor', 'undefined', 'NaN', 'global', 'globalThis'].includes(key))
    continue;
  const d = Object.getOwnPropertyDescriptor(window, key);
  if (d.value === window)
    window[key] = d.value = globalThis;
  const own = Object.getOwnPropertyDescriptor(globalThis, key)?.value;
  if (typeof key === 'symbol' || own === undefined || own !== d.value)
    Object.defineProperty(globalThis, key, {...d, configurable: true});
}
globalThis.document[PropertySymbol.defaultView] = globalThis;

// the tests import their fire/flush/resetDom from dom-shim.js; on happy-dom they get this module's
const SHIM = pathToFileURL(realpathSync(fileURLToPath(new URL('./dom-shim.js', import.meta.url)))).href;
register('data:text/javascript,' + encodeURIComponent(`
export async function resolve(specifier, context, next) {
  const r = await next(specifier, context);
  return r.url === ${JSON.stringify(SHIM)} ? {...r, url: ${JSON.stringify(import.meta.url)}} : r;
}`));

globalThis.Option = function Option(text = '', value, defaultSelected = false, selected = false) {
  const option = document.createElement('option');
  option.textContent = text;
  if (value !== undefined)
    option.value = value;
  if (defaultSelected)
    option.setAttribute('selected', '');
  if (selected)
    option.selected = true;
  return option;
};
globalThis.Option.prototype = HTMLOptionElement.prototype;

// window.event, which browsers set during a dispatch: u2's Input tells a user's edit by it
let target = Node.prototype;
while (!Object.hasOwn(target, 'dispatchEvent'))
  target = Object.getPrototypeOf(target);
const dispatch = target.dispatchEvent;
target.dispatchEvent = function(event) {
  const outer = globalThis.event;
  globalThis.event = event;
  try {
    return dispatch.call(this, event);
  }
  finally {
    globalThis.event = outer;
  }
};

// browsers refuse focus by the `disabled` attribute; happy-dom reads the property, which focusableWhenDisabled overrides
const focus = HTMLElement.prototype.focus;
HTMLElement.prototype.focus = function(...args) {
  const own = Object.getOwnPropertyDescriptor(this, 'disabled');
  if (!own)
    return focus.apply(this, args);
  delete this.disabled;
  try {
    return focus.apply(this, args);
  }
  finally {
    Object.defineProperty(this, 'disabled', own);
  }
};

// happy-dom has canvases but no 2D context class; js-api extends its prototype on load
globalThis.CanvasRenderingContext2D ??= class CanvasRenderingContext2D {};

const rows = Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'rows');
Object.defineProperty(HTMLTextAreaElement.prototype, 'rows', {...rows, get() {
  return Number(rows.get.call(this)) || 2;
}});

// node's AssertionError diffs a full inspection of both sides, which over happy-dom nodes runs until an OOM kill
const tag = (x) => x instanceof Node ? x[Symbol.for('nodejs.util.inspect.custom')]() :
  Array.isArray(x) ? `[${x.map(tag).join(', ')}]` : typeof x === 'string' ? JSON.stringify(x) : String(x);
const hasNode = (x) => x instanceof Node || (Array.isArray(x) && x.some(hasNode));
for (const target of [assert, strict]) {
  for (const [name, same, sign] of [['equal', Object.is, '!=='], ['strictEqual', Object.is, '!=='],
    ['deepEqual', isDeepStrictEqual, '!=='], ['deepStrictEqual', isDeepStrictEqual, '!=='],
    ['notEqual', (a, b) => !Object.is(a, b), '==='], ['notStrictEqual', (a, b) => !Object.is(a, b), '===']]) {
    const original = target[name];
    target[name] = function guarded(actual, expected, message) {
      if ((hasNode(actual) || hasNode(expected)) && !same(actual, expected)) {
        throw new assert.AssertionError({message: `${message ? message + ': ' : ''}${tag(actual)} ${sign} ${tag(expected)}`,
          actual, expected, operator: name, stackStartFn: guarded});
      }
      try {
        return original.apply(this, arguments);
      }
      catch (e) {
        Error.captureStackTrace(e, guarded);
        throw e;
      }
    };
  }
}

Object.defineProperty(Node.prototype, Symbol.for('nodejs.util.inspect.custom'), {
  configurable: true,
  value() {
    if (!(this instanceof Element))
      return this.nodeName;
    const cls = this.getAttribute('class');
    return `<${this.tagName.toLowerCase()}${cls ? ` class="${cls}"` : ''}>`;
  },
});

const EVENTS = {keydown: 'KeyboardEvent', keyup: 'KeyboardEvent', click: 'MouseEvent', dblclick: 'MouseEvent'};

export function fire(el, type, init = {}) {
  const pointer = type.startsWith('pointer') || type === 'mousedown' || type === 'contextmenu';
  const event = new globalThis[pointer ? 'PointerEvent' : EVENTS[type] ?? 'Event'](type,
    {bubbles: true, cancelable: true, ...init});
  for (const [key, value] of Object.entries(init)) {
    if (key !== 'target' && event[key] !== value)
      Object.defineProperty(event, key, {value, configurable: true, enumerable: true});
  }
  return el.dispatchEvent(event);
}

// happy-dom runs animation frames as immediates, which a timer alone may overtake
export function flush() {
  return new Promise((resolve) => setTimeout(() => setImmediate(resolve), 1));
}

export function resetDom() {
  document.body.replaceChildren();
  document.activeElement?.blur?.();
  if (location.search)
    history.replaceState(null, '', location.pathname);
}

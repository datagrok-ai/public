var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
  function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
  return new (P || (P = Promise))(function (resolve, reject) {
    function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
    function rejected(value) { try { step(generator['throw'](value)); } catch (e) { reject(e); } }
    function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
    step((generator = generator.apply(thisArg, _arguments || [])).next());
  });
};
import { Balloon } from './widgets';
import { toJs } from './wrappers';
let api = window;
export class Utils {
  /** @param {Iterable} iterable*/
  static firstOrNull(iterable) {
    let first = iterable[Symbol.iterator]().next();
    return first.done ? null : first.value;
  }
  /** Returns an 'identity' array where the element in idx-th position is equals to idx. */
  static identity(length) {
    let res = new Uint32Array(length);
    for (let i = 0; i < length; i++)
      res[i] = i;
    return res;
  }
}
// export class PropProxy {
//     constructor(d) {
//         this.d = d;
//         return new Proxy({}, {
//             get: function(target, prop) { return DG.toJs(api.grok_PropMixin_Get(d, prop)); },
//             set: function(target, prop, value) { api.grok_PropMixin_Set(d, prop, DG.toDart(value)); }
//         })
//     }
// }
export function _toIterable(o) {
  let iterable = {};
  // @ts-ignore
  iterable[Symbol.iterator] = () => _getIterator(o);
  // @ts-ignore
  return iterable;
}
export function _getIterator(d) {
  let iterator = api.grok_Iterable_Get_Iterator(d);
  return {
    next: function () {
      return api.grok_Iterator_MoveNext(iterator) ?
        { value: toJs(api.grok_Iterator_Current(iterator)), done: false } :
        { done: true };
    }
  };
}
export function _isDartium() {
  return Array
    .from(document.getElementsByTagName('script'))
    .some((s) => {
      let a = s.getAttribute('src');
      if (a == null)
        return null;
      return a.includes('dart.js');
    });
}
export function _jsThen(promise, f) {
  promise.then(f);
}
export function _toJson(x) {
  return x === null ? null : JSON.stringify(x);
}
export function* range(length) {
  for (let i = 0; i < length; i++)
    yield i;
}
/** Returns an 'identity' array where the element in idx-th position is equals to idx. */
export function identity(length) {
  let res = new Uint32Array(length);
  for (let i = 0; i < length; i++)
    res[i] = i;
  return res;
}
/*window.onerror = function (message, url, lineNumber, columnNumber, errorObject) {
    return api.grok_Error(message, url, lineNumber, columnNumber, errorObject);
};*/
/** Times the execution of function f
 * @param {string} name - a label for the execution time to display
 * @param {Function} f - function with no parameters that will get measured
 * @returns {value} - a value which f returns
 * */
export function time(name, f) {
  let start = new Date();
  let result = f();
  let stop = new Date();
  // @ts-ignore
  console.log(`${name}: ${stop - start} ms`);
  // @ts-ignore
  new Balloon().info(`${name}: ${stop - start} ms`);
  return result;
}
/** Times the execution of asyncronous function f
 * @async
 * @param {string} name - a label for the execution time to display
 * @param {Function} f - async function with no parameters that will get measured
 * @returns {Promise<value>} - a promise for the value which f returns
 * */
export function timeAsync(name, f) {
  return __awaiter(this, void 0, void 0, function* () {
    let start = new Date();
    let result = yield f();
    let stop = new Date();
    // @ts-ignore
    console.log(`${name}: ${stop - start} ms`);
    // @ts-ignore
    new Balloon().info(`${name}: ${stop - start} ms`);
    return result;
  });
}
export function _identityInt32(length) {
  let values = new Int32Array(length);
  for (let i = 0; i < length; i++)
    values[i] = i;
  return values;
}
/**
 * Inspired by https://github.com/Yomguithereal/mnemonist/blob/master/lru-cache.js
 * */
export class LruCache {
  constructor() {
    this.capacity = 100;
    this.forward = new Uint16Array(this.capacity);
    this.backward = new Uint16Array(this.capacity);
    this.V = new Array(this.capacity);
    this.K = new Array(this.capacity);
    this.size = 0;
    this.head = 0;
    this.tail = 0;
    this.items = {};
    this.onItemEvicted = null;
  }
  /**
     * Splays a value on top.
     * @param {number} pointer - Pointer of the value to splay on top.
     * @return {LruCache}
     */
  splayOnTop(pointer) {
    let oldHead = this.head;
    if (this.head === pointer)
      return this;
    let previous = this.backward[pointer];
    let next = this.forward[pointer];
    if (this.tail === pointer)
      this.tail = previous;
    else
      this.backward[next] = previous;
    this.forward[previous] = next;
    this.backward[oldHead] = pointer;
    this.head = pointer;
    this.forward[pointer] = oldHead;
    return this;
  }
  /**
     * Checks whether the key exists in the cache.
     *
     * @param  {any} key   - Key.
     * @return {boolean}
     */
  has(key) {
    return key in this.items;
  }
  /**
     * Sets the value for the given key in the cache.
     *
     * @param  {any} key   - Key.
     * @param  {any} value - Value.
     * @return {undefined}
     */
  set(key, value) {
    // The key already exists, we just need to update the value and splay on top
    // @ts-ignore
    let pointer = this.items[key];
    if (typeof pointer !== 'undefined') {
      this.splayOnTop(pointer);
      this.V[pointer] = value;
      return;
    }
    // The cache is not yet full
    if (this.size < this.capacity) {
      pointer = this.size++;
    }
    // Cache is full, we need to drop the last value
    else {
      pointer = this.tail;
      this.tail = this.backward[pointer];
      if (this.onItemEvicted != null)
        this.onItemEvicted(this.V[pointer]);
      // @ts-ignore
      delete this.items[this.K[pointer]];
    }
    // Storing key & value
    // @ts-ignore
    this.items[key] = pointer;
    this.K[pointer] = key;
    this.V[pointer] = value;
    // Moving the item at the front of the list
    this.forward[pointer] = this.head;
    this.backward[this.head] = pointer;
    this.head = pointer;
  }
  /**
     * Gets the value attached to the given key, and makes it the most recently used item.
     *
     * @param  {any} key   - Key.
     * @return {any}
     */
  get(key) {
    // @ts-ignore
    let pointer = this.items[key];
    if (typeof pointer === 'undefined')
      return;
    this.splayOnTop(pointer);
    return this.V[pointer];
  }
  /**
     * Returns the value with the specified key, if it already exists in the cache,
     * or creates a new one by calling the provided function.
     *
     * @param  {any} key   - Key.
     * @param  {Function} createFromKey - Function to create a new item.
     * @return {any}
     */
  getOrCreate(key, createFromKey) {
    let value = this.get(key);
    if (typeof value !== 'undefined')
      return value;
    else {
      let item = createFromKey(key);
      this.set(key, item);
      return item;
    }
  }
}
/**
 * @param {HTMLElement} element
 * @param {string | ElementOptions | null} options
 * @returns {HTMLElement}
 * */
export function _options(element, options) {
  if (options == null)
    return element;
  if (typeof options === 'string')
    element.className += ` ${options.replace(/,/g, ' ')}`;
  if (options.id != null)
    element.id = options.id;
  if (options.classes != null)
    element.className += ` ${options.classes}`;
  if (options.style != null)
    Object.assign(element.style, options.style);
  return element;
}

import {Balloon} from "./widgets";
import * as rxjs from 'rxjs';
import {toJs} from "./wrappers";

export function _toIterable(o) {
    let iterable = {};
    iterable[Symbol.iterator] = () => _getIterator(o);
    return iterable;
}

export function _getIterator(d) {
    let iterator = grok_Iterable_Get_Iterator(d);
    return {
        next: function() {
            return grok_Iterator_MoveNext(iterator) ?
                {value: toJs(grok_Iterator_Current(iterator)), done: false} :
                {done: true};
        }
    }
}

export function _isDartium() {
    return Array
        .from(document.getElementsByTagName("script"))
        .some((s) => s.getAttribute('src').includes('dart.js'));
}

export function _jsThen(promise, f) {
    promise.then(f);
}

export function _toJson(x) {
    return x === null ? null : JSON.stringify(x);
}

export function *range(length) {
    for (let i = 0; i < length; i++)
        yield i;
}

/** Returns an 'identity' array where the element in idx-th position is equals to idx. */
export function identity(length) {
    let res = new Array(length);
    for (let i = 0; i < length; i++)
        res[i] = i;
    return res;
}

/*window.onerror = function (message, url, lineNumber, columnNumber, errorObject) {
    return grok_Error(message, url, lineNumber, columnNumber, errorObject);
};*/

/** Times the execution of function f
 * @param {string} name
 * @param {Function} f - function with no parameters that will get measured
 * @returns {number} - milliseconds elapsed
 * */
export function time(name, f) {
    let start = new Date();
    let result = f();
    let stop = new Date();
    console.log(`${name}: ${stop - start}ms`);
    new Balloon().info(`${name}: ${stop - start}ms`);
    return result;
}

/** @returns {rxjs.Observable} */
export function _onSizeChanged(element) {
    if (_isDartium())
        return rxjs.empty();

    return rxjs.Observable.create(function(observer) {
        const resizeObserver = new ResizeObserver(observerEntries => {
            // trigger a new item on the stream when resizes happen
            for (const entry of observerEntries) {
                observer.next(entry);
            }
        });

        // start listening for resize events
        resizeObserver.observe(element);

        // cancel resize observer on cancelation
        return () => resizeObserver.disconnect();
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
    };


    /**
     * Checks whether the key exists in the cache.
     *
     * @param  {any} key   - Key.
     * @return {boolean}
     */
    has(key) {
        return key in this.items;
    };

    /**
     * Sets the value for the given key in the cache.
     *
     * @param  {any} key   - Key.
     * @param  {any} value - Value.
     * @return {undefined}
     */
    set(key, value) {

        // The key already exists, we just need to update the value and splay on top
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
            delete this.items[this.K[pointer]];
        }

        // Storing key & value
        this.items[key] = pointer;
        this.K[pointer] = key;
        this.V[pointer] = value;

        // Moving the item at the front of the list
        this.forward[pointer] = this.head;
        this.backward[this.head] = pointer;
        this.head = pointer;
    };

    /**
     * Gets the value attached to the given key, and makes it the most recently used item.
     *
     * @param  {any} key   - Key.
     * @return {any}
     */
    get(key) {
        let pointer = this.items[key];

        if (typeof pointer === 'undefined')
            return;

        this.splayOnTop(pointer);

        return this.V[pointer];
    };

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
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
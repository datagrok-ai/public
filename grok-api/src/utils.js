import {Balloon} from "./ui_classes";

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

let time = function(s, f) {
    let start = new Date();
    let result = f();
    let stop = new Date();
    console.log(`${s}: ${stop - start}ms`);
    Balloon.info(`${s}: ${stop - start}ms`);
    return result;
};

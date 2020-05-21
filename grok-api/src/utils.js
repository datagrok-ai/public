import {TYPE, TYPES_SCALAR} from "./const";
import * as ui from "./../ui.js";

export function _jsThen(promise, f) {
    promise.then(f);
}

export function _toJson(x) {
    return x === null ? null : JSON.stringify(x);
}

export function paramsToJs(params) {
    let result = [];
    for (let i = 0; i < params.length; i++) {
        let type = grok_GetType(params[i]);
        if (type !== null && !TYPES_SCALAR.has(type))
            result.push(_wrap(params[i]));
        else
            result.push(params[i]);
    }

    return result;
}

/*window.onerror = function (message, url, lineNumber, columnNumber, errorObject) {
    return grok_Error(message, url, lineNumber, columnNumber, errorObject);
};*/

let time = function(s, f) {
    let start = new Date();
    let result = f();
    let stop = new Date();
    console.log(`${s}: ${stop - start}ms`);
    ui.Balloon.info(`${s}: ${stop - start}ms`);
    return result;
};

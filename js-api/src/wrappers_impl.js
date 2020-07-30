import {Property} from "./entities";
import {TYPE, TYPES_SCALAR} from "./const";

/** Converts list of Dart objects to JavaScript objects by calling {@link toJs}
 * @param {object[]} params
 * @returns {object[]} - list of JavaScript objects */
export function paramsToJs(params) {
    let result = [];
    for (let i = 0; i < params.length; i++) {
        let type = grok_GetType(params[i]);
        if (type !== null && (!TYPES_SCALAR.has(type) ||
            type === TYPE.COLUMN_LIST ||
            type === TYPE.LIST ||
            type === TYPE.DATA_FRAME_LIST))
            result.push(toJs(params[i]));
        else
            result.push(params[i]);
    }

    return result;
}


/**
 * Instantiates the corresponding JS handler for the Dart object [d]. See also {@link toDart}
 * @param d - Dart handle
 * @param {boolean} check -
 * @returns JavaScript wrapper for the Dart object
 * */
export function toJs(d, check = false) {
    let wrapper = grok_GetWrapper(d);
    if (wrapper != null)
        return wrapper;
    let type = grok_GetType(d);
    if (type === TYPE.COLUMN_LIST || type === TYPE.DATA_FRAME_LIST || type === TYPE.LIST) {
        let p = [];
        for(let ii = 0; ii < d.length; ii++)
            p.push(toJs(d[ii]));
        return p;
    } else if (type === TYPE.PROPERTY)
        return new Property(d);
    if (check)
        throw `Not supported type: ${type}`;

    return d;
}

/** Extracts a Dart handle from the JavaScript wrapper. See also {@link toJs} */
export function toDart(x) {
    if (x === undefined || x === null)
        return x;
    return (typeof x.d !== 'undefined') ? x.d : x;
}



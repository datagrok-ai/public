/**
 * Internal conversion utilities for Dart-JavaScript interop.
 *
 * **Note on underscore prefix**: Functions in this module are prefixed with `_` to indicate
 * they are internal utilities intended for use within the API implementation.
 * While they are exported (and available to package authors who need low-level access),
 * they are not part of the stable public API and may change without notice.
 *
 * For most use cases, prefer the higher-level APIs in the main modules.
 *
 * @module utils_convert
 * @internal
 */

import {toJs} from "./wrappers";
import {IDartApi} from "./api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/**
 * Converts a value to JSON string, handling null.
 * @internal
 */
export function _toJson(x: any) {
    return x === null ? null : JSON.stringify(x);
}

/**
 * Wraps a Dart iterable as a JavaScript Iterable.
 * Allows using Dart collections in for...of loops.
 * @param dart - Dart iterable object
 * @returns JavaScript Iterable wrapper
 * @internal
 */
export function _toIterable(dart: any): Iterable<any> {
    let iterable = {};
    // @ts-ignore
    iterable[Symbol.iterator] = () => _getIterator(dart);
    // @ts-ignore
    return iterable;
}

/**
 * Creates a JavaScript iterator from a Dart iterator.
 * @param dart - Dart iterator object
 * @returns JavaScript iterator with next() method
 * @internal
 */
export function _getIterator(dart: any) {
    let iterator = api.grok_Iterable_Get_Iterator(dart);
    return {
        next: function () {
            return api.grok_Iterator_MoveNext(iterator) ?
                {value: toJs(api.grok_Iterator_Current(iterator)), done: false} :
                {done: true};
        }
    };
}


/**
 * Converts entity properties between JavaScript and Dart.
 * See also: {@link include}
 */
export function _propsToDart(s: string, cls: string): string {
    const jsToDart: { [indes:string] : {[index: string]: string} } = {
        'Group' : {
            'adminMemberships': 'parents.parent',
            'memberships': 'parents.parent',
            'members': 'children.child',
            'adminMembers': 'children.child',
        },
        'Project': {
            'children': 'relations.entity',
            'links': 'relations.entity'
        },
        'Function': {
            'inputs': 'params',
            'outputs': 'params'
        }
    };

    let propsMap = jsToDart[cls];
    if (!propsMap)
        return s;
    let res = '';
    if (s === res) return res;
    let ents = s.split(',');
    for (let ent of ents) {
        let props = ent.trim();

        while (props) {
            let idx = props.indexOf('.');
            let match = propsMap[props];
            if (match) res += match;
            else {
                let p = (idx === -1) ? props : props.slice(0, idx);
                res += propsMap[p] || p;
            }
            if (idx === -1) props = '';
            else {
                props = props.slice(idx + 1);
                res += '.';
            }
        }

        res += ',';
    }
    return res;
}
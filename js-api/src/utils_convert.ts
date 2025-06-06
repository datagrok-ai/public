import {toJs} from "./wrappers";
import {IDartApi} from "./api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export function _toJson(x: any) {
    return x === null ? null : JSON.stringify(x);
}

export function _toIterable(dart: any): Iterable<any> {
    let iterable = {};
    // @ts-ignore
    iterable[Symbol.iterator] = () => _getIterator(dart);
    // @ts-ignore
    return iterable;
}

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
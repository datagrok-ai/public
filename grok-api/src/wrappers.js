
export function paramsToJs(params) { return impl._paramsToJs(params);}

export function toJs(d) { return impl._wrap(d, false); }

export function toDart(x) { return impl._toDart(x);}

export let impl = {
    _toDart: () => {},
    _paramsToJs: () => {},
    _wrap: () => {}
}
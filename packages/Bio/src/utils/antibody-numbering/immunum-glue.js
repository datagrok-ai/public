/* eslint-disable */
// @ts-nocheck
//
// Browser/worker-friendly port of node_modules/immunum/immunum.js (v1.1.0).
// The upstream file uses require('fs').readFileSync to load the WASM at top
// level, which breaks in the browser. This file is byte-identical up to the
// wasm-bindgen glue; the only change is that the WASM bytes are supplied via
// an explicit initImmunum() call instead of being read from disk at import
// time. See README.md in the immunum package for the upstream source.

let wasm = null;

export class Annotator {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        AnnotatorFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_annotator_free(ptr, 0);
    }
    constructor(chains, scheme, min_confidence) {
        const ptr0 = passArrayJsValueToWasm0(chains, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(scheme, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.annotator_new(ptr0, len0, ptr1, len1,
            isLikeNone(min_confidence) ? 0x100000001 : Math.fround(min_confidence));
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0] >>> 0;
        AnnotatorFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    number(sequence) {
        const ptr0 = passStringToWasm0(sequence, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        return wasm.annotator_number(this.__wbg_ptr, ptr0, len0);
    }
    segment(sequence) {
        const ptr0 = passStringToWasm0(sequence, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        return wasm.annotator_segment(this.__wbg_ptr, ptr0, len0);
    }
}
if (typeof Symbol !== 'undefined' && Symbol.dispose) Annotator.prototype[Symbol.dispose] = Annotator.prototype.free;

const AnnotatorFinalization = (typeof FinalizationRegistry === 'undefined')
    ? {register: () => {}, unregister: () => {}}
    : new FinalizationRegistry((ptr) => wasm.__wbg_annotator_free(ptr >>> 0, 1));

function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_debug_string_5398f5bb970e0daa: function(arg0, arg1) {
            const ret = debugString(arg1);
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg___wbindgen_string_get_395e606bd0ee4427: function(arg0, arg1) {
            const obj = arg1;
            const ret = typeof(obj) === 'string' ? obj : undefined;
            const ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg___wbindgen_throw_6ddd609b62940d55: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbg_new_49d5571bd3f0c4d4: function() {
            return new Map();
        },
        __wbg_new_ab79df5bd7c26067: function() {
            return new Object();
        },
        __wbg_set_7eaa4f96924fd6b3: function() {
            return handleError(function(arg0, arg1, arg2) {
                return Reflect.set(arg0, arg1, arg2);
            }, arguments);
        },
        __wbg_set_bf7251625df30a02: function(arg0, arg1, arg2) {
            return arg0.set(arg1, arg2);
        },
        __wbindgen_cast_0000000000000001: function(arg0) {
            return arg0;
        },
        __wbindgen_cast_0000000000000002: function(arg0, arg1) {
            return getStringFromWasm0(arg0, arg1);
        },
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        './immunum_bg.js': import0,
    };
}

function addToExternrefTable0(obj) {
    const idx = wasm.__externref_table_alloc();
    wasm.__wbindgen_externrefs.set(idx, obj);
    return idx;
}

function debugString(val) {
    const type = typeof val;
    if (type == 'number' || type == 'boolean' || val == null) return `${val}`;
    if (type == 'string') return `"${val}"`;
    if (type == 'symbol') {
        const description = val.description;
        return description == null ? 'Symbol' : `Symbol(${description})`;
    }
    if (type == 'function') {
        const name = val.name;
        return (typeof name == 'string' && name.length > 0) ? `Function(${name})` : 'Function';
    }
    if (Array.isArray(val)) {
        const length = val.length;
        let debug = '[';
        if (length > 0) debug += debugString(val[0]);
        for (let i = 1; i < length; i++) debug += ', ' + debugString(val[i]);
        debug += ']';
        return debug;
    }
    const builtInMatches = /\[object ([^\]]+)\]/.exec(toString.call(val));
    let className;
    if (builtInMatches && builtInMatches.length > 1) {
        className = builtInMatches[1];
    } else {
        return toString.call(val);
    }
    if (className == 'Object') {
        try {
            return 'Object(' + JSON.stringify(val) + ')';
        } catch (_) {
            return 'Object';
        }
    }
    if (val instanceof Error) return `${val.name}: ${val.message}\n${val.stack}`;
    return className;
}

let cachedDataViewMemory0 = null;
function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null ||
        cachedDataViewMemory0.buffer.detached === true ||
        (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

function getStringFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return decodeText(ptr, len);
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0)
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    return cachedUint8ArrayMemory0;
}

function handleError(f, args) {
    try {
        return f.apply(this, args);
    } catch (e) {
        const idx = addToExternrefTable0(e);
        wasm.__wbindgen_exn_store(idx);
    }
}

function isLikeNone(x) {
    return x === undefined || x === null;
}

function passArrayJsValueToWasm0(array, malloc) {
    const ptr = malloc(array.length * 4, 4) >>> 0;
    for (let i = 0; i < array.length; i++) {
        const add = addToExternrefTable0(array[i]);
        getDataViewMemory0().setUint32(ptr + 4 * i, add, true);
    }
    WASM_VECTOR_LEN = array.length;
    return ptr;
}

function passStringToWasm0(arg, malloc, realloc) {
    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;
    const mem = getUint8ArrayMemory0();

    let offset = 0;
    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }
    if (offset !== len) {
        if (offset !== 0) arg = arg.slice(offset);
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = cachedTextEncoder.encodeInto(arg, view);
        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_externrefs.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}

const cachedTextDecoder = new TextDecoder('utf-8', {ignoreBOM: true, fatal: true});
cachedTextDecoder.decode();
function decodeText(ptr, len) {
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

const cachedTextEncoder = new TextEncoder();
if (!('encodeInto' in cachedTextEncoder)) {
    cachedTextEncoder.encodeInto = function(arg, view) {
        const buf = cachedTextEncoder.encode(arg);
        view.set(buf);
        return {read: arg.length, written: buf.length};
    };
}

let WASM_VECTOR_LEN = 0;

/** Instantiate the immunum WASM module from raw bytes or an ArrayBuffer.
 *  Must be called once before constructing any Annotator. Safe to call multiple
 *  times — subsequent calls are no-ops. */
export async function initImmunum(wasmSource) {
    if (wasm) return;
    const imports = __wbg_get_imports();
    let inst;
    if (wasmSource instanceof WebAssembly.Module) {
        inst = await WebAssembly.instantiate(wasmSource, imports);
    } else {
        const result = await WebAssembly.instantiate(wasmSource, imports);
        inst = result.instance;
    }
    wasm = inst.exports;
    wasm.__wbindgen_start();
}

export function isImmunumReady() {
    return wasm !== null;
}

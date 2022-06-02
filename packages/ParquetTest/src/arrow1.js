
let wasm;

const cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });

cachedTextDecoder.decode();

let cachegetUint8Memory0 = null;
function getUint8Memory0() {
    if (cachegetUint8Memory0 === null || cachegetUint8Memory0.buffer !== wasm.memory.buffer) {
        cachegetUint8Memory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachegetUint8Memory0;
}

function getStringFromWasm0(ptr, len) {
    return cachedTextDecoder.decode(getUint8Memory0().subarray(ptr, ptr + len));
}

const heap = new Array(32).fill(undefined);

heap.push(undefined, null, true, false);

let heap_next = heap.length;

function addHeapObject(obj) {
    if (heap_next === heap.length) heap.push(heap.length + 1);
    const idx = heap_next;
    heap_next = heap[idx];

    heap[idx] = obj;
    return idx;
}

function getObject(idx) { return heap[idx]; }

function dropObject(idx) {
    if (idx < 36) return;
    heap[idx] = heap_next;
    heap_next = idx;
}

function takeObject(idx) {
    const ret = getObject(idx);
    dropObject(idx);
    return ret;
}

let WASM_VECTOR_LEN = 0;

const cachedTextEncoder = new TextEncoder('utf-8');

const encodeString = (typeof cachedTextEncoder.encodeInto === 'function'
    ? function (arg, view) {
    return cachedTextEncoder.encodeInto(arg, view);
}
    : function (arg, view) {
    const buf = cachedTextEncoder.encode(arg);
    view.set(buf);
    return {
        read: arg.length,
        written: buf.length
    };
});

function passStringToWasm0(arg, malloc, realloc) {

    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length);
        getUint8Memory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len);

    const mem = getUint8Memory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }

    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3);
        const view = getUint8Memory0().subarray(ptr + offset, ptr + len);
        const ret = encodeString(arg, view);

        offset += ret.written;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function passArray8ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 1);
    getUint8Memory0().set(arg, ptr / 1);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

let cachegetInt32Memory0 = null;
function getInt32Memory0() {
    if (cachegetInt32Memory0 === null || cachegetInt32Memory0.buffer !== wasm.memory.buffer) {
        cachegetInt32Memory0 = new Int32Array(wasm.memory.buffer);
    }
    return cachegetInt32Memory0;
}
/**
* Read a Parquet file into Arrow data using the [`arrow`](https://crates.io/crates/arrow) and
* [`parquet`](https://crates.io/crates/parquet) Rust crates.
*
* Example:
*
* ```js
* import { tableFromIPC } from "apache-arrow";
* // Edit the `parquet-wasm` import as necessary
* import { readParquet } from "parquet-wasm/node";
*
* const resp = await fetch("https://example.com/file.parquet");
* const parquetUint8Array = new Uint8Array(await resp.arrayBuffer());
* const arrowUint8Array = readParquet(parquetUint8Array);
* const arrowTable = tableFromIPC(arrowUint8Array);
* ```
*
* @param parquet_file Uint8Array containing Parquet data
* @returns Uint8Array containing Arrow data in [IPC Stream format](https://arrow.apache.org/docs/format/Columnar.html#ipc-streaming-format). To parse this into an Arrow table, pass to `tableFromIPC` in the Arrow JS bindings.
* @param {Uint8Array} parquet_file
* @returns {Uint8Array}
*/
export function readParquet(parquet_file) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(parquet_file, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.readParquet(retptr, ptr0, len0);
        var r0 = getInt32Memory0()[retptr / 4 + 0];
        var r1 = getInt32Memory0()[retptr / 4 + 1];
        var r2 = getInt32Memory0()[retptr / 4 + 2];
        if (r2) {
            throw takeObject(r1);
        }
        return takeObject(r0);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

function _assertClass(instance, klass) {
    if (!(instance instanceof klass)) {
        throw new Error(`expected instance of ${klass.name}`);
    }
    return instance.ptr;
}
/**
* Write Arrow data to a Parquet file using the [`arrow`](https://crates.io/crates/arrow) and
* [`parquet`](https://crates.io/crates/parquet) Rust crates.
*
* For example, to create a Parquet file with Snappy compression:
*
* ```js
* import { tableToIPC } from "apache-arrow";
* // Edit the `parquet-wasm` import as necessary
* import { WriterPropertiesBuilder, Compression, writeParquet } from "parquet-wasm/node";
*
* // Given an existing arrow table under `table`
* const arrowUint8Array = tableToIPC(table, "file");
* const writerProperties = new WriterPropertiesBuilder()
*   .setCompression(Compression.SNAPPY)
*   .build();
* const parquetUint8Array = writeParquet(arrowUint8Array, writerProperties);
* ```
*
* @param arrow_file Uint8Array containing Arrow data in [IPC Stream format](https://arrow.apache.org/docs/format/Columnar.html#ipc-streaming-format). If you have an Arrow table in JS, call `tableToIPC(table)` in the JS bindings and pass the result here.
* @param writer_properties Configuration for writing to Parquet. Use the {@linkcode WriterPropertiesBuilder} to build a writing configuration, then call `.build()` to create an immutable writer properties to pass in here.
* @returns Uint8Array containing written Parquet data.
* @param {Uint8Array} arrow_file
* @param {WriterProperties} writer_properties
* @returns {Uint8Array}
*/
export function writeParquet(arrow_file, writer_properties) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(arrow_file, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        _assertClass(writer_properties, WriterProperties);
        var ptr1 = writer_properties.ptr;
        writer_properties.ptr = 0;
        wasm.writeParquet(retptr, ptr0, len0, ptr1);
        var r0 = getInt32Memory0()[retptr / 4 + 0];
        var r1 = getInt32Memory0()[retptr / 4 + 1];
        var r2 = getInt32Memory0()[retptr / 4 + 2];
        if (r2) {
            throw takeObject(r1);
        }
        return takeObject(r0);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
* Supported compression algorithms.
*/
export const Compression = Object.freeze({ UNCOMPRESSED:0,"0":"UNCOMPRESSED",SNAPPY:1,"1":"SNAPPY",GZIP:2,"2":"GZIP",BROTLI:3,"3":"BROTLI",LZ4:4,"4":"LZ4",ZSTD:5,"5":"ZSTD", });
/**
* Encodings supported by Parquet.
* Not all encodings are valid for all types. These enums are also used to specify the
* encoding of definition and repetition levels.
*/
export const Encoding = Object.freeze({
/**
* Default byte encoding.
* - BOOLEAN - 1 bit per value, 0 is false; 1 is true.
* - INT32 - 4 bytes per value, stored as little-endian.
* - INT64 - 8 bytes per value, stored as little-endian.
* - FLOAT - 4 bytes per value, stored as little-endian.
* - DOUBLE - 8 bytes per value, stored as little-endian.
* - BYTE_ARRAY - 4 byte length stored as little endian, followed by bytes.
* - FIXED_LEN_BYTE_ARRAY - just the bytes are stored.
*/
PLAIN:0,"0":"PLAIN",
/**
* **Deprecated** dictionary encoding.
*
* The values in the dictionary are encoded using PLAIN encoding.
* Since it is deprecated, RLE_DICTIONARY encoding is used for a data page, and
* PLAIN encoding is used for dictionary page.
*/
PLAIN_DICTIONARY:1,"1":"PLAIN_DICTIONARY",
/**
* Group packed run length encoding.
*
* Usable for definition/repetition levels encoding and boolean values.
*/
RLE:2,"2":"RLE",
/**
* Bit packed encoding.
*
* This can only be used if the data has a known max width.
* Usable for definition/repetition levels encoding.
*/
BIT_PACKED:3,"3":"BIT_PACKED",
/**
* Delta encoding for integers, either INT32 or INT64.
*
* Works best on sorted data.
*/
DELTA_BINARY_PACKED:4,"4":"DELTA_BINARY_PACKED",
/**
* Encoding for byte arrays to separate the length values and the data.
*
* The lengths are encoded using DELTA_BINARY_PACKED encoding.
*/
DELTA_LENGTH_BYTE_ARRAY:5,"5":"DELTA_LENGTH_BYTE_ARRAY",
/**
* Incremental encoding for byte arrays.
*
* Prefix lengths are encoded using DELTA_BINARY_PACKED encoding.
* Suffixes are stored using DELTA_LENGTH_BYTE_ARRAY encoding.
*/
DELTA_BYTE_ARRAY:6,"6":"DELTA_BYTE_ARRAY",
/**
* Dictionary encoding.
*
* The ids are encoded using the RLE encoding.
*/
RLE_DICTIONARY:7,"7":"RLE_DICTIONARY",
/**
* Encoding for floating-point data.
*
* K byte-streams are created where K is the size in bytes of the data type.
* The individual bytes of an FP value are scattered to the corresponding stream and
* the streams are concatenated.
* This itself does not reduce the size of the data but can lead to better compression
* afterwards.
*/
BYTE_STREAM_SPLIT:8,"8":"BYTE_STREAM_SPLIT", });
/**
* The Parquet version to use when writing
*/
export const WriterVersion = Object.freeze({ V1:0,"0":"V1",V2:1,"1":"V2", });
/**
* Immutable struct to hold writing configuration for `writeParquet`.
*
* Use {@linkcode WriterPropertiesBuilder} to create a configuration, then call {@linkcode
* WriterPropertiesBuilder.build} to create an instance of `WriterProperties`.
*/
export class WriterProperties {

    static __wrap(ptr) {
        const obj = Object.create(WriterProperties.prototype);
        obj.ptr = ptr;

        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.ptr;
        this.ptr = 0;

        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_writerproperties_free(ptr);
    }
}
/**
* Builder to create a writing configuration for `writeParquet`
*
* Call {@linkcode build} on the finished builder to create an immputable {@linkcode WriterProperties} to pass to `writeParquet`
*/
export class WriterPropertiesBuilder {

    static __wrap(ptr) {
        const obj = Object.create(WriterPropertiesBuilder.prototype);
        obj.ptr = ptr;

        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.ptr;
        this.ptr = 0;

        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_writerpropertiesbuilder_free(ptr);
    }
    /**
    * Returns default state of the builder.
    */
    constructor() {
        const ret = wasm.writerpropertiesbuilder_new();
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Finalizes the configuration and returns immutable writer properties struct.
    * @returns {WriterProperties}
    */
    build() {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_build(ptr);
        return WriterProperties.__wrap(ret);
    }
    /**
    * Sets writer version.
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setWriterVersion(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setWriterVersion(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets data page size limit.
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setDataPagesizeLimit(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setDataPagesizeLimit(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets dictionary page size limit.
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setDictionaryPagesizeLimit(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setDictionaryPagesizeLimit(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets write batch size.
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setWriteBatchSize(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setWriteBatchSize(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets maximum number of rows in a row group.
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setMaxRowGroupSize(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setMaxRowGroupSize(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets "created by" property.
    * @param {string} value
    * @returns {WriterPropertiesBuilder}
    */
    setCreatedBy(value) {
        const ptr = this.__destroy_into_raw();
        const ptr0 = passStringToWasm0(value, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.writerpropertiesbuilder_setCreatedBy(ptr, ptr0, len0);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets encoding for any column.
    *
    * If dictionary is not enabled, this is treated as a primary encoding for all
    * columns. In case when dictionary is enabled for any column, this value is
    * considered to be a fallback encoding for that column.
    *
    * Panics if user tries to set dictionary encoding here, regardless of dictionary
    * encoding flag being set.
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setEncoding(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setEncoding(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets compression codec for any column.
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setCompression(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setCompression(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets flag to enable/disable dictionary encoding for any column.
    *
    * Use this method to set dictionary encoding, instead of explicitly specifying
    * encoding in `set_encoding` method.
    * @param {boolean} value
    * @returns {WriterPropertiesBuilder}
    */
    setDictionaryEnabled(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setDictionaryEnabled(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets flag to enable/disable statistics for any column.
    * @param {boolean} value
    * @returns {WriterPropertiesBuilder}
    */
    setStatisticsEnabled(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setStatisticsEnabled(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets max statistics size for any column.
    * Applicable only if statistics are enabled.
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setMaxStatisticsSize(value) {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.writerpropertiesbuilder_setMaxStatisticsSize(ptr, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets encoding for a column.
    * Takes precedence over globally defined settings.
    *
    * If dictionary is not enabled, this is treated as a primary encoding for this
    * column. In case when dictionary is enabled for this column, either through
    * global defaults or explicitly, this value is considered to be a fallback
    * encoding for this column.
    *
    * Panics if user tries to set dictionary encoding here, regardless of dictionary
    * encoding flag being set.
    * @param {string} col
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setColumnEncoding(col, value) {
        const ptr = this.__destroy_into_raw();
        const ptr0 = passStringToWasm0(col, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.writerpropertiesbuilder_setColumnEncoding(ptr, ptr0, len0, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets compression codec for a column.
    * Takes precedence over globally defined settings.
    * @param {string} col
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setColumnCompression(col, value) {
        const ptr = this.__destroy_into_raw();
        const ptr0 = passStringToWasm0(col, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.writerpropertiesbuilder_setColumnCompression(ptr, ptr0, len0, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets flag to enable/disable dictionary encoding for a column.
    * Takes precedence over globally defined settings.
    * @param {string} col
    * @param {boolean} value
    * @returns {WriterPropertiesBuilder}
    */
    setColumnDictionaryEnabled(col, value) {
        const ptr = this.__destroy_into_raw();
        const ptr0 = passStringToWasm0(col, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.writerpropertiesbuilder_setColumnDictionaryEnabled(ptr, ptr0, len0, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets flag to enable/disable statistics for a column.
    * Takes precedence over globally defined settings.
    * @param {string} col
    * @param {boolean} value
    * @returns {WriterPropertiesBuilder}
    */
    setColumnStatisticsEnabled(col, value) {
        const ptr = this.__destroy_into_raw();
        const ptr0 = passStringToWasm0(col, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.writerpropertiesbuilder_setColumnStatisticsEnabled(ptr, ptr0, len0, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
    /**
    * Sets max size for statistics for a column.
    * Takes precedence over globally defined settings.
    * @param {string} col
    * @param {number} value
    * @returns {WriterPropertiesBuilder}
    */
    setColumnMaxStatisticsSize(col, value) {
        const ptr = this.__destroy_into_raw();
        const ptr0 = passStringToWasm0(col, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.writerpropertiesbuilder_setColumnMaxStatisticsSize(ptr, ptr0, len0, value);
        return WriterPropertiesBuilder.__wrap(ret);
    }
}

async function load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);

            } catch (e) {
                if (module.headers.get('Content-Type') != 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else {
                    throw e;
                }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);

    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };

        } else {
            return instance;
        }
    }
}

async function init(input) {
    if (typeof input === 'undefined') {
        // input = new URL('arrow1_bg.wasm', import.meta.url);
    }
    const imports = {};
    imports.wbg = {};
    imports.wbg.__wbindgen_string_new = function(arg0, arg1) {
        const ret = getStringFromWasm0(arg0, arg1);
        return addHeapObject(ret);
    };
    imports.wbg.__wbindgen_object_drop_ref = function(arg0) {
        takeObject(arg0);
    };
    imports.wbg.__wbg_buffer_7af23f65f6c64548 = function(arg0) {
        const ret = getObject(arg0).buffer;
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_newwithbyteoffsetandlength_ce1e75f0ce5f7974 = function(arg0, arg1, arg2) {
        const ret = new Uint8Array(getObject(arg0), arg1 >>> 0, arg2 >>> 0);
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_set_f25e869e4565d2a2 = function(arg0, arg1, arg2) {
        getObject(arg0).set(getObject(arg1), arg2 >>> 0);
    };
    imports.wbg.__wbg_length_0acb1cf9bbaf8519 = function(arg0) {
        const ret = getObject(arg0).length;
        return ret;
    };
    imports.wbg.__wbg_newwithlength_8f0657faca9f1422 = function(arg0) {
        const ret = new Uint8Array(arg0 >>> 0);
        return addHeapObject(ret);
    };
    imports.wbg.__wbindgen_throw = function(arg0, arg1) {
        throw new Error(getStringFromWasm0(arg0, arg1));
    };
    imports.wbg.__wbindgen_memory = function() {
        const ret = wasm.memory;
        return addHeapObject(ret);
    };

    if (typeof input === 'string' || (typeof Request === 'function' && input instanceof Request) || (typeof URL === 'function' && input instanceof URL)) {
        input = fetch(input);
    }



    const { instance, module } = await load(await input, imports);

    wasm = instance.exports;
    init.__wbindgen_wasm_module = module;

    return wasm;
}

export default init;


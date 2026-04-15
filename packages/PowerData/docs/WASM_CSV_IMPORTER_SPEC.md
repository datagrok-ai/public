# Arrow WASM CSV Importer for Datagrok

## Mission

Build a high-performance, streaming CSV importer that runs Apache Arrow's C++ CSV reader
compiled to WebAssembly inside a Web Worker. It parses CSV files (targeting 1–5 GB) in a
separate memory space, converts the output to Datagrok's `.d42` binary columnar format
(with dictionary encoding for strings), and transfers the result to the main thread via
`postMessage` + `Transferable` ArrayBuffers — achieving near-zero-copy handoff.

---

## 1. Architecture Overview

```
Main Thread                         Web Worker
───────────────────────────────────────────────────────────────
                                    ┌──────────────────────────┐
  File input (Blob)                 │  1. Receive File handle  │
       │                            │  2. ReadableStream chunks│
       │  postMessage(file)         │  3. Feed Arrow CSV reader│
       ├───────────────────────────►│     (WASM, streaming)    │
       │                            │  4. Arrow RecordBatch    │
       │                            │     → d42 chunk encoder  │
       │  postMessage(chunk,        │  5. Transfer chunk       │
       │    [transferables])        │                          │
       │◄───────────────────────────┤                          │
       │  ... repeat per chunk ...  │                          │
       │                            │                          │
       │  postMessage({done, meta}) │  6. Send final metadata  │
       │◄───────────────────────────┤     (column types, etc)  │
       │                            └──────────────────────────┘
       ▼
  Assemble d42 DataFrame
  (zero-copy attach buffers)
```

### Key principles

- **Streaming**: Never hold the entire file in memory. Read in chunks, parse in chunks,
  emit d42 chunks. Peak worker memory ≈ 2–3× chunk size, not file size.
- **Isolated memory**: WASM linear memory lives in the worker. Main thread heap untouched
  during parse.
- **Worker lifecycle**: Terminate the worker after import completes. WASM linear memory
  can only grow, never shrink — reuse would leak.
- **Progressive UI**: Each chunk triggers a progress update to the main thread so the UI
  can show a progress bar and row count.

---

## 2. Repository Structure

```
packages/csv-wasm-importer/
├── README.md
├── package.json
├── tsconfig.json
│
├── cpp/                          # Arrow C++ subset + glue
│   ├── CMakeLists.txt            # Emscripten build
│   ├── arrow-subset/             # Vendored minimal Arrow source
│   │   ├── csv/                  # arrow/csv/ (reader, parser, options)
│   │   ├── io/                   # arrow/io/ (InputStream, BufferReader)
│   │   ├── compute/              # Only: cast, dictionary_encode kernels
│   │   ├── array/                # Core array types
│   │   ├── buffer/               # Buffer, ResizableBuffer
│   │   ├── memory_pool/          # MemoryPool (default allocator)
│   │   ├── type/                 # DataType, Schema, Field
│   │   ├── record_batch/         # RecordBatch
│   │   └── util/                 # string_view, Status, Result, logging stubs
│   │
│   ├── csv_reader_binding.cpp    # Embind/Emscripten glue: JS ↔ Arrow CSV
│   └── d42_encoder.cpp           # Arrow RecordBatch → d42 binary chunk
│
├── src/                          # TypeScript (runs in worker + main)
│   ├── worker/
│   │   ├── csv-import.worker.ts  # Web Worker entry point
│   │   ├── wasm-loader.ts        # Load + instantiate .wasm
│   │   └── stream-feeder.ts      # ReadableStream → WASM input buffer
│   │
│   ├── main/
│   │   ├── csv-importer.ts       # Public API: importCsv(file, options) → DataFrame
│   │   ├── d42-assembler.ts      # Collect d42 chunks → Datagrok DataFrame
│   │   └── progress.ts           # Progress event types
│   │
│   ├── d42/
│   │   ├── format.ts             # d42 binary format constants & types
│   │   ├── encoder.ts            # TS reference encoder (for tests / fallback)
│   │   └── decoder.ts            # d42 chunk → typed arrays (main thread)
│   │
│   └── types.ts                  # Shared types, column type enum, options
│
├── build/                        # Build outputs
│   ├── csv_importer.wasm
│   ├── csv_importer.js           # Emscripten JS glue
│   └── csv-import.worker.js      # Bundled worker
│
└── test/
    ├── fixtures/                 # Small + medium CSV files
    ├── csv-reader.test.ts        # Unit: WASM CSV reader
    ├── d42-encoder.test.ts       # Unit: d42 encoding correctness
    ├── streaming.test.ts         # Integration: full pipeline, chunk by chunk
    └── large-file.bench.ts       # Benchmark: 100MB, 500MB synthetic CSVs
```

---

## 3. C++ / WASM Layer

### 3.1 Vendoring Arrow (Minimal Subset)

**Goal**: Keep the WASM binary ≤ 3 MB gzipped. Full Arrow C++ is ~50 MB of source; we
need a tiny fraction.

**What to include:**

| Arrow module          | Why                                            |
| --------------------- | ---------------------------------------------- |
| `arrow/csv/`          | CSV reader, column type inference, parser       |
| `arrow/io/`           | `InputStream`, `BufferedInputStream`, `RandomAccessFile` stubs |
| `arrow/array/`        | `ArrayData`, `StringArray`, `NumericArray`, `DictionaryArray` |
| `arrow/buffer/`       | `Buffer`, `ResizableBuffer`                     |
| `arrow/memory_pool/`  | Default memory pool (malloc-backed)             |
| `arrow/record_batch/` | `RecordBatch`                                   |
| `arrow/type/`         | `DataType`, `Field`, `Schema`                   |
| `arrow/compute/`      | **Only** `dictionary_encode` kernel             |
| `arrow/util/`         | `Status`, `Result`, `string_view`, iterator utils |

**What to exclude:**

- Parquet, ORC, Flight, IPC, JSON reader, Gandiva, Substrait, Dataset
- Filesystem abstractions (S3, GCS, HDFS)
- Full compute kernel library (aggregation, arithmetic, joins)
- Networking, compression codecs (unless snappy for d42)
- SIMD dispatch (Emscripten WASM SIMD is opt-in; start without, add later)

**Vendoring strategy:**

```bash
# From arrow/cpp/src/, copy only:
arrow/csv/
arrow/io/interfaces.{h,cc}
arrow/io/memory.{h,cc}
arrow/io/buffered.{h,cc}
arrow/array/
arrow/buffer/
arrow/memory_pool/
arrow/record_batch/
arrow/type/
arrow/compute/kernels/scalar_cast*.{h,cc}    # type casting
arrow/compute/kernels/hash_aggregate.cc       # dictionary_encode only
arrow/util/
arrow/result.h
arrow/status.{h,cc}
```

Then stub out anything that pulls in unwanted deps (filesystem, compression, etc.)
with no-op implementations or `#ifdef` guards.

### 3.2 CMakeLists.txt (Emscripten)

```cmake
cmake_minimum_required(VERSION 3.20)
project(csv_importer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# Emscripten-specific flags
set(WASM_FLAGS
  -O3                              # Full optimization
  -flto                            # Link-time optimization
  -s WASM=1
  -s ALLOW_MEMORY_GROWTH=1         # Critical: files can be large
  -s MAXIMUM_MEMORY=4GB            # Cap at 4GB (WASM32 limit)
  -s INITIAL_MEMORY=64MB           # Start small, grow as needed
  -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString']
  -s MODULARIZE=1
  -s EXPORT_ES6=1
  -s EXPORT_NAME=CsvImporterModule
  -s ENVIRONMENT=worker            # Only runs in Web Worker
  -s FILESYSTEM=0                  # No virtual FS needed
  -s NO_EXIT_RUNTIME=1
  -s DISABLE_EXCEPTION_CATCHING=0  # Arrow uses exceptions
  --bind                           # Embind for C++ ↔ JS
)

# Optional: WASM SIMD (enable after baseline works)
# -msimd128
# -s SIMD=1

add_executable(csv_importer
  cpp/csv_reader_binding.cpp
  cpp/d42_encoder.cpp
  # ... all vendored Arrow .cc files
)

target_compile_options(csv_importer PRIVATE ${WASM_FLAGS})
target_link_options(csv_importer PRIVATE ${WASM_FLAGS})
target_include_directories(csv_importer PRIVATE cpp/arrow-subset)
```

### 3.3 csv_reader_binding.cpp (Embind Glue)

This is the bridge between JS and Arrow's C++ CSV reader. Design it as a
streaming session object.

```cpp
// Pseudocode — actual implementation will use Arrow's API

#include <emscripten/bind.h>
#include <arrow/csv/api.h>
#include <arrow/io/memory.h>
#include "d42_encoder.h"

class CsvImportSession {
public:
    /// Configure the CSV reader.
    /// Called once before feeding data.
    void configure(
        int block_size,          // bytes per parse block (default 16MB)
        bool auto_detect_types,  // use Arrow type inference
        char delimiter,          // ',' by default
        char quote_char,         // '"' by default
        bool has_header          // first row is header
    );

    /// Feed a chunk of raw CSV bytes.
    /// Returns number of complete RecordBatches ready to consume.
    /// Arrow's streaming CSV reader buffers internally; a single feed()
    /// may produce 0 or N batches depending on block boundaries.
    int feed(uintptr_t data_ptr, int data_len);

    /// Signal EOF. Flushes the final partial block.
    /// Returns number of final RecordBatches produced.
    int finish();

    /// Get the next RecordBatch, already encoded as d42.
    /// Returns a pointer + length into WASM linear memory.
    /// JS side copies this out as an ArrayBuffer before calling next again.
    D42Chunk getNextChunk();

    /// Get the schema (column names + types) as JSON.
    /// Available after the first batch is produced.
    std::string getSchemaJson();

    /// Total rows parsed so far.
    int64_t rowCount();

private:
    std::shared_ptr<arrow::csv::StreamingReader> reader_;
    std::shared_ptr<arrow::io::BufferReader> buffer_reader_;
    D42Encoder encoder_;
    std::queue<D42Chunk> ready_chunks_;
    int64_t total_rows_ = 0;
};

EMSCRIPTEN_BINDINGS(csv_importer) {
    emscripten::class_<CsvImportSession>("CsvImportSession")
        .constructor<>()
        .function("configure", &CsvImportSession::configure)
        .function("feed", &CsvImportSession::feed)
        .function("finish", &CsvImportSession::finish)
        .function("getNextChunk", &CsvImportSession::getNextChunk)
        .function("getSchemaJson", &CsvImportSession::getSchemaJson)
        .function("rowCount", &CsvImportSession::rowCount);
}
```

**Critical implementation notes:**

1. **Arrow's StreamingReader** is the right API. It reads blocks incrementally
   and yields `RecordBatch` objects as they become available. Do NOT use
   `arrow::csv::TableReader` — that reads the entire file into a single Table.

2. **Buffer management**: `feed()` copies the JS chunk into a C++ buffer that
   Arrow's `BufferedInputStream` reads from. Use a ring buffer or double-buffer
   pattern so Arrow can parse block N while JS feeds block N+1.

3. **Type inference**: Arrow's CSV reader does multi-pass type inference on the
   first block. After the first block, types are locked. This is fine — the
   first chunk may be slower; subsequent chunks are type-stable.

4. **Dictionary encoding**: After Arrow produces a `RecordBatch`, run
   `arrow::compute::DictionaryEncode` on each string column before passing to
   the d42 encoder. This is where you save memory — a column with 10M rows
   but only 5K unique strings becomes a 5K-entry dictionary + int32 indices.

### 3.4 d42_encoder.cpp

Encodes an Arrow `RecordBatch` into Datagrok's d42 binary chunk format.

```
d42 Chunk Layout (per batch):
┌──────────────────────────────────────────┐
│ Header (fixed)                           │
│   magic:        4 bytes  "D42C"          │
│   version:      2 bytes  uint16          │
│   num_columns:  2 bytes  uint16          │
│   num_rows:     4 bytes  uint32          │
│   flags:        4 bytes  (bitfield)      │
├──────────────────────────────────────────┤
│ Column Directory (per column)            │
│   name_offset:  4 bytes  uint32          │
│   data_offset:  4 bytes  uint32          │
│   data_length:  4 bytes  uint32          │
│   type:         1 byte   enum            │
│   encoding:     1 byte   (raw|dict|rle)  │
│   null_bitmap:  4 bytes  offset to bitmap│
│   padding:      2 bytes                  │
├──────────────────────────────────────────┤
│ String Table (column names, dict values) │
│   length-prefixed UTF-8 strings          │
├──────────────────────────────────────────┤
│ Column Data Buffers                      │
│   For numeric: raw typed array           │
│   For dict-encoded string:               │
│     dict_count:   4 bytes uint32         │
│     dict_offsets:  (dict_count+1) uint32 │
│     dict_chars:    UTF-8 blob            │
│     indices:       int32[] per row       │
│   For null bitmap: packed bits           │
└──────────────────────────────────────────┘
```

**Type enum mapping:**

| d42 type | Arrow source type          | Buffer format        |
| -------- | -------------------------- | -------------------- |
| 0x01     | Int32                      | int32[]              |
| 0x02     | Int64                      | BigInt64Array        |
| 0x03     | Float32                    | float32[]            |
| 0x04     | Float64                    | float64[]            |
| 0x05     | String (dict-encoded)      | dict + int32 indices |
| 0x06     | Bool                       | packed bit array     |
| 0x07     | Timestamp (ms since epoch) | float64[]            |
| 0x08     | Date (days since epoch)    | int32[]              |

**IMPORTANT**: Align all buffer starts to 8-byte boundaries. This allows
the main thread to create TypedArray views directly over the transferred
ArrayBuffer without copying.

---

## 4. TypeScript / Worker Layer

### 4.1 csv-import.worker.ts

```typescript
// Web Worker entry point

import { CsvImporterModule } from '../build/csv_importer.js';

let wasmModule: any = null;

self.onmessage = async (e: MessageEvent) => {
  const { type, payload } = e.data;

  switch (type) {
    case 'init':
      // Load WASM module once
      wasmModule = await CsvImporterModule();
      self.postMessage({ type: 'ready' });
      break;

    case 'import':
      await handleImport(payload);
      break;
  }
};

async function handleImport(payload: {
  file: File;
  options: CsvImportOptions;
}) {
  const { file, options } = payload;
  const session = new wasmModule.CsvImportSession();

  session.configure(
    options.blockSize ?? 16 * 1024 * 1024,  // 16 MB default
    options.autoDetectTypes ?? true,
    (options.delimiter ?? ',').charCodeAt(0),
    (options.quoteChar ?? '"').charCodeAt(0),
    options.hasHeader ?? true
  );

  const STREAM_CHUNK_SIZE = 8 * 1024 * 1024; // 8 MB read chunks
  const reader = file.stream().getReader();

  let bytesRead = 0;

  // --- Streaming loop ---
  while (true) {
    const { done, value } = await reader.read();

    if (done) {
      const finalCount = session.finish();
      drainChunks(session, finalCount);
      break;
    }

    // Copy JS Uint8Array into WASM linear memory
    const ptr = wasmModule._malloc(value.byteLength);
    wasmModule.HEAPU8.set(value, ptr);

    const batchCount = session.feed(ptr, value.byteLength);
    wasmModule._free(ptr);

    bytesRead += value.byteLength;

    // Drain any ready d42 chunks
    drainChunks(session, batchCount);

    // Progress update
    self.postMessage({
      type: 'progress',
      bytesRead,
      totalBytes: file.size,
      rowCount: session.rowCount(),
    });
  }

  // Send schema + completion
  self.postMessage({
    type: 'done',
    schema: JSON.parse(session.getSchemaJson()),
    totalRows: session.rowCount(),
  });

  // Clean up — session destructor frees Arrow memory
  session.delete();
}

function drainChunks(session: any, count: number) {
  for (let i = 0; i < count; i++) {
    const chunk = session.getNextChunk();

    // Copy from WASM linear memory into a JS ArrayBuffer
    // that we can transfer to main thread
    const buffer = new ArrayBuffer(chunk.length);
    new Uint8Array(buffer).set(
      wasmModule.HEAPU8.subarray(chunk.ptr, chunk.ptr + chunk.length)
    );

    // Transfer — not copy — to main thread
    self.postMessage(
      { type: 'chunk', buffer, rowCount: chunk.rowCount },
      [buffer]  // Transferable list
    );
  }
}
```

### 4.2 csv-importer.ts (Main Thread Public API)

```typescript
export interface CsvImportOptions {
  delimiter?: string;
  quoteChar?: string;
  hasHeader?: boolean;
  autoDetectTypes?: boolean;
  blockSize?: number;           // bytes per WASM parse block
  onProgress?: (progress: ImportProgress) => void;
}

export interface ImportProgress {
  bytesRead: number;
  totalBytes: number;
  rowCount: number;
  percent: number;
}

/**
 * Import a CSV file into a Datagrok DataFrame.
 *
 * Spawns a Web Worker with Arrow WASM, streams the file through it,
 * and assembles the result as a DataFrame on the main thread.
 */
export async function importCsv(
  file: File,
  options: CsvImportOptions = {}
): Promise<DG.DataFrame> {

  const worker = new Worker(
    new URL('./worker/csv-import.worker.ts', import.meta.url),
    { type: 'module' }
  );

  return new Promise((resolve, reject) => {
    const assembler = new D42Assembler();

    worker.onmessage = (e: MessageEvent) => {
      const { type } = e.data;

      switch (type) {
        case 'ready':
          // Worker + WASM initialized, start import
          worker.postMessage({
            type: 'import',
            payload: { file, options }
          });
          break;

        case 'chunk':
          // Received a d42-encoded chunk (ArrayBuffer, already transferred)
          assembler.addChunk(e.data.buffer, e.data.rowCount);
          break;

        case 'progress':
          options.onProgress?.({
            bytesRead: e.data.bytesRead,
            totalBytes: e.data.totalBytes,
            rowCount: e.data.rowCount,
            percent: (e.data.bytesRead / e.data.totalBytes) * 100,
          });
          break;

        case 'done':
          // Assemble final DataFrame
          const df = assembler.toDataFrame(e.data.schema);
          worker.terminate();  // Kill worker to reclaim WASM memory
          resolve(df);
          break;

        case 'error':
          worker.terminate();
          reject(new Error(e.data.message));
          break;
      }
    };

    worker.onerror = (err) => {
      worker.terminate();
      reject(err);
    };

    // Initialize WASM in the worker
    worker.postMessage({ type: 'init' });
  });
}
```

### 4.3 d42-assembler.ts

```typescript
/**
 * Collects d42 binary chunks from the worker and assembles them
 * into a Datagrok DataFrame.
 *
 * Strategy: keep chunks as-is until done, then create column-level
 * TypedArray views that span across chunks (or concatenate if needed).
 */
export class D42Assembler {
  private chunks: ArrayBuffer[] = [];
  private rowCounts: number[] = [];
  private totalRows = 0;

  addChunk(buffer: ArrayBuffer, rowCount: number) {
    this.chunks.push(buffer);
    this.rowCounts.push(rowCount);
    this.totalRows += rowCount;
  }

  toDataFrame(schema: ColumnSchema[]): DG.DataFrame {
    // 1. Create empty DataFrame with schema
    const df = DG.DataFrame.create(this.totalRows);

    for (const col of schema) {
      const dgCol = df.columns.addNew(col.name, toDgType(col.type));

      // 2. For each column, walk chunks and copy data into the DG column buffer
      //    DG columns expose .getRawData() returning a typed array
      //    For dict-encoded strings, use DG's internal dictionary API
      let rowOffset = 0;

      for (let i = 0; i < this.chunks.length; i++) {
        const decoded = decodeD42Column(
          this.chunks[i], col.index, this.rowCounts[i]
        );

        if (col.encoding === 'dict') {
          // String column: merge dictionaries and remap indices
          copyDictEncodedStrings(dgCol, decoded, rowOffset);
        } else {
          // Numeric column: direct typed array copy
          dgCol.getRawData().set(decoded.values, rowOffset);
        }

        // Copy null bitmap
        if (decoded.nullBitmap) {
          copyNullBitmap(dgCol, decoded.nullBitmap, rowOffset, this.rowCounts[i]);
        }

        rowOffset += this.rowCounts[i];
      }
    }

    return df;
  }
}
```

---

## 5. Memory Management

### 5.1 Memory Budget

For a 2 GB CSV file with 16 MB parse blocks:

| Component                | Peak memory           |
| ------------------------ | --------------------- |
| File stream buffer (JS)  | ~8 MB (one read chunk)|
| WASM input buffer        | ~16 MB (one block)    |
| Arrow parse buffers      | ~48 MB (3× block)     |
| d42 encoded chunk        | ~16 MB (one output)   |
| **Worker total**         | **~90 MB peak**       |
| Main thread: d42 chunks  | Grows as import runs  |
| Main thread: DataFrame   | ~size of columnar data|

The main thread's d42 chunks are the big item. For a 2 GB CSV that compresses
to 800 MB columnar (typical with dict encoding), you'll hold ~800 MB in chunk
buffers before assembly. After assembly into the DataFrame, the chunks can be
GC'd.

### 5.2 Memory Optimization Strategies

1. **Stream chunks to DataFrame incrementally**: Instead of collecting all
   chunks then assembling, pre-allocate the DataFrame columns (requires knowing
   total row count — do a fast line-count pass first, or grow dynamically) and
   copy each chunk directly into the column buffers as it arrives. This halves
   main-thread peak memory.

2. **Two-pass approach for very large files**:
   - Pass 1 (fast): Count lines + detect types. Pure WASM, no d42 output.
     ~2 seconds for 2GB at streaming speed.
   - Pass 2: Parse with known row count → pre-allocated buffers.

3. **WASM memory ceiling**: Set `MAXIMUM_MEMORY=4GB` in Emscripten flags.
   Monitor `wasmModule.HEAP8.buffer.byteLength` and abort gracefully if
   approaching the limit.

4. **Worker termination**: Always `worker.terminate()` after import. The WASM
   linear memory (potentially hundreds of MB) is only freed when the worker
   dies.

---

## 6. Build System

### 6.1 Prerequisites

- Emscripten SDK (emsdk) ≥ 3.1.50
- CMake ≥ 3.20
- Node.js ≥ 18 (for TS build)

### 6.2 Build Steps

```bash
# 1. Install emsdk (one-time)
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk && ./emsdk install latest && ./emsdk activate latest
source ./emsdk_env.sh

# 2. Build WASM
cd packages/csv-wasm-importer
mkdir -p build && cd build
emcmake cmake ../cpp -DCMAKE_BUILD_TYPE=Release
emmake make -j$(nproc)
# Output: csv_importer.wasm, csv_importer.js

# 3. Build TypeScript
cd ..
npm run build
# Output: bundled worker + main thread API

# 4. Copy WASM artifacts into dist
cp build/csv_importer.{wasm,js} dist/
```

### 6.3 Integration with Datagrok Build

The WASM binary and JS glue should be served as static assets. In the Datagrok
webpack/vite config:

```javascript
// vite.config.ts
{
  worker: {
    format: 'es',  // ES module workers
  },
  assetsInclude: ['**/*.wasm'],
  optimizeDeps: {
    exclude: ['csv_importer.js']  // Don't try to bundle WASM glue
  }
}
```

The `.wasm` file should be loaded via `fetch()` in the worker, not inlined.
This allows the browser to compile it in streaming mode (`WebAssembly.compileStreaming`),
which is significantly faster for large WASM binaries.

---

## 7. Error Handling

### 7.1 CSV Parse Errors

Arrow's CSV reader has configurable error handling:

```cpp
auto read_options = arrow::csv::ReadOptions::Defaults();
auto parse_options = arrow::csv::ParseOptions::Defaults();
auto convert_options = arrow::csv::ConvertOptions::Defaults();

// On invalid rows: skip, error, or fill with null
convert_options.strings_can_be_null = true;

// For quoted fields with embedded newlines — Arrow handles these correctly
// but they affect chunk boundaries. Arrow's block-based parser re-syncs
// on quote boundaries automatically.
```

### 7.2 Error Propagation

```typescript
// In worker: catch and forward
try {
  // ... parse logic
} catch (err) {
  self.postMessage({
    type: 'error',
    message: err instanceof Error ? err.message : String(err),
    rowNumber: session.rowCount(),  // approximate location
  });
}
```

### 7.3 Abort / Cancel

Support cancellation via `AbortController`:

```typescript
// Main thread
const controller = new AbortController();

// Pass signal to worker via a shared flag or periodic polling
worker.postMessage({ type: 'import', payload: { file, options, signal: true } });

// To cancel:
controller.signal.addEventListener('abort', () => {
  worker.postMessage({ type: 'abort' });
});

// Worker handles abort by stopping the stream reader and cleaning up
```

---

## 8. Datagrok Integration Points

### 8.1 File Handler Registration

Register the WASM importer as a file handler in Datagrok's plugin system:

```typescript
// package.ts
export class CsvWasmPackage extends DG.Package {
  async init() {
    // Register for .csv, .tsv, .txt (tab-delimited)
    grok.functions.register(
      DG.Func.create({
        name: 'importCsvWasm',
        tags: ['file-handler'],
        inputs: [DG.FuncParam.create('file', DG.TYPE.FILE)],
        outputs: [DG.FuncParam.create('result', DG.TYPE.DATA_FRAME)],
        options: { extensions: 'csv,tsv,txt' },
      }),
      async (file: DG.FileInfo) => {
        const blob = await file.readAsBlob();
        return importCsv(new File([blob], file.name), {
          delimiter: file.extension === 'tsv' ? '\t' : ',',
          onProgress: (p) => grok.shell.setProgress(p.percent / 100),
        });
      }
    );
  }
}
```

### 8.2 Size Threshold

Use the WASM importer only for files above a size threshold (e.g., 50 MB).
For smaller files, Datagrok's built-in JS parser is fast enough and avoids
the WASM initialization overhead (~200ms).

```typescript
const WASM_THRESHOLD = 50 * 1024 * 1024; // 50 MB

async function importCsvSmart(file: File, options: CsvImportOptions) {
  if (file.size < WASM_THRESHOLD) {
    return DG.DataFrame.fromCsv(await file.text());
  }
  return importCsv(file, options);
}
```

---

## 9. Testing Strategy

### 9.1 Unit Tests

| Test                          | What it validates                                  |
| ----------------------------- | -------------------------------------------------- |
| `csv-reader.test.ts`          | WASM CSV reader produces correct RecordBatches      |
| `d42-encoder.test.ts`         | d42 binary format round-trips correctly             |
| `type-detection.test.ts`      | Int, float, string, date, bool detected accurately  |
| `dict-encoding.test.ts`       | String columns dictionary-encoded, indices correct  |
| `null-handling.test.ts`       | Empty cells → null bitmaps set correctly            |
| `delimiter.test.ts`           | Tab, pipe, semicolon delimiters work                |
| `quoted-fields.test.ts`       | Embedded commas, newlines, quotes handled           |

### 9.2 Integration Tests

| Test                          | What it validates                                  |
| ----------------------------- | -------------------------------------------------- |
| `streaming.test.ts`           | Full pipeline: File → Worker → d42 → DataFrame     |
| `large-file.test.ts`          | 100 MB synthetic CSV matches reference parse        |
| `memory-pressure.test.ts`     | Worker memory stays within budget during import     |
| `abort.test.ts`               | Cancel mid-import cleans up without leaks           |
| `progress.test.ts`            | Progress callbacks fire, monotonically increase     |

### 9.3 Benchmarks

```bash
# Generate synthetic CSVs
node test/generate-csv.js --rows=10000000 --cols=20 --output=test/fixtures/10M.csv

# Run benchmarks
npx vitest bench test/large-file.bench.ts
```

**Target performance:**

| File size | Parse time | Peak worker memory |
| --------- | ---------- | ------------------ |
| 100 MB    | < 0.5s     | < 150 MB           |
| 500 MB    | < 2s       | < 200 MB           |
| 1 GB      | < 4s       | < 250 MB           |
| 2 GB      | < 8s       | < 300 MB           |

---

## 10. Future Enhancements (Out of Scope for v1)

- **WASM SIMD**: Enable `-msimd128` for Arrow's SIMD-optimized string parsing.
  Expect ~30% speedup on Chrome/Edge (V8 SIMD support is mature).
- **SharedArrayBuffer**: If cross-origin isolation headers are configured,
  use SharedArrayBuffer to avoid the copy step between WASM memory and
  the transferable ArrayBuffer. This saves ~15% on large files.
- **Parallel column encoding**: Use multiple workers, each encoding a subset
  of columns. Only viable for very wide tables (100+ columns).
- **Compressed input**: Support gzipped CSV (`.csv.gz`) by adding a
  DecompressionStream stage before the WASM parser.
- **SDF/SD file support**: Extend the WASM layer to handle MDL SDF files
  (molecule + tag parsing). The d42 encoder already handles string columns;
  molecule columns would be dict-encoded SMILES/CTAB.
- **Write path**: Export DataFrame → CSV via Arrow's CSV writer in WASM.
  Same architecture in reverse.

---

## 11. Reference Implementations to Study

- **Perspective** (github.com/finos/perspective): JPMorgan's WASM data engine.
  Compiles Arrow C++ to WASM, runs in workers. Most mature reference for
  the Arrow-in-WASM pattern.
- **DuckDB-WASM** (github.com/duckdb/duckdb-wasm): Full SQL engine in WASM.
  Overkill for CSV import but excellent reference for Emscripten build config,
  memory management, and worker communication patterns.
- **Apache Arrow JS** (github.com/apache/arrow/js): Pure JS Arrow implementation.
  Useful for understanding the Arrow IPC format if you want d42 to be
  Arrow-IPC-compatible.
- **hyper-csv** (Rust): Fast Rust CSV parser. If you later want a Rust
  alternative to the C++ Arrow approach, this + wasm-bindgen is the path.

---

## 12. Implementation Order

**Phase 1: Foundation (1–2 weeks)**
1. Vendor Arrow C++ subset, get it compiling with Emscripten
2. Write `csv_reader_binding.cpp` — basic CSV → RecordBatch in WASM
3. Verify in a standalone HTML page (no Datagrok yet)

**Phase 2: d42 Encoding (1 week)**
4. Define d42 chunk binary format (finalize with Datagrok core team)
5. Implement `d42_encoder.cpp`
6. Implement `d42/decoder.ts` (TypeScript, main thread)
7. Round-trip tests: CSV → Arrow → d42 → TypedArrays → verify values

**Phase 3: Streaming Pipeline (1 week)**
8. Implement Web Worker (`csv-import.worker.ts`)
9. Implement stream feeder (ReadableStream → WASM chunks)
10. Implement main-thread API (`csv-importer.ts` + `d42-assembler.ts`)
11. End-to-end test with 100 MB file

**Phase 4: Datagrok Integration (1 week)**
12. Register as file handler in Datagrok package system
13. Add size-threshold routing (JS parser for small, WASM for large)
14. Wire up progress bar to Datagrok shell
15. Test with real pharma CSV datasets (SDF exports, assay plates, etc.)

**Phase 5: Optimization (ongoing)**
16. Profile with Chrome DevTools (WASM profiling)
17. Enable SIMD if benchmarks warrant
18. Tune block size and stream chunk size for optimal throughput
19. Benchmark against Datagrok's current JS parser — publish comparison

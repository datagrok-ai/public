# 10 — Performance

## Scale Considerations

| Scenario | File size | Events | Parameters | Challenge |
|---|---|---|---|---|
| Small conventional panel | 1–5 MB | 10K–50K | 6–10 | None — trivial |
| Standard clinical panel | 5–50 MB | 50K–200K | 10–20 | Parser speed, rendering |
| Large pharma HTS | 50–200 MB | 200K–1M | 10–20 | Parsing in background, subsampled render |
| Spectral cytometry | 200MB–1GB | 500K–2M | 30–50 | Large DATA segment, WebGL critical |
| CyTOF / mass cytometry | 100–500 MB | 100K–500K | 39–60 | High-dim rendering, server-side DR |
| Batch plate (96 files) | 96x file | 96x events | Same | Parallelism, memory management |

---

## Strategy 1: Web Workers for Parsing and Heavy Compute

**Problem:** Parsing a 200MB FCS file on the main thread freezes the browser UI for 2–5 seconds.

**Solution:** Offload to a Web Worker.

```typescript
// fcs-worker.ts (runs in Web Worker)
self.onmessage = (e: MessageEvent) => {
  const { arrayBuffer, fileId } = e.data;
  
  // Parse is synchronous but off main thread
  const parsed = parseFcs(arrayBuffer);
  
  // Post back structured data (transferable)
  self.postMessage({
    fileId,
    keywords: Object.fromEntries(parsed.keywords),
    // Transfer the Float32Arrays without copying (zero-copy)
    channelData: parsed.channelArrays,
    eventCount: parsed.eventCount,
    paramCount: parsed.paramCount
  }, parsed.channelArrays.map(a => a.buffer));  // Transfer ownership
};

// In main thread:
const worker = new Worker(new URL('./fcs-worker.ts', import.meta.url));
worker.postMessage({ arrayBuffer, fileId }, [arrayBuffer]);  // Transfer ownership
worker.onmessage = (e) => {
  const df = buildDataFrameFromWorkerResult(e.data);
  // Continue on main thread
};
```

**Also offload to workers:** Compensation matrix application, logicle transform computation, gate evaluation on large datasets, UMAP (client-side for <10K events).

---

## Strategy 2: Streaming / Chunked Parsing for Very Large Files

**Problem:** A 1GB FCS file cannot be read into memory as a single `ArrayBuffer` in some browsers.

**Solution:** Use the browser Streams API to parse the DATA segment in chunks.

```typescript
async function parseStreamingFcs(file: File): Promise<AsyncGenerator<Float32Array[]>> {
  const headerBuffer = await file.slice(0, 256).arrayBuffer();
  const { textStart, textEnd } = parseHeader(new DataView(headerBuffer));
  
  const textBuffer = await file.slice(textStart, textEnd + 1).arrayBuffer();
  const keywords = parseTextSegment(textBuffer);
  
  const dataStart = parseInt(keywords.get('$BEGINDATA')!);
  const nEvents = parseInt(keywords.get('$TOT')!);
  const nParams = parseInt(keywords.get('$PAR')!);
  const bytesPerEvent = nParams * 4;  // assuming $DATATYPE=F
  const chunkSize = 10000;  // events per chunk
  
  return async function*() {
    for (let eventOffset = 0; eventOffset < nEvents; eventOffset += chunkSize) {
      const batchSize = Math.min(chunkSize, nEvents - eventOffset);
      const byteStart = dataStart + eventOffset * bytesPerEvent;
      const byteEnd = byteStart + batchSize * bytesPerEvent;
      const chunk = await file.slice(byteStart, byteEnd).arrayBuffer();
      yield parseDataChunk(chunk, batchSize, nParams, keywords);
    }
  }();
}
```

---

## Strategy 3: Subsampling for Visualization

**Problem:** Rendering 1M events in a scatter plot is overkill — human perception cannot distinguish individual dots at that density.

**Solution:** Display a density-aware subsample; gate statistics always use the full dataset.

```typescript
// Density-aware subsampling
function subsampleForDisplay(
  xData: Float32Array, 
  yData: Float32Array, 
  targetN: number = 50000
): Int32Array {  // Returns indices to display
  if (xData.length <= targetN) {
    return Int32Array.from({ length: xData.length }, (_, i) => i);
  }
  
  // Grid-based subsampling: divide space into grid, sample from each cell
  const gridSize = 200;
  const grid = new Map<number, number[]>();
  
  for (let i = 0; i < xData.length; i++) {
    const gx = Math.floor(xData[i] / (1/gridSize));
    const gy = Math.floor(yData[i] / (1/gridSize));
    const key = gx * 1000 + gy;
    if (!grid.has(key)) grid.set(key, []);
    grid.get(key)!.push(i);
  }
  
  // Randomly sample from each grid cell
  const selected: number[] = [];
  const samplesPerCell = Math.ceil(targetN / grid.size);
  for (const events of grid.values()) {
    const n = Math.min(samplesPerCell, events.length);
    for (let i = 0; i < n; i++) {
      selected.push(events[Math.floor(Math.random() * events.length)]);
    }
  }
  
  return new Int32Array(selected);
}
```

**Important:** Gate evaluation always runs on the **full dataset**, never the subsample.

---

## Strategy 4: WebGL for Scatter Rendering

`regl-scatterplot` handles this automatically — it uses WebGL instancing to render millions of points efficiently. Key settings:

```javascript
const scatterplot = createScatterplot({
  canvas,
  // Performance tuning:
  pointSize: 2,                    // Smaller = faster render
  pointSizeSelected: 4,            // Slightly larger for selected
  opacity: 0.5,                    // Alpha blending for overplotting
  colorBy: 'valueA',               // Color by density (pre-computed)
  showReticle: true,
  lassoMinDelay: 10,               // ms between lasso updates
  lassoMinDist: 3,                 // px minimum drag to start lasso
});
```

For density coloring, pre-compute a 2D kernel density estimate and assign each event a density score. Map density to a blue-yellow-red colormap.

---

## Strategy 5: Typed Array Operations for Compensation

Compensation matrix application (multiply every event by a matrix) is done over millions of floats. Use `Float32Array` operations — the JIT compiler will vectorize these.

```typescript
function applyCompensation(
  rawData: Float32Array[],   // Array of columns (one per channel)
  compensationMatrix: number[][],  // n x n
  fluoroIndices: number[]   // Which column indices are fluorescence (not FSC/SSC)
): Float32Array[] {
  const nEvents = rawData[0].length;
  const n = fluoroIndices.length;
  
  // Build interleaved buffer for SIMD-friendly access
  const interleaved = new Float32Array(nEvents * n);
  for (let i = 0; i < nEvents; i++) {
    for (let j = 0; j < n; j++) {
      interleaved[i * n + j] = rawData[fluoroIndices[j]][i];
    }
  }
  
  // Apply matrix multiplication row by row
  const result = new Float32Array(nEvents * n);
  for (let i = 0; i < nEvents; i++) {
    for (let row = 0; row < n; row++) {
      let sum = 0;
      for (let col = 0; col < n; col++) {
        sum += compensationMatrix[row][col] * interleaved[i * n + col];
      }
      result[i * n + row] = sum;
    }
  }
  
  // De-interleave into output columns
  const output = fluoroIndices.map(() => new Float32Array(nEvents));
  for (let i = 0; i < nEvents; i++) {
    for (let j = 0; j < n; j++) {
      output[j][i] = result[i * n + j];
    }
  }
  
  return output;
}
```

For 200K events × 8 channels: ~3.2M float multiply-add operations. Expected: <50ms in a Web Worker.

---

## Strategy 6: Server-Side for Heavy Algorithms

**Client-side threshold:** UMAP, FlowSOM, PhenoGraph are NOT feasible in the browser for datasets >10K events.

| Algorithm | Client-side feasible? | Max events | Server-side tool |
|---|---|---|---|
| K-means | Yes | 500K | ml-kmeans |
| DBSCAN | Yes (small datasets) | 100K | density-clustering |
| UMAP | Yes | 10K | umap-learn (Python) |
| t-SNE | No | — | sklearn TSNE |
| FlowSOM | No | — | flowsom (Python) |
| PhenoGraph | No | — | phenograph (Python) |
| Leiden clustering | No | — | leidenalg (Python) |

**Trigger:** If dataset > 10K events AND user requests UMAP/FlowSOM, automatically route to Datagrok Python scripting. Show a progress indicator.

---

## Strategy 7: Batch Parallelism

For batch analysis of 96-well plates (96 FCS files):

```typescript
async function analyzeBatchParallel(
  files: DG.FileInfo[],
  template: GatingTree,
  concurrency: number = 4
): Promise<BatchResult[]> {
  const results: BatchResult[] = new Array(files.length);
  
  // Process in chunks of `concurrency` to avoid memory explosion
  for (let i = 0; i < files.length; i += concurrency) {
    const chunk = files.slice(i, i + concurrency);
    const chunkResults = await Promise.all(
      chunk.map(async (file, j) => {
        const df = await importFcsInWorker(file);
        const compensated = applyCompensationInWorker(df);
        return evaluateGatingTree(template, compensated);
      })
    );
    chunkResults.forEach((r, j) => { results[i + j] = r; });
    
    // GC hint: release DataFrames after stats are extracted
    // (Datagrok DataFrame is GC'd when all references dropped)
  }
  
  return results;
}
```

Expected throughput: 96 files × 50K events each = 4.8M total events. With 4 concurrent workers: ~30–60 seconds.

---

## Memory Budget

| Component | Memory per 100K events, 10 channels |
|---|---|
| Raw data (Float32) | 4 MB |
| Compensated data (Float32) | 4 MB |
| Transform-applied data (Float32) | 4 MB |
| Gate masks (1 bit/event × 20 gates) | 0.25 MB |
| Subsampled display (50K events × 2 coords) | 0.4 MB |
| **Total** | **~13 MB** |

For 96-well plate batch analysis processed sequentially: peak ~13 MB + overhead. For parallel (4 workers): ~52 MB peak. Both well within browser limits.

For 1M events (large spectral file): ~130 MB. Still manageable but use streaming parser and consider releasing intermediate buffers.

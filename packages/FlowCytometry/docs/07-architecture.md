# 07 — Architecture

## Plugin Identity

```
Package name:     FlowCytometry
Datagrok type:    Plugin (TypeScript + Webpack)
Entry point:      src/package.ts
Build:            webpack --config webpack.config.js
Datagrok API:     datagrok-api (npm)
```

## Package Directory Structure

```
packages/FlowCytometry/
├── src/
│   ├── package.ts                    # Entry: register all functions, viewers, file handlers
│   ├── fcs/
│   │   ├── fcs-parser.ts             # Binary FCS parser (ArrayBuffer + DataView)
│   │   ├── fcs-writer.ts             # FCS 3.1 export (Phase 3)
│   │   └── fcs-types.ts              # TypeScript interfaces for FCS data
│   ├── compensation/
│   │   ├── spillover.ts              # $SPILLOVER/$SPILL parsing
│   │   ├── compensation-matrix.ts    # Matrix inversion using ml-matrix
│   │   └── spectral-unmixing.ts      # WLS unmixing for spectral (Phase 3)
│   ├── transforms/
│   │   ├── logicle.ts                # Biexponential transform (port from Java)
│   │   ├── arcsinh.ts                # Arcsinh (trivial: Math.asinh(x/cofactor))
│   │   └── transform-types.ts        # Union type for all transform configs
│   ├── gating/
│   │   ├── gate-model.ts             # GateNode, GateTree data structures
│   │   ├── gates.ts                  # PolygonGate, RectangleGate, QuadrantGate, BooleanGate
│   │   ├── gate-evaluator.ts         # Point-in-polygon, point-in-ellipse algorithms
│   │   └── gatingml.ts               # GatingML 2.0 XML import/export
│   ├── workspaces/
│   │   └── wsp-parser.ts             # FlowJo .wsp XML parser (Phase 2)
│   ├── viewers/
│   │   ├── cytometry-scatter.ts      # Main dot/scatter plot viewer (JsViewer)
│   │   ├── histogram-viewer.ts       # 1D histogram viewer (JsViewer)
│   │   ├── density-viewer.ts         # Pseudo-color density viewer (JsViewer)
│   │   └── gating-tree-panel.ts      # Hierarchical gating tree panel (JsViewer)
│   ├── ui/
│   │   ├── compensation-dialog.ts    # Spillover matrix editor dialog
│   │   ├── transform-selector.ts     # Per-channel transform picker
│   │   └── statistics-table.ts       # Population statistics summary table
│   ├── statistics.ts                 # Population stats: count, %parent, MFI, CV
│   ├── batch.ts                      # Multi-file analysis, gating template application
│   └── scripts/
│       ├── flowkit-analysis.py       # FlowSOM, PhenoGraph, UMAP via FlowKit
│       ├── flowsom-cluster.py        # Standalone FlowSOM wrapper
│       └── flowcore-analysis.r       # flowCore/openCyto integration (alternative)
├── css/
│   └── flow-cytometry.css
├── detectors.js                      # Semantic type detectors for flow data columns
├── package.json
├── webpack.config.js
└── package.json (Datagrok manifest)
```

## Data Model: FCS Events as Datagrok DataFrame

### Mapping

```
FCS concept          →  Datagrok concept
───────────────────────────────────────────
One FCS file         →  One DG.DataFrame
Event (one cell)     →  One row
Parameter/channel    →  One column (DG.Column)
$PnN value           →  Column name (e.g., "FSC-A", "FITC-A")
$PnS value           →  Column tag: 'marker' (e.g., "CD4")
FCS TEXT keywords    →  DataFrame tags (DG.DataFrame.tags)
$SPILLOVER matrix    →  DataFrame tag: 'spillover' (JSON-serialized)
Compensated data     →  Additional columns, prefix "comp_"
```

### Construction Pattern

```typescript
// In fcs-parser.ts — after parsing binary DATA segment
function fcsToDataFrame(parsed: FcsParsed): DG.DataFrame {
  const df = DG.DataFrame.create(parsed.eventCount);

  for (let p = 0; p < parsed.paramCount; p++) {
    const col = DG.Column.fromFloat32Array(
      parsed.shortNames[p],   // e.g., "FSC-A"
      parsed.data[p]          // Float32Array of event values
    );
    col.setTag('marker', parsed.longNames[p] ?? '');  // e.g., "CD4"
    col.setTag('fcs_param_index', String(p + 1));
    df.columns.add(col);
  }

  // Store all FCS keywords as DataFrame tags
  for (const [key, value] of parsed.keywords) {
    df.setTag(key, value);
  }

  // Store serialized spillover matrix if present
  if (parsed.spillover) {
    df.setTag('spillover', JSON.stringify(parsed.spillover));
  }

  df.name = parsed.keywords.get('$FIL') ?? 'Unknown';
  return df;
}
```

### Semantic Type Detection

Register in `detectors.js`:

```javascript
//name: FlowCytometryDetector
//tags: semTypeDetector
//input: column col
//output: string semType
export function detectFlowCytometryColumn(col) {
  // Detect FSC/SSC/fluorescence channels by naming convention
  const name = col.name.toUpperCase();
  const fcsPattern = /^(FSC|SSC|FL\d|FITC|PE|APC|BV\d+|PACIFIC|PERCP|7AAD|DAPI)/;
  if (fcsPattern.test(name) && col.type === DG.COLUMN_TYPE.FLOAT) {
    return 'FlowCytometryChannel';
  }
  return null;
}
```

## Datagrok API Integration Points

### 1. File Importer Registration

```typescript
// In package.ts
//name: FcsFileImporter
//tags: fileImporter
//meta.ext: fcs
export async function importFcs(file: DG.FileInfo): Promise<DG.DataFrame[]> {
  const buffer = await file.readAsBytes();          // Returns Uint8Array
  const arrayBuffer = buffer.buffer;                 // Convert to ArrayBuffer
  const parsed = parseFcs(arrayBuffer);              // Your parser
  const df = fcsToDataFrame(parsed);
  return [df];
}
```

### 2. File Viewer Registration

```typescript
// In package.ts
//name: FcsFileViewer
//tags: fileViewer
//meta.ext: fcs
export async function viewFcs(file: DG.FileInfo): Promise<HTMLElement> {
  const buffer = await file.readAsBytes();
  const parsed = parseFcs(buffer.buffer);
  // Return a quick-view HTMLElement showing metadata + mini scatter plot
  return renderFcsPreview(parsed);
}
```

### 3. Custom Viewer Registration

```typescript
// CytometryScatterViewer extends DG.JsViewer
export class CytometryScatterViewer extends DG.JsViewer {
  xChannel: string;  // Property backed by Datagrok property system
  yChannel: string;
  gateTreeJson: string;

  constructor() {
    super();
    this.xChannel = this.string('xChannel', 'FSC-A');
    this.yChannel = this.string('yChannel', 'SSC-A');
  }

  onTableAttached() {
    this.subs.push(this.dataFrame.selection.onChanged.subscribe(() => this.render()));
    this.subs.push(this.dataFrame.filter.onChanged.subscribe(() => this.render()));
    this.render();
  }

  render() {
    // Get Float32Arrays for x/y channels
    const xData = this.dataFrame.getCol(this.xChannel).getRawData();
    const yData = this.dataFrame.getCol(this.yChannel).getRawData();
    const filter = this.dataFrame.filter.getRaw(0);  // BitArray

    // Render with regl-scatterplot (see 08-library-stack.md)
    this.scatter.setData({ x: xData, y: yData, filter });
  }
}

// Registration in package.ts:
//name: CytometryScatterViewer
//tags: viewer
DG.registerViewer('CytometryScatterViewer', () => new CytometryScatterViewer());
```

### 4. Python Scripting Integration

```typescript
// Calling FlowSOM from TypeScript via Datagrok scripting
async function runFlowSOM(df: DG.DataFrame, nClusters: number): Promise<DG.DataFrame> {
  const script = await grok.functions.eval('FlowCytometry:runFlowSOM');
  const result = await script.apply({
    data: df,
    n_clusters: nClusters
  });
  return result.get('result');  // DataFrame with 'cluster' column added
}
```

```python
# scripts/flowsom-cluster.py
#name: runFlowSOM
#language: python
#input: dataframe data
#input: int n_clusters = 10
#output: dataframe result
#conda: flowsom>=2.0, pandas

import flowsom as fs
import pandas as pd

# Select only fluorescence columns (exclude FSC/SSC/Time)
fluoro_cols = [c for c in data.columns if not c.startswith(('FSC','SSC','Time'))]
X = data[fluoro_cols].values

fsom = fs.FlowSOM(X, n_clusters=n_clusters)
data['cluster'] = fsom.cluster_labels_.astype(str)
result = data
```

## Gate Model Data Structure

```typescript
// gate-model.ts

interface GateNode {
  id: string;                    // Unique ID within experiment
  name: string;                  // Display name (e.g., "Lymphocytes")
  parentId: string | null;       // null = root (applied to all events)
  gateType: 'polygon' | 'rectangle' | 'ellipse' | 'quadrant' | 'boolean' | 'range';
  xChannel: string;              // $PnN of x-axis channel
  yChannel: string | null;       // null for 1D range gates
  transformX: TransformConfig;   // Transform active when gate was drawn
  transformY: TransformConfig | null;
  compensationRef: 'FCS' | 'uncompensated' | string;  // 'FCS' = apply FCS file compensation
  gateData: PolygonData | RectangleData | EllipseData | QuadrantData | BooleanData | RangeData;
  eventMask?: BitArray;          // Cached result — which events pass this gate
  statistics?: GateStatistics;  // Cached statistics
}

interface PolygonData {
  vertices: Array<[number, number]>;  // In transformed display space
}

interface RectangleData {
  xMin: number; xMax: number;
  yMin: number; yMax: number;  // null for 1D range gates
}

interface GateStatistics {
  count: number;
  percentParent: number;
  percentTotal: number;
  mfi: Record<string, number>;         // channel -> median fluorescence
  geometricMean: Record<string, number>;
  cv: Record<string, number>;
}

// The full tree for one FCS file's analysis:
interface GatingTree {
  fcsFileId: string;
  nodes: Map<string, GateNode>;
  rootIds: string[];  // Top-level gates (applied to all events)
}
```

## Gating Algorithm

```typescript
// gate-evaluator.ts

// Point-in-polygon using ray casting
function pointInPolygon(px: number, py: number, vertices: [number, number][]): boolean {
  let inside = false;
  const n = vertices.length;
  for (let i = 0, j = n - 1; i < n; j = i++) {
    const xi = vertices[i][0], yi = vertices[i][1];
    const xj = vertices[j][0], yj = vertices[j][1];
    const intersect = ((yi > py) !== (yj > py)) &&
      (px < (xj - xi) * (py - yi) / (yj - yi) + xi);
    if (intersect) inside = !inside;
  }
  return inside;
}

// Evaluate a gate against a DataFrame, return BitArray of passing events
function evaluateGate(node: GateNode, df: DG.DataFrame, parentMask: BitArray): BitArray {
  const xCol = df.getCol(node.xChannel).getRawData();  // Float32Array
  const yCol = node.yChannel ? df.getCol(node.yChannel).getRawData() : null;
  const xTransformed = applyTransform(xCol, node.transformX);
  const yTransformed = yCol ? applyTransform(yCol, node.transformY!) : null;

  const result = new BitArray(df.rowCount);
  for (let i = 0; i < df.rowCount; i++) {
    if (!parentMask.get(i)) continue;  // Parent gate must pass
    const x = xTransformed[i];
    const y = yTransformed ? yTransformed[i] : 0;
    result.set(i, testPoint(node, x, y));
  }
  return result;
}
```

## Multi-File / Batch Analysis Pattern

```typescript
// batch.ts

interface BatchResult {
  filename: string;
  wellId?: string;
  statistics: Map<string, GateStatistics>;  // gateId -> stats
}

async function applyGatingTemplate(
  fcsFiles: DG.FileInfo[],
  template: GatingTree
): Promise<BatchResult[]> {
  const results: BatchResult[] = [];

  for (const file of fcsFiles) {
    const df = await importFcs(file);
    const compensated = applyCompensation(df);
    const gatingResult = evaluateGatingTree(template, compensated);
    results.push({
      filename: file.name,
      wellId: df.getTag('$WELLID') ?? undefined,
      statistics: gatingResult.statistics
    });
  }

  return results;
}

// Convert batch results to a summary DataFrame for downstream analysis
function batchResultsToDataFrame(results: BatchResult[]): DG.DataFrame {
  // Rows = files/samples, columns = gate statistics
  // Suitable for Datagrok's dose-response, heatmap, etc.
}
```

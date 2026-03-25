# 08 — Library Stack

All dependencies needed to build the FlowCytometry plugin. Organized by concern.

## Core Dependencies

### FCS Parsing — Custom Implementation Required

**No suitable browser-native FCS parser exists.** The `fcs` npm package (`github.com/MorganConrad/fcs`) is Node.js-oriented and uses `Buffer`, not browser `ArrayBuffer`. Build a custom parser.

Reference implementations to study:
- **FlowIO (Python):** `github.com/whitews/FlowIO` — the cleanest, most spec-conformant reference; study `flowio/flowdata.py`
- **flowCore (R):** `bioconductor.org/packages/flowCore` — handles many vendor quirks
- **fcs npm (Node.js):** `npmjs.com/package/fcs` — simple, readable, good starting point for understanding structure

**Implementation:** Use browser `DataView` on `ArrayBuffer`. See `07-architecture.md` for pseudocode.

---

## Scatter Plot Rendering

### `regl-scatterplot` — PRIMARY CHOICE
```
npm install regl-scatterplot
github: github.com/flekschas/regl-scatterplot
stars: 800+
license: MIT
```
**Why:** WebGL-based, handles **20M+ points**, built-in **lasso selection** (critical for polygon gating), pan/zoom, color encoding by categorical or continuous variable. Purpose-built for large scatter plots.

**Usage:**
```javascript
import createScatterplot from 'regl-scatterplot';

const scatterplot = createScatterplot({
  canvas: canvasElement,
  width: 600,
  height: 600,
  pointSize: 2,
  pointColor: [0.0, 0.5, 1.0, 0.3],  // RGBA
  lassoColor: [1, 0, 0, 1],
  showRecticle: true,
});

// Set data: x and y must be Float32Array normalized to [-1, 1]
scatterplot.draw({ x: normalizedX, y: normalizedY });

// Listen for lasso selection
scatterplot.subscribe('select', ({ points }) => {
  // points = array of event indices inside lasso
  handleGateCreated(points);
});
```

**Gate overlay:** Draw polygon gate shapes on a Canvas overlay positioned above the WebGL canvas.

---

## Histogram Rendering

### `d3` (d3-array, d3-scale, d3-shape, d3-axis)
```
npm install d3
github: github.com/d3/d3
stars: 109K+
license: ISC
```
**Why:** Standard for axes, histogram bins, contour density estimation, gate path rendering.

**Usage for histogram:**
```javascript
import { bin, extent } from 'd3-array';
import { scaleLinear, scaleSymlog } from 'd3-scale';
import { area, line } from 'd3-shape';

const binner = bin().domain(extent(data)).thresholds(256);
const bins = binner(data);
// Render as SVG path or Canvas
```

**Usage for contour density:**
```javascript
import { contourDensity } from 'd3-contour';
const density = contourDensity()
  .x(d => xScale(d[0]))
  .y(d => yScale(d[1]))
  .size([width, height])
  .bandwidth(20)(events);
```

---

## Matrix Mathematics (Compensation)

### `ml-matrix`
```
npm install ml-matrix
github: github.com/mljs/matrix
stars: 200+
license: MIT
```
**Why:** Pure JavaScript matrix operations including LU decomposition and matrix inversion. Used to invert the spillover matrix to get the compensation matrix.

**Usage:**
```javascript
import { Matrix } from 'ml-matrix';

// Parse $SPILLOVER into a 2D array 'spilloverValues'
const S = new Matrix(spilloverValues);  // n x n

// Invert to get compensation matrix
const M = S.clone().pseudoInverse();  // or .inverse() if well-conditioned

// Apply to event: compensated = M * raw_fluorescence_vector
// For bulk application over all events:
const rawMatrix = new Matrix(eventData);  // events x n_channels
const compensated = rawMatrix.mmul(M.transpose());
```

---

## Dimensionality Reduction (Client-Side)

### `umap-js` — For UMAP in browser (small datasets only)
```
npm install umap-js
github: github.com/PAIR-code/umap-js
stars: 700+
license: Apache-2.0
```
**When to use client-side:** Datasets up to ~10,000 events. For larger datasets, use server-side Python.

**Usage:**
```javascript
import { UMAP } from 'umap-js';

const umap = new UMAP({ nComponents: 2, nNeighbors: 15, minDist: 0.1 });
const embedding = await umap.fitAsync(dataMatrix);
// embedding: Float64Array of [x0, y0, x1, y1, ...] pairs
```

### `@saehrimnir/druidjs` — PCA + 10 other DR algorithms
```
npm install @saehrimnir/druidjs
github: github.com/saehrimnir/druidjs
```
**Use for:** PCA (fast, client-side for any dataset size), MDS, t-SNE alternatives.

---

## Clustering (Client-Side)

### `ml-kmeans`
```
npm install ml-kmeans
github: github.com/mljs/kmeans
license: MIT
```
**Usage:** K-means for initial automated population identification.
```javascript
import { kmeans } from 'ml-kmeans';
const result = kmeans(data, 10, { initialization: 'kmeans++' });
// result.clusters: array of cluster assignments
```

### `density-clustering`
```
npm install density-clustering
npmjs: npmjs.com/package/density-clustering
license: MIT
```
**Usage:** DBSCAN for density-based population finding.
```javascript
import DBSCAN from 'density-clustering';
const dbscan = new DBSCAN();
const clusters = dbscan.run(data, epsilon, minPoints);
```

### `gaussian-mixture-model`
```
npm install gaussian-mixture-model
npmjs: npmjs.com/package/gaussian-mixture-model
```
**Usage:** GMM for fitting cell populations to Gaussian distributions.

---

## Server-Side (Python via Datagrok Scripting)

### FlowKit — FCS + gating + FlowSOM + UMAP
```
pip install flowkit
github: github.com/whitews/FlowKit
license: BSD-3
```
**Used for:** FlowSOM clustering, PhenoGraph, GatingML-compliant gating in Python, UMAP on large datasets.

```python
import flowkit as fk

# Load FCS file
sample = fk.Sample('path/to/file.fcs')

# Get compensated, transformed data as DataFrame
sample.apply_compensation(sample.metadata['spill'])
sample.apply_transform(fk.transforms.LogicleTransform('logicle', param_t=262144, param_w=0.5, param_m=4.5, param_a=0))
df = sample.as_dataframe(source='xform')

# Apply FlowSOM
# (FlowKit integrates FlowSOM via flowsom package)
```

### FlowSOM (Python)
```
pip install flowsom
github: github.com/saeyslab/FlowSOM_Python
```

### Phenograph
```
pip install phenograph
github: github.com/jacoblevine/PhenoGraph
```

### umap-learn
```
pip install umap-learn
```

### leidenalg (for Leiden clustering)
```
pip install leidenalg igraph
```

### flowCore (R alternative)
```r
BiocManager::install("flowCore")
BiocManager::install("openCyto")
```

---

## GatingML Parsing

### Browser-native `DOMParser`
No npm package needed. GatingML 2.0 is XML — use the browser's built-in DOMParser:

```javascript
const parser = new DOMParser();
const xmlDoc = parser.parseFromString(gatingmlString, 'text/xml');

// Query gates
const polygonGates = xmlDoc.querySelectorAll('PolygonGate');
polygonGates.forEach(gate => {
  const id = gate.getAttribute('gating:id');
  const parentId = gate.getAttribute('gating:parent_id');
  const vertices = Array.from(gate.querySelectorAll('vertex')).map(v => [
    parseFloat(v.querySelector(':first-child').getAttribute('data-type:value')),
    parseFloat(v.querySelector(':last-child').getAttribute('data-type:value'))
  ]);
});
```

---

## PDF/HTML Report Export

### `jspdf`
```
npm install jspdf
github: github.com/parallax/jsPDF
stars: 29K+
```
**Usage:** Capture canvas/SVG plots as images and embed in PDF.
```javascript
import jsPDF from 'jspdf';
const doc = new jsPDF();
const imgData = canvas.toDataURL('image/png');
doc.addImage(imgData, 'PNG', 10, 10, 180, 120);
doc.save('cytometry-report.pdf');
```

### `pdfmake`
```
npm install pdfmake
```
**Usage:** Structured reports with tables (statistics), headers, and layout. Combine with jsPDF for plot images embedded in structured reports.

---

## CSV Export

### `papaparse`
```
npm install papaparse
github: github.com/mholt/PapaParse
stars: 13K+
```
**Usage:** Export population statistics table to CSV.

---

## Complete `package.json` Dependencies

```json
{
  "dependencies": {
    "datagrok-api": "latest",
    "regl-scatterplot": "^9.0.0",
    "d3": "^7.0.0",
    "ml-matrix": "^6.10.0",
    "umap-js": "^1.3.3",
    "@saehrimnir/druidjs": "^0.5.0",
    "ml-kmeans": "^6.0.0",
    "density-clustering": "^1.3.0",
    "gaussian-mixture-model": "^0.3.0",
    "jspdf": "^2.5.1",
    "pdfmake": "^0.2.7",
    "papaparse": "^5.4.1"
  },
  "devDependencies": {
    "typescript": "^5.0.0",
    "webpack": "^5.0.0",
    "ts-loader": "^9.0.0"
  }
}
```

## Conda Environment for Python Scripts

```yaml
# conda-env.yml — used by Datagrok for Python scripting environment
name: flowcytometry
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - pip
  - pip:
    - flowkit>=1.1.0
    - flowsom>=2.0.0
    - flowio>=1.0.0
    - phenograph>=1.5.7
    - umap-learn>=0.5.3
    - leidenalg>=0.10.0
    - python-igraph>=0.10.0
    - pandas>=2.0.0
    - numpy>=1.24.0
    - scipy>=1.10.0
```

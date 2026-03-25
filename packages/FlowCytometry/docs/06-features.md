# 06 — Features

## Feature Priority Framework

Features are rated by segment importance:
- **P** = Pharma/biotech drug discovery (primary target)
- **A** = Academic researchers
- **C** = Clinical labs / CROs

Priority: **MUST** (MVP), **SHOULD** (Phase 2), **COULD** (Phase 3+)

---

## Core Features (MVP — Phase 1)

### FCS File Handling

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| FCS 3.0/3.1 binary parser | MUST | 5 | 5 | 5 | Primary format |
| FCS 2.0 backward compat | MUST | 4 | 3 | 3 | Many legacy files exist |
| Drag-and-drop import | MUST | 5 | 5 | 5 | |
| Multi-file batch import | MUST | 5 | 4 | 5 | |
| FCS metadata display (keywords) | MUST | 4 | 4 | 4 | |
| File preview in Datagrok browser | MUST | 5 | 5 | 4 | Register as `.fcs` file viewer |

### Compensation

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Parse `$SPILLOVER` / `$SPILL` from FCS | MUST | 5 | 5 | 5 | Both keywords required |
| Display spillover matrix as editable table | MUST | 5 | 5 | 4 | |
| Matrix inversion (spillover to compensation) | MUST | 5 | 5 | 5 | Use `ml-matrix` |
| Apply compensation to fluorescence channels | MUST | 5 | 5 | 5 | |
| Compensation from external matrix file (CSV) | MUST | 4 | 4 | 5 | Common in GvHD-style data |

### Transformations

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Logicle (biexponential) transform | MUST | 5 | 5 | 5 | ISAC standard; must port from Java/Python |
| Arcsinh transform | MUST | 5 | 4 | 3 | CyTOF standard; trivial to implement |
| Linear scale (scatter channels) | MUST | 5 | 5 | 5 | FSC/SSC always linear |
| Log scale (legacy) | MUST | 3 | 4 | 4 | Needed for FCS 2.0/3.0 data |

### Visualization (Plots)

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Scatter dot plot (FSC vs SSC, bivariate) | MUST | 5 | 5 | 5 | Core cytometry view |
| Pseudo-color density overlay on scatter | MUST | 5 | 5 | 5 | Blue-to-red by density |
| Histogram (1D per channel) | MUST | 5 | 5 | 5 | |
| Axis scale selection (logicle/arcsinh/log/linear) | MUST | 5 | 5 | 5 | |
| Pan and zoom on all plots | MUST | 5 | 5 | 4 | |
| Channel selector (X/Y axis) | MUST | 5 | 5 | 5 | |

### Manual Gating

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Polygon gate (draw on scatter plot) | MUST | 5 | 5 | 5 | Most common gate type |
| Rectangle gate | MUST | 5 | 5 | 5 | |
| Quadrant gate | MUST | 5 | 5 | 5 | CD4 vs CD8 paradigm |
| Range gate (1D on histogram) | MUST | 5 | 5 | 5 | |
| Gating tree panel (hierarchy display) | MUST | 5 | 5 | 5 | Shows parent/child structure |
| Gate editing (move vertices, resize) | MUST | 4 | 5 | 4 | |
| Gate deletion | MUST | 5 | 5 | 5 | |

### Statistics

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Event count per gate | MUST | 5 | 5 | 5 | |
| % Parent per gate | MUST | 5 | 5 | 5 | |
| % Total per gate | MUST | 5 | 5 | 5 | |
| MFI (median fluorescence intensity) per gate per channel | MUST | 5 | 5 | 5 | |
| Geometric mean per gate per channel | MUST | 4 | 4 | 4 | |
| CV (coefficient of variation) per gate | SHOULD | 4 | 4 | 5 | |
| Statistics export to CSV | MUST | 5 | 5 | 5 | |

---

## Differentiation Features (Phase 2)

### Batch Processing

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Apply gating template to multiple FCS files | MUST | 5 | 3 | 5 | Core pharma need |
| Cross-sample statistics table | MUST | 5 | 4 | 5 | |
| Statistical testing across groups (t-test, ANOVA) | SHOULD | 5 | 4 | 4 | |
| Dose-response curve fitting from population % | SHOULD | 5 | 2 | 2 | Key pharma feature |

### Interoperability

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| FlowJo .wsp workspace import | SHOULD | 5 | 5 | 3 | Migration from FlowJo |
| GatingML 2.0 export | SHOULD | 4 | 4 | 4 | Standards compliance |
| GatingML 2.0 import | SHOULD | 4 | 4 | 4 | From FlowKit/CellEngine |

### Advanced Visualization

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Contour plot overlay | SHOULD | 5 | 5 | 4 | |
| Heatmap (samples x channels) | SHOULD | 5 | 4 | 3 | |
| Overlay multiple populations | SHOULD | 5 | 5 | 3 | |
| Violin/box plots for population comparison | SHOULD | 5 | 4 | 4 | |
| PDF/HTML report export | SHOULD | 5 | 4 | 5 | |

### Plate Integration

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Link FCS files to plate well positions | SHOULD | 5 | 2 | 2 | Via $WELLID or metadata CSV |
| Plate heatmap visualization of statistics | SHOULD | 5 | 2 | 2 | |
| Dose-response with plate layout context | SHOULD | 5 | 1 | 1 | |

### Sample Metadata

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Upload metadata CSV (treatment, timepoint, dose) | SHOULD | 5 | 4 | 4 | |
| Group-by metadata for batch statistics | SHOULD | 5 | 4 | 4 | |
| Color scatter plots by metadata group | SHOULD | 4 | 4 | 3 | |

---

## Advanced Features (Phase 3+)

### Automated Analysis

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| FlowSOM clustering (server-side Python) | COULD | 4 | 5 | 3 | Via Datagrok Python scripting |
| PhenoGraph / Leiden clustering (server-side) | COULD | 4 | 5 | 3 | Via Datagrok Python scripting |
| UMAP dimensionality reduction (server-side) | COULD | 4 | 5 | 3 | Browser-side for <10K events |
| t-SNE (server-side) | COULD | 3 | 5 | 2 | Slow; UMAP preferred |
| AI-assisted gate boundary suggestion | COULD | 5 | 3 | 4 | Key differentiator |

### Spectral Cytometry

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| Spectral unmixing (WLS) | COULD | 4 | 3 | 3 | For Cytek Aurora, Sony ID7000 |
| Reference spectra library management | COULD | 4 | 3 | 3 | |
| Unmixing quality metrics (N x N plots) | COULD | 3 | 3 | 3 | |

### Compliance and QC

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| 21 CFR Part 11 audit trail | COULD | 5 | 1 | 5 | Leverage Datagrok infrastructure |
| Electronic signatures on analysis | COULD | 3 | 1 | 5 | |
| Levey-Jennings QC tracking | COULD | 3 | 2 | 5 | Instrument drift monitoring |
| Time gate (acquisition stability QC) | SHOULD | 5 | 4 | 5 | Standard first gate in workflow |

### Format Support

| Feature | Priority | P | A | C | Notes |
|---|---|---|---|---|---|
| FCS 3.2 support ($PnDATATYPE) | COULD | 3 | 3 | 3 | Few files exist yet |
| Write/export FCS 3.1 | COULD | 3 | 4 | 3 | Re-export gated populations |
| ACS container support | COULD | 2 | 3 | 3 | ZIP with FCS + GatingML |

---

## Feature Interaction Notes

- **Compensation MUST precede gating** — the plugin architecture must enforce this sequence
- **Transforms are applied for display only** — gating can be done in either raw or transformed space; GatingML must record which
- **All gate coordinates must store which transform was active** — critical for GatingML export correctness
- **Statistics always computed in compensated, untransformed space** — MFI is computed on compensated (not raw, not display-transformed) data

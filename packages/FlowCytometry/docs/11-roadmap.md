# 11 — Development Roadmap

## Phase Overview

| Phase | Name | Duration | Deliverable |
|---|---|---|---|
| **Phase 1** | MVP — Single File Analysis | 3–4 months | Functional parity with basic FlowJo for one file |
| **Phase 2** | Batch + Interop | 2–3 months | Plate-based analysis, FlowJo migration |
| **Phase 3** | Advanced Analytics | 3–4 months | AI/ML, spectral, compliance |
| **Phase 4** | Ecosystem | Ongoing | Clinical panels, omics integration |

---

## Phase 1: MVP — Single File Analysis

**Goal:** A scientist can open an FCS file, apply compensation, gate their populations manually, and export statistics. Functional parity with FlowJo for a single-file workflow.

### P1 Deliverables

#### Parser (Sprint 1–2)
- [ ] FCS 3.1 binary parser (TEXT + DATA segments, `DataView` on `ArrayBuffer`)
- [ ] FCS 3.0 and 2.0 backward compatibility
- [ ] Parse all required TEXT keywords (`$TOT`, `$PAR`, `$DATATYPE`, `$BYTEORD`, etc.)
- [ ] Parse `$SPILLOVER` and `$SPILL` keywords
- [ ] Parse `$PnN` (detector names) and `$PnS` (marker names)
- [ ] Map FCS events to `DG.DataFrame` (events = rows, channels = columns)
- [ ] Store all FCS keywords as DataFrame tags
- [ ] Register as `.fcs` file importer in Datagrok
- [ ] Register as `.fcs` file viewer (metadata preview)
- [ ] Web Worker offloading for files > 10MB

**Test:** Successfully parse all 8 canonical test files from `03-example-files.md`. DataFrame column count = `$PAR`, row count = `$TOT`.

#### Compensation (Sprint 2)
- [ ] Extract and display spillover matrix from `$SPILLOVER`/`$SPILL`
- [ ] Editable spillover matrix dialog (change individual coefficients)
- [ ] Matrix inversion using `ml-matrix`
- [ ] Apply compensation to fluorescence columns (add `comp_` prefix columns)
- [ ] Load compensation matrix from external CSV file
- [ ] Compensation undo/redo

**Test:** Apply compensation to FlowKit 8-color dataset. Compensated values match FlowKit Python output within floating point tolerance.

#### Transforms (Sprint 2)
- [ ] Logicle transform implementation (port from Java or FlowKit Python)
- [ ] Arcsinh transform (`Math.asinh(x / cofactor)`)
- [ ] Linear transform (passthrough)
- [ ] Log transform (for legacy data)
- [ ] Per-channel transform configuration
- [ ] Transform parameters UI (logicle T/M/W/A)

**Test:** Logicle output for known input values matches Wayne Moore's reference implementation and FlowKit Python output.

#### Scatter Plot Viewer (Sprint 3)
- [ ] `CytometryScatterViewer` extending `DG.JsViewer`
- [ ] `regl-scatterplot` WebGL rendering integrated
- [ ] X/Y channel selectors
- [ ] Scale mode selector (Logicle/Arcsinh/Log/Linear) per axis
- [ ] Density pseudo-color overlay (blue-to-red)
- [ ] Pan and zoom
- [ ] Datagrok filter/selection sync (filtered events dimmed)
- [ ] Subsampled display for datasets > 50K events

#### Histogram Viewer (Sprint 3)
- [ ] 1D histogram for any channel
- [ ] Scale mode (same as scatter)
- [ ] Overlay multiple populations (different colors)
- [ ] Gate range markers displayed on histogram

#### Gating (Sprint 4)
- [ ] Gating tree panel (`DG.JsViewer`)
- [ ] Polygon gate tool (click to add vertices, click-to-close)
- [ ] Rectangle gate tool (click-drag)
- [ ] Quadrant gate tool (two perpendicular lines)
- [ ] Range gate on histogram (1D)
- [ ] Gate evaluation: point-in-polygon, point-in-rectangle
- [ ] Hierarchical gating (child gates apply to parent events only)
- [ ] Gate editing: move vertices, drag gate
- [ ] Gate deletion
- [ ] Active population selection (clicks in tree update scatter plot view)

#### Statistics (Sprint 4)
- [ ] Per-gate: event count, %Parent, %Total
- [ ] Per-gate per-channel: MFI (median), Geometric Mean, CV
- [ ] Statistics table panel (live update on gate change)
- [ ] Export statistics to CSV

#### Phase 1 Acceptance Criteria
1. Open `101_DEN084Y5_15_E01_008_clean.fcs` (FlowKit 8-color)
2. Apply $SPILLOVER compensation
3. Set Logicle transform on all fluorescence channels
4. Gate: All Events → Time gate → Singlets (FSC-A vs FSC-H) → Lymphocytes (FSC vs SSC) → CD4+ (CD4 vs CD8)
5. Statistics for CD4+ gate: count, %parent, MFI of each channel
6. Export statistics CSV
7. All steps complete in < 5 seconds for a 50K-event file

---

## Phase 2: Batch + Interop

**Goal:** Apply gating templates across multiple files (plate-based experiments). Import from FlowJo. Export to GatingML.

### P2 Deliverables

#### Batch Processing (Sprint 5)
- [ ] Multi-file import (drag multiple .fcs files)
- [ ] Gating template save/load (JSON format)
- [ ] Apply template to all files in project
- [ ] Cross-sample statistics table (rows = files, columns = gate stats)
- [ ] Progress indicator for batch analysis
- [ ] Sample metadata CSV import (merge with statistics)
- [ ] Group-by-metadata for batch statistics

#### Plate Integration (Sprint 5–6)
- [ ] Detect `$WELLID`/`$PLATEID` keywords and auto-organize by plate
- [ ] Manual plate layout from metadata CSV
- [ ] Plate heatmap visualization (wrap Datagrok's existing plate viewer)
- [ ] Color wells by any gate statistic

#### Interoperability (Sprint 6)
- [ ] FlowJo .wsp XML parser
  - [ ] Extract FCS file references (sample paths)
  - [ ] Extract compensation matrices
  - [ ] Extract biexponential transform parameters per channel
  - [ ] Extract gate hierarchies (polygon, rectangle, quadrant)
  - [ ] Reconstruct gating tree in Datagrok
- [ ] GatingML 2.0 export
  - [ ] Export gate hierarchy as valid GatingML 2.0 XML
  - [ ] Include transform definitions (logicle, arcsinh)
  - [ ] Include compensation reference
- [ ] GatingML 2.0 import

#### Advanced Visualization (Sprint 7)
- [ ] Contour/density overlay on scatter plot
- [ ] Heatmap: samples vs channels (population MFI)
- [ ] Overlay multiple FCS files on same scatter plot (different colors)
- [ ] Violin/box plots for population comparison across samples
- [ ] Dose-response curve fitting from population % (IC50/EC50)
- [ ] PDF/HTML report export

#### Phase 2 Acceptance Criteria
1. Import 12 GvHD FCS files, apply external compensation, apply gating template
2. Cross-sample statistics table with %Lymphocytes, %CD4+, %CD8+, %CD19+
3. Import FlowJo .wsp and reproduce gates on same FCS file — gate coordinates within 2%
4. Export GatingML, reimport in FlowKit Python, verify gate membership matches

---

## Phase 3: Advanced Analytics

**Goal:** AI-assisted gating, spectral unmixing, 21 CFR Part 11 compliance.

### P3 Deliverables

#### Automated Analysis (Sprint 8–9)
- [ ] FlowSOM clustering via Datagrok Python scripting
- [ ] PhenoGraph / Leiden clustering
- [ ] UMAP dimensionality reduction (server-side for >10K events)
- [ ] AI-assisted gate boundary suggestion (ML model suggests gate, user confirms)
- [ ] Automated time gate (QC: exclude unstable acquisition periods)
- [ ] Automated singlet gate (FSC-A vs FSC-H ratio)

#### Spectral Cytometry (Sprint 9–10)
- [ ] Spectral unmixing (WLS implementation in TypeScript)
- [ ] Reference spectra library management (upload from CSV/FCS)
- [ ] Unmixing quality metrics (residuals, N×N plots)
- [ ] Autofluorescence as unmixing endmember
- [ ] Cytek Aurora file structure support

#### Compliance (Sprint 10)
- [ ] 21 CFR Part 11 audit trail (log all gating operations with user + timestamp)
- [ ] Electronic signature workflow for analysis approval
- [ ] Levey-Jennings QC chart for instrument monitoring
- [ ] Analysis version control (save multiple analysis versions per experiment)

#### Format Support (Sprint 10)
- [ ] FCS 3.2 support (`$PnDATATYPE` per-parameter types, new detector keywords)
- [ ] Export compensated/gated population as new FCS 3.1 file

---

## Phase 4: Ecosystem

**Goal:** Clinical panel templates, single-cell omics integration, advanced platform integrations.

### P4 Deliverables (Ongoing)

- [ ] Pre-configured clinical panel templates (lymphocyte subsets, T cell panel, B cell panel)
- [ ] CITE-seq protein data visualization (FCS-like from 10x Genomics)
- [ ] Integration with Datagrok Biological Registration (cell line metadata)
- [ ] Integration with Datagrok Compound Registration (compound-linked assays)
- [ ] CyTOF-specific optimizations (arcsinh-default UI, mass channel naming)
- [ ] Mass cytometry debarcoding
- [ ] FlowJo plugin migration tools (import FlowJo plugin results)

---

## Testing Strategy

### Unit Tests (vitest)
- FCS parser: each test file from `03-example-files.md`
- Compensation: known spillover matrix → expected compensated values
- Logicle: known input → expected output (compare to FlowKit reference)
- Gate evaluation: polygon, rectangle, ellipse with known point sets
- GatingML: round-trip export → import → same gate membership

### Integration Tests
- Full workflow test: import → compensate → gate → statistics → CSV export
- FlowJo .wsp import → gate comparison vs reference analysis
- Batch test: 12 GvHD files, verify statistics match flowCore R reference

### Performance Tests
- 100K events: end-to-end workflow < 5 seconds
- 1M events: parsing < 10 seconds (in worker)
- 96-file batch: < 5 minutes

### Test Data Location
```
packages/FlowCytometry/
├── test/
│   ├── data/                    # Symlinks or copies of canonical test files
│   │   ├── 0877408774.B08
│   │   ├── 8_color_sample.fcs   # FlowKit file
│   │   ├── GvHD_001.fcs
│   │   └── plate_A01.fcs        # flowTime well file
│   └── unit/
│       ├── fcs-parser.test.ts
│       ├── compensation.test.ts
│       ├── logicle.test.ts
│       ├── gate-evaluator.test.ts
│       └── gatingml.test.ts
```

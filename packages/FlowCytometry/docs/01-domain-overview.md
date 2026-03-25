# 01 — Domain Overview

## What Is Flow Cytometry?

Flow cytometry is a laser-based analytical technique that interrogates individual cells (or particles) as they flow single-file through one or more laser beams. Detectors capture:
- **Scattered light** — indicating cell size (FSC, Forward Scatter) and internal complexity/granularity (SSC, Side Scatter)
- **Fluorescence emission** — from dye-conjugated antibodies bound to specific cell-surface or intracellular markers

Modern instruments measure **18–50+ parameters simultaneously** per cell, producing datasets of thousands to millions of events per sample.

## Biological Questions Answered

| Question | Example assay |
|---|---|
| What immune cell subsets are present? | Immunophenotyping: CD3/CD4/CD8/CD19/CD56 panel |
| What fraction of cells are alive/dead? | Viability: Annexin V + PI |
| Are cells proliferating? | CFSE/CellTrace dye dilution |
| What phase of cell cycle? | DNA content: PI or DAPI |
| Is a drug hitting its target? | Receptor occupancy / target engagement |
| Are cells expressing a transgene? | Reporter: GFP, mCherry |
| Is a drug inducing apoptosis? | Caspase activation, cytochrome c release |
| What is the signaling state? | Phospho-flow: pSTAT3, pERK, pAKT |

## Typical Users

| Segment | Primary workflows |
|---|---|
| **Pharma/biotech drug discovery** | Target engagement, mechanism of action, immunomodulation, HTS dose-response |
| **Immunologists** | T/B/NK subset identification, activation, cytokine profiling |
| **Cell biologists** | Cell cycle, apoptosis, transfection efficiency, reporter expression |
| **Clinical researchers** | MRD (minimal residual disease), immune monitoring in clinical trials |
| **CROs** | GLP-compliant immunophenotyping, biomarker assay validation at scale |
| **Hospitals/clinical labs** | Diagnostic immunophenotyping (leukemia, HIV CD4, stem cell enumeration) |

> **Primary Datagrok target:** Pharma/biotech drug discovery and CROs — these organizations already use or evaluate Datagrok for compound registration, plate management, and analytics.

## Standard Analytical Workflow

```
1. SAMPLE PREPARATION
   └─ Single-cell suspension, antibody staining, controls (unstained, FMO, single-stain)

2. INSTRUMENT SETUP
   └─ QC beads, PMT voltage optimization, threshold/trigger configuration

3. DATA ACQUISITION
   └─ Run compensation controls → acquire samples as FCS files
   └─ Output: one .fcs file per sample/well, 10K–2M events each

4. COMPENSATION / UNMIXING
   ├─ Conventional: calculate spillover matrix from single-stain controls
   │   → invert spillover matrix → apply to fluorescence channels
   └─ Spectral: weighted least squares unmixing against reference spectra library

5. HIERARCHICAL GATING
   └─ Time gate → Scatter gate (FSC vs SSC) → Singlet gate (FSC-A vs FSC-H)
       → Viability gate → Lineage markers → Subset identification

6. STATISTICS
   └─ %Parent, %Grandparent, %Total, count, MFI, geometric mean, CV

7. REPORTING
   └─ Export statistics CSV, generate figures, save workspace/project
```

## Key Terminology Reference

| Term | Definition |
|---|---|
| **Event** | A single cell measurement; one row in the data matrix |
| **Parameter / Channel** | A measured dimension; one column (e.g., FITC-A, FSC-H) |
| **$PnN** | FCS keyword for the short channel name (detector name, e.g., "FITC-A") |
| **$PnS** | FCS keyword for the long channel name (marker name, e.g., "CD4") |
| **Scatter** | FSC (size), SSC (granularity) — not fluorescence |
| **Fluorochrome** | Fluorescent dye conjugated to antibody (FITC, PE, APC, BV421…) |
| **Spillover** | Signal from fluorochrome X leaking into detector Y |
| **Compensation** | Mathematical correction to remove spillover; requires matrix inversion |
| **Gate** | A selection rule defining a cell population |
| **Parent population** | The population a gate is applied to |
| **%Parent** | Fraction of parent population events passing the child gate |
| **MFI** | Mean (or Median) Fluorescence Intensity of a population in a channel |
| **FMO** | Fluorescence Minus One — control to set gate boundaries |
| **Singlets** | Events with one cell (vs doublets/aggregates) |
| **Logicle** | Biexponential transform mapping flow data to a displayable scale |
| **Arcsinh** | Simpler alternative transform; standard for CyTOF data |
| **GatingML** | XML standard for encoding gate hierarchies portably |
| **FlowSOM** | Self-organizing map clustering algorithm for automated gating |
| **UMAP / t-SNE** | Dimensionality reduction for visualizing high-dimensional data in 2D |

## Instrument Vendors and Their Ecosystems

| Vendor | Key Instruments | Bundled Software | Notes |
|---|---|---|---|
| BD Biosciences | FACSCanto, FACSAria, FACSLyric, FACSDiscover S8 | FACSDiva, FACSuite | Market leader; FlowJo is their analysis software |
| Beckman Coulter | CytoFLEX, DxFLEX, MoFlo, CytExpert | CytExpert, Kaluza | Strong clinical presence |
| Cytek Biosciences | Aurora, Northern Lights | SpectroFlo | Spectral cytometry leader |
| Sony | ID7000, SP6800 | Cell Sorter Software | Full-spectrum spectral |
| Miltenyi | MACSQuant | MACSQuantify | European academic strong |
| Luminex/Agilent | NovoCyte, NovoCyte Opteon | NovoExpress | Spectral + conventional |

> **Important for parser development:** Each vendor's instrument exports FCS with slightly different keyword sets and quirks. The parser must be lenient on optional keywords and strict only on required ones.

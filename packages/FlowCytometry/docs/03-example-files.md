# 03 — Example FCS Files for Testing

A curated set of canonical, publicly available FCS files covering the key parsing and analysis scenarios. Use these as the plugin's test fixtures.

## Recommended Download Order

Start with files 1–3 for basic parser development, then add the rest as each feature area is implemented.

---

## 1. `0877408774.B08` — The FCS "Hello World"

| Property | Value |
|---|---|
| **Source** | Bioconductor `flowCore` package |
| **FCS version** | 2.0 / 3.0 |
| **Parameters** | 6 (FSC-H, SSC-H, FL1-H, FL2-H, FL3-H, FL4-H) |
| **Events** | ~9,000 |
| **$SPILLOVER?** | No — compensation is external |
| **Instrument** | BD FACSCalibur |
| **What it tests** | Basic binary parser (TEXT + DATA), minimal parameter set, no inline compensation |

**Download:**
```
https://github.com/Bioconductor/flowCore/raw/master/inst/extdata/0877408774.B08
```

Notes: The canonical test file used in every flowCore vignette. Passes through any conformant FCS parser without edge cases. Start here.

---

## 2. FlowKit 8-Color Set — Inline `$SPILLOVER`, FCS 3.1

| Property | Value |
|---|---|
| **Source** | FlowKit Python library test data |
| **FCS version** | 3.1 |
| **Parameters** | 8 fluorescence + scatter (~10 total) |
| **Events** | ~50,000–100,000 per file |
| **$SPILLOVER?** | YES — inline in TEXT segment |
| **Instrument** | BD FACSCanto II |
| **What it tests** | FCS 3.1 parsing, $SPILLOVER extraction and matrix inversion, compensation |

**Download (single file):**
```
https://github.com/whitews/FlowKit/raw/master/data/8_color_data_set/fcs_files/101_DEN084Y5_15_E01_008_clean.fcs
```

**Browse all files:**
```
https://github.com/whitews/FlowKit/tree/master/data/8_color_data_set/fcs_files
```

Companion files in same repo: GatingML 2.0 gating strategy + FlowJo .wsp workspace.
This is the best all-in-one test dataset.

---

## 3. GvHD Dataset (12 files) — Multi-sample, External Compensation

| Property | Value |
|---|---|
| **Source** | FlowRepository FR-FCM-ZZY2 |
| **FCS version** | 2.0/3.0 |
| **Parameters** | 7 (FSC-H, SSC-H, FL1-H/CD15, FL2-H/CD45, FL3-H/CD14, FL4-H/CD33, FL2-A) |
| **Events** | 3,000–6,000 per file |
| **$SPILLOVER?** | No — external file |
| **Instrument** | BD FACSCalibur |
| **What it tests** | Multi-file batch loading, external compensation matrix, clinical immunophenotyping |

**Download:**
```
https://flowrepository.org/id/FR-FCM-ZZY2
```
(Free account required)

**flowGate package (first 3 files, no login):**
```
https://github.com/DillonHammill/flowGate/tree/master/inst/extdata
```

---

## 4. flowTime Plate Dataset — WELLID / PLATEID Keywords

| Property | Value |
|---|---|
| **Source** | Bioconductor `flowTime` package |
| **FCS version** | 3.1 |
| **Files** | 32 files named 1_A01.fcs through 1_H04.fcs |
| **$WELLID?** | YES |
| **$PLATEID?** | YES |
| **What it tests** | Plate layout reconstruction, batch-by-well analysis |

**Download:**
```
https://github.com/Bioconductor/flowTime/tree/master/inst/extdata/ss_example
```

---

## 5. FlowCAP-II AML Dataset — Multi-tube Clinical Panel

| Property | Value |
|---|---|
| **Source** | FlowRepository FR-FCM-ZZYA |
| **FCS version** | 3.0/3.1 |
| **Files** | 8 FCS files per patient (8 tubes), multiple patients |
| **Markers** | CD45, CD34, CD117, HLA-DR, CD13, CD33, CD7, CD19 |
| **What it tests** | Multi-tube experiment structure, clinical leukemia panel |

**Download:**
```
https://flowrepository.org/id/FR-FCM-ZZYA
```

---

## 6. Samusik CyTOF Dataset — Mass Cytometry, 39 Parameters

| Property | Value |
|---|---|
| **Source** | FlowRepository FR-FCM-ZZPH |
| **FCS version** | 3.1 |
| **Parameters** | 39 (metal isotope channels: `89Y_CD45`, `115In_CD11b`, etc.) |
| **Events** | ~170,000 |
| **$SPILLOVER?** | No — CyTOF has no fluorescence spillover |
| **Transform** | Arcsinh (cofactor=5) only |
| **What it tests** | High-dimensional rendering (39 params), arcsinh-only transform, absence of $SPILLOVER is valid |

**Download:**
```
https://flowrepository.org/id/FR-FCM-ZZPH
```

---

## 7. Spectral Flow Dataset (FR-FCM-Z4KT) — 31-Color

| Property | Value |
|---|---|
| **Source** | FlowRepository FR-FCM-Z4KT |
| **FCS version** | 3.1 |
| **Parameters** | 31 fluorescence + scatter |
| **Sample** | Human PBMCs, T cell panel |
| **Instrument** | Cytek Aurora |
| **What it tests** | Spectral unmixing (WLS), high-parameter rendering, Cytek file conventions |

**Download:**
```
https://flowrepository.org/id/FR-FCM-Z4KT
```

---

## 8. Single-Stain Compensation Controls

| Property | Value |
|---|---|
| **Source** | Bioconductor flowCore `compdata` |
| **Files** | One unstained + one per fluorochrome |
| **What it tests** | Spillover CALCULATION from controls (not just application) |

**Download:**
```
https://github.com/Bioconductor/flowCore/tree/master/inst/extdata/compdata
```

---

## Test Coverage Map

| Feature | Test File |
|---|---|
| Basic FCS parsing | File 1 |
| FCS 3.1 keyword parsing | File 2 |
| $SPILLOVER extraction | File 2 |
| Compensation matrix inversion | File 2 |
| External compensation matrix | File 3 |
| Multi-file batch loading | File 3 (12 files) |
| $WELLID / $PLATEID | File 4 |
| Multi-tube experiment | File 5 |
| High-dimensional (39+ params) | File 6 |
| Arcsinh-only (no compensation) | File 6 |
| Spectral unmixing (WLS) | File 7 |
| Spillover calculation from controls | File 8 |
| FlowJo .wsp import | File 2 (FlowKit repo has .wsp) |
| GatingML 2.0 import/export | File 2 (FlowKit repo has GatingML) |

## FlowRepository Access

FlowRepository datasets require a free account at flowrepository.org. Files are then downloadable via the web interface or API:
```
https://flowrepository.org/experiments/EXPERIMENT_ID/fcs_files/FILE_ID/download
```

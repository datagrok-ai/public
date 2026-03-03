# Sequence Annotation System — Design & Implementation Plan

## 1. Context and Motivation

### The Problem

When working with antibody or protein sequences in Datagrok, a user today can see the raw sequence rendered in the grid, scroll through aligned positions, view a WebLogo, and extract sub-regions by manually typing position ranges. But there is **no structured way to say "these positions are CDR1"**, no way to automatically flag sequence liabilities (deamidation hotspots, oxidation-prone residues), and no way to visually see region boundaries while scrolling through a long aligned sequence.

For antibody engineering — the primary use-case — this is a significant gap. A scientist looking at a panel of antibody candidates needs to:
- Know where the variable regions (VH/VL) are, and within them the framework (FR1–FR4) and complementarity-determining regions (CDR1–CDR3).
- Spot sequence liabilities (NG deamidation, DG isomerization, Met oxidation, N-linked glycosylation motifs) that may cause manufacturing or stability issues.
- Act on these annotations: extract a CDR, filter candidates by liability count, compare CDR3 loops across variants.

### What Exists Today

The Bio package already has building blocks that we can extend:

| Capability | Where | Status |
|---|---|---|
| **Position naming** — arbitrary string labels per position | Column tag `.positionNames`, read by `ISeqHandler.posList` in [seq-handler.ts](libraries/bio/src/utils/macromolecule/seq-handler.ts) | Working, but only populated manually or by external tools |
| **Region definitions** — named start/end spans | Column tag `.regions` (JSON array of `SeqRegion`), read by [get-region-func-editor.ts](src/utils/get-region-func-editor.ts) | Working, but only basic {name, description, start, end} |
| **Region extraction** — cut a sub-sequence by position names | `ISeqHandler.getRegion()` in [seq-handler.ts](libraries/bio/src/utils/macromolecule/seq-handler.ts), UI in [get-region.ts](src/utils/get-region.ts) | Working |
| **VD Regions viewer** — grid of WebLogos per FR/CDR | [vd-regions-viewer.ts](src/viewers/vd-regions-viewer.ts) using `VdRegion` interface from [vd-regions.ts](libraries/bio/src/viewers/vd-regions.ts) | Working, but requires manual VdRegion data; IMGT region defs are commented out (lines 27-48) |
| **MSA scrolling header** — WebLogo + conservation tracks above grid | [sequence-scrolling-widget.ts](src/widgets/sequence-scrolling-widget.ts), track base class `MSAHeaderTrack` in [sequence-position-scroller.ts](libraries/bio/src/utils/sequence-position-scroller.ts) | Working; extensible via new track subclasses |
| **Cell renderer** — monomer-by-monomer rendering with color coding | `MonomerPlacer` in [cell-renderer-monomer-placer.ts](libraries/bio/src/utils/cell-renderer-monomer-placer.ts), drawing via `printLeftOrCentered()` in [cell-renderer.ts](libraries/bio/src/utils/cell-renderer.ts) | Working; could accept per-position background colors |
| **Column property panel** — font, color, reference seq settings | [representations.ts](src/widgets/representations.ts) | Working; extensible with new accordion panes |

### What's Missing

1. **No antibody numbering scheme support** — no way to assign IMGT/Kabat/Chothia/AHo numbers to positions automatically.
2. **No liability scanning** — no built-in rules for detecting deamidation, oxidation, isomerization, glycosylation motifs.
3. **No visual annotation in the renderer** — no colored region bands, no liability markers in the grid cells.
4. **No annotation track in the MSA header** — no way to see region boundaries while scrolling.
5. **No annotation management UI** — no way to add, edit, remove, or toggle annotations.
6. **No annotation-driven actions** — no "extract this CDR" button from an annotation, no "filter by liability count".

---

## 2. Data Model

### 2.1 Design Principles

- **Column-level annotations** for things that apply uniformly to all rows (e.g., "CDR1 is positions 27-38 in IMGT numbering"). Stored as **column tags** so they serialize with projects automatically.
- **Row-level annotations** for per-sequence findings (e.g., "row 7 has an NG deamidation motif at position 52"). Stored as a **parallel string column** containing JSON.
- **Backward compatible** with the existing `.regions` tag — old data still works, new features are additive.

### 2.2 Core Types

New file: `libraries/bio/src/utils/macromolecule/annotations.ts`

```typescript
/** Visual mark type */
export enum AnnotationVisualType {
  Region = 'region',   // Contiguous span (e.g., CDR1)
  Point = 'point',     // Single position (e.g., oxidation-prone Met)
  Motif = 'motif',     // 2-3 consecutive positions (e.g., NG deamidation)
}

/** Grouping category */
export enum AnnotationCategory {
  Structure = 'structure',   // FR/CDR regions
  Liability = 'liability',   // Sequence liabilities
  PTM = 'ptm',               // Post-translational modifications
  Custom = 'custom',         // User-defined
}

/** Severity for liabilities */
export enum LiabilitySeverity {
  High = 'high',
  Medium = 'medium',
  Low = 'low',
  Info = 'info',
}

/**
 * A single annotation definition attached to a column.
 * Extends the existing SeqRegion concept with visual/categorical metadata.
 */
export interface SeqAnnotation {
  id: string;                              // Unique within column
  name: string;                            // "CDR1", "Deamidation (NG)"
  description?: string;
  start: string | null;                    // Position name (null for row-level-only)
  end: string | null;                      // Position name (null for point/row-level)
  visualType: AnnotationVisualType;
  category: AnnotationCategory;
  color?: string;                          // CSS color; defaults derived from category
  severity?: LiabilitySeverity;
  motifPattern?: string;                   // Regex for liability scanning
  sourceScheme?: string;                   // "IMGT", "Kabat", etc.
  autoGenerated?: boolean;
}

/** A single hit in a specific row */
export interface SeqAnnotationHit {
  annotationId: string;                    // References SeqAnnotation.id
  positionIndex: number;                   // 0-based index in the row's sequence
  positionName?: string;                   // From positionNames (for display)
  matchedMonomers: string;                 // The actual residues matched
  score?: number;
}

/** Stored per-row as JSON array in the companion annotation column */
export type RowAnnotationData = SeqAnnotationHit[];
```

### 2.3 Storage

| What | Where | Tag/Column Name | Format |
|---|---|---|---|
| Column-level annotations | Column tag | `.annotations` | JSON array of `SeqAnnotation` |
| Row-level annotation hits | Companion string column | `${colName}_annotations` | JSON array of `SeqAnnotationHit` per row |
| Numbering scheme name | Column tag | `.numberingScheme` | String: `"IMGT"`, `"Kabat"`, etc. |
| Position names (from numbering) | Column tag | `.positionNames` | Comma-separated strings (existing tag) |
| Back-reference to companion column | Column tag | `.annotationColumnName` | Column name string |

**Backward compatibility**: When reading annotations, check `.annotations` first; if absent, fall back to `.regions` and convert. When writing `.annotations`, also update `.regions` with the structure-category entries (so `GetRegionFuncEditor` continues to work).

### 2.4 New Constants

Add to `TAGS` enum in `libraries/bio/src/utils/macromolecule/consts.ts`:

```typescript
annotations = '.annotations',
numberingScheme = '.numberingScheme',
annotationColumnName = '.annotationColumnName',
```

---

## 3. Antibody Numbering Schemes

### 3.1 The Problem

Antibody variable region sequences don't have a natural position numbering. Different numbering schemes (IMGT, Kabat, Chothia, AHo, Martin) assign different numbers to the same physical residue. Each scheme also defines different CDR boundaries. Accurate numbering requires alignment against germline databases.

### 3.2 Approach: TypeScript-First with Python Fallback

**Primary approach — TypeScript heuristic numbering**:

Implement a lightweight antibody numbering engine in TypeScript that works client-side without any server dependencies. The approach:

1. **Chain type detection**: Identify VH vs VL (kappa/lambda) using conserved residue signatures:
   - Conserved Cys at positions ~23 and ~104 (IMGT) — present in all chains
   - Conserved Trp at position ~41 (IMGT) — present in all chains
   - VH-specific: conserved Trp at ~118 (IMGT FR4)
   - VK vs VL: distinguished by conserved residue patterns at specific positions

2. **Position assignment**: Use pairwise alignment (the existing `SequenceAlignment` class in [seq_align.ts](src/seq_align.ts) with BLOSUM matrices) to align the query sequence against curated germline reference sequences. The alignment maps query positions to scheme-numbered positions.

3. **Germline references**: Bundle a small set of representative germline sequences (IGHV, IGKV, IGLV families) with pre-assigned IMGT numbers. This is ~50-100 sequences, stored as a JSON data file.

4. **CDR boundary detection**: After IMGT numbering is assigned, CDR boundaries follow directly from the static region definitions (Section 3.3). For other schemes (Kabat, Chothia), apply well-known position mapping rules between IMGT and the target scheme.

**Secondary approach — Python script fallback** (for higher accuracy):

For users needing production-grade numbering (e.g., handling unusual frameworks, non-human species), provide a Python script that calls ANARCI via Datagrok's scripting infrastructure. This runs server-side but doesn't require Docker — just `pip install anarci` on the Datagrok server.

New file: `scripts/number_antibody.py`

```python
#name: Number Antibody Sequences
#language: python
#input: dataframe df
#input: column seqCol {semType: Macromolecule}
#input: string scheme {choices: ["IMGT", "Kabat", "Chothia", "AHo"]}
#output: dataframe result
from anarci import anarci
# ... process sequences, return numbered positions
```

### 3.3 Static Region Definitions

New file: `libraries/bio/src/utils/macromolecule/numbering-schemes.ts`

Contains the standard FR/CDR boundaries for each scheme and chain type. This is the data currently commented out in [vd-regions-viewer.ts:27-48](src/viewers/vd-regions-viewer.ts):

```typescript
export enum NumberingScheme { IMGT = 'IMGT', Kabat = 'Kabat', Chothia = 'Chothia', AHo = 'AHo' }
export enum ChainType { Heavy = 'Heavy', Light_Kappa = 'Light_Kappa', Light_Lambda = 'Light_Lambda' }

export interface SchemeRegionDef {
  name: string;           // "FR1", "CDR1", etc.
  type: 'FR' | 'CDR';
  startPosition: string;  // Position name in the scheme
  endPosition: string;
  chainType: ChainType;
}

export const IMGT_REGIONS: Record<ChainType, SchemeRegionDef[]> = {
  [ChainType.Heavy]: [
    { name: 'FR1',  type: 'FR',  startPosition: '1',   endPosition: '26',  chainType: ChainType.Heavy },
    { name: 'CDR1', type: 'CDR', startPosition: '27',  endPosition: '38',  chainType: ChainType.Heavy },
    { name: 'FR2',  type: 'FR',  startPosition: '39',  endPosition: '55',  chainType: ChainType.Heavy },
    { name: 'CDR2', type: 'CDR', startPosition: '56',  endPosition: '65',  chainType: ChainType.Heavy },
    { name: 'FR3',  type: 'FR',  startPosition: '66',  endPosition: '104', chainType: ChainType.Heavy },
    { name: 'CDR3', type: 'CDR', startPosition: '105', endPosition: '117', chainType: ChainType.Heavy },
    { name: 'FR4',  type: 'FR',  startPosition: '118', endPosition: '128', chainType: ChainType.Heavy },
  ],
  // ... Light_Kappa, Light_Lambda with their respective boundaries
};
// Similar for KABAT_REGIONS, CHOTHIA_REGIONS, AHO_REGIONS
```

### 3.4 Integration Flow

```
User: "Apply Numbering Scheme" dialog
  → Select column + scheme (or "Auto-detect") + engine (TypeScript / Python)
  → TypeScript path:
      → Align query sequences against bundled germline references
      → Determine chain type from conserved residue patterns
      → Map alignment positions to IMGT numbers
      → (If non-IMGT scheme selected, apply IMGT → target scheme mapping)
  → Python path (fallback):
      → Call number_antibody.py script via grok.functions.call()
      → Script runs ANARCI server-side
  → Set column tag .positionNames with scheme-specific position names
  → Set column tag .numberingScheme
  → Look up static region defs for the detected scheme + chain type
  → Convert to SeqAnnotation[] and store in .annotations tag
  → VdRegionsViewer auto-refreshes from .annotations
  → Cell renderer picks up region colors
```

---

## 4. Liability Detection

### 4.1 Approach: Client-Side Rule Engine

Liability patterns are simple regex motifs over single-letter amino acid sequences. No need for server-side computation. A TypeScript scanner iterates over rows, applies regex patterns, and records hits.

### 4.2 Built-in Rules

New file: `src/utils/annotations/liability-scanner.ts`

| Rule ID | Name | Pattern | Length | Severity | Category |
|---|---|---|---|---|---|
| `deamid-ng` | Deamidation (NG) | `NG` | 2 | High | Deamidation |
| `deamid-ns` | Deamidation (NS) | `NS` | 2 | Medium | Deamidation |
| `deamid-na` | Deamidation (NA) | `NA` | 2 | Low | Deamidation |
| `deamid-nd` | Deamidation (ND) | `ND` | 2 | Low | Deamidation |
| `deamid-nt` | Deamidation (NT) | `NT` | 2 | Low | Deamidation |
| `isom-dg` | Isomerization (DG) | `DG` | 2 | High | Isomerization |
| `isom-ds` | Isomerization (DS) | `DS` | 2 | Medium | Isomerization |
| `oxid-m` | Oxidation (Met) | `M` | 1 | Medium | Oxidation |
| `oxid-w` | Oxidation (Trp) | `W` | 1 | Low | Oxidation |
| `glyco-nxst` | N-glycosylation | `N[^P][ST]` | 3 | High | Glycosylation |
| `free-cys` | Free Cysteine | `C` | 1 | Info | Free Cysteine |

### 4.3 Scanner Logic

```typescript
function scanLiabilities(col, seqHandler, rules): { annotations: SeqAnnotation[], rowData: RowAnnotationData[] }
```

For each row:
1. Get the canonical single-letter sequence via `seqHandler.getSplitted(rowIdx)`
2. For each rule, run the regex against the canonical string
3. Record each match as a `SeqAnnotationHit` with the position index and matched monomers
4. After scanning all rows, create a `SeqAnnotation` entry for each rule that had at least one hit

Output:
- Column-level `SeqAnnotation[]` → stored in `.annotations` tag
- Row-level `RowAnnotationData[]` → stored in a new companion column `${colName}_annotations`

### 4.4 Custom Rules

Users can add custom pattern rules via the scanner dialog. Custom rules follow the same `LiabilityRule` structure and are stored as part of the `.annotations` tag with `autoGenerated: false`.

---

## 5. Rendering Integration

### 5.1 Cell Renderer Changes

**Where**: `MonomerPlacer.render()` in [cell-renderer-monomer-placer.ts](libraries/bio/src/utils/cell-renderer-monomer-placer.ts)

**Currently**, the rendering loop for single-line mode (around line 620) does:
```
for each visible position:
    printLeftOrCentered(g, monomer, x, y, w, h, opts)
```

**After the change**, a new `AnnotationRenderer` object is consulted before each monomer:
```
for each visible position:
    annotationRenderer.drawPositionBackground(g, x, y, w, h, posIdx, rowIdx)
    printLeftOrCentered(g, monomer, x, y, w, h, opts)
```

The same applies to the multiline rendering path (around line 568).

### 5.2 AnnotationRenderer

New file: `libraries/bio/src/utils/cell-renderer-annotations.ts`

This class:
1. Reads `.annotations` from the column tag and builds a position → color map for region annotations.
2. Reads the companion annotation column for row-level liability hits.
3. `drawPositionBackground()`: draws a translucent colored rectangle (opacity ~0.15) for region annotations, and a 2px colored underline for liability hits at that position.
4. Caches parsed annotation data; invalidates on tag change.

**Visual layering** (bottom to top):
1. **Region background** — translucent color band spanning the full cell height (e.g., light blue for CDR1, light green for FR2)
2. **Monomer text** — existing color-coded monomer rendering (unchanged)
3. **Liability underline** — 2px colored line at the bottom of the cell for positions with liability hits

### 5.3 MSA Header Annotation Track

**Where**: Extend `MSAHeaderTrack` base class from [sequence-position-scroller.ts:110](libraries/bio/src/utils/sequence-position-scroller.ts)

New file: `libraries/bio/src/utils/annotation-track.ts`

The annotation track renders colored horizontal bands with region labels in the MSA scrolling header:

```
┌──────────────────────────────────────────────────────────────┐
│ Column Name                                      [▼] [+]    │
├──────────────────────────────────────────────────────────────┤
│  FR1        │ CDR1  │ FR2       │ CDR2 │ FR3          │CDR3 │  ← NEW: Annotation Track
├──────────────────────────────────────────────────────────────┤
│ ████ ██ █ █████ ███ █ █ ██ ███████ █ █ ██ █████ █ █ █████  │  ← WebLogo Track (existing)
├──────────────────────────────────────────────────────────────┤
│ ▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮  │  ← Conservation Track (existing)
├──────────────────────────────────────────────────────────────┤
│ · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · │  ← Position Dots / Slider
└──────────────────────────────────────────────────────────────┘
```

- Uses the same viewport-aware lazy loading pattern as `LazyWebLogoTrack` in [sequence-scrolling-widget.ts](src/widgets/sequence-scrolling-widget.ts).
- Track height: ~20px.
- Draws colored rectangles per region, with the region name as centered text.
- Added to `initializeHeaders()` in [sequence-scrolling-widget.ts](src/widgets/sequence-scrolling-widget.ts) alongside existing tracks.

### 5.4 Tooltip Enrichment

In `MonomerPlacer.onMouseMove()` (in [cell-renderer-monomer-placer.ts](libraries/bio/src/utils/cell-renderer-monomer-placer.ts), around line 700), after the existing monomer tooltip, append annotation information:

```
Hovering over position 52:
┌─────────────────────────────┐
│ [A] Alanine                 │  ← Existing monomer tooltip
│ ─────────────────────────── │
│ Region: CDR2 (IMGT 56-65)  │  ← NEW: from column-level annotation
│ ⚠ Deamidation (NG) - High  │  ← NEW: from row-level annotation
└─────────────────────────────┘
```

### 5.5 Color Scheme

Default colors by category:
- **FR regions**: shades of blue (`#E3F2FD`, `#BBDEFB`, `#90CAF9`, `#64B5F6`)
- **CDR regions**: shades of orange/red (`#FFF3E0`, `#FFE0B2`, `#FFCC80`)
- **Liabilities**: category-specific (red for deamidation, orange for isomerization, purple for oxidation, green for glycosylation)
- All region backgrounds at 15% opacity to not obscure monomer colors.

---

## 6. User Interface

### 6.1 Top Menu

Add `Bio | Annotate` submenu in [package.ts](src/package.ts):

```
Bio
 ├── Annotate
 │    ├── Apply Numbering Scheme...     → NumberingSchemeDialog
 │    ├── Scan Liabilities...           → LiabilityScannerDialog
 │    ├── Add Region Annotation...      → ManualRegionDialog
 │    └── Manage Annotations...         → AnnotationManagerView
 ├── Analyze (existing)
 ├── Transform (existing)
 └── ...
```

### 6.2 "Apply Numbering Scheme" Dialog

```
┌─────────────────────────────────────────────┐
│ Apply Antibody Numbering                    │
│─────────────────────────────────────────────│
│ Table:      [Current Table            ▼]    │
│ Sequence:   [Heavy chain sequence     ▼]    │
│ Scheme:     [IMGT                     ▼]    │
│             [ ] Auto-detect scheme          │
│             [x] Populate FR/CDR regions     │
│             [x] Open VD Regions viewer      │
│─────────────────────────────────────────────│
│                       [Cancel]  [Apply]     │
└─────────────────────────────────────────────┘
```

### 6.3 "Scan Liabilities" Dialog

```
┌─────────────────────────────────────────────┐
│ Scan Sequence Liabilities                   │
│─────────────────────────────────────────────│
│ Table:      [Current Table            ▼]    │
│ Sequence:   [Heavy chain sequence     ▼]    │
│─────────────────────────────────────────────│
│ Rules:                                      │
│ [x] Deamidation (NG, NS, NA, ND, NT) High  │
│ [x] Isomerization (DG, DS)           High  │
│ [x] Oxidation (Met, Trp)           Medium  │
│ [x] N-glycosylation (NxS/T)          High  │
│ [ ] Free Cysteine                     Info  │
│ [+ Add Custom Rule]                        │
│─────────────────────────────────────────────│
│ Output:                                     │
│ [x] Highlight in cell renderer              │
│ [x] Create annotation column               │
│ [ ] Create summary count column             │
│─────────────────────────────────────────────│
│                       [Cancel]  [Scan]      │
└─────────────────────────────────────────────┘
```

### 6.4 Column Property Panel Extension

Extend `getMacromoleculeColumnPropertyPanel()` in [representations.ts](src/widgets/representations.ts) with a new "Annotations" accordion:

```
┌─────────────────────────────────────────────┐
│ ▸ Renderer Settings (existing)              │
│─────────────────────────────────────────────│
│ ▾ Annotations                               │
│   Numbering: IMGT                           │
│   Regions: FR1, CDR1, FR2, CDR2, FR3, ...  │
│   Liabilities: 3 rules (12 hits total)      │
│   [Show Legend] [Manage] [Clear All]        │
│─────────────────────────────────────────────│
│ ▸ Get Region (existing)                     │
└─────────────────────────────────────────────┘
```

### 6.5 Context Menu Extension

Extend `addCopyMenuUI()` in [context-menu.ts](src/utils/context-menu.ts):

```
Right-click on a cell at position 52:
  ───────────────────────────
  Annotations at position 52:
    CDR2 (IMGT 56-65)
    Deamidation (NG) - High
  ───────────────────────────
  Extract CDR2 as Column
  Filter: Rows with Liability Hits
  Select Rows with NG Motif
```

---

## 7. Actions on Annotations

| Action | Implementation | Reuses |
|---|---|---|
| **Extract region** | Call existing `getRegionDo()` with start/end from annotation | [get-region.ts](src/utils/get-region.ts) |
| **Filter by liability** | Create BitSet from row annotation column, apply to DataFrame | Standard DG filter |
| **Select annotated positions** | Set `.selectedPosition` tag on column | Existing renderer highlight |
| **Liability summary column** | New int column counting hits per row | — |
| **Compare annotations** | Property panel widget showing annotation matrix for selected rows | — |

---

## 8. Implementation Phases

### Phase 1: Data Model & Core Infrastructure
**Goal**: Types, storage, read/write utilities, backward compat with `.regions`.

**New files**:
- `libraries/bio/src/utils/macromolecule/annotations.ts` — type definitions
- `libraries/bio/src/utils/macromolecule/numbering-schemes.ts` — static region defs per scheme
- `src/utils/annotations/annotation-manager.ts` — CRUD for column annotations

**Modified files**:
- `libraries/bio/src/utils/macromolecule/consts.ts` — add new TAGS
- `src/utils/get-region-func-editor.ts` — read `.annotations` in addition to `.regions`

**Tests**: Serialization round-trip, backward compat with old `.regions` data.

### Phase 2: Liability Scanning
**Goal**: Working scanner with dialog and result columns.

**New files**:
- `src/utils/annotations/liability-scanner.ts` — rule engine + built-in rules
- `src/utils/annotations/liability-scanner-ui.ts` — dialog

**Modified files**:
- `src/package.ts` — register `scanLiabilities` function + `Bio | Annotate | Scan Liabilities` menu

**Tests**: Scan known antibody sequences (trastuzumab, adalimumab) and verify expected motifs.

### Phase 3: Rendering Integration
**Goal**: Annotations visible in the grid.

**New files**:
- `libraries/bio/src/utils/cell-renderer-annotations.ts` — `AnnotationRenderer` class

**Modified files**:
- `libraries/bio/src/utils/cell-renderer-monomer-placer.ts` — hook annotation drawing into both single-line (line ~620) and multiline (line ~568) render paths
- `libraries/bio/src/utils/cell-renderer-monomer-placer.ts` — enrich tooltip in `onMouseMove()` (line ~700)
- `src/widgets/representations.ts` — add Annotations accordion to property panel
- `src/utils/context-menu.ts` — add annotation context menu items

**Tests**: Visual verification with annotated sequences in single-line, multiline, and MSA modes.

### Phase 4: MSA Header Annotation Track
**Goal**: Annotation band in the scrolling header.

**New files**:
- `libraries/bio/src/utils/annotation-track.ts` — `AnnotationTrack` extending `MSAHeaderTrack`

**Modified files**:
- `src/widgets/sequence-scrolling-widget.ts` — add `AnnotationTrack` to header initialization

**Tests**: MSA sequences with IMGT-numbered positions; verify track alignment.

### Phase 5: Antibody Numbering
**Goal**: TypeScript-based numbering engine with Python fallback and UI.

**New files**:
- `libraries/bio/src/utils/macromolecule/antibody-numbering.ts` — TypeScript numbering engine (chain detection, germline alignment, position mapping)
- `libraries/bio/src/data/germline-references.json` — bundled representative germline sequences with IMGT numbers
- `src/utils/annotations/numbering-ui.ts` — "Apply Numbering Scheme" dialog
- `scripts/number_antibody.py` — Python/ANARCI fallback script

**Modified files**:
- `src/package.ts` — register numbering function + menu entry
- `src/viewers/vd-regions-viewer.ts` — auto-refresh from `.annotations` tag

**Tests**: Number known antibody sequences (trastuzumab VH/VL, adalimumab), verify against published IMGT reference numbering.

### Phase 6: Actions & Polish
**Goal**: Extract, filter, select from annotations; legend widget; manage dialog.

**New files**:
- `src/utils/annotations/annotation-actions.ts` — extract, filter, select, summary
- `src/widgets/annotation-legend-widget.ts` — color legend
- `src/utils/annotations/annotation-manager-ui.ts` — manage dialog
- `src/tests/annotation-tests.ts` — comprehensive test suite

---

## 9. Architecture Decisions

| Decision | Choice | Rationale |
|---|---|---|
| **Storage for column-level annotations** | Column tags (`.annotations`) | Serializes with projects automatically. Consistent with existing `.regions` tag pattern. |
| **Storage for row-level annotations** | Parallel string column with JSON | Standard Datagrok pattern for per-row metadata. Survives DataFrame serialization. Respects filtering. |
| **Numbering engine** | TypeScript heuristic (primary) + Python/ANARCI script (fallback) | TypeScript runs client-side with zero server deps. Uses existing `SequenceAlignment` + BLOSUM matrices from [seq_align.ts](src/seq_align.ts). Python fallback for higher accuracy when server available. |
| **Liability scanning** | Client-side TypeScript regex | Simple patterns, no server needed. Instant results. Works offline. Fast enough for large datasets. |
| **Rendering approach** | Transparent background overlay + underline markers | Preserves existing monomer coloring. Layered approach encodes multiple info channels without conflict. |
| **Annotation track in header** | New `MSAHeaderTrack` subclass | Extensible architecture already exists. Follows same viewport-aware lazy loading pattern. |
| **Backward compat with `.regions`** | Read both tags; write both | `GetRegionFuncEditor` and `getRegionDo()` continue to work unchanged for basic use-cases. |

---

## 10. Key Reuse Points

| Existing Component | File | How Reused |
|---|---|---|
| `SeqRegion` type | [get-region-func-editor.ts:18-23](src/utils/get-region-func-editor.ts) | Extended into `SeqAnnotation` |
| `.regions` tag | [consts.ts:34](libraries/bio/src/utils/macromolecule/consts.ts) | Read for backward compat |
| `positionNames` / `positionLabels` tags | [consts.ts:32-33](libraries/bio/src/utils/macromolecule/consts.ts) | Numbering output populates these |
| `VdRegion` interface | [vd-regions.ts:19-34](libraries/bio/src/viewers/vd-regions.ts) | Annotations convert to VdRegions |
| `VdRegionsViewer` | [vd-regions-viewer.ts](src/viewers/vd-regions-viewer.ts) | Receives annotation data |
| `getRegionDo()` | [get-region.ts](src/utils/get-region.ts) | "Extract Region" action |
| `ISeqHandler.getRegion()` | [seq-handler.ts](libraries/bio/src/utils/macromolecule/seq-handler.ts) | Core extraction, unchanged |
| `ISeqHandler.posList` | [seq-handler.ts](libraries/bio/src/utils/macromolecule/seq-handler.ts) | Position name mapping |
| `MonomerPlacer` | [cell-renderer-monomer-placer.ts](libraries/bio/src/utils/cell-renderer-monomer-placer.ts) | Extended with annotation hooks |
| `MSAHeaderTrack` | [sequence-position-scroller.ts:110](libraries/bio/src/utils/sequence-position-scroller.ts) | Base class for AnnotationTrack |
| `SequenceAlignment` (Needleman-Wunsch) | [seq_align.ts](src/seq_align.ts) | Align query vs germline references for numbering |
| Column property panel | [representations.ts](src/widgets/representations.ts) | Extended with Annotations pane |
| Context menu | [context-menu.ts](src/utils/context-menu.ts) | Extended with annotation actions |
| `selectedPosition` tag | [consts.ts:36](libraries/bio/src/utils/macromolecule/consts.ts) | Highlight-on-click |

---

## 11. Potential Challenges

1. **Performance of row-level JSON parsing**: Parsing JSON per cell on every render frame would be slow. Mitigation: `AnnotationRenderer` caches parsed data, keyed by annotation column version. Only re-parses on change.

2. **Annotation overlay vs. existing diff highlighting**: Reference sequence comparison uses transparency to dim matching monomers. Annotation backgrounds must use very low opacity (0.10-0.15) to avoid visual conflict. Underline markers provide a non-overlapping secondary channel.

3. **Position mapping for unaligned sequences**: Row-level liability positions are relative to each row's own sequence. For MSA, the renderer translates through `getSplitted()` per row.

4. **Numbering accuracy**: The TypeScript heuristic approach works well for common human antibody frameworks but may struggle with unusual germlines (camelid, shark) or heavily engineered sequences. The Python/ANARCI fallback covers these cases when server scripting is available. The UI should clearly indicate which engine was used.

5. **Backward compat with `.regions` consumers**: Keep `.regions` in sync when writing `.annotations`, so `GetRegionFuncEditor` and `getRegionDo()` continue to work.

---

## 12. Verification

After implementation, verify with this end-to-end test:

1. Load a table with aligned antibody heavy chain sequences (FASTA, MSA).
2. `Bio | Annotate | Apply Numbering Scheme` → select IMGT → verify `.positionNames` and `.annotations` tags are set, VdRegionsViewer shows FR/CDR regions.
3. `Bio | Annotate | Scan Liabilities` → enable all rules → verify companion annotation column is created, cell renderer shows region backgrounds and liability underlines.
4. Hover over an annotated position → verify tooltip shows region name + liability info.
5. Scroll through MSA header → verify annotation track shows region bands aligned with WebLogo.
6. Right-click on a CDR position → "Extract CDR2 as Column" → verify new column contains the correct sub-sequences.
7. Save and reload project → verify all annotations persist.
8. Run `npm run test` to verify no regressions.

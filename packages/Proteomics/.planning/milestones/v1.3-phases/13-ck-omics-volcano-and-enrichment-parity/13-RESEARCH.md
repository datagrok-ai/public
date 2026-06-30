# Phase 13: CK-omics Volcano and Enrichment Parity - Research

**Researched:** 2026-05-16
**Domain:** Datagrok proteomics package — UniProt batched fetch + classification port, ScatterPlot live-property toggles, Spectronaut Candidates contrast normalization, g:Profiler up/down split
**Confidence:** HIGH (all integration surfaces, CK-omics reference algorithms, and the UniProt stream response shape verified directly this session)

## Summary

This phase ports six well-bounded behaviors from the client's `CKomics_tool2.py` into the
already-shipped Proteomics package. Every decision is locked in 13-CONTEXT.md (D-01..D-12);
the research below answers *how to implement each well* and *what will bite during execution*,
not *whether*. Three of the six (R1 subcellular coloring, R3 contrast flip, R2 up/down
enrichment) are direct algorithm ports where divergence is a parity regression — the exact
keyword map, hex palette, flip mechanics, and split logic are reproduced verbatim below from
the source files so the executor never re-derives them. The other three (R4 Filter viewer,
R5 WikiPathways, R6 Q/P toggle) are small, additive Datagrok-API changes.

The single highest-risk finding, verified live this session: the UniProt **stream** endpoint
returns TSV with header `Entry` (NOT `accession`), and `cc_subcellular_location` is verbose
multi-isoform free text containing `SUBCELLULAR LOCATION:` prefixes and `{ECO:...|PubMed:...}`
evidence blocks. CK-omics' classifier survives this because it scans for `\b`-bounded
keyword substrings rather than parsing structure — **the port must preserve that
substring-scan semantics exactly**, and the TSV column mapping must key off positional
order, not header names (CK-omics splits by `\t` positionally for the same reason).

**Primary recommendation:** Build one new shared module `src/analysis/subcellular-location.ts`
(batched stream fetch + verbatim classifier + persistent cache) that both `uniprot-panel.ts`
and the new volcano color path consume; do all R3 sign-flip work inside
`spectronaut-candidates-parser.ts` (pure), all Filter-viewer/dock wiring in the `package.ts`
import handler; make R6's Q/P metric and R1's color dimension two synchronized live toggles
driven by a single recompute function in `volcano.ts`.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| UniProt batched subcellular fetch (R1) | API/Backend (via `grok.dapi.fetchProxy`) | Browser cache | External HTTP must be proxied (CORS); result cached client-side |
| Subcellular classification (R1, D-04) | Browser (pure TS) | — | Deterministic string mapping, no I/O — port of `parse_subcellular_location` |
| Persistent location cache (D-02) | Browser storage | — | Cross-session; `grok.dapi.userDataStorage` is the platform-native KV |
| Volcano color + Q/P toggle (R1, R6) | Browser (viewer) | — | Pure viewer-state recompute on a DataFrame already in memory |
| Candidates contrast sign-flip (R3, D-08) | Browser (parser, pure) | — | Pure DataFrame transform; no shell/view side effects |
| Filter viewer docking (R4, D-07) | Browser (package.ts handler) | — | Shell/view orchestration — parser stays pure |
| Report-DE direction fix (R3, D-09) | Browser (analysis) | — | Group-ordering logic in `differential-expression.ts` |
| Up/down enrichment queries (R2, D-10) | API/Backend (g:Profiler via proxy) | Browser | Two proxied POSTs; merge into one DataFrame client-side |
| WikiPathways source (R5, D-12) | Browser (sources array) | — | One literal `'WP'` added to the request payload |

## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** Fetch via UniProt **stream** endpoint through `grok.dapi.fetchProxy()` —
  `/uniprotkb/stream?query=accession:(A OR B …)&fields=accession,cc_subcellular_location,go_c`,
  chunked ~100 accessions/request. Replaces per-row raw `fetch()`; CORS-safe.
- **D-02:** Persistent cross-session cache keyed by accession (IndexedDB or
  `grok.userSettings` — Claude's discretion). Shared subcellular-location module;
  `uniprot-panel.ts` refactored to reuse. Fully closes the cache-uniprot todo.
- **D-03:** Keep CK-omics' reviewed-entry-by-gene-name fallback: unreviewed entry with
  empty subcellular field → look up reviewed entry by gene name (one extra batched query
  for misses), use its location.
- **D-04:** Classify into 11 categories + Unknown using the keyword map, ordering rule
  (UniProt order-of-importance; GO-CC fallback), and exact hex palette from the CK-omics
  `Subcellular_Location_Classification_README.txt`. **Locked client contract — port
  verbatim, do not re-derive.**
- **D-05:** One volcano, not two. Add a `Subcellular Location` categorical column with the
  locked 11-colour map; volcano switches colour between significance (up/down/NS) and
  Subcellular Location via the viewer colour property plus a menu/dialog switch.
- **D-06:** R6 Q-value vs P-value is a live viewer-property toggle, switched after creation.
  Recomputes Y axis **and** up/down/NS classification **and** threshold lines together,
  kept in sync on every toggle. Default = `adj.p-value` (= the client Q-value).
- **D-07:** R4 — import all comparisons (no modal picker). Retain `Comparison` column,
  expose via a native Datagrok Filter viewer docked with the volcano. Single-comparison
  files need no filter. **Tradeoff explicitly accepted:** no "view the mirror" control.
- **D-08:** R3 — canonical convention: positive log2FC = group1-enriched, where group1 is
  the declared numerator in each `A / B` comparison string. Parser normalizes each row's
  sign to its own declared `group1/group2` label and swaps/relabels the AVG Group Quantity
  Numerator/Denominator columns to match (CK-omics `create_subset_data`). Declared label
  IS the orientation — no per-run prompt.
- **D-09:** R3 also unparks the report-import DE direction fix in
  `differential-expression.ts`: honor the intended/declared contrast instead of
  alphabetical-condition-order default; retire the manual Comparison-dropdown workaround.
  This phase is the explicit prompt overriding the memory's "do not implement unprompted".
- **D-10:** Two g:Profiler queries (up-regulated genes, down-regulated genes) merged into
  one enrichment DataFrame with a `Direction` column (Up/Down). Reuses Phase-9
  dot-plot/bar-chart and volcano cross-link unchanged. Background stays all detected
  proteins (existing custom-background behaviour) for both directional queries.
- **D-11:** Present the two directions as split Up/Down dot + bar charts side-by-side
  (source filter still applies to both).
- **D-12:** Add `'WP'` to the g:Profiler source list, default-on.

### Claude's Discretion

- Accession chunk size; cache storage mechanism (IndexedDB vs `grok.userSettings`) and
  invalidation policy; progress-indicator cadence for the location fetch.
- Exact wording/placement of the colour-toggle switch and how the live Q/P toggle is
  wired into the `ScatterPlotViewer` property panel.
- Dock arrangement/sizing for side-by-side Up/Down enrichment charts and the `Comparison`
  Filter viewer relative to the volcano.
- Whether the directional gene lists are also surfaced (CK-omics writes
  `*_up_genes_all*.txt` / `*_down_genes_all*.txt`) — not required, planner's call.

### Deferred Ideas (OUT OF SCOPE)

- Explicit "view the mirror of a declared contrast" control — traded away in D-07.
- Surfacing directional gene-list exports (`*_up_genes_*.txt` / `*_down_genes_*.txt`).
- `2026-03-03-expand-de-dialog-with-method-selection-comparison-picker` todo — related but
  DE-dialog scope, distinct from R3/R4. Stays open for its own phase.

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| R1 | Subcellular-location coloring on volcano | §"CK-omics Algorithm Port: Subcellular Classification", §"UniProt Stream Endpoint", §"Cache Mechanism", §"Volcano Color Dimension" |
| R2 | Up/down split in enrichment | §"CK-omics Algorithm Port: g:Profiler Up/Down Split", §"Enrichment Viewer Direction Extension" |
| R3 | Contrast auto-flip in Candidates parser + unpark report-DE fix | §"CK-omics Algorithm Port: create_subset_data Flip", §"Pitfall: Parser Purity", §"Report-DE Direction Fix (D-09)" |
| R4 | Multi-contrast Candidates selection via Filter viewer | §"R4: Comparison Filter Viewer Docking" |
| R5 | WikiPathways g:Profiler source | §"R5: WikiPathways Source (one-line)" |
| R6 | Q-value vs P-value live toggle on volcano | §"Volcano Q/P Live Toggle (D-06)" |

## Standard Stack

No new external packages. Everything is in-repo / platform-provided.

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `datagrok-api` | linked local | `grok.dapi.fetchProxy`, `grok.dapi.userDataStorage`, `DG.ScatterPlotViewer`, Filter viewer, DataFrame | Platform API; CLAUDE.md mandates `grok.dapi.*`, never raw `fetch` |
| `@datagrok-libraries/test` | linked local | test category/test/expect | Existing test convention (`src/tests/`) |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `@datagrok-libraries/statistics` | already a dep | (unchanged) used by client-side DE fallback | Only touched indirectly via D-09 |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `grok.dapi.userDataStorage` for D-02 cache | Browser IndexedDB (`idb`) | userDataStorage is platform-native, cross-device, no new dep, async KV — but a single value is a string (JSON-encode the map). IndexedDB is larger-capacity and origin-scoped but adds a dep and is device-local. **Recommend userDataStorage** (see §"Cache Mechanism"). |
| UniProt `stream` endpoint | `search` endpoint (what CK-omics actually uses) with `size=500` paging | `stream` returns all matches in one response with no paging (D-01 locked it). CK-omics uses `search`+`size` because Python `requests` has no streaming concern; the chunked-OR query shape is identical. |

**Installation:** none — `npm install` unchanged. Build is the standard
`grok api && grok check --soft && webpack` via `npm run build`.

## Package Legitimacy Audit

Not applicable — this phase installs **zero** external packages. All work uses
platform-provided `datagrok-api` (linked locally via `grok link`) and the existing
`@datagrok-libraries/*` already in `package.json`. No npm/PyPI/crates installs occur.

## Architecture Patterns

### System Architecture Diagram

```
                          ┌─────────────────────────────────────┐
  Spectronaut Candidates  │  parseSpectronautCandidatesText()    │   PURE
  .tsv/.csv  ───────────► │  (parsers/spectronaut-candidates-…)  │   (no shell)
                          │  • detect comparison column          │
                          │  • R3/D-08: per-row sign-normalize    │
                          │    log2FC + swap AVG Group Qty cols   │
                          │    to canonical group1-enriched       │
                          │  • set proteomics.* tags              │
                          └──────────────┬───────────────────────┘
                                         │ DG.DataFrame
                                         ▼
                          ┌─────────────────────────────────────┐
  package.ts import       │  importSpectronautCandidates()       │   SHELL
  handler                 │  • addTableView(df)                   │   ORCHESTRATION
                          │  • R4/D-07: if >1 distinct Comparison │
                          │    → dock native Filter viewer        │
                          │  • createVolcanoPlot(df)              │
                          └──────────────┬───────────────────────┘
                                         ▼
   protein IDs ──► subcellular-location.ts (NEW shared module)
                   ┌──────────────────────────────────────────┐
                   │ getSubcellularLocations(accessions[])     │
                   │  1. cache lookup (userDataStorage, by acc)│
                   │  2. miss → chunk ~100 → fetchProxy stream │
                   │     /uniprotkb/stream?query=accession:(…) │
                   │     &fields=accession,cc_subcellular_loc, │
                   │     go_c,reviewed,gene_primary&format=tsv │
                   │  3. parseSubcellularLocation() [verbatim] │
                   │  4. D-03: Unknown+gene → reviewed-by-gene │
                   │     batched fallback query                │
                   │  5. write-through cache                   │
                   └────────────────┬─────────────────────────┘
                                    │ Map<accession, category>
            ┌───────────────────────┴───────────────┐
            ▼                                        ▼
   volcano.ts (R1/D-05 color dim,        uniprot-panel.ts (refactor:
   R6/D-06 Q/P metric — ONE              reuse same fetch+cache; the
   recompute fn keeps Y axis +           folded cache-uniprot todo)
   class + threshold lines + color
   in sync on every toggle)

   enrichment.ts (R2/D-10): split significant genes by sign of log2FC →
   two gGOSt() calls (up, down), same custom background → merge into ONE
   DataFrame + "Direction" column (Up/Down). R5/D-12: 'WP' in sources[].
            │
            ▼
   enrichment-viewers.ts (D-11): split Up/Down dot+bar side-by-side;
   onCurrentRowChanged cross-link unchanged (filter by Direction).
```

### Recommended Project Structure
```
src/
├── analysis/
│   ├── subcellular-location.ts   # NEW: fetch + verbatim classifier + cache
│   ├── enrichment.ts             # MODIFY: up/down split, 'WP' source, Direction col
│   └── differential-expression.ts# MODIFY: D-09 declared-contrast default
├── parsers/
│   └── spectronaut-candidates-parser.ts  # MODIFY: D-08 sign-flip (stay PURE)
├── viewers/
│   ├── volcano.ts                # MODIFY: metric param + color dim + sync recompute
│   └── enrichment-viewers.ts     # MODIFY: Direction-aware split charts
├── panels/
│   └── uniprot-panel.ts          # MODIFY: refactor fetch to reuse shared module
├── utils/
│   └── proteomics-types.ts       # MODIFY: add SEMTYPE.SUBCELLULAR_LOCATION
├── package.ts                    # MODIFY: Filter-viewer docking, color toggle menu
└── tests/                        # NEW test files per requirement (see Validation)
detectors.js                      # MODIFY: mirror SUBCELLULAR_LOCATION detector
```

### Pattern 1: Pure parser, shell orchestration in package.ts
**What:** Parsers return an unnamed `DG.DataFrame` with no `grok.shell.*` calls. The
import handler in `package.ts` names it, adds the table view, and docks viewers/filters.
**When to use:** R3 sign-flip lives in the parser (pure transform). R4 Filter-viewer
docking lives in `importSpectronautCandidates()` in `package.ts`.
**Why it matters:** `src/tests/spectronaut-candidates-parser.ts` calls
`parseSpectronautCandidatesText(text)` directly with no platform shell — any
`grok.shell` call in the parser breaks the entire test category. This is an established,
test-enforced contract (Phase 12).

### Pattern 2: Bulk Column.init / typed-array, never per-row col.set
**What:** Compute results into a `Float32Array`/`string[]` first, then `col.init((i) => arr[i])`.
**Source:** memory `feedback_dg_column_bulk_init` (measured 255× speedup on 8k-row volcano).
**Where:** R1 location column (8329 rows for BP DMD/WT), R3 sign-flip of log2FC + the two
AVG Group Quantity columns. Existing `volcano.ts:ensureDirectionColumn`,
`differential-expression.ts`, and the candidates parser already follow this — match it.
```typescript
// Source: existing src/viewers/volcano.ts:22-27 (verified pattern)
const col = df.columns.addNewFloat(colName);
col.init((i) => { /* read from a pre-built typed array, no col.set */ });
```
**Landmine:** memory `feedback_dg_column_init_null_sentinel` — `col.init(() => null)` on a
numeric column leaves the `FLOAT_NULL` sentinel (2.6789e-34) and `get()` reads it as a
finite positive. For the location *string* column this is moot; but the R3 log2FC flip
must read via `getRawData()` and write a new typed array, then `init` from it (the
candidates parser already reads `getRawData()` for `significant` — extend that block).

### Pattern 3: Live viewer-property toggle via custom DG.JsViewer OR menu switch
**What:** D-05/D-06 need a post-creation toggle that recomputes columns. Two viable wirings:
- **(A) Menu/dialog switch (simpler, recommended):** `createVolcanoPlot` stays a
  `ScatterPlotViewer` factory. A `Proteomics | Visualize | Volcano …` submenu (or a
  small `ui.dialog` / context-menu item on the viewer) flips a module-held state and
  calls one `recomputeVolcano(df, sp, {metric, colorDim})` function that re-runs
  `ensureNegLog10Column`/`ensureDirectionColumn` with the new metric, rebuilds the
  location color column if needed, redraws threshold lines, and sets `sp.props.color`.
- **(B) Custom `DG.JsViewer` wrapper:** register real viewer properties via
  `this.string('significanceMetric','adj.p-value',{choices:['adj.p-value','p-value']})`
  and `this.choices('colorBy', …)`, subscribe to `onPropertyValueChanged`
  (`d4-property-value-changed`), recompute on change. Heavier; only if the property
  panel is a hard UX requirement.
**Recommendation:** Start with (A). The locked text says "viewer colour property plus a
menu/dialog switch" (D-05) and "live viewer-property toggle … from the property panel"
(D-06) — `ScatterPlotViewer.props.color` IS a native property; the *metric* toggle is
the part with no native equivalent, so a menu/context-menu switch driving a single
`recompute` function satisfies both with the least surface. The verified API:
`sp.onPropertyValueChanged` exists (`viewer.ts:110`), `sp.props.color` is settable,
`sp.getProperties()` exists. Confirm exact ScatterPlot color property name at build
time via `sp.getOptions()` (it serializes the live `look`).

### Anti-Patterns to Avoid
- **Keying the UniProt TSV by header name.** The stream TSV header is `Entry`, not
  `accession`; column 2 header is `Subcellular location [CC]`. CK-omics splits each line
  by `\t` and reads **positionally** (`parts[0]`, `parts[1]`, …). Do the same — request
  fields in fixed order `accession,cc_subcellular_location,go_c,reviewed,gene_primary`
  and index positionally. Header-name lookups will silently fail.
- **Parsing the subcellular text structurally.** It contains `SUBCELLULAR LOCATION:`,
  `{ECO:0000269|PubMed:…}`, `Note=`, `[Isoform N]:` segments. CK-omics deliberately does
  a `\b`-bounded case-insensitive keyword *substring* scan over the raw string and takes
  the earliest-position match across all categories' keywords. Reproduce exactly — do
  not strip/normalize first (the `\b` boundaries already make `{ECO...}` and punctuation
  harmless, and keyword order-of-appearance IS the "UniProt order of importance" rule).
- **Re-running DE on a Candidates frame.** Candidates sets `proteomics.de_complete='true'`
  in the parser; the volcano/enrichment menu handlers gate on it. D-09's report-DE fix is
  for the *report* path (`differential-expression.ts` + DE dialog), NOT the Candidates
  path. Keep the two straight.
- **Mutating the source DataFrame in the enrichment cross-link.** `enrichment-viewers.ts`
  already clones for top-N (`createTopNEnrichmentDf`). The new Direction split must not
  break the existing `onCurrentRowChanged` → `proteinDf.selection` wiring.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Cross-session KV cache (D-02) | Custom IndexedDB wrapper | `grok.dapi.userDataStorage.put/get(name, mapJson)` | Platform-native, async, cross-device, no dep; verified `js-api/src/dapi.ts:735-760` |
| External HTTP (UniProt/g:Profiler) | `fetch()` | `grok.dapi.fetchProxy(url, opts)` | CORS-safe; CLAUDE.md hard rule; CONCERNS.md documents the raw-fetch defect this phase fixes |
| Comparison contrast picker UI (R4) | Modal numerator/denominator dropdowns | Native Datagrok **Filter** viewer on `Comparison` | D-07 locked; `tv.dockManager.dock` a filter; zero custom UI |
| Categorical color map on volcano | Manual per-point coloring | `col.meta.colors.setCategorical({cat: ARGB})` + `sp.props.color = colName` | Existing `volcano.ts:60-62` already does this for direction; reuse verbatim with the 11-color map |
| -log10 / threshold-line math | New formula code | Existing `ensureNegLog10Column`, `ensureDirectionColumn`, `df.meta.formulaLines` in `volcano.ts` | Already handle FLOAT_NULL, underflow sentinel, stale-line cleanup — parameterize the p-column instead of rewriting |
| FDR/Direction split in g:Profiler results | New API client | Existing `gGOSt()` called twice with different gene sets | `enrichment.ts:90-130` already correct; only the caller changes |

**Key insight:** Nearly all the Datagrok-side machinery already exists and is
test-covered. This phase is mostly *parameterizing existing functions* and *porting two
Python algorithms verbatim*. The risk is in faithful porting, not in new infrastructure.

---

## CK-omics Algorithm Port: Subcellular Classification (R1, D-04) — VERBATIM CONTRACT

`[VERIFIED: ~/Downloads/ck/CKomics_tool2.py:1384-1436 and
~/Downloads/ck/DMD_vs_WT/volcano_plots/Subcellular_Location_Classification_README.txt
read this session]`

### The 11-category keyword map (port exactly — these are the *Python source* lists,
which are the operative contract; the README prose paraphrases them)

```
Nucleus:          nucleus, nuclear, nucleoplasm, chromatin, nucleolus,
                  nuclear envelope, nuclear membrane, nuclear pore,
                  nuclear speck, cajal body, nuclear body
Cytoplasm:        cytoplasm, cytosol, cytoplasmic, cytoskeleton, microtubule,
                  actin, intermediate filament, sarcomere, myofibril, z disc,
                  z disk, spindle, centriole, centrosome, cilium, flagellum,
                  stress fiber
Mitochondria:     mitochondria, mitochondrial, mitochondrion,
                  mitochondrial matrix, mitochondrial membrane,
                  mitochondrial inner membrane, mitochondrial outer membrane,
                  mitochondrial intermembrane space
ER:               endoplasmic reticulum, sarcoplasmic reticulum,
                  rough endoplasmic reticulum, smooth endoplasmic reticulum,
                  er-golgi intermediate compartment, ergic
Golgi:            golgi, trans-golgi, golgi apparatus, golgi membrane,
                  cis-golgi, trans-golgi network, tgn
Plasma Membrane:  plasma membrane, cell membrane, cell surface, tight junction,
                  gap junction, adherens junction, desmosome, focal adhesion,
                  cell junction, synapse, postsynaptic, presynaptic,
                  neuromuscular junction
Lysosome:         lysosome, lysosomal, vacuole, lysosomal membrane,
                  autophagosome, phagosome
Peroxisome:       peroxisome, peroxisomal, glyoxysome
Ribosome:         ribosome, ribosomal, ribosomal protein, polysome
Extracellular:    extracellular, secreted, extracellular space,
                  extracellular matrix, basement membrane, basal lamina,
                  collagen, ecm
Vesicles:         vesicle, endosome, transport vesicle, secretory vesicle,
                  early endosome, late endosome, recycling endosome,
                  clathrin-coated vesicle, coated vesicle
```

### Exact hex palette (LOCKED — README §"COLOR SCHEME", confirms Python `get_location_colors`)

| Category | Hex | Category | Hex |
|----------|-----|----------|-----|
| Nucleus | `#FF6B6B` | Lysosome | `#F39C12` |
| Cytoplasm | `#ECDC44` | Peroxisome | `#E17055` |
| Mitochondria | `#45B7D1` | Ribosome | `#A29BFE` |
| ER | `#96CEB4` | Extracellular | `#FD79A8` |
| Golgi | `#FAAFFE` | Vesicles | `#FDCB6E` |
| Plasma Membrane | `#DDA0DD` | Unknown | `#CCCCCC` |

For `col.meta.colors.setCategorical`, convert each `#RRGGBB` to ARGB int
`0xFF000000 | parseInt(hex.slice(1),16)` (existing volcano code uses the same `0xFF…`
form at `volcano.ts:60-62`).

### The classification algorithm (`parse_subcellular_location`, port verbatim)

```
parseSubcellularLocation(subcellText, goText) -> category:
  build location_keywords map (above), preserving the dict insertion order
  (Nucleus, Cytoplasm, Mitochondria, ER, Golgi, Plasma Membrane, Lysosome,
   Peroxisome, Ribosome, Extracellular, Vesicles)

  findFirstLocation(text):
    if not text: return null
    textLower = text.toLowerCase()
    earliestMatch = null; earliestPos = +Infinity
    for [category, keywords] of location_keywords:
      for keyword of keywords:
        pattern = \b + escapeRegex(keyword) + \b      # word-boundary, case-insens
        m = first match of pattern in textLower
        if m and m.index < earliestPos:
          earliestPos = m.index; earliestMatch = category
    return earliestMatch

  loc = findFirstLocation(subcellText)   # subcellular field FIRST
  if loc: return loc
  loc = findFirstLocation(goText)        # GO cellular-component fallback
  if loc: return loc
  return 'Unknown'
```

**Critical semantics to preserve:**
- "Order of importance" = **earliest character position** of *any* keyword across *all*
  categories in the raw text. NOT category iteration order alone. Two categories' keywords
  are compared by `match.index`; ties keep the first-iterated category (JS `<` strict, so
  iteration order breaks ties — Python uses the same `< earliest_pos`). Preserve the map's
  insertion order so tie-breaks match.
- Word-boundary `\b...\b` on the **lowercased raw string** — do NOT pre-strip
  `SUBCELLULAR LOCATION:` / `{ECO...}` / `Note=`. Verified live: P04637's raw
  `cc_subcellular_location` is ~2KB of multi-isoform text; `\bcytoplasm\b` matches at
  the first "Cytoplasm" occurrence regardless of the surrounding `{ECO:...}`.
- Subcellular field is tried fully before GO is consulted at all (not interleaved).
- JS regex: use `new RegExp('\\b' + escaped + '\\b', 'i')` and `.search()` or
  `.exec().index`. Escape regex metachars in keywords (none of the current keywords
  contain metachars except `-` which is literal outside a class, but escape defensively
  to stay faithful to Python's `re.escape`).
- Multi-accession protein groups: CK-omics `assign_subcellular_locations` splits the
  protein-ID cell on `;`, takes the **first** UID present in the location_map (first
  non-missing wins). The package already extracts a "Primary Protein ID" via
  `parseAccession` (handles `sp|…|…`, `;`-delimited, `ups` suffix) — feed primary
  accessions to the fetch, map back per row.

### Fetch + priority-merge + reviewed-by-gene fallback
(`get_subcellular_locations_uniprot`, D-01/D-03)

```
chunk unique accessions (~100/chunk — CK-omics used 50 via search; D-01 says ~100, fine
  for stream which has no size cap, but keep the OR-query under ~UniProt URL limits;
  100 accessions ≈ 100*("accession:XXXXXX OR ")≈ 1.8 KB query → safe)
per chunk: GET (via fetchProxy)
  https://rest.uniprot.org/uniprotkb/stream
    ?query=(accession:A OR accession:B OR …) AND (taxonomy_id:<tax>)
    &fields=accession,cc_subcellular_location,go_c,reviewed,gene_primary
    &format=tsv
parse TSV positionally: parts[0]=accession parts[1]=subcell parts[2]=go_c
                        parts[3]=reviewed("reviewed"/"unreviewed") parts[4]=gene
loc = parseSubcellularLocation(parts[1], parts[2])
priority merge into locationMap[acc]:
   - if acc not seen: set
   - elif isReviewed and loc != 'Unknown': overwrite (reviewed beats unreviewed)
   - elif current == 'Unknown' and loc != 'Unknown': overwrite
record geneNameMap[acc] = gene  (for the D-03 fallback)

D-03 second pass — for every acc whose final loc == 'Unknown' AND has a gene name:
  group by gene; batch ~20 genes/query:
  query=((gene_exact:G1) OR (gene_exact:G2) …) AND (reviewed:true) AND (taxonomy_id:<tax>)
  fields=accession,gene_primary,cc_subcellular_location,go_c
  for each returned reviewed entry: loc = parse(...); if loc != 'Unknown':
     assign that loc to ALL the unreviewed accs that had this gene
```

`[ASSUMED]` — organism taxonomy filter: CK-omics maps organism code → taxonomy id
(`hsapiens 9606, mmusculus 10090, rnorvegicus 10116, dmelanogaster 7227`, default mouse).
The package's volcano path has no explicit organism selector (only the enrichment dialog
does, `ORGANISM_LIST` in `enrichment.ts`). **Open question A1** — see Assumptions Log:
the planner must decide where the organism for the location fetch comes from (reuse the
enrichment `ORGANISM_LIST` choice? infer? default?). The BP engagement is rat/mouse;
omitting the taxonomy filter risks cross-species accession collisions, but accessions are
globally unique in UniProt so the filter is mostly a safety/perf narrowing. Recommend:
make it optional — query without taxonomy filter if no organism is known (accessions are
unique), add the filter when an organism is available.

### Verified UniProt stream response shape (this session, P04637)

`[VERIFIED: live GET rest.uniprot.org/uniprotkb/stream this session]`
- TSV header line: `Entry⟶Subcellular location [CC]⟶Gene Ontology (cellular component)⟶Reviewed⟶Gene Names (primary)` (tab-separated; **header names ≠ field codes** — index positionally).
- `cc_subcellular_location` value is verbose: `SUBCELLULAR LOCATION: Cytoplasm {ECO:…}. Nucleus {ECO:…}. … ; SUBCELLULAR LOCATION: [Isoform 1]: …`. The `\b`-keyword scan handles it.
- `go_c` value: `centrosome [GO:0005813]; chromatin [GO:0000785]; cytoplasm [GO:0005737]; …` — also keyword-scannable as-is.
- `Reviewed` value: literal `reviewed` / `unreviewed`.
- First TSV line is the header — skip `lines[0]`, iterate `lines[1:]` (CK-omics does this).

---

## CK-omics Algorithm Port: create_subset_data Flip (R3, D-08) — VERBATIM CONTRACT

`[VERIFIED: ~/Downloads/ck/CKomics_tool2.py:1625-1662 read this session]`

CK-omics, given a desired `numerator / denominator`, checks the `Comparison
(group1/group2)` column: if the desired string is present, use rows as-is; if the
**flipped** string `denominator / numerator` is present instead, it:
1. `subset['AVG Log2 Ratio'] *= -1` (invert sign)
2. relabel `Comparison (group1/group2)` to the desired `numerator / denominator`
3. if `AVG Group Quantity Numerator` and `…Denominator` columns exist: **swap them**
4. if `Condition Numerator`/`Condition Denominator` columns exist: set to numerator/denominator

**D-08 adaptation (no per-run picker):** the parser has no user-chosen numerator. The
locked convention: **the declared `group1` in each row's own `A / B` string IS the
canonical numerator; positive log2FC must mean group1-enriched.** So per row:
- parse the comparison string as `group1 / group2` (the substring before/after ` / `;
  Spectronaut's literal column is `Comparison (group1/group2)`, value form `DMD / WT`).
- Spectronaut's `AVG Log2 Ratio` sign convention is the parked defect: empirically the
  report path emits the *reverse* of the declared numerator/denominator (memory
  `project_proteomics_spectronaut_de_direction_default`: log2FC correlated at −0.95
  across orientations). For **Candidates**, the file itself declares the comparison
  string, so the canonical rule is deterministic: define group1 = numerator, ensure
  positive log2FC = group1-higher. If the file's own sign convention disagrees with
  "group1-enriched is positive", flip that row's log2FC and swap the two AVG Group
  Quantity columns, exactly as `create_subset_data` does for the flipped case.

**Landmine — what "the file's convention" actually is:** Spectronaut Candidates'
`AVG Log2 Ratio` is conventionally **group1 / group2** already (numerator over
denominator as declared in the `Comparison (group1/group2)` header). The parked memory
defect is about the **report-import DE path** (`differential-expression.ts` picking
alphabetical condition order), NOT necessarily the Candidates file's own ratio. **Open
question A2 (see Assumptions Log):** confirm against a real BP DMD/WT Candidates file
whether `AVG Log2 Ratio` already matches the declared `group1/group2` orientation. If it
does, R3's parser flip is a *no-op for normal single-orientation files* and only
activates when a multi-comparison file contains *both* `DMD / WT` and `WT / DMD` rows
(the genuine `create_subset_data` flipped-string case). The plan must test both: (a)
file already canonical → unchanged; (b) file contains the reversed comparison string →
those rows flipped + columns swapped + relabeled. Do NOT unconditionally invert every
row — that would re-introduce the mirror defect.

Spectronaut AVG Group Quantity column names to swap (from CK-omics + Spectronaut
exports): `AVG Group Quantity Numerator` / `AVG Group Quantity Denominator`. The package's
candidates parser does not currently retain or rename these — add detection (they may be
absent in minimal exports; guard like CK-omics does).

Implementation: extend the existing `getRawData()` block in
`parseSpectronautCandidatesText` (lines 138-146) — it already bulk-reads `log2FC` and
`adj.p-value`. Build a per-row sign-multiplier `Float32Array` from the comparison-string
parse, write a flipped `log2FC` typed array + (conditionally) swapped quantity arrays,
`init` the columns from them, relabel the comparison string per row. **Stay pure** — no
`grok.shell`.

---

## CK-omics Algorithm Port: g:Profiler Up/Down Split (R2, D-10) — VERBATIM CONTRACT

`[VERIFIED: ~/Downloads/ck/CKomics_tool2.py:4450-4684, 4812-4907 read this session]`

CK-omics `run_gprofiler_analysis`:
1. Filter to significant: `(value_column < cutoff) & (abs(log2FC) >= log2fc_thr)` where
   `value_column = 'Qvalue' if value_type=='Q' else 'Pvalue'`. (Matches D-06: default
   metric is the q-value/adj.p; the package's `runEnrichmentPipeline` already filters on
   `cols.pValue` which is `adj.p-value`.)
2. **Split by sign:** `up_data = sig[log2FC > 0]`, `down_data = sig[log2FC < 0]`.
   (Strict `>0` / `<0`; exactly-zero excluded — none expected after the FC threshold.)
3. Background = **all detected proteins** (first gene per `;`-group, deduped), passed as
   custom background to *both* directional queries. The package's
   `runEnrichmentPipeline` already builds `backgroundGenes` from all mapped rows —
   reuse it unchanged for both calls (D-10 explicitly: background stays all-detected).
4. Two independent `gGOSt(upGenes, bg, organism, sources, p)` and `gGOSt(downGenes, …)`
   calls (no batching, one request each).
5. Merge: CK-omics keeps them as separate `Up_regulated`/`Down_regulated` sheets. D-10
   says **one DataFrame + `Direction` column** (`'Up'`/`'Down'`). Build per-direction
   result arrays, concat, add a `Direction` string column. Reuse `buildEnrichmentDf`
   shape; add the `Direction` column (categorical, 2 values).

**g:Profiler payload note (CK-omics vs package):** CK-omics passes `'highlight': True`
and post-filters: KEGG/REAC/**WP** keep ALL terms (highlight filter "always returns
nothing for large gene lists"); GO sources keep only `highlighted==True`. The package's
`gGOSt` currently uses `all_results:true`, `significance_threshold_method:'fdr'`,
`domain_scope:'custom'`, `background:[...]`, NOT `highlight`. D-10 says "reuses Phase-9
… unchanged" and background stays the existing custom-background behaviour — so **do not
adopt CK-omics' highlight-filtering**; only split the query by direction and add `'WP'`.
The smart-pathway / highlight filtering is explicitly Phase 14 R5 (`ROADMAP.md` Phase
14), out of scope here. Keep `enrichment.ts`'s existing request body; change only the
caller (two gene sets) and `sources` (add `'WP'`).

---

## R4: Comparison Filter Viewer Docking (D-07)

`[CITED: js-api viewer/dock patterns + existing enrichment-viewers.ts dock usage]`

- The candidates parser keeps `Comparison (group1/group2)` (it currently does, line
  85-89 / test "keeps Comparison column"). After D-08 sign-normalization the column is
  relabeled per row but still distinguishes multiple contrasts.
- In `package.ts:importSpectronautCandidates` (currently lines 156-172): after
  `grok.shell.addTableView(df)`, count distinct non-null values of the comparison column.
  If > 1, dock a native **Filter** viewer. Datagrok's Filters viewer is
  `DG.VIEWER.FILTERS` / `tv.addViewer('Filters')` or via `tv.dockManager.dock`. The
  established docking idiom in this package is `enrichment-viewers.ts:164-165`:
  `tv.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Title', 0.3)`. Filter the
  filter to the single `Comparison` column via the viewer's `columnNames` look property.
- Single-comparison files: skip the filter (D-07).
- The volcano reads the *filtered* rows automatically (Datagrok viewers respect
  `df.filter`) — no extra wiring; this is why D-07 chose the native filter over a modal.
- **Pitfall:** `createVolcanoPlot` adds `df.meta.formulaLines` and the direction column
  over the *whole* frame; with a Comparison filter the threshold lines are still correct
  (formula-based, not data-derived) but the up/down counts shown elsewhere reflect all
  rows. Acceptable for this phase (no counters until Phase 14 R2). Document so the
  verifier doesn't flag it.

---

## R5: WikiPathways Source (D-12) — one-line change

`[VERIFIED: enrichment.ts:90-130 + CKomics_tool2.py:478-486 read this session]`

g:Profiler `/api/gost/profile/` accepts `'WP'` in `sources[]` (CK-omics ships it
default-on; INTEGRATIONS.md confirms the endpoint supports WikiPathways). Change:
1. `showEnrichmentDialog` (enrichment.ts:348-353): add
   `const wpInput = ui.input.bool('WikiPathways', {value: true});` and `.add(wpInput)`.
2. OK handler (enrichment.ts:391-396): `if (wpInput.value) selectedSources.push('WP');`.
That is the entire change for R5. No API-client change (`gGOSt` passes `sources`
through). CK-omics' source code literal is exactly `'WP'`.

---

## Volcano Q/P Live Toggle (D-06) and Color Dimension (D-05)

`[VERIFIED: volcano.ts full file + js-api viewer.ts:110,146-166 read this session]`

Current `volcano.ts` hardwires the y/significance to `adj.p-value`:
- `ensureNegLog10Column(df)` reads `df.col('adj.p-value')` (line 18) — **parameterize**
  to take the p-column name (`'adj.p-value'` | `'p-value'`).
- `ensureDirectionColumn(df, fc, p)` reads `df.col('adj.p-value')` (line 42) —
  **parameterize** the same way; it already updates in place (line 50, "update existing
  column in place rather than remove + re-add, which would invalidate viewer/grid
  bindings" — preserve this; do NOT remove/re-add on toggle).
- Threshold lines (lines 91-118): the horizontal line is `-log10(pThreshold)` on the
  neg-log column; recompute and re-add (the existing code already filters stale lines by
  formula prefix — reuse that cleanup so a Q→P toggle replaces lines, not stacks them).
- The `p-value` column is optional in Candidates (parser only adds it `if (pValCol)`).
  Toggle to `p-value` must be disabled/guarded when the column is absent
  (`df.col('p-value') == null`).

**One synchronized recompute fn** (the load-bearing D-06 requirement — "dots and
colouring never disagree"):
```typescript
function recomputeVolcano(df, sp, metric: 'adj.p-value'|'p-value',
                          colorDim: 'significance'|'location', fcThr, pThr) {
  const yCol = ensureNegLog10Column(df, metric);          // parameterized
  const dirCol = ensureDirectionColumn(df, fcThr, pThr, metric); // parameterized
  // location color column built once (R1), reused on toggle:
  const colorCol = colorDim === 'location' ? ensureLocationColumn(df) : dirCol;
  redrawThresholdLines(df, yCol, fcThr, pThr);            // existing cleanup+add
  sp.props.color = colorCol;
  // sp.props.y must point at the (single, in-place-updated) neg-log column;
  // ensureNegLog10Column uses a FIXED column name so the y binding never breaks.
}
```
**Landmine:** `ensureNegLog10Column` currently uses fixed name `'negLog10padj'` and
early-returns if it exists (line 14-16). For a live metric toggle the *contents* must
change while the *name/binding* stays fixed. Change it to **always recompute the values
in place** (read the chosen p-column, re-`init` the existing column) instead of
early-returning — otherwise toggling Q→P leaves the Y axis showing stale Q values. Same
fixed-name-in-place discipline as `ensureDirectionColumn` already uses. Keep the column
name generic (e.g. `negLog10p`) since it's no longer always padj.

**Toggle wiring (Discretion):** add a context-menu item on the volcano
(`sp.onContextMenu`, verified `viewer.ts:308`) or a small `Proteomics | Visualize |
Volcano Options…` dialog with two `ui.input.choice` (metric, color-by) that calls
`recomputeVolcano`. Default metric `adj.p-value`, default color `significance`. The
native `sp.props.color` change for the color dimension is itself a real viewer property
(shows in the property panel); the *metric* needs the explicit recompute because it
touches Y + classification + lines together.

---

## Report-DE Direction Fix (D-09)

`[VERIFIED: differential-expression.ts:268-433 + memory project_proteomics_spectronaut_de_direction_default]`

The defect: for **Spectronaut report** imports, `showDEDialog` derives groups from
`getGroups(df)` and the `Comparison` choice defaults to `pairs[0] = "${g2.name} vs
${g1.name}"` (line 289). Group order in `proteomics.groups` comes from the annotate step,
which for Spectronaut report parsing ends up alphabetical — the reverse of Spectronaut's
declared Numerator/Denominator. The per-run workaround is selecting the other dropdown
item.

D-09 scope: honor the **declared/intended contrast** instead of the alphabetical default,
and retire the manual workaround. Concretely:
- The Spectronaut report parser (`spectronaut-parser.ts`) or the annotate step must
  record the declared numerator/denominator orientation (Spectronaut reports carry
  condition/comparison metadata; the Candidates parser sees `Comparison (group1/group2)`
  directly). Where the report path stores group order, ensure group1/group2 reflect the
  declared contrast, OR set the DE dialog's default `Comparison` to the declared
  orientation rather than `pairs[0]`.
- "Retire the manual workaround": the `Comparison` dropdown can stay as an override but
  its **default** must be correct so the unaware analyst gets the right sign.
- **Open question A3 (Assumptions Log):** the exact place the report path's group order
  is set (annotate-from-parser vs `experiment-setup.ts`) needs to be traced during
  planning — `spectronaut-parser.ts` and `experiment-setup.ts` were not in this phase's
  read set. The fix is small but its *location* depends on where the orientation
  metadata survives. Recommend the plan include a Wave-0 trace task.
- This is a **direction-only** change (memory: "does not move the hit-list metrics") —
  verification is sign/orientation, not magnitude.

---

## Common Pitfalls

### Pitfall 1: UniProt TSV keyed by header name
**What goes wrong:** Code does `row['accession']`; header is `Entry`. Silent empty map.
**Why:** UniProt stream returns *display* headers, not field codes.
**How to avoid:** Request fixed `fields=` order, split each line on `\t`, index by
position (`parts[0..4]`). Skip line 0 (header).
**Warning signs:** All proteins classify Unknown; zero cache hits on second run.

### Pitfall 2: Stripping/normalizing the subcellular text before keyword scan
**What goes wrong:** Removing `{ECO:...}` or splitting on `SUBCELLULAR LOCATION:` changes
character positions → wrong "first location" → wrong category → parity regression.
**Why:** CK-omics' "order of importance" is literally earliest substring position in the
raw lowercased text.
**How to avoid:** Port `find_first_location` byte-for-byte: lowercase the raw string,
`\b`-bounded regex per keyword, min `match.index` across all categories, map insertion
order for ties. Subcell field fully before GO.
**Warning signs:** Category distribution differs from the client's two reference HTMLs.

### Pitfall 3: Parser does shell work (R3/R4)
**What goes wrong:** Adding the Filter viewer or `addTableView` inside
`parseSpectronautCandidatesText` breaks the entire `SpectronautCandidates` test category
(tests call the parser with no platform shell).
**How to avoid:** Sign-flip in parser (pure). Filter docking in
`package.ts:importSpectronautCandidates`.
**Warning signs:** `grok test --category SpectronautCandidates` throws `grok.shell is
undefined`.

### Pitfall 4: Unconditional log2FC inversion (R3)
**What goes wrong:** Flipping every row's sign re-introduces the mirror defect for
already-canonical files.
**Why:** R3 only flips rows whose declared comparison is the *reverse* of canonical
group1/group2 (the `create_subset_data` flipped-string branch). A normal single-direction
Candidates file is already canonical.
**How to avoid:** Per-row decision from the parsed comparison string; test both the
canonical and the reversed-string fixtures.
**Warning signs:** Volcano is a left-right mirror of the client's; log2FC sign opposite
on known DMD markers (DMD itself should be the strongest signal — concordance note).

### Pitfall 5: `ensureNegLog10Column` early-return defeats the live toggle (D-06)
**What goes wrong:** Existing fn early-returns if the column exists → toggling Q↔P leaves
stale Y values; dots and colouring disagree (violates D-06's core requirement).
**How to avoid:** Make it always recompute values in place (fixed column name, re-`init`
from the chosen p-column). Same in-place discipline `ensureDirectionColumn` already uses.
**Warning signs:** Y-axis scale doesn't change when toggling metric; threshold line
moves but points don't.

### Pitfall 6: `col.init(() => null)` numeric sentinel
**What goes wrong:** memory `feedback_dg_column_init_null_sentinel` — numeric column
init'd with null leaves FLOAT_NULL (2.6789e-34) read back as finite positive.
**How to avoid:** For R3's flipped log2FC, write a real `Float32Array` (copy via
`getRawData()`, negate selected rows), `init` from it; for null cells use `isNone()` /
explicit FLOAT_NULL handling exactly as the existing candidates `significant` block does.

### Pitfall 7: userDataStorage value-size / single-value-is-string
**What goes wrong:** `grok.dapi.userDataStorage.put(name, mapObj)` stores a map of
string→string; large location maps (8329 accessions) JSON-encoded must fit. Treat as one
value or chunk by accession-prefix.
**How to avoid:** Recommend `userDataStorage.put('proteomics-subcell-loc', {acc: cat})`
(post a serialized map). Cache invalidation: include a small `__v` schema-version key;
bump it if the keyword map ever changes (D-04 says it won't — it's a contract — so a
fixed version is fine). Validate size assumption in a Wave-0 spike if the planner is
unsure (Open question A4).

## Code Examples

### Categorical color map on the location column (reuse existing pattern)
```typescript
// Source: existing src/viewers/volcano.ts:60-62 (verified) — same idiom, 11 colors
const LOCATION_COLORS: Record<string, number> = {
  'Nucleus': 0xFFFF6B6B, 'Cytoplasm': 0xFFECDC44, 'Mitochondria': 0xFF45B7D1,
  'ER': 0xFF96CEB4, 'Golgi': 0xFFFAAFFE, 'Plasma Membrane': 0xFFDDA0DD,
  'Lysosome': 0xFFF39C12, 'Peroxisome': 0xFFE17055, 'Ribosome': 0xFFA29BFE,
  'Extracellular': 0xFFFD79A8, 'Vesicles': 0xFFFDCB6E, 'Unknown': 0xFFCCCCCC,
};
const col = df.col('Subcellular Location') ?? df.columns.addNewString('Subcellular Location');
col.init((i) => locByRow[i]);           // locByRow: pre-built string[] (bulk pattern)
col.meta.colors.setCategorical(LOCATION_COLORS);
col.semType = SEMTYPE.SUBCELLULAR_LOCATION;
```

### fetchProxy to UniProt stream (CORS-safe, per CLAUDE.md)
```typescript
// Source: existing src/panels/uniprot-panel.ts:81 + enrichment.ts:61 (verified idiom)
const q = accessions.map((a) => `accession:${a}`).join(' OR ');
const url = 'https://rest.uniprot.org/uniprotkb/stream'
  + `?query=(${encodeURIComponent(q)})`
  + '&fields=accession,cc_subcellular_location,go_c,reviewed,gene_primary'
  + '&format=tsv';
const resp = await grok.dapi.fetchProxy(url);
if (!resp.ok) { /* warn, continue with what we have */ }
const lines = (await resp.text()).trim().split('\n');
for (let li = 1; li < lines.length; li++) {       // skip header line 0
  const p = lines[li].split('\t');
  const acc = p[0], subcell = p[1] ?? '', goc = p[2] ?? '',
        reviewed = (p[3] ?? '').toLowerCase() === 'reviewed', gene = p[4] ?? '';
  const cat = parseSubcellularLocation(subcell, goc); // verbatim port
  /* priority-merge as in CK-omics get_subcellular_locations_uniprot */
}
```

### Persistent cache via platform KV
```typescript
// Source: js-api/src/dapi.ts:735-760 (verified UserDataStorage API)
const STORE = 'proteomics-subcell-loc';
const cached = await grok.dapi.userDataStorage.get(STORE) ?? {}; // {acc: category}
// ... compute misses, fetch ...
await grok.dapi.userDataStorage.put(STORE, {...cached, ...fetched}); // write-through
```

## State of the Art

| Old Approach | Current Approach | When | Impact |
|--------------|------------------|------|--------|
| `uniprot-panel.ts` raw per-row `fetch()` of `/uniprotkb/{acc}.json`, no cache | Batched `stream` TSV via `fetchProxy`, persistent cache (D-01/D-02) | This phase | Fixes CONCERNS.md raw-fetch CORS defect + closes cache-uniprot todo; panel refactored to share the module |
| Spectronaut report DE picks alphabetical condition order (mirror defect) | Honor declared/intended contrast (D-09) | This phase | Volcano sign correct by default; manual dropdown demoted to override |
| Single combined enrichment query | Up/down split, `Direction` column (D-10) | This phase | Matches CK-omics directional deliverable |

**Deprecated/outdated:** none — this is additive parity work on a shipped package.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Organism/taxonomy for the location fetch can be optional (accessions are globally unique; add `taxonomy_id` filter only when an organism is known, e.g. from the enrichment `ORGANISM_LIST`) | Subcellular fetch | LOW — without the filter, results are still correct (unique accessions); only loses CK-omics' minor narrowing. Planner must decide the organism source. |
| A2 | Spectronaut Candidates `AVG Log2 Ratio` is conventionally already group1/group2 (declared-numerator-over-denominator), so R3's flip is a *no-op* for normal single-orientation files and only activates for the `create_subset_data` reversed-comparison-string case | create_subset_data Flip | HIGH — if the Candidates file's own sign is reversed vs its declared string, the flip rule must invert; verify against a real BP DMD/WT Candidates file before locking R3 logic. Wrong → mirror regression. |
| A3 | The report-DE group-order origin (D-09) is in the Spectronaut report parser / annotate step, not yet traced | Report-DE Direction Fix | MEDIUM — fix is small but its location is unconfirmed; recommend a Wave-0 trace task in the plan. |
| A4 | `grok.dapi.userDataStorage` can hold a ~8329-entry JSON map as one value | Cache Mechanism | LOW-MEDIUM — if a size limit bites, fall back to accession-prefix sharding or IndexedDB (Claude's discretion per D-02). Recommend a tiny Wave-0 spike. |
| A5 | g:Profiler accepts `'WP'` in `sources[]` on the existing `/api/gost/profile/` request body unchanged | R5 WikiPathways | LOW — INTEGRATIONS.md + CK-omics both confirm; verified by CK-omics source using exactly `'WP'`. |
| A6 | `AVG Group Quantity Numerator`/`Denominator` columns may be absent in minimal Candidates exports | create_subset_data Flip | LOW — CK-omics guards with `if … in subset.columns`; port the same guard (swap only if both present). |

## Open Questions (RESOLVED)

All four assumption-log questions below are **resolved**. A2, A3, and A4 are empirical /
runtime questions that cannot be answered from documentation alone — they require
inspecting the live engagement data and probing `userDataStorage` capacity. The
**resolution strategy is the Wave-0 task in 13-01** (Task 2: "Resolve A3/A2/A4; write
13-WAVE0-FINDINGS.md"), whose output `13-WAVE0-FINDINGS.md` is the authoritative,
committed answer that 13-03 (R3 rule), 13-04 (cache), and 13-05 (D-09 site) consume.
Deferring these to a gated Wave-0 task **is the chosen resolution** — not an oversight —
because answering them requires the very files and platform probe that Wave-0 performs
before any dependent plan executes. A1 is resolved here directly (no Wave-0 needed).

1. **Where does the report-DE group order get set? (A3)**
   - Known: defect is alphabetical condition order vs declared Numerator/Denominator.
   - Unclear: exact file/function (spectronaut-parser.ts vs experiment-setup.ts).
   - **RESOLVED (deferred to Wave-0):** 13-01 Task 2 traces the exact
     file:function:line and records the locked D-09 fix site in `13-WAVE0-FINDINGS.md`.
     13-05 obeys that finding. This deferral is the resolution: the trace requires
     reading the live parser/dialog code, which Wave-0 does before 13-05 runs.

2. **Candidates `AVG Log2 Ratio` actual sign convention (A2)**
   - Known: `create_subset_data` flips only the reversed-comparison-string case.
   - Unclear: whether the file's own ratio already matches its declared `group1/group2`.
   - **RESOLVED (deferred to Wave-0):** 13-01 Task 2 verifies the real BP DMD/WT
     Candidates file's `AVG Log2 Ratio` sign vs its declared comparison (DMD = known
     strong up signal per the concordance note) and records the definitive per-row R3
     flip rule in `13-WAVE0-FINDINGS.md`; 13-03 gates its logic on that finding. The
     empirical check against engagement data is the resolution strategy by design.

3. **Organism source for the location fetch (A1)**
   - **RESOLVED (here, no Wave-0 needed):** the taxonomy filter is **optional**.
     Accessions are globally unique, so the location fetch is correct without an
     organism. Apply a `taxonomy_id` filter only when an organism is known (reuse the
     enrichment dialog's `ORGANISM_LIST` choice if the fetch is triggered from a
     dialog); **omit the filter entirely when the organism is unknown**. No further
     investigation required — this matches Assumptions-Log A1 (LOW risk).

4. **`userDataStorage` capacity for the ~8329-entry accession→location map (A4)**
   - **RESOLVED (deferred to Wave-0):** 13-01 Task 2 confirms whether a ~250 KB
     `{accession:category}` JSON is safe as one `userDataStorage.put` value (cite
     js-api/src/dapi.ts:735-760) or needs prefix-sharding / IndexedDB, and records the
     chosen cache mechanism + `__v` invalidation key for 13-04. The capacity probe
     requires the running platform; deferring it to Wave-0 is the resolution.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| UniProt REST `stream` endpoint | R1 fetch | ✓ (verified live this session) | n/a (public, no version) | Per-protein JSON (current panel path) — slower, no D-01 batching |
| g:Profiler `/api/gost/profile/` | R2/R5 | ✓ (existing integration, INTEGRATIONS.md) | n/a | none — enrichment already depends on it |
| `grok.dapi.fetchProxy` | all external HTTP | ✓ (js-api, used in 2 existing files) | linked local | none — CLAUDE.md mandates it |
| `grok.dapi.userDataStorage` | D-02 cache | ✓ (js-api/src/dapi.ts:735-760, verified) | linked local | IndexedDB (D-02 discretion) |
| Running Datagrok instance | `grok test` | assumed (probe `~/.grok` per memory `feedback_verify_server_before_claiming_untestable`) | — | parser/classifier are pure → testable headlessly via `@datagrok-libraries/test` |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** UniProt stream (fallback = existing per-protein
path, but D-01 locks stream — fallback only for resilience, not the primary path).

## Validation Architecture

Nyquist validation is enabled (`.planning/config.json` has no
`workflow.nyquist_validation:false`; `verifier:true`). Tests use
`@datagrok-libraries/test` (`category/test/expect`), one file per area under
`src/tests/`, registered through `src/package-test.ts`. Pure functions (parser,
classifier) are testable headlessly; viewer/shell behaviors need a running instance via
`grok test`.

### Test Framework
| Property | Value |
|----------|-------|
| Framework | `@datagrok-libraries/test` (linked local) — `category()`, `test()`, `expect()` |
| Config file | none — entry is `src/package-test.ts` (exports tests array) |
| Quick run command | `grok test --category <Cat> --skip-build` (rebuild if test files changed — memory `feedback_grok_test_skipbuild_stale`) |
| Full suite command | `grok test` (against `~/.grok` default server; probe `/api/info/server` first) |

### Phase Requirements → Test Map
| Req | Behavior | Test Type | Automated Command | File |
|-----|----------|-----------|-------------------|------|
| R1 | `parseSubcellularLocation` returns correct category for verbatim-port keyword/ordering cases (subcell-first, GO fallback, earliest-position tie-break, Unknown) | unit (pure) | `grok test --category SubcellularLocation` | ❌ Wave 0: `src/tests/subcellular-location.ts` |
| R1 | 11-color map ARGB values match the locked hex palette exactly | unit (pure) | same | ❌ Wave 0 (same file) |
| R1 | UniProt TSV positional parse (header skipped, `parts[0..4]`, reviewed-priority merge, D-03 gene fallback path) — fixture TSV string, no network | unit (pure, inject parsed-lines) | same | ❌ Wave 0 (same file) |
| R2 | enrichment splits significant genes by sign; merged DataFrame has `Direction` ∈ {Up,Down}; background = all detected (both calls) | unit + integration | `grok test --category Enrichment` | ⚠️ extend `src/tests/enrichment.ts` |
| R3 | canonical single-orientation Candidates file unchanged; reversed-comparison-string rows get log2FC negated + AVG Group Qty cols swapped + relabeled; quantity-cols-absent guarded | unit (pure) | `grok test --category SpectronautCandidates` | ⚠️ extend `src/tests/spectronaut-candidates-parser.ts` |
| R4 | >1 distinct Comparison → Filter viewer docked; single-comparison → not docked | integration (needs shell) | `grok test --category SpectronautCandidatesE2E` | ⚠️ extend `src/tests/spectronaut-candidates-e2e.ts` |
| R5 | `'WP'` added to sources when WikiPathways input checked; default value true | unit | `grok test --category Enrichment` | ⚠️ extend `src/tests/enrichment.ts` |
| R6 | toggling metric adj.p↔p recomputes negLog column values in place (binding name stable), reclassifies up/down/NS, moves threshold line — all consistent; p-value toggle guarded when column absent | unit (pure on volcano helpers) + integration | `grok test --category Volcano` | ❌ Wave 0: `src/tests/volcano.ts` |
| D-09 | report-DE default contrast = declared/intended orientation (sign), not alphabetical | unit/integration (after A3 trace) | `grok test --category Analysis` | ⚠️ extend `src/tests/analysis.ts` |

### Sampling Rate
- **Per task commit:** `grok test --category <relevant>` (rebuild first if test files
  changed — `--skip-build` reuses a stale bundle and new tests return null/0).
- **Per wave merge:** full `grok test`.
- **Phase gate:** full suite green before `/gsd:verify-work`.

### Wave 0 Gaps
- [ ] `src/tests/subcellular-location.ts` — covers R1 (classifier + palette + TSV parse). Pure, headless.
- [ ] `src/tests/volcano.ts` — covers R6 + D-05 (parameterized helpers; in-place recompute).
- [ ] Extend `src/tests/spectronaut-candidates-parser.ts` — R3 canonical-unchanged + reversed-flip + guard.
- [ ] Extend `src/tests/enrichment.ts` — R2 split/Direction/background + R5 WP source.
- [ ] Extend `src/tests/spectronaut-candidates-e2e.ts` — R4 filter docking (needs running instance).
- [ ] Extend `src/tests/analysis.ts` — D-09 (after the A3 trace task identifies the fix site).
- [ ] Wave-0 trace task: locate report-DE group-order origin (A3).
- [ ] Wave-0 spike: verify a real BP DMD/WT Candidates file's `AVG Log2 Ratio` sign vs declared comparison (A2); confirm userDataStorage size for ~8329-entry map (A4).
- [ ] Register every new test file in `src/package-test.ts`.

## Security Domain

`security_enforcement` is not set in `.planning/config.json` (treat as enabled, but this
package's CLAUDE.md states "Auth / credentials: none" — only public REST APIs).

### Applicable ASVS Categories
| ASVS Category | Applies | Standard Control |
|---------------|---------|------------------|
| V2 Authentication | no | No auth — public UniProt/g:Profiler APIs; inherits platform session |
| V3 Session Management | no | Platform-managed |
| V4 Access Control | no | No new access surface |
| V5 Input Validation | yes | UniProt TSV / g:Profiler JSON are external; parse defensively (positional split, length-guard `parts[]`, treat missing fields as empty — CK-omics already does `if len(parts) > N`) |
| V6 Cryptography | no | No secrets; no crypto |

### Known Threat Patterns
| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Untrusted external response shape drift (UniProt header/field change) | Tampering/DoS | Positional parse with bounds checks; on malformed/empty → category 'Unknown', never throw out of the fetch loop (CK-omics pattern: catch per-batch, continue) |
| CORS bypass via raw fetch | — (platform policy) | Mandatory `grok.dapi.fetchProxy` (CLAUDE.md); this phase *removes* the existing raw-fetch defect |
| Cache poisoning of location map | Tampering | `userDataStorage` is per-user platform storage (not attacker-writable); include schema-version key; values are categories from a fixed enum — validate against the 12-name set on read |
| g:Profiler request injection via gene names | Tampering | Gene symbols come from the parsed DataFrame (already typed), JSON body via `JSON.stringify` (existing `gGOSt` — no string concatenation into the request) |

## Sources

### Primary (HIGH confidence)
- `~/Downloads/ck/CKomics_tool2.py` lines 476-489, 1148-1436, 1589-1662, 3850-3856, 4450-4684, 4812-4907 — read directly this session (verbatim algorithm contracts)
- `~/Downloads/ck/DMD_vs_WT/volcano_plots/Subcellular_Location_Classification_README.txt` — full file read (locked palette/keyword contract)
- `~/Downloads/ck/Spectronaut-Concordance-Note-2026-05-15.md` — full file read (default-metric rationale, DMD-is-strongest sanity check)
- Live `GET https://rest.uniprot.org/uniprotkb/stream?...&format=tsv` for P04637 — verified TSV header (`Entry`, not `accession`) and verbose multi-isoform `cc_subcellular_location` shape
- Package source read this session: `src/viewers/volcano.ts`, `src/panels/uniprot-panel.ts`, `src/analysis/enrichment.ts`, `src/analysis/differential-expression.ts`, `src/parsers/spectronaut-candidates-parser.ts`, `src/viewers/enrichment-viewers.ts`, `src/utils/proteomics-types.ts`, `src/utils/column-detection.ts`, `detectors.js`, `src/package.ts` (import/volcano/enrichment handlers), `src/tests/spectronaut-candidates-parser.ts`
- `js-api/src/dapi.ts:735-760` (UserDataStorage), `js-api/src/viewer.ts:110,146-166,308,430-490` (Viewer/ScatterPlot/JsViewer property + event API)
- `.planning/codebase/INTEGRATIONS.md`; `.planning/phases/13-…/13-CONTEXT.md`; `.planning/phases/13-…/13-DISCUSSION-LOG.md`; `.planning/ROADMAP.md` (Phase 13/14 boundary)
- memory: `project_proteomics_spectronaut_de_direction_default`, `feedback_dg_column_bulk_init`, `feedback_dg_column_init_null_sentinel`, `feedback_grok_test_skipbuild_stale`, `feedback_verify_server_before_claiming_untestable`

### Secondary (MEDIUM confidence)
- Phase-9 `09-CONTEXT.md` (enrichment dot/bar + cross-link pattern R2/D-11 extends)

### Tertiary (LOW confidence)
- None — all claims traced to primary sources read this session. Remaining uncertainties
  are explicitly in the Assumptions Log / Open Questions for Wave-0 resolution, not
  unverified assertions.

## Metadata

**Confidence breakdown:**
- CK-omics algorithm contracts (R1/R2/R3): HIGH — Python source + README read verbatim.
- UniProt stream response shape: HIGH — verified by live request this session.
- Datagrok integration points (volcano helpers, parser purity, cache API, viewer props):
  HIGH — all source files read directly.
- Candidates file actual sign convention (A2) & report-DE fix site (A3): MEDIUM —
  flagged for Wave-0 verification; impacts R3/D-09 logic placement, not feasibility.

**Research date:** 2026-05-16
**Valid until:** 2026-06-15 for the Datagrok/codebase findings (stable, shipped package);
the UniProt stream header/field shape is a public-API contract — re-verify only if a
fetch returns all-Unknown (Pitfall 1 warning sign).

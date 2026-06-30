# Phase 4: Annotation and Alternatives - Research

**Researched:** 2026-03-01
**Domain:** UniProt info panels, DEqMS differential expression, Datagrok panel system
**Confidence:** HIGH

## Summary

Phase 4 has two independent workstreams: (1) a UniProt info panel that appears in the context panel when a user clicks a protein ID cell, and (2) a DEqMS R script for peptide-count-weighted differential expression as an alternative to the existing limma DE.

The UniProt info panel follows a well-established Datagrok pattern used extensively in the Chem and Bio packages -- register a function with `@grok.decorators.panel` and `meta: {role: 'widgets'}`, accepting a string parameter constrained to the semantic type `Proteomics-ProteinId`. The panel fetches data from the UniProt REST API (`rest.uniprot.org/uniprotkb/{accession}.json`) and renders protein name, function, and GO terms in a widget. The DEqMS R script extends the existing limma pattern (already in `scripts/limma_de.R`) by adding a `spectraCounteBayes()` call after `eBayes()`, using the "Peptides" or "Unique peptides" column from MaxQuant as the count input.

**Primary recommendation:** Implement the UniProt panel as an async `DG.Widget` using `@grok.decorators.panel` with semantic type binding, and create `scripts/deqms_de.R` following the identical input/output contract as `limma_de.R` but with an additional `peptideCount` integer vector input.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| ANNOT-01 | User can click a protein ID in the grid and see UniProt details (name, function, GO terms) in the context panel | Panel decorator pattern from Chem/Bio packages; UniProt REST API JSON structure verified; `ui.tableFromMap` + `ui.wait` for async rendering |
| ANNOT-02 | UniProt info panel triggers automatically for columns with Proteomics:ProteinId semantic type | Semantic type binding via `@grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}})` -- same pattern as Chem molecule panels |
| ANLY-04 | User can run DEqMS differential expression (peptide-count-weighted) via R script | DEqMS Bioconductor package extends limma with `spectraCounteBayes()`; existing limma_de.R provides the template; MaxQuant demo data already has "Peptides" column |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| UniProt REST API | v2 (rest.uniprot.org) | Protein annotation lookup | Official UniProt programmatic access; JSON format |
| DEqMS | 1.28.0+ (Bioconductor) | Peptide-count-weighted DE | Standard proteomics DE method; extends limma |
| limma | 3.58+ (Bioconductor) | Linear model fitting (DEqMS prerequisite) | Already used in limma_de.R |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| datagrok-api | ^1.25.0 | Panel registration, UI widgets, DataFrame | Core platform API (already in package.json) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| UniProt REST API | Proteins API (ebi.ac.uk) | UniProt REST is authoritative; Proteins API is a wrapper |
| DEqMS | MSqRob2, proDA | DEqMS is simpler (extends limma directly), well-cited, matches existing limma pattern |

## Architecture Patterns

### Recommended Project Structure
```
packages/Proteomics/
  src/
    panels/
      uniprot-panel.ts      # UniProt info panel widget
    analysis/
      differential-expression.ts  # Add DEqMS dialog option
  scripts/
    limma_de.R              # Existing limma script
    deqms_de.R              # New DEqMS script (same I/O contract + peptide counts)
```

### Pattern 1: Semantic Type Info Panel (Cell-Level)
**What:** A function decorated with `@grok.decorators.panel` that receives a cell value (string) when the user clicks a cell with a matching semantic type. The platform automatically shows the returned `DG.Widget` in the context panel.
**When to use:** Any time you want to show details for a cell value based on its semantic type.
**Example:**
```typescript
// Source: Verified from Chem package (packages/Chem/src/package.ts lines 1325-1332)
// and Bio package (packages/Bio/src/package.ts lines 719-728)

// In package.ts:
@grok.decorators.panel({
  name: 'Proteomics | UniProt',
  description: 'UniProt protein details',
  meta: {role: 'widgets', domain: 'bio'},
})
static uniprotPanel(
  @grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}}) proteinId: string,
): DG.Widget {
  return new DG.Widget(ui.wait(async () => {
    const data = await fetchUniProtData(proteinId);
    return renderUniProtWidget(data);
  }));
}
```

```typescript
// Corresponding entry in package.g.ts (must be manually added):
//name: Proteomics | UniProt
//description: UniProt protein details
//input: string proteinId { semType: Proteomics-ProteinId }
//output: widget result
//meta.role: panel
//meta.domain: bio
export function uniprotPanel(proteinId: string): DG.Widget {
  return PackageFunctions.uniprotPanel(proteinId);
}
```

### Pattern 2: Async Widget with Loading State
**What:** Use `ui.wait()` inside a `DG.Widget` to show a spinner while fetching external data.
**When to use:** Any panel that makes HTTP requests.
**Example:**
```typescript
// Source: Verified from Chem identifiers widget (packages/Chem/src/widgets/identifiers.ts line 138-141)
export async function identifiersWidget(molfile: string): Promise<DG.Widget> {
  const map = await getIdentifiersSingle(molfile);
  return map ? new DG.Widget(ui.tableFromMap(map)) : new DG.Widget(ui.divText('Malformed molecule'));
}

// Alternative pattern using ui.wait for synchronous return:
export function myPanel(value: string): DG.Widget {
  return new DG.Widget(ui.wait(async () => {
    const data = await fetchData(value);
    return ui.tableFromMap(data);
  }));
}
```

### Pattern 3: R Script with grok.functions.call
**What:** R scripts in `scripts/` directory with metadata comments are registered as platform functions. TypeScript calls them via `grok.functions.call('PackageName:scriptName', params)`.
**When to use:** Server-side statistical computation (DEqMS).
**Example:**
```typescript
// Source: Existing pattern in differential-expression.ts (lines 133-140)
const result: DG.DataFrame = await grok.functions.call('Proteomics:deqmsDE', {
  exprDf: exprDf,
  nGroup1: group1Cols.length,
  peptideCounts: peptideCountCol,  // NEW: integer column with peptide counts
  fcThreshold: fcThreshold,
  pThreshold: pThreshold,
});
```

### Pattern 4: UniProt REST API Fetch
**What:** Direct fetch to `rest.uniprot.org/uniprotkb/{accession}.json` with field filtering.
**When to use:** Fetching protein annotation for a single accession.
**Example:**
```typescript
// Source: Verified via direct API testing (rest.uniprot.org/uniprotkb/P04637.json)
async function fetchUniProtData(accession: string): Promise<UniProtEntry | null> {
  const url = `https://rest.uniprot.org/uniprotkb/${encodeURIComponent(accession)}.json` +
    '?fields=accession,protein_name,gene_names,organism_name,cc_function,go';
  try {
    const resp = await fetch(url);
    if (!resp.ok) return null;
    return await resp.json();
  } catch {
    return null;
  }
}
```

**UniProt JSON structure (verified):**
```
proteinDescription.recommendedName.fullName.value  -> "Cellular tumor antigen p53"
genes[0].geneName.value                            -> "TP53"
organism.scientificName                            -> "Homo sapiens"
comments[].texts[].value  (where commentType = "FUNCTION") -> function description
uniProtKBCrossReferences[] (where database = "GO"):
  .id                                              -> "GO:0005634"
  .properties[0].value                             -> "C:nucleus" / "F:DNA binding" / "P:apoptosis"
```

### Anti-Patterns to Avoid
- **Bundling UniProt data locally:** Always fetch live from API -- data changes with each UniProt release.
- **Blocking the UI thread:** Use `ui.wait()` or return `Promise<DG.Widget>` -- never synchronously fetch.
- **Passing complex objects to R scripts:** Use simple column types (dataframe, int, double). The existing `buildExpressionDf` pattern (renaming columns to s1, s2, ...) avoids encoding issues.
- **Creating a new panel registration system:** Use the exact same decorator pattern from Chem/Bio packages.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Protein annotation lookup | Custom database or scraper | UniProt REST API | Authoritative, maintained, versioned data |
| Peptide-count-weighted statistics | Custom variance estimation | DEqMS `spectraCounteBayes()` | Validated statistical method with publications |
| Info panel registration | Custom event listeners on grid clicks | `@grok.decorators.panel` + semantic type binding | Platform handles panel lifecycle, caching, display |
| Async widget loading spinner | Custom loading indicator | `ui.wait()` | Built-in platform pattern with consistent UX |
| FDR correction in DEqMS | Manual BH adjustment | DEqMS `outputResult()` with `sca.adj.pval` | Already computed by the package |

**Key insight:** The Datagrok platform already provides the entire info panel infrastructure. The UniProt panel is purely: fetch JSON, extract fields, render in `ui.tableFromMap`. No custom UI framework needed.

## Common Pitfalls

### Pitfall 1: Protein ID Format Mismatch
**What goes wrong:** MaxQuant protein IDs can be semicolon-delimited groups (e.g., "P04637;Q9H3D4"), but UniProt API expects a single accession.
**Why it happens:** The "Protein IDs" column contains protein groups, not individual accessions.
**How to avoid:** Use the "Primary Protein ID" column (already created by the parser, extracts first accession) for the panel. The semantic type `Proteomics-ProteinId` is assigned to both the original and primary columns, so the panel will trigger on either. When the value contains semicolons, extract the first accession before querying.
**Warning signs:** 404 responses from UniProt API.

### Pitfall 2: package.g.ts Not Updated for Panels
**What goes wrong:** The panel function exists in package.ts but is not registered with the platform because package.g.ts is missing the corresponding metadata comment entry.
**Why it happens:** `grok api` may not fully regenerate decorator-based registrations (noted in STATE.md decision [03-02]).
**How to avoid:** Manually add the panel entry to package.g.ts with the correct `//meta.role: panel` and `//input: string proteinId { semType: Proteomics-ProteinId }` metadata.
**Warning signs:** Panel does not appear in context panel when clicking a protein ID cell.

### Pitfall 3: DEqMS Package Not Installed on Server
**What goes wrong:** The R script fails because DEqMS is not available in the server's R environment.
**Why it happens:** DEqMS is a Bioconductor package, not CRAN. Requires `BiocManager::install("DEqMS")`.
**How to avoid:** Follow the same fallback pattern as limma_de.R: `hasDeqms <- suppressWarnings(require(DEqMS, quietly = TRUE))`. If DEqMS is unavailable, fall back to limma (which itself falls back to t-test).
**Warning signs:** R script error about missing package.

### Pitfall 4: DEqMS Output Column Name Mismatch
**What goes wrong:** The TypeScript code expects columns named `log2FC`, `p.value`, `adj.p.value` but DEqMS `outputResult()` produces `logFC`, `sca.P.Value`, `sca.adj.pval`.
**Why it happens:** DEqMS has its own column naming convention different from limma's `topTable()`.
**How to avoid:** In the R script, rename DEqMS output columns to match the existing contract: `log2FC`, `p.value`, `adj.p.value`, `significant`.
**Warning signs:** Null columns when copying results back to the DataFrame.

### Pitfall 5: UniProt API Unavailability
**What goes wrong:** Network errors or UniProt downtime causes the panel to show an error.
**Why it happens:** External API dependency.
**How to avoid:** Wrap fetch in try/catch, return a graceful "Unable to fetch UniProt data" message in the widget. Use `DG.Widget` with `ui.divText()` for error states.
**Warning signs:** Console errors on every cell click.

### Pitfall 6: Peptide Count Column Not Available
**What goes wrong:** DEqMS requires peptide counts but the DataFrame may not have a "Peptides" or "Unique peptides" column (e.g., if data was imported differently).
**Why it happens:** Not all proteomics data sources include peptide count columns.
**How to avoid:** In the DE dialog, add a column picker for peptide counts (defaulting to "Peptides" or "Unique peptides" if found). If no peptide count column exists, disable the DEqMS option and only offer limma.
**Warning signs:** R script fails with missing column error.

## Code Examples

### UniProt Panel Widget
```typescript
// Source: Pattern verified from Chem package identifiers widget + UniProt API testing

interface UniProtEntry {
  proteinDescription?: {
    recommendedName?: { fullName?: { value?: string } };
  };
  genes?: Array<{ geneName?: { value?: string } }>;
  organism?: { scientificName?: string; commonName?: string };
  comments?: Array<{
    commentType?: string;
    texts?: Array<{ value?: string }>;
  }>;
  uniProtKBCrossReferences?: Array<{
    database?: string;
    id?: string;
    properties?: Array<{ key?: string; value?: string }>;
  }>;
}

function extractGoTerms(entry: UniProtEntry): {mf: string[]; bp: string[]; cc: string[]} {
  const mf: string[] = [];
  const bp: string[] = [];
  const cc: string[] = [];
  for (const xref of entry.uniProtKBCrossReferences ?? []) {
    if (xref.database !== 'GO') continue;
    const term = xref.properties?.find(p => p.key === 'GoTerm')?.value ?? '';
    if (term.startsWith('F:')) mf.push(term.substring(2));
    else if (term.startsWith('P:')) bp.push(term.substring(2));
    else if (term.startsWith('C:')) cc.push(term.substring(2));
  }
  return {mf, bp, cc};
}

function renderUniProtWidget(entry: UniProtEntry): HTMLElement {
  const map: Record<string, any> = {};
  map['Protein'] = entry.proteinDescription?.recommendedName?.fullName?.value ?? 'Unknown';
  map['Gene'] = entry.genes?.[0]?.geneName?.value ?? 'Unknown';
  map['Organism'] = entry.organism?.scientificName ?? 'Unknown';

  const funcComment = entry.comments?.find(c => c.commentType === 'FUNCTION');
  if (funcComment?.texts?.[0]?.value)
    map['Function'] = funcComment.texts[0].value;

  const go = extractGoTerms(entry);
  if (go.mf.length) map['Molecular Function'] = go.mf.slice(0, 5).join(', ');
  if (go.bp.length) map['Biological Process'] = go.bp.slice(0, 5).join(', ');
  if (go.cc.length) map['Cellular Component'] = go.cc.slice(0, 5).join(', ');

  return ui.tableFromMap(map);
}
```

### DEqMS R Script Structure
```r
# Source: DEqMS Bioconductor vignette + existing limma_de.R pattern
#name: deqmsDE
#description: DEqMS peptide-count-weighted differential expression
#language: r
#input: dataframe exprDf
#input: int nGroup1
#input: column_list peptideCounts {type: numerical}
#input: double fcThreshold = 1.0
#input: double pThreshold = 0.05
#output: dataframe result

library(limma)
library(DEqMS)

exprMat <- as.matrix(exprDf[, -ncol(exprDf)])  # exclude peptide count col
counts <- exprDf[, ncol(exprDf)]                # last column is peptide counts

group <- factor(c(rep("ctrl", nGroup1), rep("treat", ncol(exprMat) - nGroup1)),
                levels = c("ctrl", "treat"))
design <- model.matrix(~ group)

fit <- lmFit(exprMat, design)
fit <- eBayes(fit)

# DEqMS: assign peptide counts and run spectra-count-adjusted eBayes
fit$count <- counts
fit <- spectraCounteBayes(fit)

# Extract DEqMS results
tt <- outputResult(fit, coef = 2)

result <- data.frame(
  log2FC = tt$logFC,
  p.value = tt$sca.P.Value,
  adj.p.value = tt$sca.adj.pval,
  significant = (abs(tt$logFC) >= fcThreshold) & (tt$sca.adj.pval <= pThreshold),
  check.names = FALSE
)
```

### Modifying DE Dialog for Method Selection
```typescript
// Source: Pattern from existing showDEDialog in differential-expression.ts

// Add method selector to dialog:
const methodInput = ui.input.choice('Method', {
  value: 'limma',
  items: ['limma', 'DEqMS'],
  nullable: false,
});

// Add peptide count column picker (shown only when DEqMS selected):
const peptideColInput = ui.input.column('Peptide count column', {
  table: df,
  value: df.col('Peptides') ?? df.col('Unique peptides') ?? null,
  filter: (col: DG.Column) => col.type === DG.COLUMN_TYPE.INT || col.type === DG.COLUMN_TYPE.FLOAT,
});
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Client-side t-test | limma moderated t-test (server-side R) | Phase 2->3 of this project | Better statistical power via empirical Bayes |
| limma (uniform variance) | DEqMS (peptide-count-weighted variance) | DEqMS published 2020 | More accurate for proteins with varying peptide counts |
| UniProt JAPI / XML API | UniProt REST API v2 (rest.uniprot.org) | 2022 | JSON format, field filtering, modern REST interface |

**Deprecated/outdated:**
- UniProt XML API (www.uniprot.org/uniprot/P04637.xml): Still works but REST JSON is preferred
- DEP R package: Similar goals to DEqMS but less widely adopted

## Open Questions

1. **DEqMS R environment availability**
   - What we know: limma_de.R uses `suppressWarnings(require(limma))` for graceful fallback. DEqMS depends on limma.
   - What's unclear: Whether the Datagrok server's R environment has BiocManager and can install DEqMS. Same concern flagged in STATE.md for limma.
   - Recommendation: Use the same `require(DEqMS, quietly = TRUE)` fallback pattern. If DEqMS unavailable, fall back to limma. If limma unavailable, fall back to t-test.

2. **UniProt API rate limiting**
   - What we know: The API works for single-entry lookups. No specific rate limit documented for cell-click panels (one request at a time).
   - What's unclear: Exact requests/second limit. The EBI Proteins API allows 200/sec.
   - Recommendation: Not a concern for cell-click panels (one request per user click). Add simple error handling for 429 responses.

3. **Peptide count column naming across MaxQuant versions**
   - What we know: Demo data has "Peptides" and "Unique peptides" columns. DEqMS recommends using "Unique peptides" or PSM counts.
   - What's unclear: Whether all MaxQuant versions use the same column names.
   - Recommendation: Auto-detect both "Peptides" and "Unique peptides", prefer "Unique peptides" as default (more specific), allow user to pick any numerical column.

## Sources

### Primary (HIGH confidence)
- Chem package source code (packages/Chem/src/package.ts) - panel decorator pattern, identifiers widget
- Bio package source code (packages/Bio/src/package.ts) - panel decorator pattern with semantic types
- UniProt REST API direct testing (rest.uniprot.org/uniprotkb/P04637.json) - JSON structure verified
- Existing limma_de.R and differential-expression.ts - R script calling pattern

### Secondary (MEDIUM confidence)
- [DEqMS Bioconductor page](https://bioconductor.org/packages/release/bioc/html/DEqMS.html) - package info
- [DEqMS vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html) - workflow
- [spectraCounteBayes docs](https://rdrr.io/bioc/DEqMS/man/spectraCounteBayes.html) - function signature
- [outputResult docs](https://rdrr.io/bioc/DEqMS/man/outputResult.html) - output column names (logFC, sca.P.Value, sca.adj.pval)

### Tertiary (LOW confidence)
- UniProt API rate limits - not found in official docs; estimated based on EBI Proteins API (200/sec)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - UniProt REST API verified via direct testing; DEqMS well-documented on Bioconductor
- Architecture: HIGH - Panel pattern verified from two existing Datagrok packages (Chem, Bio) with working code
- Pitfalls: HIGH - Derived from actual codebase analysis and known issues (STATE.md blockers)

**Research date:** 2026-03-01
**Valid until:** 2026-03-31 (stable APIs and libraries, no fast-moving dependencies)

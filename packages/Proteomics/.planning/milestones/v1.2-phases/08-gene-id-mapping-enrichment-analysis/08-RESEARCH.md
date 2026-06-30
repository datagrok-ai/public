# Phase 8: Gene ID Mapping & Enrichment Analysis - Research

**Researched:** 2026-03-06
**Domain:** Proteomics enrichment analysis via g:Profiler REST API
**Confidence:** HIGH

## Summary

This phase adds gene ID mapping (UniProt to gene symbol) and functional enrichment analysis (GO, KEGG, Reactome) to the Proteomics package. The user has locked the decision to use g:Profiler's REST API for both operations -- g:Convert for ID mapping and g:GOSt for enrichment. This is a pure TypeScript implementation with no R scripts needed.

The g:Profiler API is well-documented, free, and handles the statistical heavy lifting (hypergeometric test, multiple testing correction). The implementation follows existing patterns: a menu entry in `PackageFunctions`, a config dialog, an async API call with progress indicator, and results displayed in a new Datagrok table view (matching the PCA pattern for different-shaped data).

**Primary recommendation:** Implement as two modules: `analysis/enrichment.ts` (API client + orchestration + dialog) and the menu wiring in `package.ts`. The g:Profiler API accepts mixed identifier types and returns structured JSON, making the integration straightforward.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Menu item: `Proteomics | Enrichment Analysis...` opens a config dialog
- Dialog fields: significance thresholds (FC, p-value), enrichment sources (GO BP/MF/CC, KEGG, Reactome), organism dropdown
- Auto-fill thresholds from existing DE columns (log2FC, adj.p-value) with defaults FC>=1, padj<=0.05
- All three enrichment sources (GO, KEGG, Reactome) enabled by default -- user can uncheck any
- Live count of significant proteins updates as thresholds change: "42 of 1,200 proteins are significant with these thresholds"
- Progress dialog blocks while API call runs: "Running enrichment analysis..."
- Results open as a new Datagrok table view (follows PCA pattern for different-shaped data)
- Full schema: Source (GO:BP/GO:MF/GO:CC/KEGG/REAC), Term ID, Term Name, P-value, FDR, Gene Count, Gene Ratio, Intersection (member genes)
- Show all results, no pre-filtering -- user can sort/filter manually
- Significant terms (FDR < 0.05) highlighted but not hidden
- Gene mapping is automatic and transparent -- triggered as part of enrichment, not a separate menu item
- If a gene name column with SEMTYPE.GENE_SYMBOL already exists, use it directly -- skip API mapping
- If no gene names exist, map UniProt accessions to gene symbols via g:Profiler g:Convert
- Add a 'Gene Symbol (mapped)' column to the main table after mapping (only when no gene name column exists)
- Show mapping statistics in progress dialog: "850/1,200 proteins mapped (70.8%). 350 unmapped."
- Organism dropdown in the enrichment dialog, defaulting to Homo sapiens
- Short curated list of 8-10 common proteomics organisms: human, mouse, rat, yeast, E. coli, zebrafish, Drosophila, Arabidopsis, C. elegans
- g:Profiler REST API handles both ID mapping (g:Convert) and enrichment (g:GOSt) -- pure TypeScript, no R scripts needed
- Use `parseAccession()` from uniprot-panel.ts to extract clean accessions before sending to g:Profiler
- Background set for enrichment must be all quantified proteins (not whole genome) -- use g:Profiler's `domain_scope: "custom"` parameter

### Claude's Discretion
- Exact highlighting mechanism for significant terms (color coding, row tags, etc.)
- Error handling for API failures (timeout, rate limit, no results)
- Exact wording of mapping statistics message
- How to handle proteins that map to multiple gene symbols

### Deferred Ideas (OUT OF SCOPE)
- Enrichment visualization (dot plots, bar charts) -- Phase 9 (VIZ-02, VIZ-03)
- Volcano plot integration (select GO term -> highlight proteins) -- Phase 9 (ENRICH-04)
- GSEA (rank-based enrichment) -- tracked as ENRICH-05 in future requirements
- Functional annotation clustering -- tracked as ENRICH-06 in future requirements
- Standalone ID mapping menu item -- add if users request it
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| MAP-01 | User can map UniProt accessions to gene symbols and Entrez IDs for downstream enrichment | g:Convert API endpoint (`POST /api/convert/convert/`) returns `name` field (gene symbol) for each input; `parseAccession()` extracts clean UniProt IDs; batch operation in single request |
| MAP-02 | User can select the organism (default: human) for ID mapping | g:Convert and g:GOSt both accept `organism` param using binomial code (e.g., "hsapiens", "mmusculus"); curated organism list of 9 common species |
| ENRICH-01 | User can run GO overrepresentation analysis (BP, MF, CC) on DE protein list with proteomics-specific background set | g:GOSt `sources: ["GO:BP","GO:MF","GO:CC"]` with `domain_scope: "custom"` and `background: [all quantified gene symbols]` |
| ENRICH-02 | User can run KEGG pathway enrichment on DE protein list | g:GOSt source `"KEGG"` in same API call; FDR via `significance_threshold_method: "fdr"` |
| ENRICH-03 | User can run Reactome pathway enrichment on DE protein list | g:GOSt source `"REAC"` in same API call |
| VIZ-01 | User can view enrichment results in an interactive, sortable, filterable table with term ID, description, p-value, FDR, gene count, and member genes | Results returned as DG.DataFrame with standard schema; opened via `grok.shell.addTableView()` -- Datagrok grid provides native sort/filter/highlight |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| g:Profiler REST API | current (biit.cs.ut.ee) | ID mapping + enrichment analysis | Free, CORS-enabled, handles GO/KEGG/Reactome in single API; supports custom backgrounds; standard in proteomics/genomics |
| datagrok-api | workspace version | UI, DataFrame, viewers | Platform API -- already in use |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Native `fetch()` | browser built-in | HTTP POST to g:Profiler | Already used in uniprot-panel.ts pattern |
| `parseAccession()` | existing in package | Extract clean UniProt accessions | Before sending IDs to g:Convert |
| `findProteomicsColumns()` | existing in package | Detect protein ID, gene symbol, log2FC, p-value columns | Dialog auto-fill and pipeline orchestration |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| g:Profiler API | clusterProfiler (R) | Requires R scripting infrastructure; more complex; offers GSEA (deferred to future) |
| g:Profiler API | Enrichr REST API | Enrichr does not support custom background sets (critical for proteomics); limited organisms |
| g:Profiler API | UniProt ID Mapping + separate enrichment API | Two APIs to maintain; UniProt API requires polling for batch results |

**No additional npm packages needed.** This phase uses only browser-native `fetch()` and existing Datagrok APIs.

## Architecture Patterns

### Recommended Project Structure
```
src/
├── analysis/
│   └── enrichment.ts        # g:Profiler API client, ID mapping, enrichment orchestration, dialog, result builder
├── utils/
│   ├── proteomics-types.ts  # Existing SEMTYPE constants (no changes needed)
│   └── column-detection.ts  # Existing findColumn/findProteomicsColumns (no changes needed)
├── panels/
│   └── uniprot-panel.ts     # Existing parseAccession() -- reused, not modified
└── package.ts               # New menu entry: Proteomics | Enrichment Analysis...
```

### Pattern 1: g:Convert API Call (ID Mapping)
**What:** Map UniProt accessions to gene symbols via POST request
**When to use:** When no SEMTYPE.GENE_SYMBOL column exists in the DataFrame
**Example:**
```typescript
// Source: gprofiler2 R package source + REST API docs
const GPROFILER_BASE = 'https://biit.cs.ut.ee/gprofiler';

interface ConvertResult {
  incoming: string;      // input ID (e.g. 'P04637')
  converted: string;     // target namespace ID (e.g. 'ENSG00000141510')
  name: string;          // gene symbol (e.g. 'TP53') -- always returned regardless of target
  description: string;   // gene description
  namespaces: string;    // detected input namespace
  n_incoming: number;    // position in input
  n_converted: number;   // conversion order
}

async function gConvert(
  accessions: string[],
  organism: string = 'hsapiens',
): Promise<ConvertResult[]> {
  const response = await fetch(`${GPROFILER_BASE}/api/convert/convert/`, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({
      organism: organism,
      query: accessions,
      target: 'ENSG',
      numeric_ns: '',
      output: 'json',
    }),
  });
  if (!response.ok) throw new Error(`g:Convert failed: ${response.status}`);
  const data = await response.json();
  return data.result ?? [];
}
```

### Pattern 2: g:GOSt API Call (Enrichment)
**What:** Run overrepresentation analysis on significant gene list with custom background
**When to use:** After ID mapping, for all selected enrichment sources
**Example:**
```typescript
// Source: gprofiler2 R package source + REST API docs
interface GostResult {
  native: string;              // term ID (e.g. 'GO:0006915', 'KEGG:04210')
  name: string;                // term name
  source: string;              // 'GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC'
  p_value: number;             // adjusted p-value (already FDR-corrected)
  significant: boolean;
  term_size: number;           // total genes annotated to this term
  query_size: number;          // input gene list size (after filtering by annotations)
  intersection_size: number;   // overlapping genes count
  effective_domain_size: number;
  precision: number;           // intersection_size / query_size
  recall: number;              // intersection_size / term_size
  intersections: string[][];   // positional per-query-gene array (when no_evidences=false)
}

async function gGOSt(
  queryGenes: string[],
  backgroundGenes: string[],
  organism: string,
  sources: string[],
  threshold: number = 0.05,
): Promise<GostResult[]> {
  const response = await fetch(`${GPROFILER_BASE}/api/gost/profile/`, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({
      organism: organism,
      query: queryGenes,
      sources: sources,
      user_threshold: threshold,
      significance_threshold_method: 'fdr',
      domain_scope: 'custom',
      background: backgroundGenes,
      all_results: true,         // return all -- user filters in Datagrok grid
      ordered: false,
      no_evidences: false,       // include intersection gene info
      combined: false,
      measure_underrepresentation: false,
      no_iea: false,
      numeric_ns: '',
      output: 'json',
    }),
  });
  if (!response.ok) throw new Error(`g:GOSt failed: ${response.status}`);
  const data = await response.json();
  return data.result ?? [];
}
```

### Pattern 3: Enrichment Dialog (follows DE Dialog pattern)
**What:** Config dialog with live preview count, similar to `showDEDialog()`
**When to use:** Entry point from menu
**Example:**
```typescript
// Follows existing dialog pattern from differential-expression.ts
const ORGANISM_LIST = [
  {display: 'Homo sapiens (Human)', code: 'hsapiens'},
  {display: 'Mus musculus (Mouse)', code: 'mmusculus'},
  {display: 'Rattus norvegicus (Rat)', code: 'rnorvegicus'},
  {display: 'Saccharomyces cerevisiae (Yeast)', code: 'scerevisiae'},
  {display: 'Escherichia coli (K12)', code: 'ecoli'},
  {display: 'Danio rerio (Zebrafish)', code: 'drerio'},
  {display: 'Drosophila melanogaster (Fruit fly)', code: 'dmelanogaster'},
  {display: 'Arabidopsis thaliana', code: 'athaliana'},
  {display: 'Caenorhabditis elegans', code: 'celegans'},
] as const;

function showEnrichmentDialog(df: DG.DataFrame): void {
  if (df.getTag('proteomics.de_complete') !== 'true') {
    grok.shell.warning('Run Differential Expression first');
    return;
  }

  const fcInput = ui.input.float('|log2FC| threshold', {value: 1.0});
  const pInput = ui.input.float('Adj. p-value threshold', {value: 0.05});

  const organismInput = ui.input.choice('Organism', {
    value: 'Homo sapiens (Human)',
    items: ORGANISM_LIST.map(o => o.display),
    nullable: false,
  });

  // Enrichment source checkboxes (all enabled by default)
  const goBpInput = ui.input.bool('GO: Biological Process', {value: true});
  const goMfInput = ui.input.bool('GO: Molecular Function', {value: true});
  const goCcInput = ui.input.bool('GO: Cellular Component', {value: true});
  const keggInput = ui.input.bool('KEGG Pathways', {value: true});
  const reactomeInput = ui.input.bool('Reactome Pathways', {value: true});

  // Live count of significant proteins
  const countDiv = ui.divText('');
  const updateCount = () => { /* count proteins matching thresholds, update countDiv */ };
  updateCount();
  fcInput.onChanged.subscribe(updateCount);
  pInput.onChanged.subscribe(updateCount);

  ui.dialog('Enrichment Analysis')
    .add(countDiv)
    .add(fcInput).add(pInput).add(organismInput)
    .add(goBpInput).add(goMfInput).add(goCcInput)
    .add(keggInput).add(reactomeInput)
    .onOK(async () => {
      const pi = DG.TaskBarProgressIndicator.create('Running enrichment analysis...');
      try {
        // 1. Map IDs (if needed)
        // 2. Run g:GOSt enrichment
        // 3. Build result DataFrame
        // 4. grok.shell.addTableView(enrichDf)
      } finally { pi.close(); }
    })
    .show();
}
```

### Pattern 4: Result DataFrame Construction
**What:** Build enrichment results as a new DataFrame with the locked schema
**When to use:** After g:GOSt returns results
**Example:**
```typescript
function buildEnrichmentDf(
  results: GostResult[],
  queryGenes: string[],
): DG.DataFrame {
  const n = results.length;
  const df = DG.DataFrame.create(n);
  df.name = 'Enrichment Results';

  const sourceCol = df.columns.addNewString('Source');
  const termIdCol = df.columns.addNewString('Term ID');
  const termNameCol = df.columns.addNewString('Term Name');
  const pValCol = df.columns.addNewFloat('P-value');
  const fdrCol = df.columns.addNewFloat('FDR');
  const geneCountCol = df.columns.addNewInt('Gene Count');
  const geneRatioCol = df.columns.addNewFloat('Gene Ratio');
  const intersectionCol = df.columns.addNewString('Intersection');

  for (let i = 0; i < n; i++) {
    const r = results[i];
    sourceCol.set(i, r.source);
    termIdCol.set(i, r.native);
    termNameCol.set(i, r.name);
    pValCol.set(i, r.p_value);
    fdrCol.set(i, r.p_value);  // g:GOSt with fdr method returns adjusted p-values
    geneCountCol.set(i, r.intersection_size);
    geneRatioCol.set(i, r.precision);  // precision = intersection_size / query_size

    // Derive member genes from positional intersections array
    const memberGenes: string[] = [];
    if (r.intersections) {
      for (let j = 0; j < queryGenes.length && j < r.intersections.length; j++) {
        if (r.intersections[j] && r.intersections[j].length > 0)
          memberGenes.push(queryGenes[j]);
      }
    }
    intersectionCol.set(i, memberGenes.join(', '));
  }

  df.setTag('proteomics.enrichment', 'true');
  return df;
}
```

### Pattern 5: Menu Integration
**What:** Single menu entry following existing package.ts patterns
**When to use:** Wiring enrichment to the Proteomics menu
**Example:**
```typescript
// In PackageFunctions class in package.ts
@grok.decorators.func({'top-menu': 'Proteomics | Enrichment Analysis...'})
static async enrichmentAnalysis(): Promise<void> {
  const df = grok.shell.tv?.dataFrame;
  if (!df) { grok.shell.warning('No table open'); return; }
  showEnrichmentDialog(df);
}
```

### Anti-Patterns to Avoid
- **Sending raw protein IDs to g:Profiler:** Always use `parseAccession()` first to strip sp|/tr| prefixes and semicolons.
- **Using whole genome as background:** Always use `domain_scope: "custom"` with all quantified proteins as background.
- **Sending accessions directly to g:GOSt:** Map to gene symbols first via g:Convert, then send gene symbols to g:GOSt for better intersection reporting.
- **Blocking the UI during API calls:** Use `DG.TaskBarProgressIndicator` (established pattern).
- **Storing enrichment results in the protein DataFrame:** Different row semantics -- always create a separate DataFrame.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| GO/KEGG/Reactome enrichment statistics | Custom hypergeometric test + FDR | g:Profiler g:GOSt API | Handles annotation databases, multiple testing correction, term hierarchy |
| UniProt-to-gene-symbol mapping | Custom UniProt REST API with polling | g:Profiler g:Convert API | Synchronous single-POST; handles >100 ID namespaces |
| Multiple testing correction for enrichment | Custom BH/FDR correction | g:Profiler `significance_threshold_method: "fdr"` | Already built into API response |
| Protein accession parsing | New parsing logic | Existing `parseAccession()` in uniprot-panel.ts | Already handles sp|, tr|, semicolons, CPTAC formats |
| Column detection (protein ID, gene name) | Manual column scanning | Existing `findProteomicsColumns()` | Already finds by semType or name hints |
| Gene annotation databases | Local GO/KEGG/Reactome files | g:Profiler (updated quarterly via Ensembl) | Database maintenance is enormous |

**Key insight:** g:Profiler handles the entire enrichment pipeline (ID resolution, annotation lookup, statistical test, multiple testing correction, background handling) in a single API call. Building any component individually would require maintaining gene annotation databases.

## Common Pitfalls

### Pitfall 1: Wrong Background Set
**What goes wrong:** Using default genome-wide background inflates significance. In proteomics, only ~2,000-6,000 proteins are typically quantified, not 20,000+ coding genes.
**Why it happens:** Most enrichment tutorials use transcriptomics examples where whole-genome background is more appropriate.
**How to avoid:** ALWAYS pass `domain_scope: "custom"` and `background: [all quantified gene symbols]` to g:GOSt. Background = ALL proteins in the dataset, not just significant ones.
**Warning signs:** If nearly every GO term is significant, the background is likely wrong.

### Pitfall 2: Accession Format Mismatch
**What goes wrong:** g:Convert returns no matches because accessions still have sp|/tr| prefixes.
**Why it happens:** Forgetting to run `parseAccession()` on raw protein IDs.
**How to avoid:** Always clean accessions through `parseAccession()` before any API call. Filter out nulls.
**Warning signs:** Mapping success rate near 0%.

### Pitfall 3: Missing Intersection Data
**What goes wrong:** The intersection/member genes column is empty in results.
**Why it happens:** Not setting `no_evidences: false` in the g:GOSt request. The API field name is confusing -- it's a double negative.
**How to avoid:** Explicitly set `no_evidences: false` to get gene intersection lists. Note: The R client uses `evcodes = TRUE` but the REST API uses the inverse `no_evidences: false`.
**Warning signs:** `intersections` field missing or empty in API response.

### Pitfall 4: Unmapped Proteins Silently Excluded
**What goes wrong:** Proteins that fail ID mapping are silently dropped, biasing the background set and confusing the user.
**Why it happens:** Only successfully mapped genes are sent to g:GOSt as query/background.
**How to avoid:** Report mapping statistics clearly in progress dialog. Background should include all mappable proteins, query includes only significant AND mappable proteins.
**Warning signs:** Background set much smaller than total protein count.

### Pitfall 5: One-to-Many Mapping
**What goes wrong:** A UniProt accession maps to multiple gene symbols, inflating gene counts.
**Why it happens:** g:Convert can return multiple targets for one input. Protein groups (P12345;Q67890) map to different genes.
**How to avoid:** Use first match per input accession (g:Convert returns entries ordered by relevance). For protein groups, `parseAccession()` already takes the first entry. Deduplicate gene symbols before sending to g:GOSt.
**Warning signs:** Gene count exceeds protein count.

### Pitfall 6: Organism Mismatch
**What goes wrong:** User data is from mouse but enrichment runs against human annotations. Results are meaningless.
**Why it happens:** Organism defaults to human; user forgets to change.
**How to avoid:** Prominent organism selector in dialog. Consider checking if DataFrame has an organism tag from import.
**Warning signs:** Very low mapping rate can indicate organism mismatch.

## Code Examples

### Organism List (Curated for Proteomics)
```typescript
const ORGANISM_LIST = [
  {display: 'Homo sapiens (Human)', code: 'hsapiens'},
  {display: 'Mus musculus (Mouse)', code: 'mmusculus'},
  {display: 'Rattus norvegicus (Rat)', code: 'rnorvegicus'},
  {display: 'Saccharomyces cerevisiae (Yeast)', code: 'scerevisiae'},
  {display: 'Escherichia coli (K12)', code: 'ecoli'},
  {display: 'Danio rerio (Zebrafish)', code: 'drerio'},
  {display: 'Drosophila melanogaster (Fruit fly)', code: 'dmelanogaster'},
  {display: 'Arabidopsis thaliana', code: 'athaliana'},
  {display: 'Caenorhabditis elegans', code: 'celegans'},
] as const;
```

### Complete Enrichment Pipeline Flow
```typescript
async function runEnrichmentPipeline(
  df: DG.DataFrame,
  fcThreshold: number,
  pThreshold: number,
  organismCode: string,
  sources: string[],
): Promise<DG.DataFrame> {
  const cols = findProteomicsColumns(df);
  if (!cols.proteinId) throw new Error('No protein ID column found');

  // Step 1: Get gene symbols (existing or mapped)
  let geneForRow: Map<number, string>;  // row index -> gene symbol
  const existingGeneCol = cols.geneName;

  if (existingGeneCol) {
    geneForRow = new Map();
    for (let i = 0; i < df.rowCount; i++) {
      const val = existingGeneCol.get(i);
      if (val && !existingGeneCol.isNone(i)) geneForRow.set(i, val as string);
    }
  } else {
    // Extract clean accessions and map via g:Convert
    const accToRows = new Map<string, number[]>();
    for (let i = 0; i < df.rowCount; i++) {
      const raw = cols.proteinId!.get(i) as string;
      const acc = parseAccession(raw);
      if (acc) {
        if (!accToRows.has(acc)) accToRows.set(acc, []);
        accToRows.get(acc)!.push(i);
      }
    }

    const uniqueAccessions = [...accToRows.keys()];
    const convertResults = await gConvert(uniqueAccessions, organismCode);

    // Build accession -> gene symbol map (first match per accession)
    const accToGene = new Map<string, string>();
    for (const r of convertResults) {
      if (r.incoming && r.name && r.name !== 'N/A' && !accToGene.has(r.incoming))
        accToGene.set(r.incoming, r.name);
    }

    // Map rows to gene symbols
    geneForRow = new Map();
    for (const [acc, rows] of accToRows) {
      const gene = accToGene.get(acc);
      if (gene) { for (const row of rows) geneForRow.set(row, gene); }
    }

    // Add gene symbol column to main table
    const geneCol = df.columns.addNewString('Gene Symbol (mapped)');
    geneCol.semType = SEMTYPE.GENE_SYMBOL;
    for (const [row, gene] of geneForRow) geneCol.set(row, gene);

    // Report mapping stats
    const mapped = accToGene.size;
    const total = uniqueAccessions.length;
    const unmapped = total - mapped;
    // Show in progress: `${mapped}/${total} proteins mapped (${(mapped/total*100).toFixed(1)}%). ${unmapped} unmapped.`
  }

  // Step 2: Build significant and background gene sets
  const significantGenes: Set<string> = new Set();
  const backgroundGenes: Set<string> = new Set();

  for (const [row, gene] of geneForRow) {
    backgroundGenes.add(gene);
    const adjP = cols.pValue && !cols.pValue.isNone(row) ? cols.pValue.get(row) as number : null;
    const fc = cols.log2fc && !cols.log2fc.isNone(row) ? Math.abs(cols.log2fc.get(row) as number) : null;
    if (adjP !== null && fc !== null && adjP <= pThreshold && fc >= fcThreshold)
      significantGenes.add(gene);
  }

  const queryArray = [...significantGenes];
  const bgArray = [...backgroundGenes];

  // Step 3: Call g:GOSt enrichment
  const results = await gGOSt(queryArray, bgArray, organismCode, sources);

  // Step 4: Build result DataFrame
  return buildEnrichmentDf(results, queryArray);
}
```

### Highlighting Significant Terms (Claude's Discretion)
```typescript
// Recommended: Add a boolean 'Significant' column and use color coding on FDR
function highlightSignificantTerms(df: DG.DataFrame): void {
  const fdrCol = df.col('FDR');
  if (!fdrCol) return;

  // Boolean column for easy filtering
  const sigCol = df.columns.addNewBool('Significant');
  for (let i = 0; i < df.rowCount; i++)
    sigCol.set(i, !fdrCol.isNone(i) && fdrCol.get(i) < 0.05);

  // Color coding on FDR column (green = significant, red = not)
  fdrCol.setTag('color-coding-type', 'Linear');
  fdrCol.setTag('color-coding-linear', '{"0":"#2ecc71","0.05":"#f39c12","1":"#e74c3c"}');
}
```

### Error Handling (Claude's Discretion)
```typescript
// Recommended: wrap API calls with timeout and informative error messages
async function fetchWithTimeout(url: string, options: RequestInit, timeoutMs: number = 30000): Promise<Response> {
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs);
  try {
    const resp = await fetch(url, {...options, signal: controller.signal});
    return resp;
  } catch (e: any) {
    if (e.name === 'AbortError')
      throw new Error('g:Profiler API request timed out. The service may be temporarily unavailable.');
    throw new Error(`g:Profiler API error: ${e.message}`);
  } finally {
    clearTimeout(timer);
  }
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| DAVID (web form, outdated DBs) | g:Profiler REST API (free, current DBs) | ~2018 | DAVID databases last updated ~2016; no longer recommended |
| clusterProfiler R only | g:Profiler REST API (language-agnostic) | 2019 (JSON API) | No R dependency needed for enrichment |
| Genome-wide background | Custom proteomics background | Best practice always | Critical for proteomics; prevents false enrichment |
| g:SCS correction (g:Profiler default) | BH FDR (standard in proteomics) | N/A | Use `significance_threshold_method: "fdr"` |

**Deprecated/outdated:**
- DAVID: Databases last updated ~2016. Do not use.
- gProfileR (old R package, capital R): Replaced by gprofiler2.
- GOrilla: Limited to human/mouse, no Reactome/KEGG.

## Open Questions

1. **E. coli organism code**
   - What we know: g:Profiler uses binomial naming (first letter + family). Most organisms are straightforward.
   - What's unclear: E. coli K12 code might be "ecoli" or "ecolik12".
   - Recommendation: Test g:Convert with "ecoli" at implementation time. If it fails, try "ecolik12". The CONTEXT.md lists "E. coli" in the curated organism list.

2. **Intersection gene extraction accuracy**
   - What we know: The `intersections` array in g:GOSt response is positional against the query gene list. Non-empty entries indicate membership.
   - What's unclear: Exact structure when `no_evidences: false` -- each entry may contain evidence code arrays.
   - Recommendation: Verify by testing with a small known gene set. The positional extraction pattern is well-established in the R client.

3. **P-value vs FDR column distinction**
   - What we know: g:GOSt with `significance_threshold_method: "fdr"` returns adjusted p-values in the `p_value` field.
   - What's unclear: Whether the returned value is the raw p-value or the FDR-adjusted value. The R documentation says "adjusted p-value."
   - Recommendation: The CONTEXT.md schema calls for both "P-value" and "FDR" columns. If the API returns only the adjusted value, set both P-value and FDR to the same value, or run a separate call with no correction for raw p-values. Most users only care about FDR anyway.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | @datagrok-libraries/test (Puppeteer-based) |
| Config file | packages/Proteomics/src/package-test.ts |
| Quick run command | `grok test --host localhost --category "Enrichment"` |
| Full suite command | `grok test --host localhost` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| MAP-01 | parseAccession extracts clean accessions; gConvert maps to gene symbols; gene column added to DataFrame | unit + integration | `grok test --host localhost --test "ID mapping"` | No -- Wave 0 |
| MAP-02 | Organism parameter passed to g:Convert correctly | unit | `grok test --host localhost --test "organism selection"` | No -- Wave 0 |
| ENRICH-01 | GO enrichment with custom background returns terms with correct schema | integration | `grok test --host localhost --test "GO enrichment"` | No -- Wave 0 |
| ENRICH-02 | KEGG enrichment returns pathway terms | integration | `grok test --host localhost --test "KEGG enrichment"` | No -- Wave 0 |
| ENRICH-03 | Reactome enrichment returns pathway terms | integration | `grok test --host localhost --test "Reactome enrichment"` | No -- Wave 0 |
| VIZ-01 | Result DataFrame has correct schema (Source, Term ID, Term Name, P-value, FDR, Gene Count, Gene Ratio, Intersection) | unit | `grok test --host localhost --test "enrichment result schema"` | No -- Wave 0 |

### Sampling Rate
- **Per task commit:** `grok test --host localhost --category "Enrichment"`
- **Per wave merge:** `grok test --host localhost`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `src/tests/enrichment.ts` -- new test file covering MAP-01, MAP-02, ENRICH-01, ENRICH-02, ENRICH-03, VIZ-01
- [ ] Add `import './tests/enrichment';` to `package-test.ts`
- [ ] No new framework install needed -- @datagrok-libraries/test already configured

### Test Strategy Notes
- **Unit tests (no network):** Test `buildEnrichmentDf()` with mock g:GOSt response data, test significant protein counting logic with known thresholds, test organism code lookup, test gene symbol deduplication.
- **Integration tests (network required):** Test actual g:Convert call with small known gene set (e.g., P04637/TP53, Q09472/EP300), test actual g:GOSt call with small gene list and custom background. These tests hit the live g:Profiler API and may be flaky -- use generous timeout (60s) and mark as integration tests.
- **Manual verification:** Dialog layout, live count updates, progress indicator behavior, highlighting. These are UI-level and cannot be fully automated with the current test framework.

## Sources

### Primary (HIGH confidence)
- [g:Profiler official documentation](https://biit.cs.ut.ee/gprofiler/page/docs) -- API capabilities, supported sources, organism coverage
- [g:Profiler API page](https://biit.cs.ut.ee/gprofiler/page/apis) -- REST endpoint reference
- [gprofiler2 R package source](https://rdrr.io/cran/gprofiler2/src/R/gprofiler2.R) -- exact REST endpoint URLs and JSON request/response structures extracted from R client code
- [gost() function documentation](https://rdrr.io/cran/gprofiler2/man/gost.html) -- result columns, domain_scope options, correction methods
- [gconvert() function documentation](https://rdrr.io/cran/gprofiler2/man/gconvert.html) -- target namespaces, response fields
- [g:Profiler 2023 update paper (NAR)](https://academic.oup.com/nar/article/51/W1/W207/7152869) -- API architecture, quarterly Ensembl data updates
- Existing Proteomics package code -- `parseAccession()`, `findProteomicsColumns()`, `SEMTYPE`, dialog patterns, menu wiring patterns

### Secondary (MEDIUM confidence)
- [gprofiler2 vignette](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html) -- usage examples, custom background workflow
- [gprofiler2 CRAN reference](https://cran.r-project.org/web/packages/gprofiler2/gprofiler2.pdf) -- parameter documentation

### Tertiary (LOW confidence)
- E. coli organism code -- needs runtime verification against g:Profiler organism list

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- g:Profiler is the locked decision; API endpoints verified from R package source
- Architecture: HIGH -- follows established patterns from DE analysis, PCA, and UniProt panel in this codebase
- Pitfalls: HIGH -- background set issue is well-documented in proteomics literature; API format verified from R package source
- API response format: MEDIUM -- derived from R package source code and documentation, not direct REST API testing

**Research date:** 2026-03-06
**Valid until:** 2026-04-06 (30 days -- g:Profiler API is stable with quarterly data updates)

# External Integrations

**Analysis Date:** 2026-05-11

This package integrates with three external HTTP services and two platform-internal packages.

## External HTTP APIs

### g:Profiler — Enrichment Analysis

- **Service:** g:Profiler (`https://biit.cs.ut.ee/gprofiler`) — University of Tartu's gene-function profiling toolkit
- **Endpoints used:**
  - `/api/convert/convert/` — gene/protein ID mapping (UniProt → gene symbol, etc.)
  - `/api/gost/profile/` — over-representation analysis against GO, KEGG, Reactome, WikiPathways, miRBase, TF, HPA, CORUM, HP
- **Called from:** `src/analysis/enrichment.ts` (`gConvert` at line ~77, `gGOSt` at line ~99)
- **Auth:** None — public API
- **Transport:** Raw `fetch()` with a 30-second timeout via `AbortController` in `fetchWithTimeout` (`src/analysis/enrichment.ts:54-71`)
- **Failure mode:** Returns descriptive error string; UI shows `grok.shell.error()` and the enrichment DataFrame is not produced
- **Concern:** Raw `fetch()` does not use `grok.dapi.fetchProxy()` → CORS risk in production environments. See CONCERNS.md "Raw `fetch()` to g:Profiler Without Proxy".

### UniProt — Protein Metadata

- **Service:** UniProt REST API (`https://rest.uniprot.org/uniprotkb/{accession}.json`)
- **Endpoint used:** Per-protein JSON metadata (gene names, descriptions, organism, function comments, GO terms)
- **Called from:** `src/panels/uniprot-panel.ts:74-83` (`fetchUniProtData` function)
- **Auth:** None — public API
- **Transport:** Raw `fetch()` (no proxy, no caching, no debounce)
- **Trigger:** `@grok.decorators.panel` registers `uniprotPanel` to fire on `currentRow` change for any DataFrame row with a `Proteomics-ProteinId` semantic-typed cell
- **Failure mode:** Panel shows "Unable to fetch UniProt data for X" with no specific error
- **Concerns:** Same CORS risk as g:Profiler; no caching → one HTTP request per row navigation. Filed as todo `2026-03-03-cache-uniprot-protein-content-to-avoid-repeated-api-calls.md`.

### Datagrok Server R Scripts

- **Service:** Internal Datagrok scripting infrastructure (R worker pool on the server)
- **Called from:** `src/analysis/differential-expression.ts`, `src/analysis/normalization.ts`
- **Transport:** `grok.functions.call('Proteomics:LimmaDE' | 'Proteomics:DeqmsDE' | 'Proteomics:VsnNormalize', {params})`
- **Cold-start:** ~3–10 seconds per call (fresh R session, library load); see CONCERNS.md for details and proposed mitigation
- **Failure mode:** All three calls have JS-side `try/catch` with client-side fallbacks (Welch's t-test, quantile normalization)

## Platform-Internal Package Dependencies

### Dendrogram Package (Datagrok Bio)

- **Purpose:** Hierarchical clustering for the heatmap viewer
- **Called from:** `src/viewers/heatmap.ts:138-145` — `DG.Func.find({package: 'Dendrogram', name: 'hierarchicalClustering'})`
- **Coupling:** Loose — discovered at runtime via `DG.Func.find`. If Dendrogram is not installed, the heatmap silently falls back to significance-based row sort.
- **Version pinning:** None — only an arity check (`inputs.length !== 4`) before calling
- **Concern:** No user-visible indication when the fallback fires; user thinks they have a "clustered" heatmap but they don't.

### `@datagrok-libraries/bio` Tree Helpers

- **Purpose:** Tree/dendrogram utilities (currently unused in active code)
- **Imported in:** `src/viewers/heatmap.ts:5-7` — imports `getTreeHelper`, `getDendrogramService`, `DistanceMetric` from `@datagrok-libraries/bio/src/trees/...`
- **Status:** Imported but not called — the code path that used them is commented out at `src/viewers/heatmap.ts:151-169`
- **Cleanup needed:** Remove the unused imports and dead code block. Filed implicitly in CONCERNS.md "Heatmap — Commented-Out Code from Prior Debug".

### `@datagrok-libraries/statistics`

- **Purpose:** Welch's t-test and Benjamini-Hochberg FDR correction
- **Called from:** `src/analysis/differential-expression.ts:runDifferentialExpression` (client-side fallback DE)
- **Functions used:** `tTest`, `fdrcorrection`
- **Coupling:** Tight — npm dependency listed in `package.json`

## Cross-Package Function Discovery Pattern

The package uses `DG.Func.find` to locate platform-wide functions at runtime rather than importing them. Two consequences:

1. **Optional dependencies work at runtime:** If Dendrogram is missing, the catch block handles it; the package still loads.
2. **No compile-time guarantees:** A typo in the function name (`'hierarchicalClustering'` vs `'hierarchical_clustering'`) only fails at runtime.

## File-Share / App Data

Demo data files (`files/demo/*.txt`, `*.tsv`) are bundled with the package and accessed via:

```ts
const text = await _package.files.readAsText('demo/proteinGroups.txt');
```

The platform exposes these at `${pkg.webRoot}/files/demo/<filename>` for HTTP access — useful for the `proteomicsDemo()` function and for `grok.data.loadTable(url)`-style integrations.

## Webhooks / Event Producers / Message Queues

None. The package is request-response only — no background tasks, no platform event listeners beyond the per-DataFrame `onCurrentRowChanged` and `onSelectionChanged` subscriptions used in `src/viewers/enrichment-viewers.ts`.

## Auth / Identity

None directly. The package inherits the Datagrok platform's authenticated session and does not manage credentials of its own. No client secrets are stored or exchanged.

---

*Integrations analysis: 2026-05-11*

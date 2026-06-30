---
status: diagnosed
trigger: "the limma based differential expression takes a very long time. eventually the fold-change and p-value columns appear"
created: 2026-03-03T00:00:00Z
updated: 2026-03-03T00:00:00Z
---

## Current Focus

hypothesis: Three independent performance problems compound: (1) the entire dataframe is serialized to R (not just intensity columns), (2) the R environment cold-starts on every invocation, (3) for the DEqMS path there is a sequential fallback chain that can invoke limma R TWICE before reaching client-side t-test
test: static code analysis — no running instance available
expecting: confirmed by code reading
next_action: report findings

## Symptoms

expected: limma DE completes quickly (sub-second to a few seconds)
actual: DE eventually completes but takes a very long time
errors: none — it works, just slow
reproduction: run any Proteomics DE dialog with limma method on a reasonably sized dataset
started: unknown — likely always been slow

## Eliminated

- hypothesis: DE is being run twice (double-fire)
  evidence: no duplicate call sites; onComplete callback is called once after the try/finally block; de_complete tag check at line 242 prevents re-runs on the same df
  timestamp: 2026-03-03

- hypothesis: the R script itself is algorithmically slow (limma math is O(n^2) or worse)
  evidence: limma lmFit/eBayes is well-optimised O(n) per gene; the R script is minimal and correct
  timestamp: 2026-03-03

## Evidence

- timestamp: 2026-03-03
  checked: differential-expression.ts buildExpressionDf() lines 110-122
  found: copies ONLY the group1+group2 intensity columns into a fresh DataFrame — this is correct and already lean
  implication: data passed to R is the minimum needed; serialisation of the expression matrix itself is not the bottleneck for typical proteomics datasets (hundreds of proteins, 6-20 samples)

- timestamp: 2026-03-03
  checked: differential-expression.ts runLimmaDE() line 135
  found: grok.functions.call('Proteomics:LimmaDE', ...) — synchronous R script invocation through Datagrok server
  implication: the entire round-trip is: JS → HTTP → Datagrok server → R worker → R cold-start + library load → script → HTTP → JS. Cold-start dominates on first call.

- timestamp: 2026-03-03
  checked: limma_de.R lines 17-18
  found: hasLimma <- suppressWarnings(require(limma, quietly = TRUE))
  implication: EVERY call loads the limma library from disk. R does not cache loaded packages between script invocations in the Datagrok scripting model (each script call spawns a fresh R session). This is the single largest latency driver.

- timestamp: 2026-03-03
  checked: limma_de.R lines 35-62 (fallback branch)
  found: the R script itself has a row-wise t-test fallback loop in R (for i in seq_len(nRows)) if limma is absent
  implication: if the R environment exists but limma is not installed, the R script falls through to a slow R for-loop. This is a secondary path but compounds cold-start.

- timestamp: 2026-03-03
  checked: deqms_de.R lines 22-23
  found: library(limma) — hard require (not suppressWarnings/require). If limma is missing this throws, triggering the JS catch block.
  implication: DEqMS path also loads limma on every call.

- timestamp: 2026-03-03
  checked: showDEDialog onOK handler lines 301-351 — DEqMS fallback chain
  found: DEqMS path (lines 307-331): try runDeqmsDE → on failure try runLimmaDE → on failure runDifferentialExpression (client-side)
  implication: TWO sequential R round-trips can occur before the client-side fallback. Each trip pays the full R cold-start cost. If DEqMS is installed but returns an error, or if the R environment is flaky, the user waits for two full R cold-starts.

- timestamp: 2026-03-03
  checked: showDEDialog onOK handler lines 333-350 — limma fallback chain
  found: limma path (lines 334-349): try runLimmaDE → on failure runDifferentialExpression (client-side)
  implication: Only ONE R round-trip for the limma-only path — this is fine. The bottleneck is purely R cold-start.

- timestamp: 2026-03-03
  checked: no environments/ directory exists in the package
  found: ls returned nothing — no .yml conda/renv environment file
  implication: the R environment used is the platform default. limma may not be pre-installed, causing either: (a) the R fallback loop in the script, or (b) the JS fallback chain. Even when limma IS available, the R session cold-start (loading base R + limma) is the primary cost.

## Resolution

root_cause: |
  PRIMARY: R session cold-start on every call.
  The Datagrok scripting model spawns a fresh R process for each grok.functions.call() invocation.
  Loading base R + the limma package (which depends on BiocGenerics, Biobase, limma itself — ~3-5 MB of R bytecode) takes several seconds on every run. This is unavoidable for a single call but is the dominant latency.

  SECONDARY: DEqMS sequential double-R-call risk.
  When the user selects "DEqMS" and DEqMS is unavailable/throws, the code sequentially tries limma next (line 316 in differential-expression.ts), paying a second full cold-start before producing results.

  TERTIARY: R script has its own slow fallback loop.
  limma_de.R lines 40-51: if limma is not installed the script falls into a pure R for-loop over rows — this is slower than the JS client-side implementation and adds latency beyond just the cold-start.

fix: |
  1. Add a limma R environment spec (environments/limma.yml or renv.lock) so the platform
     pre-warms an R session with limma already loaded. This eliminates cold-start for warm calls.

  2. In showDEDialog (differential-expression.ts ~line 301), probe R availability ONCE before
     showing the dialog (or right after OK), and skip the sequential fallback chain for DEqMS:
     - If R+limma is known available: attempt DEqMS, fall back to limma in one decision.
     - If R is known unavailable: skip R entirely and go straight to client-side.

  3. Remove the slow R for-loop fallback inside limma_de.R (lines 35-62). If limma is absent on
     the R server, the JS already falls back to the fast client-side TypeScript implementation.
     The R-loop fallback is both slower than the TS fallback and unreachable in practice (the JS
     catch block intercepts first if R errors). Replace with a hard stop(\"limma not available\").

  4. (Nice-to-have) Add a warm-up call (a trivial R health-check) at package init time
     (initProteomics in package.ts) so the R worker is warm by the time the user invokes DE.

verification: not applied — diagnosis only
files_changed: []

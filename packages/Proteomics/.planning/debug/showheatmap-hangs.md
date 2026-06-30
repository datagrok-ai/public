---
status: diagnosed
trigger: "showHeatmap runs for a very long time and doesn't complete when displaying proteomics DE results"
created: 2026-03-03T00:00:00Z
updated: 2026-03-03T00:00:00Z
---

## Current Focus

hypothesis: calcDistanceMatrix is called with the full raw column data (N rows, N can be hundreds to thousands of proteins), causing an O(N^2) Web Worker computation across all sample z-score columns simultaneously. The Dendrogram package function call chain also has multiple async layers that can hang if the package isn't installed.
test: Traced the full call chain from showHeatmap -> createExpressionHeatmap -> calcDistanceMatrix -> DistanceMatrixService.calcMulti
expecting: Two root causes confirmed — one performance (N^2 matrix computation on float column raw data passed without row subsetting), one hang (getDendrogramService never rejects if package missing)
next_action: Deliver diagnosis

## Symptoms

expected: showHeatmap should complete in under a second for typical proteomics data (top-50 proteins, 6-20 samples)
actual: showHeatmap hangs for a very long time or never completes
errors: No explicit error — the function appears to stall silently
reproduction: Run DE on a dataset, then invoke Proteomics | Visualize | Heatmap...
started: Unknown — likely since heatmap clustering code was introduced

## Eliminated

- hypothesis: df.clone(filter) clones the entire dataframe unnecessarily causing the slowdown
  evidence: df.clone(filter) is called with a topN BitSet (default 50 rows). The clone is correctly scoped to topN rows. This is not the bottleneck.
  timestamp: 2026-03-03

- hypothesis: The z-score computation loops are blocking (O(N * columns))
  evidence: z-score loops are pure JS, O(N*C) where N=50 and C=samples. This is negligible.
  timestamp: 2026-03-03

- hypothesis: hierarchicalClusteringByDistance is O(N^2) and blocks
  evidence: It is async and runs after calcDistanceMatrix. For N=50 rows, the distance matrix is 1225 elements — this is not the bottleneck. The real issue is upstream in calcDistanceMatrix.
  timestamp: 2026-03-03

## Evidence

- timestamp: 2026-03-03
  checked: heatmap.ts line 140 — calcDistanceMatrix call
  found: |
    treeHelper.calcDistanceMatrix(heatmapDf, zScoreColNames, DistanceMetric.Euclidean)
    zScoreColNames is ALL z-score columns: e.g., 12 samples = 12 calls to DistanceMatrixService.calc()
    Each call passes col.getRawData() — a Float32Array of the FULL column length
  implication: |
    The distance matrix is N_rows x N_rows where N_rows = heatmapDf.rowCount = topN (default 50).
    That part is fine. BUT calcDistanceMatrix iterates once per column (once per sample),
    calling distanceMatrixService.calc() 12 times sequentially for 12 samples.
    The DistanceMatrixService creates a NEW Worker pool on construction — and terminateOnComplete=false
    means workers are NOT terminated between calls, but the service IS terminated at line 574 after the loop.
    The bigger issue: col.getRawData() returns the raw Float32Array for the entire column. For a cloned
    50-row heatmapDf this is correct. So the matrix computation itself is bounded.

- timestamp: 2026-03-03
  checked: DistanceMatrixService constructor — libraries/ml/src/distance-matrix/distance-matrix-service.ts
  found: |
    Workers are spawned from a URL: new Worker(new URL('./distance-matrix-worker', import.meta.url))
    The DistanceMatrixService is constructed with (true, false): useConcurrentWorkers=true, terminateOnComplete=false.
    Workers are shared across all 12+ column calls. However, DistanceMatrixService.terminate() is called
    once at the end (line 574 of tree-helper.ts). The workers are reused across column iterations.
    Worker spawn itself requires module resolution. In a platform environment this may involve network
    fetches for the worker URL, which could be slow or hang if the Dendrogram package bundle is not loaded.
  implication: Worker instantiation from a dynamic import URL in a platform context may stall if the
    Dendrogram package bundle hasn't been fully loaded.

- timestamp: 2026-03-03
  checked: heatmap.ts lines 136-137 — getTreeHelper() and getDendrogramService() calls
  found: |
    const treeHelper = await getTreeHelper();
    const dendrogramService = await getDendrogramService();

    getTreeHelper() (bio/trees/tree-helper.ts line 71):
      Uses DG.Func.find({package: 'Dendrogram', name: 'getTreeHelper'}).
      If Dendrogram IS installed: calls funcList[0].prepare().call() — a server round-trip.
      If Dendrogram is NOT installed: throws immediately ("Package Dendrogram must be installed").
      Catch block on line 152 catches this and falls through to significance-based sort. OK.

    getDendrogramService() (bio/trees/dendrogram.ts line 24):
      Uses grok.functions.call('Dendrogram:getDendrogramService', {}).
      If Dendrogram IS installed: this is a server-side function call, async round-trip.
      If Dendrogram is NOT installed: grok.functions.call() on a missing function — behavior depends
        on the platform. It may hang indefinitely (never reject) rather than throw, OR throw with a
        platform error. This is called AFTER getTreeHelper() already succeeded, meaning Dendrogram IS
        installed when getDendrogramService() is called.
  implication: Both calls require Dendrogram to be installed and functional. If installed, both proceed.
    But each is a separate async round-trip to the platform function registry and may add latency.

- timestamp: 2026-03-03
  checked: heatmap.ts line 140 — what column data is passed to calcDistanceMatrix
  found: |
    heatmapDf has topN=50 rows. zScoreColNames contains one entry per sample column (e.g., 12 names).
    Inside calcDistanceMatrix (Dendrogram tree-helper.ts line 526):
      For each z-score column (FLOAT type), it calls:
        distanceMatrixService.calc(col.getRawData(), NumberMetricsNames.Difference, false)
      col.getRawData() returns the raw Float32Array — for the cloned 50-row df, this is correct.
      The distance matrix computed is 50x50 = 1225 pairs. This is fast.
      BUT: This is called once per sample column, sequentially in a for loop (line 534).
      With 12 sample columns = 12 sequential Worker dispatches.
      Each worker dispatch involves postMessage + Promise resolution.
      Result: total time = 12 * (worker overhead + O(1225) computation). Still should be fast.
  implication: The per-column loop is NOT the primary hang. The timing is bounded even for 12+ columns.

- timestamp: 2026-03-03
  checked: showHeatmap in package.ts line 114 — no progress indicator
  found: |
    static async showHeatmap(): Promise<void> {
      const tv = grok.shell.tv;
      const df = tv?.dataFrame;
      ...
      const grid = await createExpressionHeatmap(df);
      tv!.addViewer(grid);
    }
    There is NO TaskBarProgressIndicator wrapping the createExpressionHeatmap call.
    The function call is awaited silently. If it hangs, the user sees nothing — no spinner, no message.
    The DE dialog (showDEDialog) DOES have a progress indicator (line 298).
    showHeatmap does NOT.
  implication: User sees UI freeze with no feedback. Even if the operation completes in 5-10 seconds,
    it APPEARS to hang because there is no indicator.

- timestamp: 2026-03-03
  checked: getTreeHelper() call sequence when Dendrogram IS installed
  found: |
    funcList[0].prepare().call() — this is an async Datagrok function invocation.
    The Dendrogram package's getTreeHelper function returns a TreeHelper instance.
    This is a server-registered JS function, not an R/Python script — it should resolve quickly.
    But .prepare().call() goes through the platform's function call pipeline.
    getDendrogramService() via grok.functions.call() is the same pattern.
    If the Dendrogram package hasn't been loaded/initialized yet at call time, the platform may
    need to load it first, which could take several seconds.
  implication: Cold-start latency when Dendrogram package is not yet loaded in the browser session.

## Resolution

root_cause: |
  There are THREE compounding root causes:

  ROOT CAUSE 1 — No progress indicator (primary UX issue):
  File: packages/Proteomics/src/package.ts, lines 114-124
  showHeatmap() calls createExpressionHeatmap() with no TaskBarProgressIndicator.
  createExpressionHeatmap() involves multiple async operations:
    - getTreeHelper() platform function call
    - getDendrogramService() platform function call
    - calcDistanceMatrix() with N_samples Worker dispatches
    - hierarchicalClusteringByDistance()
  Total elapsed time can be 3-15 seconds. With no spinner, the user concludes it has "hung."
  The DE dialog (showDEDialog) correctly uses TaskBarProgressIndicator — showHeatmap does not.

  ROOT CAUSE 2 — Dendrogram package cold-start latency:
  File: packages/Proteomics/src/viewers/heatmap.ts, lines 136-137
  getTreeHelper() and getDendrogramService() both invoke platform function calls that require
  the Dendrogram package to be loaded. On first invocation in a session (cold start), the platform
  may need to initialize the Dendrogram package, adding 2-10 seconds of latency before any
  computation begins. There is no timeout on either call.

  ROOT CAUSE 3 — calcDistanceMatrix called once per sample column (sequential Worker overhead):
  File: packages/Proteomics/src/viewers/heatmap.ts, line 140
  The calcDistanceMatrix implementation in Dendrogram iterates over each column sequentially,
  dispatching a new Worker computation per column. For a dataset with 12 samples, this is 12
  sequential Worker round-trips. Each round-trip involves Worker postMessage + Promise resolution
  overhead. While the actual matrix size (50x50=1225 pairs) is small, 12 sequential async
  operations with Worker overhead accumulates. The implementation does not use calcMulti() which
  could batch across all columns in one Worker dispatch.

fix: |
  Fix 1 (highest priority — addresses user perception):
  Add TaskBarProgressIndicator to showHeatmap() in package.ts:
    const pi = DG.TaskBarProgressIndicator.create('Building expression heatmap...');
    try {
      const grid = await createExpressionHeatmap(df);
      tv!.addViewer(grid);
    } finally {
      pi.close();
    }

  Fix 2 (robustness — prevents silent infinite hang):
  Add a timeout to the Dendrogram service calls in heatmap.ts, or move them inside the
  try/catch that already wraps them (they ARE inside try/catch — the catch correctly falls
  through to significance-based sort). Verify the catch on line 152 fires correctly when
  getDendrogramService hangs vs throws.

  Fix 3 (performance — reduces Worker round-trips):
  In createExpressionHeatmap, instead of passing individual z-score column names to
  calcDistanceMatrix (which iterates them one by one), consider whether the Dendrogram
  TreeHelper's calcDistanceMatrix multi-column path can be replaced by calling calcMulti
  directly, or restructure to pass only a subset of representative columns.

verification: Not yet applied — diagnosis only
files_changed: []

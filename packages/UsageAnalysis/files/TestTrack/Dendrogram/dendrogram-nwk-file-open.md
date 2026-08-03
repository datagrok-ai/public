---
feature: dendrogram
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [dendrogram.cp.newick-file-open-via-files-browser]
realizes: [dendrogram]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - dendrogram-nwk-file-open-spec.ts
gate_verdicts:
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T17:12:11Z
    spec_runs:
      - spec: dendrogram-nwk-file-open-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 117
        failure_keys: []
---

# Dendrogram — Open a .nwk / .newick file via the Files browser

Covers opening and previewing a `.nwk` / `.newick` file from the Files
browser. Two file-role handlers are exercised end-to-end:

- `importNwk` (the `meta.role: fileHandler` for `ext: nwk, newick`,
  `public/packages/Dendrogram/src/package.ts#L274-L294`) — fires when the user
  **opens** a `.nwk` / `.newick` file; parses the newick string via
  `TreeHelper.newickToDf` and opens a `DendrogramApp` on the resulting tree
  DataFrame.
- `previewNewick` (the `meta.role: fileViewer` for `fileViewer: nwk,newick`,
  `public/packages/Dendrogram/src/package.ts#L297-L317`) — fires when the user
  **previews** a `.nwk` / `.newick` file in the Files browser (file selected
  but not opened); reads the file via `file.readAsString()`, calls
  `TreeHelper.newickToDf` on it, and renders a `PhylocanvasGL` preview view.

Scenario 1 is the smoke open path through `importNwk` → `DendrogramApp`.
Scenario 2 is the regression preview path through `previewNewick` →
`PhylocanvasGL`. Both scenarios share a single small synthetic 4-leaf newick
fixture written into the Dendrogram package's `AppData` folder during Setup so
the same file feeds both flows; this preserves the round-trip property that
`importNwk` and `previewNewick` produce the same leaf set from the same
on-disk bytes.

## Setup

- A clean Datagrok session (no preloaded tables, no other dendrograms
  attached); user is logged in with permission to write under
  `System:AppData/Dendrogram/`.
- The Dendrogram package is installed and registered. In particular the
  following two function registrations are present:
  - `Dendrogram:importNwk` — `meta.role: fileHandler`, `meta.ext: nwk, newick`
    (registered at `public/packages/Dendrogram/src/package.ts#L274`).
  - `Dendrogram:previewNewick` — `meta.role: fileViewer`,
    `meta.fileViewer: nwk,newick` (registered at
    `public/packages/Dendrogram/src/package.ts#L297`).
- A small synthetic 4-leaf newick fixture is written to the Dendrogram
  package's `AppData` files folder at `System:AppData/Dendrogram/sample.nwk`
  (the path the Files browser surfaces under `Files | App Data |
  Dendrogram | sample.nwk`). Newick content:
  `((a:1,b:1):1,(c:1,d:1):1);`. The four leaf names are `a, b, c, d`.
- The fixture is written via the Datagrok JS API
  (`grok.dapi.files.writeAsText('System:AppData/Dendrogram/sample.nwk',
  '((a:1,b:1):1,(c:1,d:1):1);')`) so the test owns its own data and does
  not depend on a pre-seeded server-side file. Cleanup at the end of the
  scenario file removes the fixture
  (`grok.dapi.files.delete('System:AppData/Dendrogram/sample.nwk')`).

## Scenarios

### Scenario 1 — Open `.nwk` from Files browser, `importNwk` opens `DendrogramApp` (smoke)

Steps:

1. Open the Files browser and navigate to `Files | App Data | Dendrogram`.
   Wait for the folder listing to populate and assert the fixture
   `sample.nwk` is visible in the listing.

   * Expected result: a row with the name `sample.nwk` is present in the
     Files browser listing for the `Dendrogram` AppData folder.
2. Open the file. The canonical open gesture in the Files browser is the
   double-click on the file row (equivalent to right-click row → `Open`);
   when neither is available to the driver, the open path can also be
   exercised by calling
   `grok.functions.call('Dendrogram:importNwk',
   {fileContent: await grok.dapi.files.readAsText('System:AppData/Dendrogram/sample.nwk')})`
   which routes to the exact same `importNewick` handler the file-open
   gesture dispatches to (`public/packages/Dendrogram/src/package.ts#L285`).

   * Expected result: the `importNewick` handler runs without throwing and
     a new view is added to the workspace. No console error is emitted
     during the open.
3. Wait for the `DendrogramApp` view to mount. The app surface is the
   `DendrogramApp` view object instantiated at
   `public/packages/Dendrogram/src/package.ts#L290` (`new DendrogramApp();
   await app.init(df)`) where `df` is the result of
   `TreeHelper.newickToDf(fileContent, '')`.

   * Expected result: a new view becomes the active view. The view hosts
     a `[name="viewer-Dendrogram"]` viewer element (the embedded
     Dendrogram viewer) and its `canvas` child has `canvas.width > 0` and
     `canvas.height > 0`.
4. Assert the parse-then-render round trip. The newick file content
   `((a:1,b:1):1,(c:1,d:1):1);` parses into a DataFrame with one row per
   tree node (leaves first, then internal nodes), per the
   `TreeHelper.newickToDf` contract at
   `public/packages/Dendrogram/src/utils/tree-helper.ts#L30`. The resulting
   DataFrame is tagged with `.newick` (the original newick string) and
   `.newickJson` (the parsed JSON form).

   * Expected result: the active view's underlying DataFrame is tagged
     with `.newick` whose value equals the file content
     `((a:1,b:1):1,(c:1,d:1):1);` exactly. The rendered Dendrogram
     viewer's `treeNewick` resolves to the same string (the `.newick`
     tag path of the priority order at
     `public/packages/Dendrogram/src/viewers/dendrogram.ts#L309` —
     property is unset on this code path because `importNwk` does not
     set the viewer's `newick` property, only the DataFrame's
     `.newick` tag via `TreeHelper.newickToDf`).
5. Assert the leaf set. Read the leaves off the rendered tree via
   `TreeHelper.getLeafList(rootNode)` (or, equivalently, off the parsed
   DataFrame's `leaf` column where the row's `leaf` value is `true`).

   * Expected result: the leaf set is `{a, b, c, d}` (set semantics,
     no duplicates). No `Non unique key tree leaf name` exception is
     raised.

### Scenario 2 — Preview `.nwk` from Files browser, `previewNewick` renders `PhylocanvasGL` (regression)

Continuing from the Scenario 1 fixture (the `sample.nwk` file already
written; re-write at the start of this scenario if cleanup ran between
scenarios). The preview path is distinct from open: clicking the file row
without opening it surfaces the registered `fileViewer` for the extension,
which for `.nwk` / `.newick` is `Dendrogram:previewNewick`
(`public/packages/Dendrogram/src/package.ts#L297`).

Steps:

1. Open the Files browser and navigate to `Files | App Data | Dendrogram`.
   Click the `sample.nwk` row exactly once (a click that selects the row
   without opening the file — the canonical preview gesture). Wait for the
   preview pane to populate.

   * Expected result: the file preview pane (the right-hand or bottom
     pane that the Files browser surfaces for the selected file) becomes
     non-empty. No console error is emitted during preview mount.
2. Assert that the preview pane hosts a `PhylocanvasGL` viewer root. The
   preview function reads the file via `file.readAsString()`, builds a
   tree DataFrame via `TreeHelper.newickToDf(newickString,
   file.fileName.slice(0, -4))`, then mounts a `PhylocanvasGL` viewer via
   `df.plot.fromType('PhylocanvasGL', {})` and wraps the resulting root
   into a `DG.View.fromRoot(viewerRoot)`
   (`public/packages/Dendrogram/src/package.ts#L306-L316`).

   * Expected result: the preview pane contains an element whose viewer
     `name` attribute is `viewer-PhylocanvasGL` (or the conventional
     `[name="viewer-PhylocanvasGL"]` selector the PhyloTreeViewer package
     surfaces). The viewer root's `width` and `height` are both set to
     `100%` per the explicit style override at
     `public/packages/Dendrogram/src/package.ts#L313-L314`.
3. Assert the parse contract. The `TreeHelper.newickToDf` call inside
   `previewNewick` runs on the exact same byte content as
   Scenario 1's `importNwk` — both flows feed the file's raw newick
   string into `newickToDf`. The resulting DataFrame's leaf set must
   match the leaf set asserted in Scenario 1.

   * Expected result: the leaves rendered by the PhylocanvasGL preview
     view are `{a, b, c, d}` — identical to Scenario 1's leaf set.
     (Where the driver cannot read leaf labels off the PhylocanvasGL
     canvas, assert through the preview view's underlying DataFrame:
     count rows where `leaf` is `true` and read off the `node`
     column for those rows. Cardinality should be 4 and the set should
     match.)
4. Cleanup. Close the preview pane (click away from the file row or
   navigate away from the Files browser). Delete the fixture file via
   `grok.dapi.files.delete('System:AppData/Dendrogram/sample.nwk')` so
   the AppData folder is restored to its pre-test state and re-running
   the scenario file is idempotent.

   * Expected result: `grok.dapi.files.list('System:AppData/Dendrogram')`
     no longer returns a `sample.nwk` entry. No exception is raised
     during delete.

## Notes

- Selector note: the DOM selectors used for the `DendrogramApp` view
  and the `PhylocanvasGL` preview view (`[name="viewer-Dendrogram"]`,
  `[name="viewer-PhylocanvasGL"]`, the Files-browser row selector, the
  preview-pane selector) were not verified against a live `.nwk`
  sample during authoring — confirm them against the running UI
  before relying on this spec as a regression baseline.

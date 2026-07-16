---
feature: powerpack
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [powerpack.cp.xlsx-file-open]
realizes: []
produced_from: atlas-driven
original_path: D:\work\datagrok\reddata\public\packages\UsageAnalysis\files\TestTrack\Powerpack\xlsx-open.md
migration_date: 2026-05-23
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-19329
realized_as:
  - xlsx-open-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-05-22-powerpack-migrate-05
    timestamp: 2026-05-23T00:00:00Z
    review_round: 1
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T14:05:00Z
    failure_keys: []
  d:
    verdict: EVIDENCE_GAP
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T17:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T13:14:00Z
    spec_runs:
      - spec: xlsx-open-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 39
        failure_keys: []
---

# PowerPack — XLSX file open across all entry paths (GROK-19329 regression)

Regression test for GROK-19329: PowerPack's XLSX file handler must
keep working across every way Datagrok lets a user open an `.xlsx`
file — Browse / My Files, Recent files, drag-and-drop into the window,
the top menu **File > Open**, and Shared with me. Each of these paths
hands the file's bytes to the same handler, which parses the workbook
and surfaces one table (DataFrame) per sheet. In the 1.27.0
regression, this broke completely: XLSX files did not open from any
of the 5 entry paths, no data was shown, and no error was surfaced.

## Setup

1. Open Datagrok with PowerPack installed (default platform load —
   `powerPackInit` runs at startup and registers `xlsxFileHandler`
   as the `@fileHandler` for the `xlsx` extension via
   `grok.functions.register(...)`).
2. Ensure an XLSX test fixture is available in two locations:
   - **Platform fixture (`System:DemoFiles/`).** Verify
     `System:DemoFiles/SPGI-linked.xlsx` (or any pre-existing
     `.xlsx` under `System:DemoFiles/`) exists and is readable.
     If no XLSX exists under `System:DemoFiles/`, upload one via
     `My Files` first (Step 2 below).
   - **User fixture (`Home/`).** Upload `Home/xlsx-open-test.xlsx`
     into the current user's `My Files` (Home directory). The
     fixture should be a small (< 1 MB) multi-sheet workbook —
     for example a 3-sheet workbook with sheets `Customers`,
     `Orders`, and `Products`, each with a header row plus 5-10
     data rows of mixed-type columns (numeric, string, date). If
     no convenient fixture exists locally, use Datagrok's
     `Customers.xlsx` sample if present under `tables/` of any
     installed package, or generate a fresh fixture with `grok s
     files save Home/xlsx-open-test.xlsx --from-file <local-path>`.
3. Open the browser DevTools Console (F12) so console errors are
   visible while opening XLSX files — the GROK-19329 regression
   surfaces as either a thrown error in the console (handler /
   Web Worker init failure) or as a silent no-op (file does not
   open, no view rendered).

## Scenarios

### Scenario 1: Open XLSX from Browse / My Files (canonical GROK-19329 reproduction)

1. **Navigate to Browse.** Click the `Browse` icon in the sidebar
   (or `Browse` in the top-left menu) to open the file browser
   view.
2. **Open the `My Files` (Home) tree.** Click `My Files` in the
   browse navigation so the current user's home directory listing
   surfaces.
3. **Locate the XLSX fixture.** The `xlsx-open-test.xlsx` fixture
   uploaded in Setup is visible in the listing.
4. **Single-click (or double-click per default Browse behavior)
   the XLSX file.** The `xlsxFileHandler` registered for the
   `xlsx` extension is invoked with the file's bytes; the
   `ExcelJSService` Web Worker parses the file.
5. **Verify the file opens.** A new TableView opens as the active
   view rendering the first sheet (e.g. `Customers`) as a Datagrok
   grid. The header row maps to grid columns; data rows surface in
   the grid body.
6. **Verify additional sheets surface.** If the fixture is
   multi-sheet, the second / third sheets surface as additional
   tabs / sibling views (one DataFrame per sheet per
   `ExcelJSService.parse`).
7. **Verify no console errors.** DevTools Console shows no
   exception (no `Cannot read property ... of undefined` from the
   handler, no Web-Worker init failure, no `ExcelJSService` null
   reference). No Datagrok red error balloon.

Expected:
- The XLSX file opens into a Datagrok TableView with data
  rendered from the first sheet — this is the canonical
  GROK-19329 invariant.
- Multi-sheet workbooks surface one view per sheet.
- No thrown exception in the console; no error balloon.

### Scenario 2: Open XLSX from Recent files

1. **Pre-condition.** Scenario 1 has already been run in this
   session so the XLSX fixture appears in the platform's Recent
   files list. If Scenario 1 was not run, open the fixture once
   via Browse / My Files first to populate Recent files.
2. **Navigate to Recent files.** Click `Recent files` in the
   sidebar (or `Browse > Recent` if the Recent files entry lives
   there).
3. **Locate the XLSX fixture** in the Recent files listing.
4. **Click the XLSX file.** The handler is invoked via the
   Recent-files navigation path.
5. **Verify the file opens.** Same expected behavior as Scenario
   1 Steps 5-7: TableView opens with data rendered; multi-sheet
   sheets surface as additional views; no console errors.

Expected:
- Opening via Recent files routes through the same registered
  `xlsxFileHandler` as Browse / My Files — both must succeed.
- No regression-specific divergence between the Recent files
  navigation path and the Browse / My Files path.

### Scenario 3: Open XLSX via drag-and-drop into the platform window

1. **Navigate to the home / Welcome view** (or any non-Browse
   view — the drag-and-drop handler is registered globally on
   the Datagrok window).
2. **Drag the XLSX file from the local filesystem** (the OS file
   manager — Windows Explorer / macOS Finder) onto the Datagrok
   browser window.
3. **Drop the file anywhere on the Datagrok window.** The
   platform's global drag-and-drop handler reads the dropped
   file's bytes and dispatches to `xlsxFileHandler` (registered
   for `xlsx` extension).
4. **Verify the file opens.** Same expected behavior as Scenario
   1 Steps 5-7: TableView opens with data rendered; multi-sheet
   sheets surface as additional views; no console errors.

Expected:
- Drag-and-drop dispatches to the same registered XLSX handler
  as Browse / My Files and Recent files — all three paths must
  succeed.
- The dropped file is parsed in-memory (no upload to `My Files`
  required for the open action to succeed).

### Scenario 4: Open XLSX via top menu `File > Open`

1. **Open the top menu `File > Open...`** (or the equivalent
   open-file action in the current Datagrok build's top menu).
2. **The OS native file picker opens** (Datagrok delegates to
   the browser's file-input element backed by the OS picker).
3. **Select the XLSX fixture** (`xlsx-open-test.xlsx`) from the
   local filesystem.
4. **Click `Open` in the picker.** The handler is invoked with
   the selected file's bytes.
5. **Verify the file opens.** Same expected behavior as Scenario
   1 Steps 5-7: TableView opens with data rendered; multi-sheet
   sheets surface as additional views; no console errors.

Expected:
- The top menu `File > Open` path routes through the same
  registered `xlsxFileHandler` as the other entry paths — all
  four exercised so far must succeed.

### Scenario 5: Open XLSX from Shared with me

1. **Pre-condition.** A second user (or the current user via
   `grok s shares add`) must have shared an XLSX file with the
   current user. If no shared XLSX is available, run
   `grok s shares add "OtherUser:xlsx-shared-test.xlsx" <current-user>
   --access View` from a precondition setup step, OR ask another
   QA user to share an XLSX file via the Datagrok share dialog
   before running this scenario.
2. **Navigate to `Shared with me`.** Click `Shared with me` in
   the sidebar / Browse navigation.
3. **Locate the shared XLSX file** in the `Shared with me`
   listing.
4. **Click the shared XLSX file.** The handler is invoked.
5. **Verify the file opens.** Same expected behavior as Scenario
   1 Steps 5-7: TableView opens with data rendered; multi-sheet
   sheets surface as additional views; no console errors.

Expected:
- The `Shared with me` path routes through the same registered
  `xlsxFileHandler` as Browse / My Files, Recent files,
  drag-and-drop, and File > Open — all 5 entry paths must
  succeed (the GROK-19329 invariant: XLSX opens from ANY entry
  path).

### Scenario 6: Multi-sheet selection via sheetName parameter (optional)

This scenario covers the optional `sheetName` parameter: if the
platform UI exposes a sheet-selector dialog on multi-sheet workbooks,
this scenario verifies that path too.

1. **Open the multi-sheet XLSX fixture via any of Scenarios 1-5's
   paths** — Browse / My Files is the simplest.
2. **If the platform prompts with a sheet-selector dialog** on
   multi-sheet workbooks, select a non-default sheet (e.g.
   `Orders` instead of `Customers`).
3. **Verify the selected sheet's data is rendered** as the
   primary TableView (header row mapped to grid columns; rows
   from the selected sheet).
4. **If no sheet-selector dialog is exposed**, the default
   behavior (one DataFrame per sheet, all surfaced) applies —
   verify all sheets are accessible per Scenario 1 Step 6.

Expected:
- If a sheet selector is exposed, selecting a non-default sheet
  renders that sheet's data.
- If no selector is exposed, all sheets surface as separate
  TableViews / tabs.
- `sheetName` parameter behavior is consistent with the atlas
  declaration (`xlsxFileHandler(bytes, sheetName?)` — the
  parameter is optional and defaults to opening all sheets).

## Notes

- **Related bug.** GROK-19329: XLSX files stopped opening from every
  entry path (Browse / My Files, Recent files, drag-and-drop, File >
  Open, Shared with me) in the 1.27.0 dev build — no data was shown
  and no error was surfaced. The handler must work from all five entry
  paths again.
- **Deferrals.**
  - A headless test that exercises the parser directly (one DataFrame
    per sheet, header de-duplication, default column/sheet names) was
    deferred — this regression was about routing files to the handler
    from each entry path, not about parsing correctness.
  - A boundary test for the handler's 80 MB file-size limit was
    deferred — this scenario's fixture is intentionally small (under
    1 MB) for speed, and the size limit isn't part of the GROK-19329
    reproduction.
  - Scenario 5 needs a second user account to share a file with the
    current user; this is the same precondition as `data-enrichment.md`'s
    cross-user sub-scenario, and can be run manually if no such
    account is available.

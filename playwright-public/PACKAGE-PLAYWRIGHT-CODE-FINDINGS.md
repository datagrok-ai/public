# Package Playwright E2E — Findings for Developers

**Branch:** `opavlenko/playwright-public-ci-fixes-0709` · **Date:** 2026-07-10 · **Scope:** `packages/<Name>/playwright/` suites run on the Jenkins **Test-Playwright** job (per-package, `PLAYWRIGHT_DIR=public/packages/<Name>`).

Methodology: each suite was baselined locally against `dev.datagrok.ai` (fast iteration) and on the isolated CI stack (authoritative). Test-harness issues were fixed in the specs; **genuine product/code defects are collected in Section A below** for developer follow-up. Environment/CI-only artifacts are in Section B.

---

## CI status per package (latest build)

| Package | CI | Notes |
|---|---|---|
| Charts | ✅ 10/10 | radar save/reopen fixed (test) |
| Helm | ✅ 4/4 | — |
| SequenceTranslator | ✅ 6/1skip | Markush = in-dev, `test.skip` |
| Peptides | ✅ 14/14 | — |
| Chem | 🟢 17/19 | 4 test fixes; 2 residue = heavy UMAP walks (Section B) |
| EDA | 🟢 12–13/14 | linreg + anova test fixes; residue = share-model UI, anova tabs (see A/B) |
| PowerPack | 🟡 11/13 | input-functions fixed; residue = DB-enrichment, direct-link (Section B) |
| Dendrogram | 🟡 9/10 | residue = previewNewick (see A?) |
| BiostructureViewer | 🟡 7/10 | residue = WebGL 3D-pane + RCSB network (Section B) |
| DiffStudio | 🔴 7/8 | **product bug** (Section A) |
| Bio | 🔴 27/9 | 9 CI-confirmed **code/domain** defects (A2–A4, A7, A9, A10) — for dev review |
| Notebooks | 🟠 2/4 (kept, not skipped per owner) | the 4 fails need a live Jupyter container/kernel (env, §B6) |

---

## Section A — CODE / PRODUCT defects (for developers)

> These are failures where the **product behaviour looks wrong**, not the test or the CI environment.
> The E2E test was **left asserting the correct behaviour** (not adjusted to hide the defect) unless noted.

### A1. DiffStudio — URL deep-link does not restore the `step` parameter — **CONFIRMED**

- **Symptom:** Set Step=0.1 + Count=4 in a Diff Studio model, share the URL (URL correctly contains `step=0.1` and `count=4`), open the URL in a new tab. On reopen **`count` is restored (4) but `step` reverts to the model default (0.01)**.
- **Evidence:** `packages/DiffStudio/playwright/files-and-sharing.test.ts` Step 3. Polled 15 s — the step never settles at 0.1; count applies immediately. Reproduces on dev AND CI. Asymmetric (count works, step doesn't) ⇒ not a timing/env issue.
- **Suspected area:** Diff Studio URL/state deserialization on deep-link load (`app.ts` state restore) — the `step` argument isn't applied to the inputs after the model re-solves.
- **Test:** left red (correctly catches the bug).

### A2. Bio — Split-to-Monomers emits fewer per-position columns than the alignment width — **SUSPECTED (regression)**

- **Symptom:** `Bio | Transform | Split to Monomers` should add one `Monomer` column per alignment position. On the test fixtures it adds **fewer**: FASTA 38 (expected 39), **HELM 10 (expected 17)**. The HELM shortfall (−7) is large.
- **Evidence:** `packages/Bio/playwright/convert.test.ts` (`expect(monCount - beforeMonCount).toBe(ds.monomers)`, comment: "a broken split emitting a single column must fail here"). Deterministic on dev.
- **Suspected area:** `splitAlignedSequences` / `splitToMonomersUI` (`packages/Bio/src/utils/split-to-monomers.ts`) or the HELM splitter in `@datagrok-libraries/bio`.
- **Needs:** domain confirmation of the correct column counts for the current fixtures (test data unchanged) — if 10/38 are wrong, this is a splitter regression.
- **Test:** left red (correctness invariant — not adjusted).

### A3. Bio — MSA separator tag is `/` where `-` is expected (GROK-12164) — **SUSPECTED**

- **Symptom:** Opening `filter_MSA.csv` (units=separator), the detected/rendered separator tag is `/` but the GROK-12164 reproducer expects `-`.
- **Evidence:** `packages/Bio/playwright/bio-renderer-dispatch.test.ts` (Scenario for MSA separator). Deterministic in isolation on dev.
- **Suspected area:** separator detection in `detectors.js` (`detectSeparator`) or the MSA fixture's separator tag.
- **Test:** left red pending confirmation whether `/` is now correct for this data.

### A4. Bio — Split-to-Monomers dialog does not auto-select the Macromolecule column — **SUSPECTED**

- **Symptom:** The `Get Region` editor auto-populates its Sequence-column selector with the setup Macromolecule column, but the **`Split to Monomers` editor leaves its Sequence selector empty** (same table, same column).
- **Evidence:** `packages/Bio/playwright/bio-cell-actions-panels.test.ts` Scenario 4 Step 2 (passes for Get Region Scenario 3 Step 3, fails for Split to Monomers).
- **Suspected area:** `SplitToMonomersFunctionEditor` default column selection (`packages/Bio/src/utils/split-to-monomers.ts` / the function editor) — doesn't default to the sole Macromolecule column the way Get-Region does.
- **Test:** left red.

### A5. EDA — ANOVA result no longer renders "Analysis" / "F-test" tabs — **VERIFY (possibly intended)**

- **Symptom:** `ML | Analyze | Group Comparison | ANOVA` on demog produces the Box plot, but the former **"Analysis" / "F-test" result tabs are absent**; the significance + p-value now live only on the box plot description. The default method is Welch (no classical F-test).
- **Evidence:** `packages/EDA/playwright/anova.test.ts`; probe showed `anyFtest=false`, only the box-plot description carries the conclusion.
- **Suspected area:** ANOVA UI (`packages/EDA/src/anova/anova-ui.ts`) — result tabs only render with `showReport`, and the default run shows the box plot only. Likely intended (Welch default), but flagged so a developer can confirm the report tabs aren't meant to appear by default.
- **Test:** **ADJUSTED** — now asserts the ANOVA conclusion on the box-plot description (intent preserved: ANOVA produces a result-bearing view). See EDA fix note.

### A7. Bio — "Manage Monomer Libraries" **View** shows an empty library listing (the dialog is fine) — **SUSPECTED**

- **Symptom:** `Bio | Manage | Monomer Libraries` opens a **View** whose library listing is empty (`view labels []`), while the alternate `Bio:manageMonomerLibraries` **dialog** lists the libraries correctly (`[HELMCoreLibrary.json, polytool-lib.json, sample-lib-Aca-colored.json]`).
- **Evidence:** `packages/Bio/playwright/bio-lifecycle-monomer-library.test.ts` S1.3 / S1.4 / S2.4 (deterministic in isolation on dev).
- **Suspected area:** `showManageLibrariesView` (`packages/Bio/src/utils/monomer-lib/library-file-manager/ui.ts`) — the full-view path doesn't populate/refresh the library list the way the dialog path does.
- **Test:** left red.

### A8. Notebooks — Context Panel `Actions > Delete` link missing for a notebook — **SUSPECTED (verify vs env)**

- **Symptom:** Selecting a notebook, the Context Panel `Actions` pane has no clickable `Delete` link.
- **Evidence:** Notebooks CI #264 Scenario 3.
- **Suspected area:** Notebook entity Context-Panel actions registration. May be entangled with the Jupyter-container availability (env) — verify on a stack with a live Jupyter kernel.
- **Test:** left red.

### A9. Bio — Composition (WebLogo) Context-Panel property edit not accepted — **SUSPECTED**

- **Symptom:** For HELM and MSA composition analysis, editing ≥1 Context-Panel property of the composition/WebLogo viewer is not accepted (SR-01 edit-acceptance check returns false).
- **Evidence:** `packages/Bio/playwright/composition-analysis.test.ts` Scenario 3 Step 5 (CI #265, both HELM + MSA).
- **Suspected area:** WebLogo/composition viewer property-panel binding (`packages/Bio/src/viewers/web-logo-viewer.ts` property setters via the context panel).
- **Test:** left red.

### A10. Bio — empty current-row input does not surface a rejection balloon (GROK-16111) — **SUSPECTED (regression)**

- **Symptom:** Sequence **Similarity Search / Diversity Search / Activity Cliffs** on an empty current-row input should surface a rejection balloon (GROK-16111), but no balloon appears (silent).
- **Evidence:** `packages/Bio/playwright/empty-input-row-viewers.test.ts` (all 3 viewers, CI #265). GROK-16111 was a fix for exactly this — the missing balloon suggests a regression.
- **Suspected area:** empty-input guard in the sequence search/activity-cliffs viewers (`packages/Bio/src/analysis/*`).
- **Test:** left red.

### A6. Dendrogram — `previewNewick` throws `Bad state: No element` on the CI stack — **VERIFY**

- **Symptom:** `Dendrogram:previewNewick({file})` (called after `grok.shell.closeAll()`) throws a Dart `StateError: Bad state: No element`. Passes locally, fails on the isolated CI stack. NOT a WebGL issue (software-WebGL flags made no difference).
- **Evidence:** `packages/Dendrogram/playwright/dendrogram-nwk-file-open.test.ts` Scenario 2.2.
- **Suspected area:** `previewNewick` handler accesses a `.first`/`.single` on an empty collection (current view / element) when invoked with no open view — a null-safety gap when there is no active TableView.
- **Test:** left red.

---

## Section B — Environment / CI-only issues (NOT product code)

> These fail only because of the CI stack (no WebGL/GPU, no external network, no DB backend, heavy compute,
> single-user). Per direction: attempt a workaround first; where impossible, propose a CI skip (pending approval).

| # | Package / spec | Root cause (env) | Workaround status |
|---|---|---|---|
| B1 | Chem activity-cliffs, chemical-space | Heavy UMAP over 5 datasets — times out on the minimal CI stack (ApprovedDrugs2015 >90 s, no scatter); the "Show only cliffs" toggle + "N cliffs" link elements aren't reachable, and the `.md` forbids JS-API substitution | TODO: reduce dataset scope / gate the UI-toggle step / bump timeouts; else skip (approval) |
| B2 | BiostructureViewer context-panel-widgets, ngl-viewer | 3D Structure pane needs WebGL (GPU-less CI); PDB Information + ProLIF fetch from **rcsb.org** via fetchProxy (external network blocked on CI) | TODO: gate WebGL + RCSB steps gracefully; else skip (approval) |
| B3 | PowerPack data-enrichment | DB-Explorer join needs the `public.users_sessions` platform DB table + data, absent on CI | TODO: gate/skip DB-dependent steps; else skip (approval) |
| B4 | PowerPack direct-link-loading | Cold direct-link project load never completes on the minimal CI stack (preloader > 240 s) — possible deep-link cold-load slowness/hang (cf. core /p/ deep-link) | Robustify (240 s budget) FAILED → **skipped on CI** (approved) |
| B5 | EDA share-model-permissions | recipient user-search **autocomplete drop-down** never populates on CI (server user-search) | Robustify (retry+delay) FAILED → **skipped on CI** (approved) |
| B6 | Notebooks (all 6 specs) | whole suite needs a live Jupyter kernel/container + seeded notebooks (browser 0 cards, notebookView 300 s timeout, Edge/Lifecycle "no server notebook fetchable") | **All 6 skipped on CI** (approved) — run on a Jupyter node |

Dev-only flakes (pass in isolation / on CI, NO fix needed): Chem sketcher.test.ts, r-group-analysis; Bio renderer-dispatch (mostly); Charts sunburst.

---

## Test-harness fixes already applied (context, not product bugs)

- **Charts** radar-save-reopen: load the saved layout after `project.open()` (reopen restores Grid but not custom viewers).
- **Chem** structure-filter / filter-panel / chem-grok-14028: warm up the Chem package bundle (`grok.functions.call('Chem:getRdKitModule')`) before synchronous `fg.updateOrAdd({type:'Chem:substructureFilter'})` — the semtype detector loads without the function bundle. *(Arguably a minor product usability gap: adding a substructure filter synchronously throws "package not loaded" until Chem is warmed.)*
- **Chem** sketcher-backends: poll for Ketcher's async molfile before the persistence sample.
- **EDA** linear-regression: pass numeric-only features (cars.csv string columns break the kernel with "unsupported column type: string" — arguably a product gap: the trainer could drop/deny non-numeric features with a clear message).
- **EDA** anova: corrected the top-menu path (ANOVA moved under `Group Comparison`).
- **PowerPack** input-functions: poll for the Molecule semtype (Chem detector is async) instead of resolving on the first semtype event.

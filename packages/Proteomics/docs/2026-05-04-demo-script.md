# Cytokinetics Demo Script — Datagrok Proteomics

**Practice goal:** be able to run this in 15 min cold, recover from any single broken step, and adapt to a 5-min lightning version if needed.

---

## Pre-flight (do this 30 minutes before the meeting)

```bash
cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics
rm -rf dist node_modules/.cache && npx webpack && grok publish local
```

Then in the browser:

1. Open `http://localhost:8080`, log in.
2. Verify the menu: top bar → **Proteomics** should be visible. If not, hard-refresh.
3. Run a dry pass through Section A below using `proteinGroups.txt` (the small synthetic file). If anything errors, adjust before the call.
4. Have a backup tab pointed at `https://datagrok.ai` for the platform pitch if your laptop misbehaves.
5. Close all open Datagrok views and refresh once before the call (avoids stale state).

---

## The Pitch (memorize — say once at the start)

> "Datagrok is a data analytics platform built specifically for scientific data. The Proteomics package brings end-to-end LC-MS analysis into that environment — from raw vendor exports through normalization, imputation, statistical testing, and biological enrichment, with everything visualized interactively against your data. The point is not to replace your statisticians' R scripts — it's to put the same statistics in their hands without leaving the browser, on data they can explore, share, and version."

Lead with the hook, then go to the demo. Don't read feature lists.

---

## Section A — Full 15-minute Demo

### A1 (2 min) — Open the platform, set the stage

- Browser at `http://localhost:8080`. Show the home page briefly, point at left sidebar (Files, Apps, Databases, Platform).
- Say: *"Datagrok is one platform — connectors, viewers, scientific compute, all sharing one in-memory data engine. The proteomics tools we'll see are a plugin, but they speak the same language as everything else here."*
- Click the **Proteomics** menu in the top bar. Show the items: Import, Annotate Experiment, Analyze, Visualize, Enrichment Analysis. Don't click yet.

### A2 (3 min) — Load a real benchmark dataset

- **Proteomics → Import → MaxQuant…**
- Pick `cptac-spike-in.txt` from `packages/Proteomics/files/demo/` (have it open in Finder ready, or use the dialog's file browser to navigate to it).
- Talking points while it loads:
  - *"This is the CPTAC spike-in benchmark — UPS proteins spiked into a constant yeast background at known ratios across two conditions. Public, peer-reviewed, ~1,500 proteins, 6 samples."*
  - *"The importer auto-detects MaxQuant column conventions, splits intensity columns by sample, and tags semantic types — protein IDs, gene names, intensities — so the rest of the platform knows what they are."*
- When the table opens, scroll right to show LFQ intensity columns. Click a Protein IDs cell and point at the **Properties panel** on the right (gene info populates from UniProt).

### A3 (1 min) — Annotate the experiment **(prerequisite for everything that follows)**

- **Proteomics → Annotate Experiment…**
- Two columns of intensity samples appear side-by-side. Drag (or shift-click) the **6A** samples into Group 1, the **6B** samples into Group 2. Name them e.g. `6A` and `6B`.
- Talking points:
  - *"This is where the platform learns what your experiment is — which samples are condition vs control. One step, no R script preamble defining a design matrix."*
  - *"Group assignment gets written back to column tags, so every downstream step — QC, normalize, DE — knows the design without you re-specifying it. Skip this step and the QC dashboard will tell you to come back here first."*
- Click OK.

### A4 (1 min) — Quick QC

- **Proteomics → Visualize → QC Dashboard…**
- *(Requires A3 — without annotated groups, the dashboard prints "Annotate experimental groups first" and exits. If you see that warning, you skipped A3.)*
- Talk while it builds:
  - *"Before any statistics, I want to know if my samples are comparable. Sample-wise distributions, missingness pattern, sample correlation, density plots — all colored by group now that we've annotated. The dashboard also computes per-group CV and an MA-trend, so you can see batch-like effects between conditions before they corrupt your DE result."*
- Point out anything visibly off (sample with shifted distribution, high CV in one group) — *"if this were a clinical study, this is where I'd flag a sample for re-injection before the stats team got hold of it."*

### A5 (2 min) — Normalize

- **Proteomics → Analyze → Normalize…**
- Show the dialog: method dropdown (Quantile / VSN / median / log2), live box-plot preview reacting to method selection.
- Talking points:
  - *"VSN is the Bioconductor standard for label-free MS — runs in R server-side. Quantile is pure JS, instant. Both methods, one dialog."*
  - *"If the data is already normalized — like Spectronaut output with PG.NormalizedQuantity — the dialog detects it and warns you."*
- Pick **VSN** (or Quantile if R env isn't reliable on your laptop today). Click OK.

### A6 (2 min) — Impute, then DE

- **Proteomics → Analyze → Impute Missing Values…** — pick **kNN** (k=5), valid-values filter at 50%. Click OK. *"Missing-not-at-random vs missing-at-random matters for proteomics — we expose the choice."*
- **Proteomics → Analyze → Differential Expression…**
- Show the dialog: method (limma / DEqMS / t-test), comparison direction picker, FC and p-value thresholds with hint text showing what they mean for your data. Group 1 / Group 2 are pre-filled from the annotation step.
- Pick **limma**. *"limma's moderated t-test is the de facto standard for proteomics DE. It runs server-side in R via Bioconductor — and we just declared the conda environment so the server pre-warms it instead of cold-starting per call."*
- Click OK. New columns appear in the dataframe (log2FC, p-value, adjusted p-value, significance flag).

### A7 (2 min) — Visualize the result

- **Proteomics → Visualize → Show All Visualizations…**
- Volcano plot, PCA, heatmap pop up side-by-side, all linked to the same DataFrame.
- *"This is the demo punchline — every viewer is reading the same in-memory column store. Click a point on the volcano, the heatmap row highlights, the PCA dot highlights. Filter on the volcano, the heatmap subset updates."*
- On the volcano: lasso a few significantly-up proteins. Show that they're recognizable UPS spike-ins (gene names visible in tooltips).

### A8 (2 min) — Biology

- **Proteomics → Enrichment Analysis…** — pass the gene-symbol column. Pick organism: Human.
- *"This calls g:Profiler — GO/KEGG/Reactome — with proteomics-aware background. Returns enriched terms with FDR-corrected p-values."*
- **Proteomics → Visualize → Enrichment Charts…** — show dot plot or bar chart.
- *"From raw vendor TSV to enrichment in the same browser tab, ~5 minutes. No data export, no separate Bioconductor session."*

### A9 (1 min) — Close

- *"That was end-to-end on a public benchmark. The whole stack is open source under the public Datagrok repo, the package itself ships with these demo datasets, and everything you saw runs against your real data the same way. If you have proteomics workflows you'd want to reproduce, that's the conversation I'd love to continue."*
- Don't oversell. Stop talking.

---

## Section B — 5-minute Lightning Version

If time is short, run only:

1. **A2** with `proteinGroups.txt` (small, fast, predictable) — 1 min
2. **A3** Annotate Experiment — Sample1–3 → Group 1, Sample4–6 → Group 2. **Don't skip this** — DE, QC, and the linked viewers all key off the group tags. — 30 sec
3. **A5** Normalize with Quantile method (no R round-trip, instant) — 45 sec
4. **A6** DE with **t-test** method (client-side, no R) — 45 sec
5. **A7** Show All Visualizations — 1 min, click a point to demonstrate viewer linking
6. Close with one sentence — 30 sec

---

## Likely Questions — Prepared Answers

**"Is this just a UI on top of R?"**
No. The data engine is native — columnar in-memory, runs the same in browser and on the server. R/Python are integrations for specific functions (limma, VSN). Quantile normalization, t-tests, viewers, filtering, enrichment API calls — all native, no R required.

**"How does this compare to Perseus / Skyline / Spectronaut Pivot?"**
Perseus is desktop, Windows-only, and you can't share a session as a URL. Spectronaut Pivot is great inside Spectronaut but doesn't go beyond DIA-MS. We're a multi-omics platform — proteomics is one package among many (cheminformatics, sequence analysis, db connectors). The unique angle is that everything is one engine you can build on.

**"Can I run my own R / Python scripts?"**
Yes. The platform exposes scripts as functions with annotated metadata (inputs, outputs, environments). Our limma DE is a 70-line R script with a `#environment:` declaration that lets the platform pre-warm a Bioconductor session. You can drop in your own scripts the same way and they get a UI for free.

**"What about DIA / Spectronaut data?"**
Yes — there's an `Import → Spectronaut…` path that handles the long-format report directly. We have a Spectronaut HYE three-species mix as a demo dataset (`spectronaut-hye-mix.tsv`) if you want to see it.

**"How do you handle large cohorts? 100+ samples?"**
The data engine is built for this — ChEMBL with 2.7 M molecules loads in the browser. For proteomics, the bottleneck moves to your DE method, not the platform. limma scales linearly. We've not pushed past a few hundred samples with this specific package yet — happy to benchmark something representative.

**"Reproducibility / audit?"**
Every transformation is a function call recorded in the project's history. Datagrok projects can be saved as URLs, versioned, and replayed. That's a platform-level capability, not specific to our package.

**"What's missing?"**
Be honest. Currently: PTM-aware DE, MS imaging, label-based quant (TMT/iTRAQ) are not yet first-class. The package is at v1.3 — focused on label-free DDA and DIA. Roadmap includes those.

**"Can I have it?"**
The plugin is open source in the Datagrok public repo. The platform itself is commercial; there's a free tier on `public.datagrok.ai` and an enterprise version they'd talk to you about.

---

## If Something Breaks

**`Proteomics` menu missing after publish:**
Hard-refresh (`Cmd+Shift+R`). If still missing, run `grok s packages list --filter "Proteomics"` in another terminal to confirm publish landed; if version `admin` shows, the bundle is on the server but the browser cached the old function registry.

**QC Dashboard prints "Annotate experimental groups first":**
You skipped A3, or the annotation didn't stick (e.g., re-imported the file mid-demo and lost the column tags). Run **Proteomics → Annotate Experiment…** and try the dashboard again. Same warning will appear from any analysis step that needs groups, so this is the canonical recovery.

**R-based normalize/DE hangs or errors:**
Switch the method to **Quantile** (normalize) and **t-test** (DE). Both are pure JavaScript and never touch R. The flow still tells the story.

**Heatmap appears to hang:**
Known issue — captured in `.planning/debug/showheatmap-hangs.md`. If it stalls > 10 sec, skip it. Volcano + PCA still demo the linked-viewer story.

**Enrichment Analysis returns nothing:**
g:Profiler requires internet. If you're on a flaky connection, skip A8 and close on A7.

**Whole platform looks broken:**
Refresh once. If still broken, switch to the backup tab at `https://datagrok.ai`, give the pitch from there, and offer to follow up with a recorded walkthrough.

---

## Practice Checklist

Run this checklist twice today before the meeting:

- [ ] `grok publish local` succeeded with no warnings about stale `dist/`
- [ ] Section A2 — CPTAC file imports cleanly, ~1500 rows visible
- [ ] Section A3 — Annotate Experiment dialog accepts the 6A / 6B split, OK closes cleanly
- [ ] Section A4 — QC Dashboard renders **without** the "Annotate experimental groups first" warning (proves A3 stuck)
- [ ] Section A5 — VSN normalize completes (or pre-decided to use Quantile)
- [ ] Section A6 — limma DE returns log2FC + p-value columns
- [ ] Section A7 — clicking a volcano point highlights heatmap row + PCA dot
- [ ] Section A8 — Enrichment returns at least one enriched GO term
- [ ] You can give the opening pitch in under 30 seconds without notes
- [ ] You know which step you'll skip first if you're running long (suggested order: drop A8 → drop A7 heatmap → skip A4 QC)

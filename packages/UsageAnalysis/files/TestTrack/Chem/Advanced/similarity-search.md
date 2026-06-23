---
feature: chem
sub_features_covered: [chem.search.similarity, chem.search.similarity.viewer, chem.search.similarity.top-menu]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/similarity-search.md
migration_date: 2026-05-11
migration_report: similarity-search-migration-report.md
related_bugs: []
---

# Similarity Search

End-to-end exercise of the Chem Similarity Search viewer property-panel
surface on `smiles.csv`. Single phase, four beats:

- **Add** the Chem Similarity Search viewer via **Top menu > Chem >
  Search > Similarity Search...** — search results display with the
  most-similar molecules listed.
- **Open** the viewer's property panel via the gear icon on the Chem
  Similarity Search panel header.
- **Exercise** the six property-panel knobs in turn — **Fingerprint**,
  **Limit**, **Distance Metric**, **Size**, **Molecule Properties**,
  **Cutoff** — verifying each modification re-queries the viewer
  without errors or crashes.

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario
(`depends_on: []`), `classification: medium`, `pyramid_layer: integration`,
`target_layer: playwright`, strategy `simple`. Realizes atlas critical
path `chem.cp.similarity-search-viewer` (p0) end-to-end across the six
property-axis dimensions. UI coverage owned
(`ui_coverage_delegated_to: null`).

## Setup

1. **Provision linked dataset.** The scenario consumes a single bundled
   Datagrok file (no external provisioning required):
   - `System:DemoFiles/chem/smiles.csv` — primary SMILES dataset
     (`canonical_smiles` molecule column, the input to the Chem
     Similarity Search viewer).
2. **Confirm Chem package is loaded** so that **Chem > Search >
   Similarity Search...** is registered as a top-menu entry (per atlas
   `chem.search.similarity.top-menu`, package source
   `Chem/src/package.ts#L543`) and the Chem Similarity Search viewer
   (`chem.search.similarity.viewer`) is registered.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. The
   dataset is opened fresh inside the scenario.

## Scenarios

### Similarity Search viewer add + property-panel six-knob walk

End-to-end walk: open dataset (via TestTrack 'star' icon "Open test
data") → add viewer via top menu → open gear / property panel → exercise
six property knobs (Fingerprint, Limit, Distance Metric, Size, Molecule
Properties, Cutoff) in turn.

1. Open the `smiles.csv` dataset by pressing the **star** icon in
   TestTrack (tooltip: **Open test data**). The
   `System:DemoFiles/chem/smiles.csv` dataset opens into a table view.
   Verify the `canonical_smiles` column auto-detects as a Molecule
   column (RDKit renderer renders cells; SMILES auto-detection
   populates units / semType / cell renderer).
2. From the top menu, run **Chem > Search > Similarity Search...**.
3. Verify: the **Chem Similarity Search** viewer is added to the
   active table view. The search results display the most-similar
   molecules from `smiles.csv` (top-N hits ordered by similarity
   score; default Fingerprint Morgan + default Distance Metric
   Tanimoto + default Limit). No console errors fire during viewer
   add or initial query.
4. Click the **gear** icon on the Chem Similarity Search panel header
   to open the viewer's property settings. Verify: the property panel
   renders on the Context Panel (viewer-properties mode) without
   errors; the six property controls — **Fingerprint**, **Limit**,
   **Distance Metric**, **Size**, **Molecule Properties**, **Cutoff**
   — are visible and interactable.
5. **Fingerprint knob.** Change the **Fingerprint** property and
   re-test the viewer for each available fingerprint type (Morgan,
   RDKit, MACCS, AtomPair, TopologicalTorsion, Pattern — per atlas
   `chem.search.similarity` description). Verify: each fingerprint
   change re-queries the viewer and renders the result set
   ("Most similar structures" content updates); no console errors
   appear and the viewer does not crash for any fingerprint type.
6. **Limit knob.** Increase and decrease the **Limit** property
   (number of hits returned). Verify: the rendered hit count updates
   to reflect the new limit; the viewer re-queries cleanly; no
   console errors.
7. **Distance Metric knob.** Change the **Distance Metric** property
   and re-test the viewer for each available metric (Tanimoto,
   Asymmetric, Cosine, Sokal — per atlas `chem.search.similarity`
   description). Verify: each metric change re-queries the viewer
   and renders updated similarity scores; no console errors; no
   viewer crashes for any metric.
8. **Size knob.** Test the **Size** property across all available
   values (**small**, **normal**, **large**). Verify: the rendered
   hit-tile size updates accordingly; viewer layout adapts without
   errors; no console errors for any size value.
9. **Molecule Properties knob.** Add a few entries from the
   **Molecule Properties** multi-select list. Verify: the selected
   properties are appended to the "Most similar structures"
   information for each hit tile (per-hit property strip rendered
   beneath the molecule); no console errors.
10. **Cutoff knob.** Set the **Cutoff** property to **1**. Verify:
    only one molecule remains in the Chem Similarity Search viewer
    when the cutoff equals 1 (the exact-match self-hit at similarity
    1.0); the viewer re-queries cleanly; no console errors.
11. **Final assertion.** Across all six property modifications
    (steps 5-10), verify that no errors or crashes occurred — the
    viewer remained interactable throughout the property-panel walk
    and re-queried correctly for each knob change.

## Notes

- **`coverage_type: regression`** — `pyramid_layer: integration` per
  chain rev 2; the scenario walks six property-axis knobs against the
  single Chem Similarity Search viewer. It is not a single-bug repro
  (no `related_bugs`), not a one-knob smoke (six knobs exercised), and
  not a perf test (no timing or volume thresholds). It IS a
  property-axis regression walk over the chem.cp.similarity-search-
  viewer p0 critical path — protecting against feature surface erosion
  as fingerprint / metric / size enums are extended or property-panel
  layout is refactored. The section-level A-STRUCT-02 edge|perf
  coverage path is satisfied chain-wide via the 10
  `bug_focused_candidates[]` in `scenario-chains/chem.yaml` rev 2 (see
  migration report SR-01).
- **Top-menu phrasing.** The canonical menu path is **Top menu > Chem
  | Search | Similarity Search...** (with trailing ellipsis), per
  atlas `chem.search.similarity.top-menu` description. The original
  scenario's "Similatity search" was a source-text typo silently
  corrected during migration.
- **Gear icon discovery.** The gear icon on the Chem Similarity
  Search panel header is the in-viewer affordance for "open property
  settings"; clicking it switches the Context Panel into
  viewer-properties mode for the Similarity Search viewer. The
  property-panel surface is the six controls enumerated under
  `ui_coverage_responsibility` in chain YAML
  (`chem-similarity-search-fingerprint-knob`,
  `chem-similarity-search-limit-knob`,
  `chem-similarity-search-distance-metric-knob`,
  `chem-similarity-search-size-knob`,
  `chem-similarity-search-molecule-properties-multi-select`,
  `chem-similarity-search-cutoff-knob`) plus the gear-open action
  itself (`chem-similarity-search-gear-property-panel`).
- **No JS API substitution.** Every step in this scenario is a
  declared UI flow per chain `ui_coverage_responsibility`. UI driving
  via Playwright is required: the top-menu add
  (`chem-add-similarity-search`), the gear-icon click
  (`chem-similarity-search-gear-property-panel`), and each property-
  panel knob exercise. Substituting `grok.dapi.*` or direct viewer
  `setOptions({...})` for the property-panel UI driving would defeat
  the regression's intent (it asserts the property-panel WIDGETS
  drive the underlying viewer state correctly, not just that the
  viewer accepts options programmatically). The existing sibling
  `Advanced/similarity-search-spec.ts` currently uses `setOptions`
  for fingerprint/limit/metric/cutoff knobs as a property-existence
  smoke; the migrated scenario tightens the contract to property-
  panel UI driving for the regression walk.
- **Fingerprint enum.** Six fingerprint types per atlas
  `chem.search.similarity` (Morgan / RDKit / MACCS / AtomPair /
  TopologicalTorsion / Pattern). The original "Change and test
  different fingerprints" is enumerated explicitly in step 5 as
  "each available fingerprint type" with the six values surfaced —
  resolution of an implicit enumeration in the source text. Spec-
  time the Automator should pin to the runtime-available enum
  (read from the property-panel dropdown or the viewer's property
  metadata) rather than hardcoding the list.
- **Distance Metric enum.** Four metrics per atlas
  `chem.search.similarity` (Tanimoto / Asymmetric / Cosine / Sokal).
  Same resolution pattern as Fingerprint — original "Change and
  test different distance metrics" enumerated explicitly in step 7;
  Automator should pin to runtime enum at spec time.
- **Size enum.** Three values per scenario text — small / normal /
  large. Step 8 walks all three.
- **Molecule Properties multi-select.** Step 9's "a few properties
  from the list" is intentionally bounded (not "all properties") —
  the assertion is that selected properties appear in the hit-tile
  information strip, not that every available property has been
  selected. Automator picks 2-3 deterministic property entries at
  spec time (e.g. MW, LogP, HBA from the standard molecule-property
  set).
- **Cutoff = 1 expected hit count.** Step 10's "only one molecule
  should remain when the cutoff equals 1" is the exact-match self-
  hit assertion — the query molecule's own row matches itself at
  similarity 1.0; all other hits are below 1.0 by definition (no
  duplicates in `smiles.csv`). If the dataset contains exact
  duplicate SMILES strings, the hit count at cutoff = 1 may exceed
  1; current `smiles.csv` is assumed deduplicated per dataset
  curation.
- **Helpers usage.** No registered helper in `helpers-registry.yaml`
  currently abstracts the "open dataset → add similarity-search
  viewer → open gear / drive property knobs" pattern. The spec
  Author (Automator) will use existing `helpers-registry.yaml`
  entries — `loginToDatagrok`, `softStep`, `closeAllViews` for the
  standard section harness — and may inline the Similarity-Search-
  specific DOM driving until a dedicated helper lands. Candidate
  helpers surfaced in migration report Decisions section.
- **Order in chain.** `order: 3` per source JSON; tie-broken
  lexicographically AFTER `calculate.md` per chain `order_from_files`
  (line 92-99). No `must_run_last` constraint.
- **Source-text defect fixes (silent).** Two occurrences of typo
  "Similatity" → "Similarity" silently corrected during migration
  (scenario title and Step 2 menu path). Top-menu phrasing aligned
  to canonical atlas form "Chem | Search | Similarity Search..."
  (with trailing ellipsis). Per chain rev 2 directive (source-text
  defects in scenario bodies — Migrator fixes during migration; see
  decision-log Olena 2026-05-11 directive).

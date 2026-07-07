---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [monomer_library_file]
realized_as:
  - st-lifecycle-monomer-library-file-spec.ts
related_bugs: []
---

# SequenceTranslator — Monomer Library File Lifecycle

Checks how SequenceTranslator's monomer library file is created, loaded, and
falls back to a bundled sample library when the configured one is missing —
covering the CSV-to-JSON library generator, lazy library loading on first
use, and the monomer-to-molecular-weight lookup.

## Setup

1. Ensure the SequenceTranslator package is loaded.
2. Open Datagrok.
3. Verify that `System:AppData/SequenceTranslator/monomers/monomer-lib.json`
   is accessible. The package falls back to
   `System:AppData/SequenceTranslator/monomers-sample/monomer-lib.json` when
   the configured library path is missing.

## Scenarios

### Scenario 1: createMonomerLibraryForPolyTool — CSV to JSON download

Steps:
1. In Datagrok Files browser, navigate to a CSV file that represents a
   monomer library (e.g., a CSV with columns for monomer symbol, SMILES,
   molecular weight).
2. Right-click the CSV file.
3. Select `SequenceTranslator | createMonomerLibraryForPolyTool`.
4. Observe the browser download prompt.

Expected:
- `createMonomerLibraryForPolyTool` reads the CSV via `file.readAsString()`,
  runs `PolyToolCsvLibHandler.getJson()` to build the monomer library JSON,
  and triggers a download of the `.json` file via `DG.Utils.download`.
- The download dialog appears and the `.json` file is saved successfully.

### Scenario 2: initLibData lazy loading and fallback

Steps:
1. Note the current `MonomersPath` package property in SequenceTranslator's
   package settings.
2. Open any Oligo Toolkit app (e.g., Oligo Translator) — this triggers
   `initLibData()` if not yet cached.
3. Confirm the app opens without errors (library loaded from the configured
   path or from the `monomers-sample` fallback).
4. (Optional) Set `MonomersPath` to a non-existent path in package settings.
5. Re-open the Translator app.

Expected:
- On first open, `initLibData()` runs once (loading JSON data files) and
  caches the promise; subsequent opens reuse the cache.
- When `MonomersPath` points to a missing location, the app still opens
  using the `monomers-sample` fallback path; a user-visible notice or
  warning indicates fallback is active.
- No crash or blank view occurs.

### Scenario 3: getCodeToWeightsMap reflects the loaded monomer library

Steps:
1. After the Oligo Toolkit app has been opened (ensuring `initLibData` ran),
   call via JS API:
   ```js
   const weights = await grok.functions.call(
     'SequenceTranslator:getCodeToWeightsMap');
   ```
2. Inspect the returned map.

Expected:
- The map is a non-empty `Record<string, number>`.
- Keys correspond to monomer codes from the loaded library
  (e.g., `A`, `C`, `G`, `U` for RNA standard monomers).
- Values are positive numbers (molecular weights in Da).
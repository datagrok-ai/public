# Phase 0 V1 — B-factor auto-coloring verification

Blueprint reference: section 5, "Phase 0 — Scaffold the new package", V1 paragraph.
Open question: OQ13.

## Why

The whole v1 renderer story is *frequency → B-factor → Mol\* colors the consensus pharmacophore*.
If Mol\* does NOT apply a B-factor color scheme automatically when given a PDB block with varied
B-factors, the visual UX is broken — frequency is encoded but invisible. We need to know this
**before** locking the renderer design in Phase 5b.

## How to run (~5 minutes in a Datagrok session)

1. `grok publish dev` from `packages/ConsensusPharmacophore/` (publishes the empty scaffold; nothing
   user-visible registers).
2. In a Datagrok browser tab, open the **Console** (Tools → Console) and paste:

   ```js
   const pdb = await fetch('/api/files/system/AppData/ConsensusPharmacophore/_reference/v1-bfactor-test.pdb').then(r => r.text());
   const df = grok.data.demo.demog(1);
   const view = grok.shell.addTableView(df);
   const viewer = await df.plot.fromType('Biostructure', {pdb: pdb});
   view.dockManager.dock(viewer, 'right', null, 'Mol*', 0.5);
   ```

   (If the test PDB is easier to paste inline, copy `v1-bfactor-test.pdb` into the JS string
   directly.)

3. **Look at the 7 atoms in the docked Mol\* viewer.** They lie on the X-axis at 3 Å spacing.
4. Decide:
   - **PASS** — atoms render in a B-factor gradient (typically blue→white→red, or similar).
     Phase 0 done; proceed to Phase 6.
   - **PASS-with-toggle** — atoms render in element CPK by default, BUT changing the Mol\*
     property panel's "Coloring" or similar to "B-factor" / "Temperature" shows the gradient.
     Phase 0 done with an extra task: surface the property toggle from the orchestrator
     (`viewer.setOptions({colorScheme: 'b-factor'})` or equivalent — find the actual key by
     inspecting the property panel in Mol\*).
   - **FAIL** — neither element CPK nor any property toggle produces a frequency gradient.
     Re-scope: drop the "frequency" visual encoding from v1 (consensus model still ships; just
     no per-atom intensity), or invest in MolViewSpec / custom Mol\* extension (out of v1 scope).

## What to write back

After running, paste the outcome into this file under a `## Result` heading and re-run
`npm run build` (no code change needed — this is a documentation artefact).

## Result

(pending — fill in after manual verification)

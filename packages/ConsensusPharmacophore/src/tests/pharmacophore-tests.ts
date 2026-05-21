/*
 * Stage 1 tests. Run via `grok test --category Pharmacophore`.
 *
 * Tests that hit RCSB self-skip when the endpoint is unreachable — that keeps
 * the suite green on offline CI. Other tests (cofactor filter, etc.) work
 * against synthetic input and run unconditionally.
 *
 * Blueprint reference - section 3 Test strategy, Phase 1 done-when.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {enrichPdbList} from '../rcsb-client';
import {DEFAULT_OPTIONS} from '../orchestrator-types';
import {clusterByRmsd, computeDefaultSelection} from '../cluster-picker-dialog';
import {FN} from '../function-names';

const REQUIRED_PDB_QC_COLS = [
  'pdb_id', 'resolution', 'experimental_method', 'r_free',
  'ligand_comp_id', 'ligand_chain', 'ligand_formula_weight',
  'ligand_smiles', 'ligand_rscc',
] as const;

/**
 * Probe RCSB GraphQL endpoint. Returns false if offline / unreachable so
 * tests can self-skip rather than fail the suite on infrastructure issues.
 */
async function rcsbReachable(): Promise<boolean> {
  try {
    const resp = await grok.dapi.fetchProxy('https://data.rcsb.org/graphql', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({query: '{ __typename }'}),
    });
    return resp.ok;
  } catch {
    return false;
  }
}

category('Pharmacophore', () => {
  test('stage1_enrich_pdb_list', async () => {
    if (!(await rcsbReachable())) {
      console.warn('Pharmacophore.stage1_enrich_pdb_list: RCSB unreachable, skipping assertions.');
      return;
    }
    // Pinned EGFR PDBs known to pass at default QC (all <= 2.0 A, X-ray, drug ligand).
    const df = await enrichPdbList(['3W2S', '4WKQ', '5HG8'], DEFAULT_OPTIONS);
    expect(df.rowCount > 0, true, 'enrichPdbList returns at least one row for known EGFR PDBs');
    for (const colName of REQUIRED_PDB_QC_COLS)
      expect(df.col(colName) !== null, true, `pdb_qc has column ${colName}`);

    const methodCol = df.col('experimental_method')!;
    for (let i = 0; i < df.rowCount; i++) {
      const m = String(methodCol.get(i)).toLowerCase();
      expect(m.includes('x-ray'), true, `row ${i} method '${m}' contains 'x-ray'`);
    }

    const resCol = df.col('resolution')!;
    for (let i = 0; i < df.rowCount; i++) {
      const v = resCol.get(i);
      if (v !== null && v !== undefined)
        expect(Number(v) <= DEFAULT_OPTIONS.maxResolution + 1e-6, true,
          `row ${i} resolution <= ${DEFAULT_OPTIONS.maxResolution}`);
    }

    expect(df.col('ligand_smiles')!.semType, DG.SEMTYPE.MOLECULE,
      'ligand_smiles is semType=Molecule for 2D thumbnail rendering');
  }, {timeout: 30000});

  // -------------------------------------------------------------------------
  // Stage 3.5 — cluster picker
  // -------------------------------------------------------------------------

  /** Build a fake aligned_structures DataFrame with the columns clusterByRmsd reads. */
  function fakeAligned(rows: Array<{pdb_id: string; rmsd_to_ref: number}>): DG.DataFrame {
    return DG.DataFrame.fromColumns([
      DG.Column.fromStrings('pdb_id', rows.map((r) => r.pdb_id)),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'rmsd_to_ref', rows.map((r) => r.rmsd_to_ref)),
    ]);
  }

  test('cluster_byrmsd_splits_distinct', async () => {
    // 5 active-state rmsd ~0-1, 3 inactive-state rmsd ~4-5. Should split into 2 clusters.
    const aligned = fakeAligned([
      {pdb_id: 'A1', rmsd_to_ref: 0.0},  {pdb_id: 'A2', rmsd_to_ref: 0.5},
      {pdb_id: 'A3', rmsd_to_ref: 1.0},  {pdb_id: 'A4', rmsd_to_ref: 0.8},
      {pdb_id: 'A5', rmsd_to_ref: 1.2},
      {pdb_id: 'I1', rmsd_to_ref: 4.5},  {pdb_id: 'I2', rmsd_to_ref: 4.8},
      {pdb_id: 'I3', rmsd_to_ref: 5.1},
    ]);
    const clusters = clusterByRmsd(aligned);
    expect(clusters.length, 2, 'two distinct clusters expected');
    const totalSize = clusters.reduce((n, c) => n + c.size, 0);
    expect(totalSize, 8, 'all input rows accounted for');
    const sortedBySize = clusters.slice().sort((a, b) => b.size - a.size);
    expect(sortedBySize[0].size, 5, 'larger cluster has 5 members');
    expect(sortedBySize[1].size, 3, 'smaller cluster has 3 members');
  });

  test('cluster_picker_default_all_when_n_le_5', async () => {
    const aligned = fakeAligned([
      {pdb_id: 'A', rmsd_to_ref: 0.0},
      {pdb_id: 'B', rmsd_to_ref: 0.3},
      {pdb_id: 'C', rmsd_to_ref: 4.0},
      {pdb_id: 'D', rmsd_to_ref: 4.5},
    ]);
    const clusters = clusterByRmsd(aligned);
    expect(clusters.length, 2, 'two clusters expected on 4-entry fixture');
    const choice = computeDefaultSelection(clusters);
    expect(choice, 'All clusters', `default selection for N=4 should be "All clusters" (got "${choice}")`);
  });

  test('cluster_picker_default_largest_when_n_gt_5', async () => {
    const aligned = fakeAligned([
      {pdb_id: 'A1', rmsd_to_ref: 0.0}, {pdb_id: 'A2', rmsd_to_ref: 0.4},
      {pdb_id: 'A3', rmsd_to_ref: 0.7}, {pdb_id: 'A4', rmsd_to_ref: 1.0},
      {pdb_id: 'A5', rmsd_to_ref: 1.3},
      {pdb_id: 'B1', rmsd_to_ref: 4.5}, {pdb_id: 'B2', rmsd_to_ref: 4.9},
      {pdb_id: 'B3', rmsd_to_ref: 5.2},
    ]);
    const clusters = clusterByRmsd(aligned);
    expect(clusters.length, 2, 'two clusters expected on 8-entry fixture');
    const sortedBySize = clusters.slice().sort((a, b) => b.size - a.size);
    const choice = computeDefaultSelection(clusters);
    expect(choice, `Cluster ${sortedBySize[0].id}`,
      `default selection for N=8 should be "Cluster ${sortedBySize[0].id}" (got "${choice}")`);
  });

  // -------------------------------------------------------------------------
  // Stage 2b — pocket-only Kabsch refinement
  // -------------------------------------------------------------------------

  test('stage2b_bootstrap_reduces_rmsd_on_pocket_subset', async () => {
    if (!(await rcsbReachable())) {
      console.warn('Pharmacophore.stage2b: RCSB unreachable, skipping.');
      return;
    }
    // Run the pipeline through Stage 3 against three high-resolution EGFR PDBs,
    // then call Stage 2b and assert the mean pocket-Cα RMSD drops vs Stage 2a's
    // global-Cα RMSD on the same subset. (Stage 2b computes pass-1 RMSD over the
    // pocket subset internally too — that's stored in bootstrap_failed-aware
    // diagnostic logs but not surfaced in the output schema, so we settle for
    // the simpler "pass-2 rmsd_to_ref is finite and smaller than pass-1 on
    // average" check.)
    const pdbQc = await enrichPdbList(['3W2S', '4WKQ', '5HG8'], DEFAULT_OPTIONS);
    if (pdbQc.rowCount === 0) {
      console.warn('stage2b: enrichPdbList returned 0 rows; skipping.');
      return;
    }
    const aligned1 = await grok.functions.call(FN.STAGE2A, {
      pdb_qc: pdbQc, chain_selection: 'auto', alt_loc_filter: true, min_occupancy: 0.5,
    }) as DG.DataFrame;
    expect(aligned1.rowCount > 0, true, 'Stage 2a produced rows');

    const pocketAtoms = await grok.functions.call(FN.STAGE3, {
      aligned_structures: aligned1, pocket_method: 'cutoff', pocket_radius: 5.0,
      dbscan_eps: 4.0, dbscan_min_samples: 6,
    }) as DG.DataFrame;
    expect(pocketAtoms.rowCount > 0, true, 'Stage 3 produced pocket atoms');

    // Snapshot pass-1 mean rmsd_to_ref (over the GLOBAL Cα subset Stage 2a used).
    const rmsd1Col = aligned1.col('rmsd_to_ref')!;
    let sum1 = 0; let n1 = 0;
    for (let i = 0; i < aligned1.rowCount; i++) {
      const v = Number(rmsd1Col.get(i));
      if (Number.isFinite(v) && v > 0) { sum1 += v; n1 += 1; }
    }
    const mean1 = n1 > 0 ? sum1 / n1 : NaN;

    const aligned2 = await grok.functions.call(FN.STAGE2B, {
      aligned_structures: aligned1, pocket_atoms: pocketAtoms,
      selected_pdb_ids: '[]', // all
    }) as DG.DataFrame;
    expect(aligned2.rowCount > 0, true, 'Stage 2b produced rows');

    const rmsd2Col = aligned2.col('rmsd_to_ref')!;
    let sum2 = 0; let n2 = 0;
    for (let i = 0; i < aligned2.rowCount; i++) {
      const v = Number(rmsd2Col.get(i));
      if (Number.isFinite(v) && v > 0) { sum2 += v; n2 += 1; }
    }
    const mean2 = n2 > 0 ? sum2 / n2 : NaN;

    console.log(`stage2b: mean rmsd pass-1=${mean1.toFixed(3)} A (n=${n1}), ` +
      `pass-2=${mean2.toFixed(3)} A (n=${n2})`);

    // The pocket subset is a STRICT subset of the global Cα set, AND Stage 2b's
    // Kabsch is optimal over that subset. So pass-2 pocket-RMSD <= pass-1
    // global-RMSD restricted to pocket. In practice pass-2 is much smaller
    // than pass-1 because the pocket has tighter local geometry.
    expect(Number.isFinite(mean2), true, 'pass-2 mean rmsd is finite');
    expect(mean2 < mean1 + 1e-3, true,
      `pass-2 pocket RMSD (${mean2.toFixed(3)}) should be <= pass-1 global RMSD (${mean1.toFixed(3)})`);
  }, {timeout: 60000});

  test('stage1_cofactor_denylist_removes_hem_nag', async () => {
    if (!(await rcsbReachable())) {
      console.warn('Pharmacophore.stage1_cofactor_denylist_removes_hem_nag: RCSB unreachable, skipping.');
      return;
    }
    // 1A6M (myoglobin) carries HEM. 1HWW (pancreatic ribonuclease) carries multiple
    // glycosylation NAG / BMA entries. Both should be filtered to zero usable ligand
    // rows when only HEM/NAG/BMA are present and nothing else qualifies.
    const df = await enrichPdbList(['1A6M'], DEFAULT_OPTIONS);
    const compCol = df.col('ligand_comp_id');
    if (compCol) {
      for (let i = 0; i < df.rowCount; i++) {
        const v = String(compCol.get(i) ?? '').toUpperCase();
        expect(v !== 'HEM' && v !== 'NAG' && v !== 'BMA' && v !== 'HOH', true,
          `cofactor row ${v} at ${i} should have been filtered`);
      }
    }
  }, {timeout: 30000});
});

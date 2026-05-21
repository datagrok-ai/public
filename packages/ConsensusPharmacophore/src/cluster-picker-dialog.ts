/*
 * Mid-pipeline UI gate. Surfaces conformational clusters detected from Stage 2a's
 * rmsd_to_ref values; lets the user pick one (typical kinase case is active vs.
 * inactive conformations) before Stage 2b runs.
 *
 * Algorithm - greedy 1D clustering on rmsd_to_ref (blueprint OQ9, threshold 1.5 A).
 * Two structures with similar rmsd_to_ref can still be in different conformations,
 * but in practice the rmsd-to-best-resolution-template scalar is a strong enough
 * signal for the active-vs-inactive split that this is designed to catch.
 *
 * Blueprint reference - section 0 Q13, Phase 3.5, OQ9.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Distance threshold in A for greedy 1D clustering of rmsd_to_ref values. */
export const RMSD_CLUSTER_THRESHOLD_A = 1.5;

export interface ClusterPick {
  clusterIds: string[];    // ['all'] when "All clusters" is selected; otherwise ['C1'] etc.
  selectedPdbIds: string[]; // flattened list of PDB IDs the user wants Stage 2b/4/5a to operate on
}

export interface ClusterCard {
  id: string;
  size: number;
  representativePdbId: string;
  avgRmsd: number;
  ligandPdbIds: string[]; // members of the cluster
}

interface RmsdEntry { pdbId: string; rmsd: number; }

/**
 * Greedy 1D clustering on `rmsd_to_ref`. Structures are sorted by rmsd ascending;
 * a new cluster is seeded whenever the gap from the running cluster centroid
 * exceeds `threshold`. The first member of each cluster (lowest-rmsd) becomes
 * the representative.
 *
 * For the EGFR demo (5 active + 2 inactive at ~2 A): the active set lands near
 * rmsd 0-1.5 A, the inactive set sits at rmsd 4-6 A, so they split cleanly.
 */
export function clusterByRmsd(
  aligned: DG.DataFrame, threshold = RMSD_CLUSTER_THRESHOLD_A,
): ClusterCard[] {
  const pdbIdCol = aligned.col('pdb_id');
  const rmsdCol = aligned.col('rmsd_to_ref');
  if (!pdbIdCol || !rmsdCol) return [];

  const entries: RmsdEntry[] = [];
  for (let i = 0; i < aligned.rowCount; i++) {
    const rmsd = Number(rmsdCol.get(i));
    if (!Number.isFinite(rmsd)) continue;
    entries.push({pdbId: String(pdbIdCol.get(i)), rmsd});
  }
  if (entries.length === 0) return [];
  entries.sort((a, b) => a.rmsd - b.rmsd);

  // Gap-based 1D clustering — new cluster when current rmsd exceeds the running
  // centroid by more than `threshold`. The centroid drifts as new members join,
  // which is the standard "greedy single-link" behaviour.
  const buckets: RmsdEntry[][] = [[entries[0]]];
  let centroid = entries[0].rmsd;
  for (let i = 1; i < entries.length; i++) {
    const e = entries[i];
    if (e.rmsd - centroid <= threshold) {
      buckets[buckets.length - 1].push(e);
      const cur = buckets[buckets.length - 1];
      centroid = cur.reduce((s, x) => s + x.rmsd, 0) / cur.length;
    } else {
      buckets.push([e]);
      centroid = e.rmsd;
    }
  }

  return buckets.map((bucket, idx) => ({
    id: `C${idx + 1}`,
    size: bucket.length,
    representativePdbId: bucket[0].pdbId,
    avgRmsd: bucket.reduce((s, x) => s + x.rmsd, 0) / bucket.length,
    ligandPdbIds: bucket.map((x) => x.pdbId),
  }));
}

/** Default selection per blueprint Phase 3.5: "All" if ≤5 entries, else largest cluster. */
export function computeDefaultSelection(clusters: ClusterCard[]): string {
  const total = clusters.reduce((n, c) => n + c.size, 0);
  if (total <= 5) return 'All clusters';
  const largest = clusters.slice().sort((a, b) => b.size - a.size)[0];
  return `Cluster ${largest.id}`;
}

/** Build a small canvas drawing of a ligand from SMILES, or null on failure. */
function makeLigandThumb(smiles: string | undefined, w = 140, h = 80): HTMLElement | null {
  if (!smiles) return null;
  try {
    const canvas = ui.canvas(w, h);
    // grok.chem.canvasMol works synchronously after Chem's RDKit module is initialized
    // — which is true here because the orchestrator has already used it in Stage 1.
    (grok as any).chem.canvasMol(0, 0, w, h, canvas, smiles, '', {});
    return canvas;
  } catch {
    return null;
  }
}

/**
 * Show the cluster-picker dialog. Returns null when the user cancels.
 *
 *  - `aligned` — Stage 2a output with at least `pdb_id` and `rmsd_to_ref`.
 *  - `ligandSmilesByPdbId` — optional; when present, the cards render a 2D
 *    thumbnail of one representative ligand per cluster.
 *
 * Auto-skips (returns `{clusterIds: ['all']}` without showing the dialog) when
 * the input clusters into a single bucket — common for homogeneous demo sets.
 */
export async function pickClusterViaDialog(
  aligned: DG.DataFrame,
  ligandSmilesByPdbId?: Map<string, string>,
): Promise<ClusterPick | null> {
  const clusters = clusterByRmsd(aligned);
  const allPdbIds = clusters.flatMap((c) => c.ligandPdbIds);
  if (clusters.length <= 1) return {clusterIds: ['all'], selectedPdbIds: allPdbIds};

  const defaultChoice = computeDefaultSelection(clusters);
  const choices = ['All clusters', ...clusters.map((c) => `Cluster ${c.id}`)];

  // One card per cluster.
  const cards = clusters.map((c) => {
    const thumb = makeLigandThumb(ligandSmilesByPdbId?.get(c.representativePdbId));
    const memberPreview = c.ligandPdbIds.slice(0, 6).join(', ') +
      (c.ligandPdbIds.length > 6 ? ` (+${c.ligandPdbIds.length - 6})` : '');
    return ui.div([
      ui.h3(`Cluster ${c.id}`),
      ui.divText(`${c.size} structure(s) · avg RMSD ${c.avgRmsd.toFixed(2)} A`),
      ui.divText(`Representative: ${c.representativePdbId}`),
      ui.divText(`Members: ${memberPreview}`),
      thumb ?? ui.divText('(no SMILES)'),
    ], 'cp-cluster-card');
  });

  const choiceInput = ui.input.choice<string>('Selection', {value: defaultChoice, items: choices});

  return new Promise<ClusterPick | null>((resolve) => {
    let resolved = false;
    const dlg = ui.dialog({title: 'Conformational clusters detected'})
      .add(ui.divText('The aligned structures fall into multiple conformational groups ' +
        '(typical when active and inactive kinase states are mixed). Pick which to carry ' +
        'into pocket-level alignment and feature extraction.'))
      .add(ui.divH(cards, 'cp-cluster-cards-row'))
      .add(choiceInput.root)
      .onOK(() => {
        resolved = true;
        if (choiceInput.value === 'All clusters') {
          resolve({clusterIds: ['all'], selectedPdbIds: allPdbIds});
        } else {
          const cid = choiceInput.value!.replace(/^Cluster\s+/, '');
          const cluster = clusters.find((c) => c.id === cid);
          resolve({clusterIds: [cid], selectedPdbIds: cluster?.ligandPdbIds ?? []});
        }
      })
      .onCancel(() => { if (!resolved) resolve(null); });
    dlg.show();
  });
}

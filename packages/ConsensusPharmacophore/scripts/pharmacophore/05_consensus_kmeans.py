#name: PharmacophoreStage5aConsensusKmeans
#description: Stage 5a - per-family k-means consensus over raw ligand_features. Keeps the top N clusters per family that are supported by at least min_cluster_size_fraction of the input ligands.
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.11, numpy, pandas, scikit-learn]
#meta.cache: none
#input: dataframe ligand_features {caption: Stage 4 output}
#input: int kq = 7 {caption: Average features per cluster (sets n_clusters = ceil(n_features / kq))}
#input: int top_cluster_number = 4 {caption: Maximum clusters retained per family}
#input: double min_cluster_size_fraction = 0.75 {caption: Min n_distinct_ligands as fraction of total}
#input: int kmeans_random_state = 42 {caption: KMeans random_state for reproducibility}
#output: dataframe consensus_model {caption: family x y z frequency n_ligands}

# Adapted from TeachOpenCADD T009 (Sydow et al.), MIT-licensed.
# Constants pinned at module top so retuning is one line.

import math

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


KMEANS_N_INIT = 10


def n_total_ligands(df):
    """Distinct (pdb_id, ligand_comp_id) pairs in the feature table."""
    if 'pdb_id' not in df.columns or 'ligand_comp_id' not in df.columns:
        return 0
    return len(df[['pdb_id', 'ligand_comp_id']].drop_duplicates())


# --- Main body -----------------------------------------------------------

if not isinstance(ligand_features, pd.DataFrame) or len(ligand_features) == 0:
    raise ValueError('Stage 5a: ligand_features is empty.')

# Drop diagnostic rows (Stage 4 emits skip_reason on failures).
feat = ligand_features.copy()
if 'skip_reason' in feat.columns:
    feat = feat[feat['skip_reason'].fillna('') == '']
if 'family' in feat.columns:
    feat = feat[feat['family'].fillna('') != '']

if len(feat) == 0:
    print('Stage 5a: no usable feature rows after filtering diagnostics.')
    consensus_model = pd.DataFrame(columns=['family', 'x', 'y', 'z', 'frequency', 'n_ligands'])
else:
    total_ligands = n_total_ligands(feat)
    min_size = max(1, math.ceil(min_cluster_size_fraction * total_ligands))
    print(f'Stage 5a: {len(feat)} features across {total_ligands} ligand(s); ' +
          f'min cluster size = {min_size} (= ceil({min_cluster_size_fraction} * {total_ligands}))')

    out_rows = []
    for family, group in feat.groupby('family'):
        coords = group[['x', 'y', 'z']].astype(float).values
        n_features = len(coords)
        if n_features == 0:
            continue
        n_clusters = max(1, math.ceil(n_features / max(1, kq)))
        n_clusters = min(n_clusters, n_features)

        km = KMeans(n_clusters=n_clusters, n_init=KMEANS_N_INIT,
                    random_state=int(kmeans_random_state))
        labels = km.fit_predict(coords)

        # Per cluster: count distinct ligands (filter by min_size), record centroid + n_ligands.
        clusters = []
        for lbl in range(n_clusters):
            mask = (labels == lbl)
            if not mask.any():
                continue
            distinct = group.loc[mask, ['pdb_id', 'ligand_comp_id']].drop_duplicates()
            n_distinct = len(distinct)
            if n_distinct < min_size:
                continue
            centroid = km.cluster_centers_[lbl]
            # Cluster spread: max distance from centroid to any feature that
            # contributed to this cluster. Renderer uses this to size the
            # consensus sphere so the visualization shows BOTH the conserved
            # position AND the dispersion across input ligands. Clamped to
            # [0.8, 3.5] Å — under 0.8 Å is essentially a point (no visible
            # spread); over 3.5 Å is a "diffuse" cluster that probably should
            # have been split (a warning rather than a meaningful single point).
            cluster_coords = coords[mask]
            distances = np.linalg.norm(cluster_coords - centroid, axis=1)
            cluster_radius_a = float(np.max(distances))
            cluster_radius_a = max(0.8, min(3.5, cluster_radius_a))
            clusters.append({
                'family': family,
                'x': float(centroid[0]),
                'y': float(centroid[1]),
                'z': float(centroid[2]),
                'frequency': n_distinct / max(1, total_ligands),
                'n_ligands': int(n_distinct),
                'cluster_radius_a': cluster_radius_a,
                '_size': int(mask.sum()),
            })

        # Keep top `top_cluster_number` by n_ligands (then by total cluster size as tiebreaker).
        clusters.sort(key=lambda c: (-c['n_ligands'], -c['_size']))
        kept = clusters[:int(top_cluster_number)]
        for c in kept:
            c.pop('_size', None)
            out_rows.append(c)
        print(f'Stage 5a: family {family} — {n_features} features, ' +
              f'{n_clusters} k-means clusters, {len(kept)} retained')

    if out_rows:
        consensus_model = pd.DataFrame(out_rows)
    else:
        consensus_model = pd.DataFrame(columns=[
            'family', 'x', 'y', 'z', 'frequency', 'n_ligands', 'cluster_radius_a'])

print(f'Stage 5a: emitting {len(consensus_model)} consensus rows across ' +
      f'{consensus_model["family"].nunique() if len(consensus_model) else 0} families')

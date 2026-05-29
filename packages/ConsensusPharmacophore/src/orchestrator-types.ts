/*
 * Pipeline-wide option types. Lives in its own module so rcsb-client.ts,
 * future Python-wrapper helpers, and tests can import without pulling in the
 * orchestrator's runtime code.
 *
 * Blueprint reference - section 0 (PipelineOptions), Appendix B.
 */

export type PocketMethod = 'cutoff' | 'dbscan';
export type PocketRep = 'spacefill' | 'gaussian-surface';

export interface PipelineOptions {
  maxResolution: number;       // A; default 2.5
  /** DEPRECATED. Use allowXray / allowNmr / allowCryoEm instead. Kept on the
   *  interface so legacy code still type-checks; the rcsb-client filter does
   *  NOT consult this field anymore. */
  requireXray: boolean;
  // Experimental-method filters. By default we accept the three traditional
  // experimental methods (X-ray, NMR, cryo-EM). AlphaFold / predicted models
  // are off by default — RCSB's main entry index doesn't include them, so
  // toggling this on has no effect unless the user pastes a CSM identifier
  // and the API path supports it.
  allowXray: boolean;          // default true
  allowNmr: boolean;           // default true
  allowCryoEm: boolean;        // default true
  allowAlphaFold: boolean;     // default false
  minLigandMw: number;         // Da; default 100
  pocketMethod: PocketMethod;  // default 'cutoff'
  pocketRadius: number;        // A; default 5.0
  pocketRep: PocketRep;        // default 'spacefill' — viz only
  refPdbId?: string;           // optional reference overlay
  // Stage 5a consensus k-means knobs. Pinned defaults match the values that
  // were hard-coded in the Python script before they were exposed to the UI.
  kq: number;                  // Avg features per cluster (n_clusters = ceil(n_features/kq)); default 7
  minClusterSizeFraction: number; // 0.0-1.0; cluster kept only if >= this fraction of ligands contribute; default 0.75
  topClusterNumber: number;    // Max clusters retained per family; default 4
  // Consensus curation (Step 4 checkboxes). PDB ids (UPPERCASE) and individual
  // interaction keys the user excluded from the consensus. Empty = use
  // everything. Both are in the fingerprint so toggling re-runs Stage 5a.
  consensusExcludedPdbs?: string[];
  consensusExcludedInteractions?: string[];
  // Step 1 "Use" checkboxes — PDB ids (UPPERCASE) the user deselected from the
  // Accepted PDBs table. STRONGER than consensusExcludedPdbs: these are dropped
  // BEFORE Stage 2 (alignment), so they never enter pocket / features /
  // consensus and never appear in any later step's 3D view. In the fingerprint
  // so toggling invalidates Stage 2..5 caches on next run.
  excludedInputPdbs?: string[];
}

export const DEFAULT_OPTIONS: PipelineOptions = {
  maxResolution: 2.5,
  requireXray: false,          // deprecated, kept for type compatibility
  allowXray: true,
  allowNmr: true,
  allowCryoEm: true,
  allowAlphaFold: false,
  minLigandMw: 100,
  pocketMethod: 'cutoff',
  pocketRadius: 5.0,
  pocketRep: 'spacefill',
  kq: 7,
  minClusterSizeFraction: 0.75,
  topClusterNumber: 4,
  consensusExcludedPdbs: [],
  consensusExcludedInteractions: [],
  excludedInputPdbs: [],
};
